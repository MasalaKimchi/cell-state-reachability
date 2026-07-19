#!/usr/bin/env python3
"""Run the source-bound Schmidt CRISPRa/CRISPRi functional-screen stress test.

The benchmark is deliberately narrow. It measures whole-universe and source-selected
rank concordance across two fixed donors, two perturbation modalities (with different
guide libraries), and two cytokine/cell-type contexts. It does not test transcriptomic
reachability, guide-held-out generalization, or donor-population effects.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import stat
import statistics
import sys
from typing import Any, Iterable
import zipfile

import numpy as np
from scipy.stats import kendalltau, spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError  # noqa: E402


DEFAULT_CONFIG = ROOT / "configs" / "schmidt_external_validation.json"
DEFAULT_OUTPUT = ROOT / "results" / "schmidt_external_validation.json"

SGRNA_HEADER = (
    "sgrna",
    "Gene",
    "control_count",
    "treatment_count",
    "control_mean",
    "treat_mean",
    "LFC",
    "control_var",
    "adj_var",
    "score",
    "p.low",
    "p.high",
    "p.twosided",
    "FDR",
    "high_in_treatment",
)
NUMERIC_FIELDS = SGRNA_HEADER[2:-1]
SCREEN_IDS = (
    "CRISPRa_IFNG",
    "CRISPRa_IL2",
    "CRISPRi_IFNG",
    "CRISPRi_IL2",
)
TRANSFER_CLASSES = (
    "same_screen_held_donor",
    "donor_plus_modality_library_same_context",
    "donor_plus_cross_context_cytokine_cell_type_same_modality",
)


def _resolve(path: str | Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else ROOT / candidate


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _selected_hash(genes: Iterable[str]) -> str:
    return _sha256_bytes(("\n".join(genes) + "\n").encode("utf-8"))


def _file_digests(path: Path) -> tuple[str, str]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5()  # noqa: S324 - archival identity, not a security primitive
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            sha256.update(block)
            md5.update(block)
    return sha256.hexdigest(), md5.hexdigest()


def _exact_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    if set(value) != expected:
        raise InputError(
            f"{label} keys differ: missing={sorted(expected - set(value))}, "
            f"extra={sorted(set(value) - expected)}"
        )


def validate_config(config: dict[str, Any]) -> None:
    _exact_keys(
        config,
        {
            "schema_version",
            "report_version",
            "benchmark_id",
            "source",
            "input",
            "archive_members",
            "author_script",
            "screens",
            "analysis",
            "claim_contract",
        },
        "config",
    )
    if config["schema_version"] != 1 or config["report_version"] != 1:
        raise InputError("unsupported Schmidt config/report version")
    if not isinstance(config["benchmark_id"], str) or not config["benchmark_id"]:
        raise InputError("benchmark_id must be non-empty")
    _exact_keys(
        config["input"], {"path", "url", "bytes", "sha256", "md5"}, "input"
    )
    if int(config["input"]["bytes"]) <= 0:
        raise InputError("input byte length must be positive")
    for field, length in (("sha256", 64), ("md5", 32)):
        value = config["input"][field]
        if not isinstance(value, str) or len(value) != length:
            raise InputError(f"input {field} must be a {length}-character hex digest")
        try:
            int(value, 16)
        except ValueError as exc:
            raise InputError(f"input {field} is not hexadecimal") from exc

    members = config["archive_members"]
    if not isinstance(members, list) or not members or len(members) != len(set(members)):
        raise InputError("archive_members must be a non-empty duplicate-free list")
    _exact_keys(
        config["author_script"],
        {"member", "bytes", "sha256", "required_text"},
        "author_script",
    )
    if config["author_script"]["member"] not in members:
        raise InputError("author script is absent from archive member allow-list")

    if tuple(config["screens"]) != SCREEN_IDS:
        raise InputError("screen order/identities differ from the frozen four-screen design")
    combinations: set[tuple[str, str]] = set()
    for screen_id, spec in config["screens"].items():
        _exact_keys(
            spec,
            {"modality", "context", "cell_system", "member", "rows", "bytes", "sha256"},
            f"screen {screen_id}",
        )
        combinations.add((spec["modality"], spec["context"]))
        if spec["member"] not in members:
            raise InputError(f"{screen_id} member is absent from archive allow-list")
        if int(spec["rows"]) <= 0 or int(spec["bytes"]) <= 0:
            raise InputError(f"{screen_id} row/byte expectations must be positive")
    if combinations != {
        ("CRISPRa", "IFNG"),
        ("CRISPRa", "IL2"),
        ("CRISPRi", "IFNG"),
        ("CRISPRi", "IL2"),
    }:
        raise InputError("screen modality/context Cartesian product differs")

    analysis = config["analysis"]
    _exact_keys(
        analysis,
        {
            "donor_suffixes",
            "orientation",
            "excluded_gene",
            "eligibility",
            "aggregation",
            "minimum_guides_per_gene_per_donor",
            "primary_minimum_guides",
            "top_k",
            "primary_top_k",
            "expected_common_genes",
            "sensitivity_status",
        },
        "analysis",
    )
    if analysis["donor_suffixes"] != {"r0": "Donor1", "r1": "Donor2"}:
        raise InputError("donor suffix contract must remain r0/r1 -> Donor1/Donor2")
    if analysis["orientation"] != {"CRISPRa": 1, "CRISPRi": -1}:
        raise InputError("orientation must remain +CRISPRa/-CRISPRi")
    thresholds = [int(value) for value in analysis["minimum_guides_per_gene_per_donor"]]
    top_k = [int(value) for value in analysis["top_k"]]
    if thresholds != sorted(set(thresholds)) or not thresholds or thresholds[0] <= 0:
        raise InputError("minimum-guide thresholds must be sorted unique positive integers")
    if top_k != sorted(set(top_k)) or not top_k or top_k[0] < 3:
        raise InputError("top_k must be sorted unique integers of at least three")
    if int(analysis["primary_minimum_guides"]) not in thresholds:
        raise InputError("primary minimum-guide threshold is outside sensitivity grid")
    if int(analysis["primary_top_k"]) not in top_k:
        raise InputError("primary top-k is outside sensitivity grid")
    if set(analysis["expected_common_genes"]) != {str(value) for value in thresholds}:
        raise InputError("expected common-gene counts differ from threshold grid")
    if analysis["sensitivity_status"] != "EXPLORATORY_POST_HOC_DESCRIPTIVE":
        raise InputError("sensitivity must remain explicitly exploratory/post-hoc")
    if config["claim_contract"].get("tier") != "STRESS":
        raise InputError("Schmidt benchmark must remain STRESS tier")


def load_config(path: Path) -> tuple[dict[str, Any], str]:
    raw = path.read_bytes()
    try:
        config = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise InputError(f"invalid config JSON: {path}") from exc
    if not isinstance(config, dict):
        raise InputError("config must be a JSON object")
    validate_config(config)
    return config, _sha256_bytes(raw)


def verify_input(spec: dict[str, Any]) -> tuple[Path, dict[str, Any]]:
    """Verify the entire archive before opening it as ZIP."""

    path = _resolve(spec["path"])
    if not path.is_file():
        raise FileNotFoundError(path)
    size = path.stat().st_size
    if size != int(spec["bytes"]):
        raise InputError(f"archive byte length differs: {size} != {spec['bytes']}")
    sha256, md5 = _file_digests(path)
    if sha256 != spec["sha256"]:
        raise InputError("archive SHA-256 differs from the frozen identity")
    if md5 != spec["md5"]:
        raise InputError("archive MD5 differs from the frozen identity")
    return path, {
        "path": spec["path"],
        "url": spec["url"],
        "bytes": size,
        "sha256_expected": spec["sha256"],
        "sha256_actual": sha256,
        "md5_expected": spec["md5"],
        "md5_actual": md5,
        "hash_verified_before_zip_parse": True,
    }


def _validate_zip_members(zf: zipfile.ZipFile, allowed: list[str]) -> None:
    infos = zf.infolist()
    names = [info.filename for info in infos]
    if len(names) != len(set(names)):
        raise InputError("archive contains duplicate member names")
    for info in infos:
        name = info.filename
        path = PurePosixPath(name)
        if (
            not name
            or "\\" in name
            or path.is_absolute()
            or ".." in path.parts
            or (path.parts and ":" in path.parts[0])
        ):
            raise InputError(f"unsafe ZIP member path: {name!r}")
        if info.flag_bits & 0x1:
            raise InputError(f"encrypted ZIP member is not allowed: {name}")
        mode = info.external_attr >> 16
        if mode and stat.S_ISLNK(mode):
            raise InputError(f"symbolic-link ZIP member is not allowed: {name}")
    if names != allowed:
        missing = sorted(set(allowed) - set(names))
        extra = sorted(set(names) - set(allowed))
        raise InputError(f"archive member allow-list differs: missing={missing}, extra={extra}")


def _verified_member(zf: zipfile.ZipFile, spec: dict[str, Any]) -> bytes:
    try:
        info = zf.getinfo(spec["member"])
    except KeyError as exc:
        raise InputError(f"missing registered member {spec['member']}") from exc
    if info.is_dir():
        raise InputError(f"registered data member is a directory: {spec['member']}")
    if info.file_size != int(spec["bytes"]):
        raise InputError(f"registered member byte length differs: {spec['member']}")
    payload = zf.read(info)
    if len(payload) != info.file_size or _sha256_bytes(payload) != spec["sha256"]:
        raise InputError(f"registered member hash differs: {spec['member']}")
    return payload


def _parse_float(value: str, field: str, row_number: int, screen_id: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise InputError(f"{screen_id} row {row_number} has invalid {field}") from exc
    if not math.isfinite(parsed):
        raise InputError(f"{screen_id} row {row_number} has non-finite {field}")
    return parsed


def parse_screen(
    payload: bytes,
    screen_id: str,
    expected_rows: int,
    donor_suffixes: dict[str, str],
) -> dict[str, Any]:
    """Parse a complete MAGeCK sgRNA summary with paired guide identities."""

    try:
        text = payload.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise InputError(f"{screen_id} is not strict UTF-8") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    if tuple(reader.fieldnames or ()) != SGRNA_HEADER:
        raise InputError(f"{screen_id} sgRNA schema differs from the frozen 15 columns")

    seen_sgrna: set[str] = set()
    values: dict[str, dict[str, dict[str, float]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    rows = 0
    suffixes = set(donor_suffixes)
    for row_number, row in enumerate(reader, start=2):
        rows += 1
        if set(row) != set(SGRNA_HEADER) or any(value is None for value in row.values()):
            raise InputError(f"{screen_id} row {row_number} has malformed fields")
        sgrna = row["sgrna"]
        gene = row["Gene"]
        if not sgrna or not gene or sgrna in seen_sgrna:
            raise InputError(f"{screen_id} row {row_number} has empty/duplicate identity")
        seen_sgrna.add(sgrna)
        try:
            base_guide, donor = sgrna.rsplit("_", 1)
        except ValueError as exc:
            raise InputError(f"{screen_id} row {row_number} lacks donor suffix") from exc
        if donor not in suffixes or not base_guide:
            raise InputError(f"{screen_id} row {row_number} has malformed donor suffix")
        for field in NUMERIC_FIELDS:
            parsed = _parse_float(row[field], field, row_number, screen_id)
            if field == "LFC":
                lfc = parsed
        if row["high_in_treatment"] not in {"True", "False"}:
            raise InputError(f"{screen_id} row {row_number} has invalid Boolean field")
        if base_guide in values[gene][donor]:
            raise InputError(f"{screen_id} repeats a guide identity within donor")
        values[gene][donor][base_guide] = lfc

    if rows != int(expected_rows):
        raise InputError(f"{screen_id} row count differs: {rows} != {expected_rows}")
    effects: dict[str, dict[str, float]] = {}
    guide_counts: dict[str, dict[str, int]] = {}
    expected_donors = set(donor_suffixes)
    for gene, donor_values in values.items():
        if set(donor_values) != expected_donors:
            raise InputError(f"{screen_id} gene {gene} lacks one donor")
        guide_sets = [set(donor_values[donor]) for donor in donor_suffixes]
        if any(guides != guide_sets[0] for guides in guide_sets[1:]):
            raise InputError(f"{screen_id} gene {gene} has donor-unpaired guides")
        effects[gene] = {
            donor: float(statistics.median(donor_values[donor].values()))
            for donor in donor_suffixes
        }
        guide_counts[gene] = {
            donor: len(donor_values[donor]) for donor in donor_suffixes
        }
    return {
        "rows": rows,
        "genes_including_no_target": len(effects),
        "effects": effects,
        "guide_counts": guide_counts,
    }


def eligible_genes(
    screens: dict[str, dict[str, Any]],
    threshold: int,
    donors: Iterable[str],
    excluded_gene: str,
) -> list[str]:
    """Build the universe from identities and guide coverage only."""

    common = set.intersection(*(set(screen["guide_counts"]) for screen in screens.values()))
    common.discard(excluded_gene)
    return sorted(
        gene
        for gene in common
        if all(
            screen["guide_counts"][gene][donor] >= threshold
            for screen in screens.values()
            for donor in donors
        )
    )


def _vector(
    screens: dict[str, dict[str, Any]],
    specs: dict[str, dict[str, Any]],
    orientation: dict[str, int],
    screen_id: str,
    donor: str,
    genes: list[str],
) -> np.ndarray:
    sign = float(orientation[specs[screen_id]["modality"]])
    return sign * np.asarray(
        [screens[screen_id]["effects"][gene][donor] for gene in genes], dtype=float
    )


def _donor_mean_vector(
    screens: dict[str, dict[str, Any]],
    specs: dict[str, dict[str, Any]],
    orientation: dict[str, int],
    screen_id: str,
    donors: list[str],
    genes: list[str],
) -> np.ndarray:
    return np.mean(
        np.vstack(
            [
                _vector(screens, specs, orientation, screen_id, donor, genes)
                for donor in donors
            ]
        ),
        axis=0,
    )


def metric_bundle(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    if x.shape != y.shape or x.ndim != 1 or x.size < 3:
        raise InputError("rank metrics require aligned vectors with at least three values")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise InputError("rank metrics received non-finite values")
    x_norm = float(np.linalg.norm(x))
    y_norm = float(np.linalg.norm(y))
    if x_norm == 0.0 or y_norm == 0.0:
        raise InputError("cosine is undefined for an all-zero screen vector")
    values = {
        "signed_spearman": float(spearmanr(x, y).statistic),
        "signed_kendall": float(kendalltau(x, y).statistic),
        "sign_agreement": float(np.mean(np.sign(x) == np.sign(y))),
        "absolute_effect_spearman": float(spearmanr(np.abs(x), np.abs(y)).statistic),
        "cosine": float(np.dot(x, y) / (x_norm * y_norm)),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise InputError("a concordance metric is undefined or non-finite")
    return values


def select_top_k(genes: list[str], effects: np.ndarray, top_k: int) -> list[str]:
    """Select by source-only absolute effect; break ties by gene identifier."""

    if len(genes) != effects.size or top_k > len(genes):
        raise InputError("top-k selection is incompatible with the frozen universe")
    pairs = sorted(zip(genes, effects, strict=True), key=lambda item: (-abs(item[1]), item[0]))
    return [gene for gene, _ in pairs[:top_k]]


def _screen_id(specs: dict[str, dict[str, Any]], modality: str, context: str) -> str:
    matches = [
        screen_id
        for screen_id, spec in specs.items()
        if spec["modality"] == modality and spec["context"] == context
    ]
    if len(matches) != 1:
        raise InputError("screen modality/context lookup is not one-to-one")
    return matches[0]


def _target_screen(
    specs: dict[str, dict[str, Any]], source_screen: str, transfer_class: str
) -> str:
    spec = specs[source_screen]
    if transfer_class == "same_screen_held_donor":
        return source_screen
    if transfer_class == "donor_plus_modality_library_same_context":
        other_modality = "CRISPRi" if spec["modality"] == "CRISPRa" else "CRISPRa"
        return _screen_id(specs, other_modality, spec["context"])
    if transfer_class == "donor_plus_cross_context_cytokine_cell_type_same_modality":
        other_context = "IL2" if spec["context"] == "IFNG" else "IFNG"
        return _screen_id(specs, spec["modality"], other_context)
    raise InputError(f"unsupported transfer class {transfer_class}")


def _whole_universe(
    screens: dict[str, dict[str, Any]],
    specs: dict[str, dict[str, Any]],
    orientation: dict[str, int],
    donors: list[str],
    genes: list[str],
) -> dict[str, list[dict[str, Any]]]:
    donor_rows = []
    for screen_id in SCREEN_IDS:
        donor_rows.append(
            {
                "screen": screen_id,
                "n_genes": len(genes),
                "donor_a": donors[0],
                "donor_b": donors[1],
                "metrics": metric_bundle(
                    _vector(screens, specs, orientation, screen_id, donors[0], genes),
                    _vector(screens, specs, orientation, screen_id, donors[1], genes),
                ),
            }
        )
    modality_rows = []
    for context in ("IFNG", "IL2"):
        screen_a = _screen_id(specs, "CRISPRa", context)
        screen_i = _screen_id(specs, "CRISPRi", context)
        modality_rows.append(
            {
                "context": context,
                "screen_a": screen_a,
                "screen_b": screen_i,
                "n_genes": len(genes),
                "metrics": metric_bundle(
                    _donor_mean_vector(screens, specs, orientation, screen_a, donors, genes),
                    _donor_mean_vector(screens, specs, orientation, screen_i, donors, genes),
                ),
            }
        )
    context_rows = []
    for modality in ("CRISPRa", "CRISPRi"):
        screen_ifng = _screen_id(specs, modality, "IFNG")
        screen_il2 = _screen_id(specs, modality, "IL2")
        context_rows.append(
            {
                "modality": modality,
                "screen_a": screen_ifng,
                "screen_b": screen_il2,
                "n_genes": len(genes),
                "metrics": metric_bundle(
                    _donor_mean_vector(screens, specs, orientation, screen_ifng, donors, genes),
                    _donor_mean_vector(screens, specs, orientation, screen_il2, donors, genes),
                ),
            }
        )
    return {
        "donor_same_reagent": donor_rows,
        "modality_plus_library_same_context": modality_rows,
        "cross_context_cytokine_plus_cell_type": context_rows,
    }


def _transfer_rows(
    screens: dict[str, dict[str, Any]],
    specs: dict[str, dict[str, Any]],
    orientation: dict[str, int],
    donors: list[str],
    universes: dict[int, list[str]],
    thresholds: list[int],
    top_values: list[int],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for threshold in thresholds:
        genes = universes[threshold]
        for top_k in top_values:
            if top_k > len(genes):
                raise InputError(f"top-{top_k} exceeds the threshold-{threshold} universe")
            for source_screen in SCREEN_IDS:
                for training_donor in donors:
                    held_donor = donors[1] if training_donor == donors[0] else donors[0]
                    source_all = _vector(
                        screens,
                        specs,
                        orientation,
                        source_screen,
                        training_donor,
                        genes,
                    )
                    selected = select_top_k(genes, source_all, top_k)
                    selected_hash = _selected_hash(selected)
                    source_selected = _vector(
                        screens,
                        specs,
                        orientation,
                        source_screen,
                        training_donor,
                        selected,
                    )
                    for transfer_class in TRANSFER_CLASSES:
                        target_screen = _target_screen(specs, source_screen, transfer_class)
                        target_selected = _vector(
                            screens,
                            specs,
                            orientation,
                            target_screen,
                            held_donor,
                            selected,
                        )
                        target_all = _vector(
                            screens,
                            specs,
                            orientation,
                            target_screen,
                            held_donor,
                            genes,
                        )
                        target_top = select_top_k(genes, target_all, top_k)
                        overlap = len(set(selected) & set(target_top))
                        rows.append(
                            {
                                "minimum_guides": threshold,
                                "top_k": top_k,
                                "transfer_class": transfer_class,
                                "source_screen": source_screen,
                                "target_screen": target_screen,
                                "training_donor": training_donor,
                                "held_donor": held_donor,
                                "universe_genes": len(genes),
                                "source_selected_gene_sha256": selected_hash,
                                "held_target_global_top_k_sha256": _selected_hash(target_top),
                                "held_target_global_top_k_overlap_count": overlap,
                                "held_target_global_top_k_overlap_fraction": overlap / top_k,
                                "metrics": metric_bundle(source_selected, target_selected),
                            }
                        )
    return rows


def _summary(records: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    metric_names = tuple(next(iter(records))["metrics"]) + (
        "held_target_global_top_k_overlap_fraction",
    )
    for name in metric_names:
        values = [
            float(row[name] if name in row else row["metrics"][name]) for row in records
        ]
        output[name] = {
            "median": float(statistics.median(values)),
            "minimum": min(values),
            "maximum": max(values),
        }
    return output


def run(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config, config_sha256 = load_config(config_path)
    archive_path, input_identity = verify_input(config["input"])
    donors = list(config["analysis"]["donor_suffixes"])
    screen_data: dict[str, dict[str, Any]] = {}
    member_identity: dict[str, dict[str, Any]] = {}

    try:
        with zipfile.ZipFile(archive_path) as zf:
            _validate_zip_members(zf, config["archive_members"])
            script_payload = _verified_member(zf, config["author_script"])
            script_text = script_payload.decode("utf-8", errors="strict")
            for required in config["author_script"]["required_text"]:
                if required not in script_text:
                    raise InputError(f"author transformation text is missing: {required!r}")
            member_identity["author_script"] = {
                "member": config["author_script"]["member"],
                "bytes": len(script_payload),
                "sha256": _sha256_bytes(script_payload),
                "required_text_verified": True,
            }
            for screen_id, spec in config["screens"].items():
                payload = _verified_member(zf, spec)
                screen_data[screen_id] = parse_screen(
                    payload,
                    screen_id,
                    int(spec["rows"]),
                    config["analysis"]["donor_suffixes"],
                )
                member_identity[screen_id] = {
                    "member": spec["member"],
                    "rows": screen_data[screen_id]["rows"],
                    "bytes": len(payload),
                    "sha256": _sha256_bytes(payload),
                    "schema_verified": list(SGRNA_HEADER),
                }
    except (zipfile.BadZipFile, UnicodeDecodeError) as exc:
        raise InputError("registered Schmidt archive cannot be parsed safely") from exc

    thresholds = [
        int(value)
        for value in config["analysis"]["minimum_guides_per_gene_per_donor"]
    ]
    universes = {
        threshold: eligible_genes(
            screen_data,
            threshold,
            donors,
            config["analysis"]["excluded_gene"],
        )
        for threshold in thresholds
    }
    for threshold, genes in universes.items():
        expected = int(config["analysis"]["expected_common_genes"][str(threshold)])
        if len(genes) != expected:
            raise InputError(
                f"threshold-{threshold} universe differs: {len(genes)} != {expected}"
            )
    counts = [len(universes[threshold]) for threshold in thresholds]
    if any(later > earlier for earlier, later in zip(counts, counts[1:])):
        raise InputError("common-gene counts increase with a stricter guide threshold")

    orientation = config["analysis"]["orientation"]
    whole_by_threshold = {
        str(threshold): _whole_universe(
            screen_data,
            config["screens"],
            orientation,
            donors,
            universes[threshold],
        )
        for threshold in thresholds
    }
    transfer_rows = _transfer_rows(
        screen_data,
        config["screens"],
        orientation,
        donors,
        universes,
        thresholds,
        [int(value) for value in config["analysis"]["top_k"]],
    )
    primary_threshold = int(config["analysis"]["primary_minimum_guides"])
    primary_top_k = int(config["analysis"]["primary_top_k"])
    primary_summary = {}
    for transfer_class in TRANSFER_CLASSES:
        selected_rows = [
            row
            for row in transfer_rows
            if row["minimum_guides"] == primary_threshold
            and row["top_k"] == primary_top_k
            and row["transfer_class"] == transfer_class
        ]
        if len(selected_rows) != 8:
            raise InputError(f"{transfer_class} does not contain eight donor directions")
        primary_summary[transfer_class] = _summary(selected_rows)

    return {
        "schema_version": 1,
        "report_version": int(config["report_version"]),
        "benchmark_id": config["benchmark_id"],
        "execution_status": "PASS",
        "config_sha256": config_sha256,
        "provenance": {
            "source": config["source"],
            "input": input_identity,
            "archive_member_allowlist_count": len(config["archive_members"]),
            "archive_member_allowlist_verified_exactly": True,
            "members": member_identity,
            "donor_suffixes": config["analysis"]["donor_suffixes"],
            "same_donors_across_modalities": True,
            "same_donors_basis": config["source"]["same_donors_basis"],
        },
        "analysis_contract": {
            "orientation": orientation,
            "aggregation": config["analysis"]["aggregation"],
            "eligibility": config["analysis"]["eligibility"],
            "excluded_gene": config["analysis"]["excluded_gene"],
            "minimum_guides_per_gene_per_donor": thresholds,
            "primary_minimum_guides": primary_threshold,
            "top_k": [int(value) for value in config["analysis"]["top_k"]],
            "primary_top_k": primary_top_k,
            "sensitivity_status": config["analysis"]["sensitivity_status"],
            "selection": "Top-k genes are selected only by absolute oriented effect in the source screen and training donor; all target screens and the held donor remain unread until the gene set is frozen.",
            "target_overlap": "Held-target top-k is selected globally over the complete frozen universe, never within the source-selected subset.",
        },
        "eligibility": {
            "common_gene_counts_by_minimum_guides": {
                str(threshold): len(universes[threshold]) for threshold in thresholds
            },
            "outcome_fields_used": [],
            "identity_and_guide_coverage_only": True,
        },
        "whole_universe_by_minimum_guides": whole_by_threshold,
        "source_selected_transfer": {
            "rows": transfer_rows,
            "primary_summary": primary_summary,
            "primary_directions_per_class": 8,
            "correlated_descriptive_challenges": True,
        },
        "claim_contract": config["claim_contract"],
        "interpretation": {
            "status": "CONDITIONAL_TOP_EFFECT_CONCORDANCE_ONLY",
            "text": "Whole-universe reproducibility is limited. Conditional on source-donor extreme-effect selection, same-screen held-donor concordance is strongest, joint donor-plus-context concordance is intermediate and confounded, and joint donor-plus-modality/library concordance is weaker. These results are descriptive and do not establish scale, guide, modality, cell-type, cytokine, donor-population, or biological-state generality.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true", help="write the canonical report")
    group.add_argument("--check", action="store_true", help="byte-check the canonical report")
    args = parser.parse_args()

    report = run(args.config)
    rendered = _canonical_json(report)
    output = _resolve(args.output)
    if args.write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
        print(f"wrote {output.relative_to(ROOT) if output.is_relative_to(ROOT) else output}")
        return 0
    if args.check:
        if not output.is_file():
            raise FileNotFoundError(output)
        if output.read_text(encoding="utf-8") != rendered:
            raise AssertionError(f"frozen report differs: {output}")
        print(f"PASS: {output.relative_to(ROOT) if output.is_relative_to(ROOT) else output}")
        return 0
    sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
