#!/usr/bin/env python3
"""Run the source-bound Goudy GSE306915 cross-experiment stress test.

This is a deliberately negative benchmark.  It asks whether one measured
FAS+RC3H1+SUV39H1 triple is represented by target-matched constituent singles
under a declared cone model.  Single/triple status is confounded with
experiment, control type, and guide burden, and the triple guide identities are
unresolved.  The output therefore cannot validate additivity or interaction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re
import statistics
import sys
from typing import Any, Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError, ProjectionResult, project_cone  # noqa: E402


DEFAULT_CONFIG = ROOT / "configs" / "goudy_combination_validation.json"


def _threshold_key(value: float) -> str:
    return f"{float(value):g}"


def _resolve(path: str | Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else ROOT / candidate


def _file_digests(path: Path) -> tuple[str, str]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5()  # noqa: S324 - archival identity, never a security primitive
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            sha256.update(block)
            md5.update(block)
    return sha256.hexdigest(), md5.hexdigest()


def verify_input(spec: dict[str, Any]) -> dict[str, Any]:
    """Verify the complete object before any decompression or parsing."""

    path = _resolve(spec["path"])
    if not path.is_file():
        raise FileNotFoundError(path)
    size = path.stat().st_size
    if size != int(spec["bytes"]):
        raise InputError(f"{spec['path']} byte length differs: {size} != {spec['bytes']}")
    sha256, md5 = _file_digests(path)
    if sha256 != spec["sha256"]:
        raise InputError(f"{spec['path']} SHA-256 differs from the frozen identity")
    if md5 != spec["md5"]:
        raise InputError(f"{spec['path']} MD5 differs from the frozen identity")
    identity = {
        "path": spec["path"],
        "file": spec.get("file", path.name),
        "url": spec["url"],
        "bytes": size,
        "sha256_expected": spec["sha256"],
        "sha256_actual": sha256,
        "md5_expected": spec["md5"],
        "md5_actual": md5,
        "hash_verified_before_parse": True,
    }
    if "commit" in spec:
        identity["commit"] = spec["commit"]
    return identity


def _iter_specs(
    config: dict[str, Any],
) -> Iterable[tuple[str, str, str, dict[str, str]]]:
    for donor, donor_spec in config["donors"].items():
        for group in ("component", "unrelated", "controls"):
            for role, spec in donor_spec[group].items():
                yield donor, group, role, spec
        yield donor, "triple", "triple", donor_spec["triple"]


def validate_config(config: dict[str, Any]) -> None:
    if set(config["inputs"]) != {"counts", "soft", "author_key"}:
        raise InputError("counts, SOFT, and immutable author key are required")
    expected = config["expected"]
    donors = tuple(config["donors"])
    if len(donors) != int(expected["donors"]) or len(set(donors)) != len(donors):
        raise InputError("donor declarations differ from the frozen contract")
    component_order = tuple(config["models"]["component_order"])
    unrelated_order = tuple(config["models"]["unrelated_order"])
    if len(component_order) != 3 or len(set(component_order)) != 3:
        raise InputError("exactly three unique component roles are required")
    if len(unrelated_order) != 3 or len(set(unrelated_order)) != 3:
        raise InputError("exactly three unique unrelated roles are required")
    controls_per_donor = int(expected["controls_per_donor"])
    control_roles = ("AAVS1_Guide1", "AAVS1_Guide2", "NTC_multiplexing")
    gsms: list[str] = []
    columns: list[str] = []
    for donor, donor_spec in config["donors"].items():
        if tuple(donor_spec["component"]) != component_order:
            raise InputError(f"{donor} component order differs from the frozen order")
        if tuple(donor_spec["unrelated"]) != unrelated_order:
            raise InputError(f"{donor} unrelated order differs from the frozen order")
        if (
            len(donor_spec["controls"]) != controls_per_donor
            or tuple(donor_spec["controls"]) != control_roles
        ):
            raise InputError(f"{donor} requires {controls_per_donor} controls")
    for _, _, _, spec in _iter_specs(config):
        required = {"gsm", "column", "title", "genotype", "time"}
        if set(spec) != required:
            raise InputError("every sample role requires exact GSM/column/title/genotype/time")
        if not re.fullmatch(r"GSM[0-9]+", spec["gsm"]):
            raise InputError("invalid GSM accession in sample contract")
        if not re.fullmatch(r"sample_[A-Za-z0-9_]+", spec["column"]):
            raise InputError("invalid count-column name in sample contract")
        gsms.append(spec["gsm"])
        columns.append(spec["column"])
    if len(gsms) != int(expected["declared_analysis_samples"]):
        raise InputError("declared analysis-sample count differs")
    if len(set(gsms)) != len(gsms) or len(set(columns)) != len(columns):
        raise InputError("GSM accessions and count columns must map one-to-one")
    target_mask = config["preprocessing"]["target_mask"]
    if tuple(item["symbol"] for item in target_mask) != component_order:
        raise InputError("target-mask symbols must equal the component order")
    if len({item["gene_id"] for item in target_mask}) != len(target_mask):
        raise InputError("target-mask gene IDs must be unique")
    thresholds = tuple(
        float(value)
        for value in config["preprocessing"][
            "control_mean_cpm_sensitivity_thresholds"
        ]
    )
    if thresholds != (0.5, 1.0, 2.0, 5.0, 10.0):
        raise InputError("control-only expression-filter grid differs")
    if float(config["preprocessing"]["control_mean_cpm_min"]) not in thresholds:
        raise InputError("canonical control threshold must occur in the sensitivity grid")
    expected_filter_counts = expected["analysis_genes_before_mask_by_threshold"]
    if set(expected_filter_counts) != {_threshold_key(value) for value in thresholds}:
        raise InputError("expected expression-filter counts differ from the threshold grid")
    if config["claim_contract"].get("constituent_relationship") != (
        "target-matched constituent singles"
    ):
        raise InputError("constituent relationship must remain target-matched only")
    if config["claim_contract"].get("triple_guide_identity") != "UNRESOLVED":
        raise InputError("triple guide identity must remain unresolved")
    conflict = config["metadata_conflict"]
    if conflict["status"] != "FLAGGED" or conflict["computational_resolution"] != "Day7":
        raise InputError("the frozen donor-1 time conflict must remain explicitly FLAGGED")


AUTHOR_KEY_HEADER = (
    "sample",
    "Experiment ID",
    "Experimental Day Collected",
    "Donor",
    "Cas",
    "Target",
    "sgRNA notes",
)


def parse_author_key(path: Path) -> dict[str, Any]:
    """Parse the immutable author crosswalk after its full-byte identity gate."""

    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != AUTHOR_KEY_HEADER:
            raise InputError("author key schema differs from the frozen seven columns")
        for row_number, row in enumerate(reader, start=2):
            if set(row) != set(AUTHOR_KEY_HEADER) or any(value is None for value in row.values()):
                raise InputError(f"author key row {row_number} has malformed fields")
            sample = row["sample"]
            if re.fullmatch(r"sample_[0-9]+", sample) is None:
                raise InputError(f"author key row {row_number} has invalid sample ID")
            if sample in rows:
                raise InputError(f"author key repeats {sample}")
            rows[sample] = {key: str(value) for key, value in row.items()}
    return {"header": AUTHOR_KEY_HEADER, "rows": rows}


def _expected_author_target(
    group: str, role: str, component_order: tuple[str, ...]
) -> set[str]:
    if group in {"component", "unrelated"}:
        return {role}
    if group == "controls" and role.startswith("AAVS1_"):
        return {"AAVS1"}
    if group == "controls" and role == "NTC_multiplexing":
        return {"NTC"}
    if group == "triple":
        return set(component_order)
    raise InputError(f"unsupported author-key role {group}:{role}")


def _target_tokens(value: str) -> set[str]:
    return {token.strip() for token in value.split(",") if token.strip()}


def validate_author_key(
    author_key: dict[str, Any], counts: dict[str, Any], config: dict[str, Any]
) -> dict[str, Any]:
    """Crosswalk declared roles without inventing metadata for absent rows."""

    rows = author_key["rows"]
    expected = config["expected"]
    if len(rows) != int(expected["author_key_samples"]):
        raise InputError("author-key sample count differs from the frozen contract")
    count_columns = set(counts["count_columns"])
    overlap = set(rows) & count_columns
    if len(overlap) != int(expected["author_key_count_overlap"]):
        raise InputError("author-key/count overlap differs from the frozen contract")
    author_only = sorted(set(rows) - count_columns)
    if author_only != sorted(expected["author_key_rows_absent_from_counts"]):
        raise InputError("author-key rows absent from counts differ")

    component_order = tuple(config["models"]["component_order"])
    declared_columns = {spec["column"] for _, _, _, spec in _iter_specs(config)}
    missing_declared = sorted(declared_columns - set(rows))
    if missing_declared != sorted(expected["declared_rows_absent_from_author_key"]):
        raise InputError("declared rows absent from author key differ")

    declared_roles: dict[str, dict[str, Any]] = {}
    single_experiment_ids: set[str] = set()
    multiplex_experiment_ids: set[str] = set()
    component_guide_notes: dict[str, set[str]] = {
        name: set() for name in component_order
    }
    unresolved_experiment_columns: list[str] = []
    for donor, group, role, spec in _iter_specs(config):
        column = spec["column"]
        key = f"{group}:{role}"
        if column not in rows:
            unresolved_experiment_columns.append(column)
            declared_roles.setdefault(donor, {})[key] = {
                "column": column,
                "author_key_status": "ABSENT_FROM_AUTHOR_KEY",
                "experiment_id": "UNRESOLVED",
                "experimental_day": "UNRESOLVED",
                "guide_identity": "UNRESOLVED",
            }
            continue
        row = rows[column]
        donor_number = donor.removeprefix("D")
        if row["Donor"] != f"Donor {donor_number}":
            raise InputError(f"{column} author-key donor differs")
        if row["Cas"] != config["claim_contract"]["modality"]:
            raise InputError(f"{column} author-key Cas modality differs")
        day_match = re.search(
            r"([0-9]+)", config["metadata_conflict"]["computational_resolution"]
        )
        if day_match is None or row["Experimental Day Collected"] != day_match.group(1):
            raise InputError(f"{column} author-key day differs")
        expected_target = _expected_author_target(group, role, component_order)
        if _target_tokens(row["Target"]) != expected_target:
            raise InputError(f"{column} author-key target differs")
        experiment_id = row["Experiment ID"]
        if not experiment_id:
            raise InputError(f"{column} author-key experiment ID is blank")
        if group in {"component", "unrelated"} or (
            group == "controls" and role.startswith("AAVS1_")
        ):
            single_experiment_ids.add(experiment_id)
        else:
            multiplex_experiment_ids.add(experiment_id)
        guide_notes = row["sgRNA notes"]
        if group == "component":
            if not guide_notes:
                raise InputError(f"{column} constituent single guide note is blank")
            component_guide_notes[role].add(guide_notes)
        declared_roles.setdefault(donor, {})[key] = {
            "column": column,
            "author_key_status": "CONFIRMED",
            "experiment_id": experiment_id,
            "experimental_day": int(row["Experimental Day Collected"]),
            "target": row["Target"],
            "guide_identity": guide_notes or "UNRESOLVED",
        }
    if len(single_experiment_ids) != 1 or len(multiplex_experiment_ids) != 1:
        raise InputError("author key does not define one experiment per declared role family")
    if single_experiment_ids == multiplex_experiment_ids:
        raise InputError("single and multiplex role families must remain cross-experiment")
    if any(len(notes) != 1 for notes in component_guide_notes.values()):
        raise InputError("constituent-single guide notes differ across donors")
    cas_row_counts = {
        name: int(sum(row["Cas"] == name for row in rows.values()))
        for name in sorted({row["Cas"] for row in rows.values()})
    }
    declared_role_count = sum(
        len(donor_roles) for donor_roles in declared_roles.values()
    )
    confirmed_role_count = declared_role_count - len(unresolved_experiment_columns)
    return {
        "author_key_row_count": len(rows),
        "cas_modality_row_counts": cas_row_counts,
        "count_matrix_overlap_count": len(overlap),
        "declared_analysis_role_count": declared_role_count,
        "author_key_confirmed_declared_role_count": confirmed_role_count,
        "declared_roles": declared_roles,
        "experiment_role_summary": {
            "constituent_singles_unrelated_singles_and_aavs1": {
                "experiment_ids": sorted(single_experiment_ids),
                "status": "AUTHOR_KEY_CONFIRMED_ALL_DONORS",
            },
            "multiplex_ntc_and_triple": {
                "experiment_ids": sorted(multiplex_experiment_ids),
                "author_key_confirmed_donors": ["D1", "D2"],
                "unresolved_donors": ["D3", "D4"],
            },
        },
        "author_key_rows_absent_from_count_matrix": author_only,
        "declared_rows_absent_from_author_key": missing_declared,
        "unresolved_experiment_id_columns": sorted(unresolved_experiment_columns),
        "constituent_single_guide_notes": {
            name: next(iter(notes)) for name, notes in component_guide_notes.items()
        },
        "triple_guide_identity": "UNRESOLVED",
        "constituent_relationship": "target-matched constituent singles",
        "same_guide_match_claimed": False,
    }


def parse_soft(path: Path) -> dict[str, Any]:
    """Parse only the bounded GEO SOFT fields needed for identity validation."""

    samples: dict[str, dict[str, Any]] = {}
    current: dict[str, Any] | None = None
    series_accession: str | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SERIES = "):
                series_accession = line.split(" = ", 1)[1]
            elif line.startswith("^SAMPLE = "):
                gsm = line.split(" = ", 1)[1]
                if gsm in samples:
                    raise InputError(f"SOFT contains duplicate sample block {gsm}")
                current = {
                    "gsm": gsm,
                    "title": None,
                    "geo_accession": None,
                    "characteristics": {},
                    "descriptions": [],
                }
                samples[gsm] = current
            elif current is not None and line.startswith("!Sample_") and " = " in line:
                key, value = line.split(" = ", 1)
                key = key.removeprefix("!Sample_")
                if key == "title":
                    current["title"] = value
                elif key == "geo_accession":
                    current["geo_accession"] = value
                elif key == "characteristics_ch1":
                    if ": " not in value:
                        raise InputError(f"{current['gsm']} has malformed characteristic")
                    field, field_value = value.split(": ", 1)
                    if field in current["characteristics"]:
                        raise InputError(
                            f"{current['gsm']} repeats characteristic {field}"
                        )
                    current["characteristics"][field] = field_value
                elif key == "description":
                    current["descriptions"].append(value)
    return {"series_accession": series_accession, "samples": samples}


def _description_count_filename(config: dict[str, Any]) -> str:
    filename = str(config["inputs"]["counts"]["file"])
    accession_prefix = f"{config['source']['accession']}_"
    if filename.startswith(accession_prefix):
        filename = filename[len(accession_prefix) :]
    return filename.removesuffix(".gz")


def validate_soft(soft: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    if soft["series_accession"] != config["source"]["accession"]:
        raise InputError("SOFT series accession differs from the frozen contract")
    samples = soft["samples"]
    if len(samples) != int(config["expected"]["soft_samples"]):
        raise InputError("SOFT sample count differs from the frozen contract")
    described_file = _description_count_filename(config)
    role_report: dict[str, dict[str, Any]] = {}
    observed_conflicts: list[str] = []
    for donor, group, role, spec in _iter_specs(config):
        gsm = spec["gsm"]
        if gsm not in samples:
            raise InputError(f"SOFT is missing declared sample {gsm}")
        sample = samples[gsm]
        if sample["geo_accession"] != gsm:
            raise InputError(f"{gsm} GEO accession field differs")
        if sample["title"] != spec["title"]:
            raise InputError(f"{gsm} title differs from the frozen sample role")
        characteristics = sample["characteristics"]
        required_characteristics = {
            "cell type": "Bulk T cells",
            "genotype": spec["genotype"],
            "treatment": "CRISPRoff",
            "time": spec["time"],
        }
        for field, expected_value in required_characteristics.items():
            if characteristics.get(field) != expected_value:
                raise InputError(f"{gsm} characteristic {field!r} differs")
        library_description = f"Library name: {spec['column'].removeprefix('sample_')}"
        if library_description not in sample["descriptions"]:
            raise InputError(f"{gsm} library suffix does not match its count column")
        if described_file not in sample["descriptions"]:
            raise InputError(f"{gsm} does not declare the frozen count file")
        title_day_match = re.search(r"Day\s*([0-9]+)", spec["title"], flags=re.I)
        characteristic_day_match = re.fullmatch(r"Day\s*([0-9]+)", spec["time"], flags=re.I)
        if (
            title_day_match is not None
            and characteristic_day_match is not None
            and title_day_match.group(1) != characteristic_day_match.group(1)
        ):
            observed_conflicts.append(gsm)
        role_report.setdefault(donor, {})[f"{group}:{role}"] = {
            "gsm": gsm,
            "column": spec["column"],
            "title": spec["title"],
            "characteristic_time": spec["time"],
        }
    expected_conflicts = sorted(config["metadata_conflict"]["affected_gsms"])
    if sorted(observed_conflicts) != expected_conflicts:
        raise InputError("observed title/characteristic metadata conflict set differs")
    return {
        "series_accession": soft["series_accession"],
        "soft_sample_count": len(samples),
        "declared_roles": role_report,
        "observed_time_conflict_gsms": sorted(observed_conflicts),
    }


def load_counts(path: Path, config: dict[str, Any]) -> dict[str, Any]:
    """Load the hash-verified gzip CSV on its deposited linear-CPM scale."""

    declared_columns = [spec["column"] for _, _, _, spec in _iter_specs(config)]
    expected_rows = int(config["expected"]["gene_rows"])
    gene_pattern = re.compile(config["preprocessing"]["gene_id_regex"])
    csv.field_size_limit(1024 * 1024)
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise InputError("count table is empty") from exc
        if len(header) != len(set(header)):
            raise InputError("count-table columns must be unique")
        if len(header) < 3 or header[0] != "gene_id" or header[-1] != "gene_symbol":
            raise InputError("count-table schema must be gene_id, sample columns, gene_symbol")
        count_columns = header[1:-1]
        if len(count_columns) != int(config["expected"]["count_samples"]):
            raise InputError("count-sample column count differs from the frozen contract")
        if any(re.fullmatch(r"sample_[A-Za-z0-9_]+", name) is None for name in count_columns):
            raise InputError("count-table sample columns have unexpected names")
        missing = sorted(set(declared_columns) - set(count_columns))
        if missing:
            raise InputError(f"count table is missing declared sample columns: {missing}")
        header_index = {name: index for index, name in enumerate(header)}
        selected_indices = np.asarray([header_index[name] for name in declared_columns], dtype=int)
        matrix = np.empty((expected_rows, len(declared_columns)), dtype=np.float64)
        gene_ids: list[str] = []
        gene_symbols: list[str] = []
        seen_genes: set[str] = set()
        row_count = 0
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise InputError(f"count row {row_number} has the wrong field count")
            if row_count >= expected_rows:
                raise InputError("count table has more gene rows than frozen")
            gene_id = row[0]
            if gene_pattern.fullmatch(gene_id) is None:
                raise InputError(f"count row {row_number} is not an Ensembl-version gene ID")
            if gene_id in seen_genes:
                raise InputError(f"count table contains duplicate gene ID {gene_id}")
            seen_genes.add(gene_id)
            try:
                all_values = np.asarray([float(value) for value in row[1:-1]], dtype=float)
            except ValueError as exc:
                raise InputError(f"count row {row_number} contains a non-numeric value") from exc
            if not np.all(np.isfinite(all_values)) or np.any(all_values < 0):
                raise InputError("count values must be finite and non-negative")
            # selected_indices refer to full CSV fields; sample values start at field 1.
            matrix[row_count, :] = all_values[selected_indices - 1]
            gene_ids.append(gene_id)
            gene_symbols.append(row[-1])
            row_count += 1
    if row_count != expected_rows:
        raise InputError(f"count gene-row count differs: {row_count} != {expected_rows}")
    library_sizes = np.sum(matrix, axis=0)
    if not np.all(np.isfinite(library_sizes)) or np.any(library_sizes <= 0):
        raise InputError("declared samples require positive finite library sums")
    # The deposited "normalized_counts" values are the author's linear CPM-like
    # normalization and have column sums near one million.  Renormalizing them
    # again changes the frozen >=1 gene universe, so the source-bound contract
    # treats these deposited values as CPM and applies only log2(CPM + 1).
    cpm = matrix.copy()
    if not np.all(np.isfinite(cpm)):
        raise InputError("CPM conversion produced non-finite values")
    return {
        "gene_ids": tuple(gene_ids),
        "gene_symbols": tuple(gene_symbols),
        "columns": tuple(declared_columns),
        "cpm": cpm,
        "library_sizes": library_sizes,
        "count_columns": tuple(count_columns),
    }


def _names_hash(names: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(names) + "\n").encode("utf-8")).hexdigest()


def _canonical_config_sha256(config: dict[str, Any]) -> str:
    payload = json.dumps(
        config, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def prepare_profiles(
    counts: dict[str, Any],
    config: dict[str, Any],
    author_contract: dict[str, Any],
) -> dict[str, Any]:
    """Build effects and all control-only masks without inventing experiment IDs."""

    cpm = np.asarray(counts["cpm"], dtype=float)
    columns = tuple(counts["columns"])
    column_index = {name: index for index, name in enumerate(columns)}
    gene_ids = np.asarray(counts["gene_ids"], dtype=object)
    gene_symbols = np.asarray(counts["gene_symbols"], dtype=object)
    log_expression = np.log2(cpm + 1.0)
    all_control_columns = [
        spec["column"]
        for donor_spec in config["donors"].values()
        for spec in donor_spec["controls"].values()
    ]
    control_indices_by_donor = {
        donor: np.asarray(
            [
                column_index[spec["column"]]
                for spec in donor_spec["controls"].values()
            ],
            dtype=int,
        )
        for donor, donor_spec in config["donors"].items()
    }
    all_control_indices = np.concatenate(tuple(control_indices_by_donor.values()))
    thresholds = tuple(
        float(value)
        for value in config["preprocessing"][
            "control_mean_cpm_sensitivity_thresholds"
        ]
    )
    global_filter_masks: dict[str, np.ndarray] = {}
    global_filter_counts: dict[str, int] = {}
    expected_counts = config["expected"]["analysis_genes_before_mask_by_threshold"]
    control_mean = np.mean(cpm[:, all_control_indices], axis=1)
    for threshold in thresholds:
        key = _threshold_key(threshold)
        mask = control_mean >= threshold
        count = int(np.sum(mask))
        if count != int(expected_counts[key]):
            raise InputError(
                f"control-only CPM gene universe at {key} differs: "
                f"{count} != {expected_counts[key]}"
            )
        global_filter_masks[key] = mask
        global_filter_counts[key] = count

    canonical_key = _threshold_key(
        float(config["preprocessing"]["control_mean_cpm_min"])
    )
    universe_mask = global_filter_masks[canonical_key]
    universe_count = int(np.sum(universe_mask))
    if universe_count != int(config["expected"]["analysis_genes_before_mask"]):
        raise InputError("canonical control-only CPM gene universe differs")

    target_positions: dict[str, int] = {}
    for target in config["preprocessing"]["target_mask"]:
        matches = np.flatnonzero(
            (gene_ids == target["gene_id"]) & (gene_symbols == target["symbol"])
        )
        if len(matches) != 1 or not universe_mask[int(matches[0])]:
            raise InputError(
                f"target mask requires one canonically retained "
                f"{target['symbol']}/{target['gene_id']} row"
            )
        target_positions[target["symbol"]] = int(matches[0])
    target_row_mask = np.zeros(len(gene_ids), dtype=bool)
    target_row_mask[np.asarray(list(target_positions.values()), dtype=int)] = True
    transcriptome_mask = universe_mask & ~target_row_mask
    if int(np.sum(transcriptome_mask)) != int(
        config["expected"]["analysis_genes_after_mask"]
    ):
        raise InputError("masked transcriptome gene count differs from the frozen contract")

    component_order = tuple(config["models"]["component_order"])
    unrelated_order = tuple(config["models"]["unrelated_order"])
    full_profiles: dict[str, Any] = {}
    profiles: dict[str, Any] = {}
    on_target: dict[str, Any] = {}
    on_target_profiles: dict[str, Any] = {}
    aavs1_controls: dict[str, np.ndarray] = {}
    pairing: dict[str, Any] = {}
    target_index_vector = np.asarray(
        [target_positions[name] for name in component_order], dtype=int
    )
    for donor, donor_spec in config["donors"].items():
        single_control_roles = ("AAVS1_Guide1", "AAVS1_Guide2")
        single_control_indices = [
            column_index[donor_spec["controls"][role]["column"]]
            for role in single_control_roles
        ]
        triple_control_index = column_index[
            donor_spec["controls"]["NTC_multiplexing"]["column"]
        ]
        single_control_baseline = np.mean(
            log_expression[:, single_control_indices], axis=1
        )
        triple_control_baseline = log_expression[:, triple_control_index]

        def single_effect(spec: dict[str, str]) -> np.ndarray:
            return (
                log_expression[:, column_index[spec["column"]]]
                - single_control_baseline
            )

        def triple_effect(spec: dict[str, str]) -> np.ndarray:
            return (
                log_expression[:, column_index[spec["column"]]]
                - triple_control_baseline
            )

        component_full = np.vstack(
            [single_effect(donor_spec["component"][name]) for name in component_order]
        )
        unrelated_full = np.vstack(
            [single_effect(donor_spec["unrelated"][name]) for name in unrelated_order]
        )
        triple_full = triple_effect(donor_spec["triple"])
        full_profiles[donor] = {
            "component": component_full,
            "unrelated": unrelated_full,
            "triple": triple_full,
        }
        profiles[donor] = {
            "component": component_full[:, transcriptome_mask],
            "unrelated": unrelated_full[:, transcriptome_mask],
            "triple": triple_full[transcriptome_mask],
        }
        aavs1_controls[donor] = log_expression[
            :, np.asarray(single_control_indices, dtype=int)
        ].T[:, transcriptome_mask]
        on_target[donor] = {
            name: {
                "gene_id": str(gene_ids[target_positions[name]]),
                "single_effect": float(component_full[index, target_positions[name]]),
                "triple_effect": float(triple_full[target_positions[name]]),
            }
            for index, name in enumerate(component_order)
        }
        on_target_profiles[donor] = {
            "component": component_full[:, target_index_vector],
            "triple": triple_full[target_index_vector],
        }
        author_roles = author_contract["declared_roles"][donor]
        single_role_keys = [
            *(f"component:{name}" for name in component_order),
            *(f"unrelated:{name}" for name in unrelated_order),
            "controls:AAVS1_Guide1",
            "controls:AAVS1_Guide2",
        ]
        single_experiment_ids = sorted(
            {author_roles[key]["experiment_id"] for key in single_role_keys}
        )
        triple_experiment_id = author_roles["triple:triple"]["experiment_id"]
        pairing[donor] = {
            "constituent_relationship": "target-matched constituent singles",
            "single_role_experiment_ids": single_experiment_ids,
            "single_control_columns": [
                donor_spec["controls"][role]["column"]
                for role in single_control_roles
            ],
            "triple_role_experiment_id": triple_experiment_id,
            "triple_control_experiment_id": author_roles[
                "controls:NTC_multiplexing"
            ]["experiment_id"],
            "triple_control_column": donor_spec["controls"]["NTC_multiplexing"][
                "column"
            ],
            "cross_experiment_confounding": True,
        }
    masked_gene_ids = tuple(str(item) for item in gene_ids[transcriptome_mask])
    masked_gene_symbols = tuple(str(item) for item in gene_symbols[transcriptome_mask])
    return {
        "profiles": profiles,
        "full_profiles": full_profiles,
        "on_target": on_target,
        "on_target_profiles": on_target_profiles,
        "aavs1_controls": aavs1_controls,
        "masked_gene_ids": masked_gene_ids,
        "masked_gene_symbols": masked_gene_symbols,
        "universe_gene_ids": tuple(str(item) for item in gene_ids[universe_mask]),
        "full_gene_ids": tuple(str(item) for item in gene_ids),
        "full_gene_symbols": tuple(str(item) for item in gene_symbols),
        "cpm": cpm,
        "target_row_mask": target_row_mask,
        "global_filter_masks": global_filter_masks,
        "global_filter_counts": global_filter_counts,
        "control_indices_by_donor": control_indices_by_donor,
        "control_columns": tuple(all_control_columns),
        "control_pairing": pairing,
        "library_sizes": {
            column: float(counts["library_sizes"][column_index[column]])
            for column in columns
        },
    }


def _metrics(prediction: np.ndarray, target: np.ndarray) -> dict[str, float]:
    prediction = np.asarray(prediction, dtype=float)
    target = np.asarray(target, dtype=float)
    if prediction.shape != target.shape or prediction.ndim != 1:
        raise InputError("prediction and target vectors do not align")
    if not np.all(np.isfinite(prediction)) or not np.all(np.isfinite(target)):
        raise InputError("prediction metrics require finite vectors")
    target_norm = float(np.linalg.norm(target))
    if target_norm == 0:
        raise InputError("measured triple is zero after target masking")
    prediction_norm = float(np.linalg.norm(prediction))
    residual_norm = float(np.linalg.norm(target - prediction))
    cosine = (
        float(np.dot(prediction, target) / (prediction_norm * target_norm))
        if prediction_norm > 0
        else 0.0
    )
    return {
        "cosine": cosine,
        "normalized_rmse": residual_norm / target_norm,
        "norm_ratio": prediction_norm / target_norm,
    }


def _vector_hash(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def _geometry(fit: ProjectionResult, genes: tuple[str, ...]) -> dict[str, Any]:
    certificate: dict[str, Any] | None
    if fit.dual_separator is None:
        certificate = None
    else:
        separator = np.asarray(fit.dual_separator, dtype=float)
        order = np.argsort(-np.abs(separator), kind="stable")[:8]
        certificate = {
            "sha256_float64_little_endian": _vector_hash(separator),
            "separation_margin": float(fit.separation_margin),
            "polarity_violation": float(fit.polarity_violation),
            "orthogonality_error": float(fit.orthogonality_error),
            "largest_absolute_coordinates": [
                {"gene_id": genes[int(index)], "value": float(separator[int(index)])}
                for index in order
            ],
            "scope": "model-relative to this measured atom cone and masked gene metric",
        }
    return {
        "strict_status": fit.geometry_status,
        "kkt_violation": float(fit.kkt_violation),
        "projection_residual_fraction": float(fit.residual_fraction),
        "certificate": certificate,
    }


def _best_single(
    effects: np.ndarray, target: np.ndarray, names: tuple[str, ...]
) -> tuple[np.ndarray, str, float]:
    norms = np.sum(effects * effects, axis=1)
    dots = effects @ target
    scales = np.divide(
        np.maximum(dots, 0.0),
        norms,
        out=np.zeros_like(dots),
        where=norms > 0,
    )
    losses = np.sum((scales[:, None] * effects - target[None, :]) ** 2, axis=1)
    best = int(np.flatnonzero(np.isclose(losses, np.min(losses), rtol=0, atol=1e-12))[0])
    return scales[best] * effects[best], names[best], float(scales[best])


def fit_lodo(
    component_by_donor: dict[str, np.ndarray],
    triple_by_donor: dict[str, np.ndarray],
    held_out_donor: str,
    component_names: tuple[str, ...] = ("FAS", "RC3H1", "SUV39H1"),
) -> dict[str, Any]:
    """Fit shared NNLS and best-single comparators without reading held target values."""

    if set(component_by_donor) != set(triple_by_donor):
        raise InputError("LODO component and triple donor sets differ")
    if held_out_donor not in component_by_donor:
        raise InputError("LODO held donor is undeclared")
    training_donors = tuple(donor for donor in component_by_donor if donor != held_out_donor)
    if not training_donors:
        raise InputError("LODO requires at least one training donor")
    shapes = {np.asarray(component_by_donor[donor]).shape for donor in component_by_donor}
    if len(shapes) != 1:
        raise InputError("LODO component matrices must share one shape")
    atom_count, gene_count = next(iter(shapes))
    if atom_count != 3 or any(
        np.asarray(triple_by_donor[donor]).shape != (gene_count,)
        for donor in triple_by_donor
    ):
        raise InputError("LODO expects three atoms and aligned donor targets")
    if len(component_names) != atom_count or len(set(component_names)) != atom_count:
        raise InputError("LODO component names must identify the three atoms")
    training_effects = np.hstack(
        [np.asarray(component_by_donor[donor], dtype=float) for donor in training_donors]
    )
    training_target = np.concatenate(
        [np.asarray(triple_by_donor[donor], dtype=float) for donor in training_donors]
    )
    training_fit = project_cone(training_effects, training_target)
    held_prediction = training_fit.coefficients @ np.asarray(
        component_by_donor[held_out_donor], dtype=float
    )
    mean_training_triple = np.mean(
        np.vstack([triple_by_donor[donor] for donor in training_donors]), axis=0
    )
    _, selected_name, selected_scale = _best_single(
        training_effects, training_target, component_names
    )
    selected_index = component_names.index(selected_name)
    held_best_single_prediction = selected_scale * np.asarray(
        component_by_donor[held_out_donor], dtype=float
    )[selected_index]
    return {
        "training_donors": training_donors,
        "coefficients": training_fit.coefficients.copy(),
        "prediction": held_prediction,
        "training_mean_triple": mean_training_triple,
        "training_fit": training_fit,
        "training_selected_best_single": {
            "selected": selected_name,
            "nonnegative_scale": selected_scale,
            "prediction": held_best_single_prediction,
        },
    }


def _on_target_summary(on_target: dict[str, Any], genes: tuple[str, ...]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for gene in genes:
        singles = [float(on_target[donor][gene]["single_effect"]) for donor in on_target]
        triples = [float(on_target[donor][gene]["triple_effect"]) for donor in on_target]
        result[gene] = {
            "single_effect_mean": float(statistics.fmean(singles)),
            "single_effect_median": float(statistics.median(singles)),
            "single_effect_range": [float(min(singles)), float(max(singles))],
            "single_negative_donors": int(sum(value < 0 for value in singles)),
            "triple_effect_mean": float(statistics.fmean(triples)),
            "triple_effect_median": float(statistics.median(triples)),
            "triple_effect_range": [float(min(triples)), float(max(triples))],
            "triple_negative_donors": int(sum(value < 0 for value in triples)),
        }
    return result


def _on_target_coordinate_analysis(
    profiles: dict[str, Any],
    component_names: tuple[str, ...],
    target_gene_ids: tuple[str, ...],
) -> dict[str, Any]:
    """Three-coordinate saturated calibration, kept outside validation metrics."""

    donor_results: dict[str, Any] = {}
    for donor, profile in profiles.items():
        component = np.asarray(profile["component"], dtype=float)
        target = np.asarray(profile["triple"], dtype=float)
        cone = project_cone(component, target)
        equal_sum = component.sum(axis=0)
        best_prediction, best_name, best_scale = _best_single(
            component, target, component_names
        )
        donor_results[donor] = {
            "models": {
                "component_cone": {
                    "coefficients": [float(value) for value in cone.coefficients],
                    "metrics": _metrics(cone.fitted, target),
                    "geometry": _geometry(cone, target_gene_ids),
                },
                "equal_sum": {
                    "coefficients": [1.0, 1.0, 1.0],
                    "metrics": _metrics(equal_sum, target),
                },
                "best_single_oracle": {
                    "selected": best_name,
                    "nonnegative_scale": best_scale,
                    "metrics": _metrics(best_prediction, target),
                    "selection_status": "IN_SAMPLE_ORACLE",
                },
            }
        }
    return {
        "scope": "three on-target coordinates fitted by three atoms",
        "warning_status": "SATURATED_3X3_CALIBRATION_DIAGNOSTIC_NOT_VALIDATION",
        "interpretation": (
            "The three-atom cone is fit on only three on-target coordinates. Exact or "
            "near-exact interpolation is a saturated calibration diagnostic, not "
            "mechanistic, combination, or additivity validation."
        ),
        "per_donor": donor_results,
        "donor_ranges": _summaries(donor_results),
    }


def _summaries(donor_results: dict[str, Any]) -> dict[str, Any]:
    model_names = tuple(next(iter(donor_results.values()))["models"])
    metrics = ("cosine", "normalized_rmse", "norm_ratio")
    output: dict[str, Any] = {}
    for model in model_names:
        output[model] = {}
        for metric in metrics:
            values = [
                float(result["models"][model]["metrics"][metric])
                for result in donor_results.values()
            ]
            output[model][metric] = {
                "mean": float(statistics.fmean(values)),
                "median": float(statistics.median(values)),
                "donor_range": [float(min(values)), float(max(values))],
            }
    return output


def _cosine(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
    return float(np.dot(first, second) / denominator) if denominator > 0 else 0.0


def _pairwise_cosines(vectors: dict[str, np.ndarray]) -> dict[str, Any]:
    pairs = []
    for donor_a, donor_b in itertools.combinations(vectors, 2):
        pairs.append(
            {
                "donor_a": donor_a,
                "donor_b": donor_b,
                "cosine": _cosine(vectors[donor_a], vectors[donor_b]),
            }
        )
    values = [item["cosine"] for item in pairs]
    return {
        "pair_count": len(pairs),
        "pairs": pairs,
        "summary": {
            "mean": float(statistics.fmean(values)),
            "median": float(statistics.median(values)),
            "range": [float(min(values)), float(max(values))],
        },
    }


def _aavs1_replicate_noise(prepared: dict[str, Any]) -> dict[str, Any]:
    per_donor: dict[str, Any] = {}
    for donor, controls in prepared["aavs1_controls"].items():
        guide_1 = np.asarray(controls[0], dtype=float)
        guide_2 = np.asarray(controls[1], dtype=float)
        metrics = _metrics(guide_1, guide_2)
        noise_norm = float(np.linalg.norm(guide_1 - guide_2))
        triple_norm = float(np.linalg.norm(prepared["profiles"][donor]["triple"]))
        component_norms = np.linalg.norm(
            prepared["profiles"][donor]["component"], axis=1
        )
        per_donor[donor] = {
            "prediction_role": "AAVS1 Guide1 expression predicts AAVS1 Guide2 expression",
            "metrics": metrics,
            "difference_l2_norm": noise_norm,
            "difference_norm_over_triple_effect_norm": (
                noise_norm / triple_norm if triple_norm > 0 else None
            ),
            "difference_norm_over_median_constituent_effect_norm": (
                noise_norm / float(np.median(component_norms))
                if float(np.median(component_norms)) > 0
                else None
            ),
        }
    return {
        "definition": (
            "Guide1-versus-Guide2 uses target-masked log2(CPM+1) expression; "
            "normalized_rmse is ||Guide1-Guide2||/||Guide2||. Difference-norm "
            "ratios compare replicate noise with measured effect magnitudes."
        ),
        "per_donor": per_donor,
        "donor_ranges": {
            metric: {
                "median": float(
                    statistics.median(
                        item["metrics"][metric] for item in per_donor.values()
                    )
                ),
                "range": [
                    float(min(item["metrics"][metric] for item in per_donor.values())),
                    float(max(item["metrics"][metric] for item in per_donor.values())),
                ],
            }
            for metric in ("cosine", "normalized_rmse")
        },
    }


def _filter_sensitivity(
    prepared: dict[str, Any], config: dict[str, Any]
) -> dict[str, Any]:
    full_profiles = prepared["full_profiles"]
    thresholds = tuple(
        float(value)
        for value in config["preprocessing"][
            "control_mean_cpm_sensitivity_thresholds"
        ]
    )
    per_threshold: dict[str, Any] = {}
    median_component_cosines: list[float] = []
    for threshold in thresholds:
        key = _threshold_key(threshold)
        mask = prepared["global_filter_masks"][key] & ~prepared["target_row_mask"]
        donor_results: dict[str, Any] = {}
        for donor, profile in full_profiles.items():
            component = np.asarray(profile["component"][:, mask], dtype=float)
            target = np.asarray(profile["triple"][mask], dtype=float)
            fit = project_cone(component, target)
            equal_sum = component.sum(axis=0)
            donor_results[donor] = {
                "models": {
                    "component_cone": {
                        "metrics": _metrics(fit.fitted, target),
                        "strict_status": fit.geometry_status,
                    },
                    "equal_sum": {"metrics": _metrics(equal_sum, target)},
                }
            }
        summaries = _summaries(donor_results)
        component_median = summaries["component_cone"]["cosine"]["median"]
        equal_sum_median = summaries["equal_sum"]["cosine"]["median"]
        component_nrmse = summaries["component_cone"]["normalized_rmse"]["median"]
        equal_sum_nrmse = summaries["equal_sum"]["normalized_rmse"]["median"]
        median_component_cosines.append(component_median)
        retained = int(np.sum(mask))
        per_threshold[key] = {
            "control_mean_cpm_min": threshold,
            "retained_genes_before_target_mask": prepared["global_filter_counts"][key],
            "retained_transcriptome_genes": retained,
            "gene_ids_sha256": _names_hash(
                np.asarray(prepared["full_gene_ids"], dtype=object)[mask]
            ),
            "per_donor": donor_results,
            "donor_ranges": summaries,
            "conclusion": (
                f"At control-mean CPM >= {key}, {retained} target-masked genes remain; "
                f"median component-cone cosine/nRMSE={component_median:.6f}/"
                f"{component_nrmse:.6f} and median equal-sum cosine/nRMSE="
                f"{equal_sum_median:.6f}/{equal_sum_nrmse:.6f}. This remains a confounded "
                "cross-experiment descriptive stress test."
            ),
        }
    cosine_span = float(max(median_component_cosines) - min(median_component_cosines))
    filter_status = (
        "FILTER_SENSITIVE"
        if cosine_span >= 0.05
        else "NO_LARGE_SHIFT_ON_DECLARED_FILTER_GRID"
    )
    return {
        "selection_contract": (
            "Each mask uses only the 12 declared controls, never perturbation outcomes."
        ),
        "per_threshold": per_threshold,
        "component_cone_median_cosine_span": cosine_span,
        "filter_sensitivity_status": filter_status,
        "conclusion": (
            "Component-cone alignment changes materially across the declared "
            "control-only expression thresholds; no threshold rescues the geometric "
            "model."
            if filter_status == "FILTER_SENSITIVE"
            else "The declared threshold grid does not materially shift median cone cosine."
        ),
    }


def _fold_gene_mask(
    prepared: dict[str, Any], training_donors: tuple[str, ...], threshold: float
) -> np.ndarray:
    """Select a LODO mask using only training-donor control columns."""

    if not training_donors or any(
        donor not in prepared["control_indices_by_donor"] for donor in training_donors
    ):
        raise InputError("LODO training-control donor set is invalid")
    training_control_indices = np.concatenate(
        [prepared["control_indices_by_donor"][donor] for donor in training_donors]
    )
    mask = (
        np.mean(prepared["cpm"][:, training_control_indices], axis=1) >= threshold
    )
    return mask & ~prepared["target_row_mask"]


def _top_residual_genes(
    residual: np.ndarray,
    gene_ids: tuple[str, ...],
    gene_symbols: tuple[str, ...],
    limit: int = 12,
) -> list[dict[str, Any]]:
    order = np.argsort(-np.abs(residual), kind="stable")[:limit]
    return [
        {
            "gene_id": gene_ids[int(index)],
            "gene_symbol": gene_symbols[int(index)],
            "residual": float(residual[int(index)]),
            "absolute_residual": float(abs(residual[int(index)])),
        }
        for index in order
    ]


def _residual_summary(
    residual: np.ndarray,
    target: np.ndarray,
    gene_ids: tuple[str, ...],
    gene_symbols: tuple[str, ...],
) -> dict[str, Any]:
    target_norm = float(np.linalg.norm(target))
    residual_norm = float(np.linalg.norm(residual))
    return {
        "residual_l2_norm": residual_norm,
        "residual_norm_over_triple_norm": residual_norm / target_norm,
        "mean_residual": float(np.mean(residual)),
        "median_residual": float(np.median(residual)),
        "mean_absolute_residual": float(np.mean(np.abs(residual))),
        "median_absolute_residual": float(np.median(np.abs(residual))),
        "positive_gene_count": int(np.sum(residual > 0)),
        "negative_gene_count": int(np.sum(residual < 0)),
        "top_absolute_residual_genes": _top_residual_genes(
            residual, gene_ids, gene_symbols
        ),
    }


def analyze(prepared: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    profiles = prepared["profiles"]
    genes = prepared["masked_gene_ids"]
    component_names = tuple(config["models"]["component_order"])
    unrelated_names = tuple(config["models"]["unrelated_order"])
    symbols = prepared["masked_gene_symbols"]
    donor_results: dict[str, Any] = {}
    additive_residuals: dict[str, np.ndarray] = {}
    for donor, profile in profiles.items():
        component = np.asarray(profile["component"], dtype=float)
        unrelated = np.asarray(profile["unrelated"], dtype=float)
        target = np.asarray(profile["triple"], dtype=float)
        component_fit = project_cone(component, target)
        unrelated_fit = project_cone(unrelated, target)
        equal_sum = np.sum(component, axis=0)
        best_prediction, best_name, best_scale = _best_single(
            component, target, component_names
        )
        on_target_component = np.asarray(
            prepared["on_target_profiles"][donor]["component"], dtype=float
        )
        on_target_triple = np.asarray(
            prepared["on_target_profiles"][donor]["triple"], dtype=float
        )
        potency_fit = project_cone(on_target_component, on_target_triple)
        potency_prediction = potency_fit.coefficients @ component
        signed_coefficients, _, _, _ = np.linalg.lstsq(component.T, target, rcond=None)
        signed_prediction = signed_coefficients @ component
        additive_residuals[donor] = target - equal_sum
        donor_results[donor] = {
            "models": {
                "component_cone": {
                    "coefficients": [float(value) for value in component_fit.coefficients],
                    "metrics": _metrics(component_fit.fitted, target),
                    "geometry": _geometry(component_fit, genes),
                    "selection_status": "IN_SAMPLE_ORACLE_FIT",
                },
                "equal_sum": {
                    "coefficients": [1.0, 1.0, 1.0],
                    "metrics": _metrics(equal_sum, target),
                },
                "best_single_oracle": {
                    "selected": best_name,
                    "nonnegative_scale": best_scale,
                    "metrics": _metrics(best_prediction, target),
                    "selection_status": "IN_SAMPLE_ORACLE",
                },
                "unrelated_cone": {
                    "atom_order": list(unrelated_names),
                    "coefficients": [float(value) for value in unrelated_fit.coefficients],
                    "metrics": _metrics(unrelated_fit.fitted, target),
                    "geometry": _geometry(unrelated_fit, genes),
                },
                "zero": {"metrics": _metrics(np.zeros_like(target), target)},
                "signed_least_squares": {
                    "coefficients": [float(value) for value in signed_coefficients],
                    "metrics": _metrics(signed_prediction, target),
                    "interpretation": "in-sample signed upper-bound diagnostic",
                },
                "on_target_calibrated_component_cone": {
                    "coefficient_order": list(component_names),
                    "coefficients": [
                        float(value) for value in potency_fit.coefficients
                    ],
                    "calibration_fit": {
                        "coordinate_count": 3,
                        "atom_count": 3,
                        "metrics": _metrics(potency_fit.fitted, on_target_triple),
                        "strict_status": potency_fit.geometry_status,
                        "warning_status": (
                            "SATURATED_3X3_CALIBRATION_DIAGNOSTIC_NOT_VALIDATION"
                        ),
                    },
                    "metrics": _metrics(potency_prediction, target),
                    "scoring_contract": (
                        "coefficients fit only on the three on-target coordinates, "
                        "then frozen and scored on the disjoint target-masked transcriptome"
                    ),
                },
            }
        }

    for held_out_donor in profiles:
        training_donors = tuple(
            donor for donor in profiles if donor != held_out_donor
        )
        fold_mask = _fold_gene_mask(
            prepared,
            training_donors,
            float(config["preprocessing"]["control_mean_cpm_min"]),
        )
        component_by_donor = {
            donor: np.asarray(profile["component"][:, fold_mask], dtype=float)
            for donor, profile in prepared["full_profiles"].items()
        }
        triple_by_donor = {
            donor: np.asarray(profile["triple"][fold_mask], dtype=float)
            for donor, profile in prepared["full_profiles"].items()
        }
        lodo = fit_lodo(
            component_by_donor,
            triple_by_donor,
            held_out_donor,
            component_names,
        )
        target = triple_by_donor[held_out_donor]
        fold_gene_ids = tuple(
            np.asarray(prepared["full_gene_ids"], dtype=object)[fold_mask]
        )
        fold_filter = {
            "selection_scope": "training-donor controls only",
            "training_donors": list(training_donors),
            "held_donor_controls_used": False,
            "control_mean_cpm_min": float(
                config["preprocessing"]["control_mean_cpm_min"]
            ),
            "retained_transcriptome_genes": int(np.sum(fold_mask)),
            "gene_ids_sha256": _names_hash(fold_gene_ids),
        }
        donor_results[held_out_donor]["models"]["lodo_component_cone"] = {
            "training_donors": list(lodo["training_donors"]),
            "coefficient_order": list(component_names),
            "coefficients": [float(value) for value in lodo["coefficients"]],
            "metrics": _metrics(lodo["prediction"], target),
            "training_fit": {
                "strict_status": lodo["training_fit"].geometry_status,
                "kkt_violation": float(lodo["training_fit"].kkt_violation),
                "projection_residual_fraction": float(
                    lodo["training_fit"].residual_fraction
                ),
            },
            "gene_filter": fold_filter,
            "design_status": "CONDITIONAL_TRANSDUCTIVE",
            "leakage_contract": (
                "fold mask uses training-donor controls only; coefficients do not read "
                "the held triple; scoring remains conditional on author-normalized "
                "matrix values and the held donor's measured constituent singles"
            ),
        }
        training_selected = lodo["training_selected_best_single"]
        donor_results[held_out_donor]["models"][
            "lodo_training_selected_best_single"
        ] = {
            "training_donors": list(lodo["training_donors"]),
            "selected": training_selected["selected"],
            "nonnegative_scale": training_selected["nonnegative_scale"],
            "metrics": _metrics(training_selected["prediction"], target),
            "gene_filter": fold_filter,
            "selection_contract": (
                "identity and scale selected on stacked training donors only and "
                "frozen before application to the held donor"
            ),
        }
        donor_results[held_out_donor]["models"]["training_donor_mean_triple"] = {
            "training_donors": list(lodo["training_donors"]),
            "metrics": _metrics(lodo["training_mean_triple"], target),
            "gene_filter": fold_filter,
        }

    pairwise = {
        "triple": _pairwise_cosines(
            {
                donor: np.asarray(profile["triple"], dtype=float)
                for donor, profile in profiles.items()
            }
        ),
        "components": {
            name: _pairwise_cosines(
                {
                    donor: np.asarray(profile["component"][index], dtype=float)
                    for donor, profile in profiles.items()
                }
            )
            for index, name in enumerate(component_names)
        },
    }
    pairwise_medians = [
        pairwise["triple"]["summary"]["median"],
        *(
            pairwise["components"][name]["summary"]["median"]
            for name in component_names
        ),
    ]
    snr_status = (
        "SNR_LIMITED_LOW_DONOR_REPRODUCIBILITY"
        if max(abs(value) for value in pairwise_medians) < 0.2
        else "NOT_CLASSIFIED_LOW_REPRODUCIBILITY"
    )
    filter_sensitivity = _filter_sensitivity(prepared, config)
    overall_reliability = (
        "SNR_LIMITED_FILTER_SENSITIVE"
        if snr_status == "SNR_LIMITED_LOW_DONOR_REPRODUCIBILITY"
        and filter_sensitivity["filter_sensitivity_status"] == "FILTER_SENSITIVE"
        else "DESCRIPTIVE_RELIABILITY_DIAGNOSTICS_ONLY"
    )
    additive_per_donor = {
        donor: _residual_summary(
            residual,
            np.asarray(profiles[donor]["triple"], dtype=float),
            genes,
            symbols,
        )
        for donor, residual in additive_residuals.items()
    }
    mean_residual = np.mean(np.vstack(tuple(additive_residuals.values())), axis=0)
    return {
        "donor_results": donor_results,
        "donor_ranges": _summaries(donor_results),
        "on_target_controls": {
            "per_donor": prepared["on_target"],
            "summary": _on_target_summary(prepared["on_target"], component_names),
            "saturated_calibration_diagnostic": _on_target_coordinate_analysis(
                prepared["on_target_profiles"],
                component_names,
                tuple(
                    item["gene_id"] for item in config["preprocessing"]["target_mask"]
                ),
            ),
            "excluded_from_transcriptome_metrics": True,
        },
        "additive_residual": {
            "definition": "measured triple effect minus the unit-coefficient sum of constituent-single effects",
            "per_donor": additive_per_donor,
            "across_donor_mean_residual": _residual_summary(
                mean_residual,
                np.mean(
                    np.vstack(
                        [profiles[donor]["triple"] for donor in profiles]
                    ),
                    axis=0,
                ),
                genes,
                symbols,
            ),
            "causal_interpretation": (
                "UNAVAILABLE: residuals cannot distinguish interaction from experiment, "
                "control-type, or guide-burden differences."
            ),
        },
        "module_error": {
            "status": "UNAVAILABLE_NO_PREREGISTERED_MODULE_SET",
            "reason": (
                "No preregistered module set was supplied; post hoc module selection "
                "would not be a valid holdout analysis."
            ),
        },
        "filter_sensitivity": filter_sensitivity,
        "reliability_diagnostics": {
            "overall_status": overall_reliability,
            "snr_status": snr_status,
            "pairwise_donor_cosines": pairwise,
            "aavs1_guide_replicate_noise": _aavs1_replicate_noise(prepared),
            "interpretation": (
                "Low cross-donor effect cosine and sensitivity to the control-only "
                "expression filter limit transcriptome-wide reliability."
            ),
        },
    }


def run(config: dict[str, Any]) -> dict[str, Any]:
    validate_config(config)
    # Every complete-byte gate intentionally precedes any input parsing.
    input_identity = {
        name: verify_input(spec) for name, spec in config["inputs"].items()
    }
    author_key = parse_author_key(_resolve(config["inputs"]["author_key"]["path"]))
    soft = parse_soft(_resolve(config["inputs"]["soft"]["path"]))
    soft_contract = validate_soft(soft, config)
    counts = load_counts(_resolve(config["inputs"]["counts"]["path"]), config)
    author_contract = validate_author_key(author_key, counts, config)
    prepared = prepare_profiles(counts, config, author_contract)
    analysis = analyze(prepared, config)
    ranges = analysis["donor_ranges"]
    component_median = ranges["component_cone"]["cosine"]["median"]
    equal_sum_cosine = ranges["equal_sum"]["cosine"]["median"]
    equal_sum_nrmse = ranges["equal_sum"]["normalized_rmse"]["median"]
    best_oracle_median = ranges["best_single_oracle"]["cosine"]["median"]
    lodo_median = ranges["lodo_component_cone"]["cosine"]["median"]
    lodo_best_median = ranges["lodo_training_selected_best_single"]["cosine"][
        "median"
    ]
    potency_median = ranges["on_target_calibrated_component_cone"]["cosine"][
        "median"
    ]
    pairwise = analysis["reliability_diagnostics"]["pairwise_donor_cosines"]
    threshold_cosines = {
        key: value["donor_ranges"]["component_cone"]["cosine"]["median"]
        for key, value in analysis["filter_sensitivity"]["per_threshold"].items()
    }
    report = {
        "report_version": 2,
        "benchmark_id": config["benchmark_id"],
        "config_canonical_json_sha256": _canonical_config_sha256(config),
        "execution_status": "PASS",
        "geometric_model_status": "FAILS_DECLARED_GEOMETRIC_MODEL",
        "biological_interpretation_status": (
            "INCONCLUSIVE_CROSS_EXPERIMENT_CONFOUNDING_LOW_RELIABILITY"
        ),
        "source": config["source"],
        "input_identity": input_identity,
        "sample_contract": {
            "geo_soft": soft_contract,
            "author_key": author_contract,
        },
        "metadata_conflict": {
            **config["metadata_conflict"],
            "observed_affected_gsms": soft_contract["observed_time_conflict_gsms"],
        },
        "preprocessing": {
            **config["preprocessing"],
            "observed_gene_rows": len(counts["gene_ids"]),
            "observed_count_samples": len(counts["count_columns"]),
            "observed_gene_universe_before_mask": len(prepared["universe_gene_ids"]),
            "observed_transcriptome_genes_after_mask": len(prepared["masked_gene_ids"]),
            "observed_gene_universe_before_mask_by_threshold": prepared[
                "global_filter_counts"
            ],
            "universe_gene_ids_sha256": _names_hash(prepared["universe_gene_ids"]),
            "masked_gene_ids_sha256": _names_hash(prepared["masked_gene_ids"]),
            "control_columns": list(prepared["control_columns"]),
            "control_pairing": prepared["control_pairing"],
            "selected_sample_library_sums": prepared["library_sizes"],
        },
        "models": config["models"],
        "claim_contract": config["claim_contract"],
        **analysis,
        "interpretation": {
            "numerical_summary": {
                "component_cone_median_cosine": component_median,
                "equal_sum_median_cosine": equal_sum_cosine,
                "equal_sum_median_normalized_rmse": equal_sum_nrmse,
                "best_single_in_sample_oracle_median_cosine": best_oracle_median,
                "lodo_component_cone_median_cosine": lodo_median,
                "lodo_training_selected_best_single_median_cosine": lodo_best_median,
                "on_target_calibrated_disjoint_transcriptome_median_cosine": (
                    potency_median
                ),
                "triple_pairwise_donor_median_cosine": pairwise["triple"][
                    "summary"
                ]["median"],
                "component_pairwise_donor_median_cosines": {
                    name: value["summary"]["median"]
                    for name, value in pairwise["components"].items()
                },
                "filter_threshold_component_cone_median_cosines": threshold_cosines,
            },
            "text": (
                f"The in-sample component cone reaches median cosine "
                f"{component_median:.3f}, only descriptively, while the equal-sum "
                f"model has median cosine {equal_sum_cosine:.3f} and normalized RMSE "
                f"{equal_sum_nrmse:.3f}. LODO median cosine is {lodo_median:.3f} "
                f"versus {lodo_best_median:.3f} for the training-selected best single. "
                f"On-target-frozen potency scoring reaches median transcriptome cosine "
                f"{potency_median:.3f}. Low donor reproducibility and expression-filter "
                "sensitivity, together with perfect single/triple experiment-control-guide "
                "burden confounding, make the biological result inconclusive."
            ),
            "claim_boundary": (
                "Execution succeeded, but the declared geometric model fails. This is "
                "a negative cross-experiment stress test, not combination/additivity "
                "validation, an interaction test, or biological reachability evidence."
            ),
        },
    }
    # Reject NaN/Infinity before the CLI can print or freeze the report.
    json.dumps(report, allow_nan=False)
    return report


def _canonical_json(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    output = parser.add_mutually_exclusive_group()
    output.add_argument("--check", type=Path, help="compare against an existing JSON report")
    output.add_argument("--write", type=Path, help="write the deterministic JSON report")
    args = parser.parse_args()

    with args.config.open(encoding="utf-8") as handle:
        config = json.load(handle)
    report = run(config)
    rendered = _canonical_json(report)
    if args.check is not None:
        with args.check.open(encoding="utf-8") as handle:
            expected = json.load(handle)
        if report != expected:
            raise SystemExit(f"Goudy report differs from {args.check}")
        print(f"Goudy report matches {args.check}")
    elif args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
        print(f"Wrote {args.write}")
    else:
        sys.stdout.write(rendered)


if __name__ == "__main__":
    main()
