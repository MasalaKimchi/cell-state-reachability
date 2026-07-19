#!/usr/bin/env python3
"""Validate Zhu screen effects against source-selected arrayed follow-up assays.

The benchmark compares nine Stim8hr genome-scale Perturb-seq profiles with arrayed
bulk RNA-seq and donor-normalized IL-10/IL-21 flow measurements. The panel was
selected upstream by the source study, so all permutation diagnostics are
conditional on the selected panel and support no held-out-discovery claim.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError  # noqa: E402


DEFAULT_CONFIG = ROOT / "configs" / "zhu_arrayed_validation.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_input(spec: dict[str, Any], *, verify_hash: bool = True) -> dict[str, Any]:
    path = ROOT / spec["path"]
    if not path.is_file():
        raise FileNotFoundError(path)
    size = path.stat().st_size
    if size != spec["bytes"]:
        raise InputError(f"{spec['path']} byte length differs: {size} != {spec['bytes']}")
    observed_hash = sha256_file(path) if verify_hash else None
    if verify_hash and observed_hash != spec["sha256"]:
        raise InputError(f"{spec['path']} SHA-256 differs from the frozen identity")
    return {
        "path": spec["path"],
        "bytes": size,
        "sha256_expected": spec["sha256"],
        "sha256_actual": observed_hash,
        "hash_verified": verify_hash,
    }


def _decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values],
        dtype=object,
    )


def _categorical(group: Any) -> np.ndarray:
    categories = _decode(group["categories"][:])
    codes = np.asarray(group["codes"][:], dtype=np.int64)
    if np.any(codes < 0) or np.any(codes >= len(categories)):
        raise InputError("H5AD categorical field contains missing or invalid codes")
    return categories[codes]


def load_screen(config: dict[str, Any]) -> dict[str, Any]:
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError("arrayed validation requires requirements-external.txt") from exc

    spec = config["screen"]
    expected = config["expected"]
    targets = tuple(config["panel"]["perturbations"])
    path = ROOT / config["inputs"]["screen"]["path"]
    with h5py.File(path, "r") as handle:
        genes = _decode(handle[f"var/{spec['gene_id_field']}"][:])
        symbols = _decode(handle[f"var/{spec['gene_symbol_field']}"][:])
        conditions = _categorical(handle[f"obs/{spec['condition_field']}"])
        perturbations = _categorical(handle[f"obs/{spec['perturbation_field']}"])
        layer = handle[f"layers/{spec['layer']}"]
        if layer.shape != (len(conditions), len(genes)):
            raise InputError("screen effect layer does not match observation/gene axes")
        if len(conditions) != expected["screen_rows"] or len(genes) != expected["screen_genes"]:
            raise InputError("screen axes differ from the frozen contract")
        if len(set(genes)) != len(genes) or len(set(symbols)) != len(symbols):
            raise InputError("screen gene IDs and symbols must each be unique")
        rows = []
        row_indices = []
        for target in targets:
            matches = np.flatnonzero(
                (conditions == spec["condition"]) & (perturbations == target)
            )
            if len(matches) != 1:
                raise InputError(f"expected one {target}/{spec['condition']} screen row")
            row_indices.append(int(matches[0]))
            values = np.asarray(layer[int(matches[0]), :], dtype=float)
            if not np.all(np.isfinite(values)):
                raise InputError(f"screen {target} profile contains non-finite values")
            rows.append(values)
    return {
        "targets": targets,
        "genes": tuple(str(item) for item in genes),
        "symbols": tuple(str(item) for item in symbols),
        "row_indices": row_indices,
        "profiles": np.vstack(rows),
    }


def _finite(value: str, *, field: str, row: int) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(f"row {row}: {field} is not numeric") from exc
    if not math.isfinite(result):
        raise InputError(f"row {row}: {field} is not finite")
    return result


def load_bulk(config: dict[str, Any]) -> dict[str, dict[str, float]]:
    targets = tuple(config["panel"]["perturbations"])
    expected = config["expected"]
    result: dict[str, dict[str, float]] = {target: {} for target in targets}
    path = ROOT / config["inputs"]["bulk_rna"]["path"]
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"variable", "log_fc", "contrast"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise InputError("bulk RNA table schema differs from the frozen contract")
        count = 0
        for count, row in enumerate(reader, start=1):
            target = row["contrast"].strip()
            gene = row["variable"].strip()
            if target not in result or not gene:
                raise InputError(f"bulk RNA row {count} has an unexpected target/gene")
            if gene in result[target]:
                raise InputError(f"bulk RNA has duplicate {target}/{gene}")
            result[target][gene] = _finite(row["log_fc"], field="log_fc", row=count)
    if count != expected["bulk_rows"]:
        raise InputError("bulk RNA row count differs from the frozen contract")
    if any(len(values) != expected["bulk_genes_per_perturbation"] for values in result.values()):
        raise InputError("bulk RNA perturbations do not share the frozen gene count")
    gene_sets = [set(values) for values in result.values()]
    if any(genes != gene_sets[0] for genes in gene_sets[1:]):
        raise InputError("bulk RNA perturbations do not share one gene namespace")
    return result


def load_flow(config: dict[str, Any]) -> list[dict[str, Any]]:
    expected = config["expected"]
    targets = set(config["panel"]["perturbations"])
    cytokine_fields = tuple(config["panel"]["cytokines"].values())
    control_prefix = config["panel"]["control_prefix"]
    path = ROOT / config["inputs"]["flow"]["path"]
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"Sample", "Donor", "Perturbation", *cytokine_fields}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise InputError("flow table schema differs from the frozen contract")
        for row_number, row in enumerate(reader, start=1):
            perturbation = row["Perturbation"].strip()
            if perturbation not in targets and not perturbation.startswith(control_prefix):
                raise InputError(f"flow row {row_number} has unexpected perturbation")
            parsed = {
                "sample": row["Sample"].strip(),
                "donor": row["Donor"].strip(),
                "perturbation": perturbation,
            }
            for field in cytokine_fields:
                value = _finite(row[field], field=field, row=row_number)
                if value <= 0 or value > 100:
                    raise InputError(
                        "flow percentages must be in (0, 100] for log-ratio analysis"
                    )
                parsed[field] = value
            rows.append(parsed)
    if len(rows) != expected["flow_rows"]:
        raise InputError("flow row count differs from the frozen contract")
    if sorted({row["donor"] for row in rows}) != sorted(expected["donors"]):
        raise InputError("flow donor labels differ from the frozen contract")
    keys = [(row["donor"], row["perturbation"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise InputError("flow table has duplicate donor/perturbation rows")
    return rows


def build_flow_effects(
    rows: list[dict[str, Any]], config: dict[str, Any]
) -> list[dict[str, Any]]:
    cytokines = config["panel"]["cytokines"]
    control_prefix = config["panel"]["control_prefix"]
    controls: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        if row["perturbation"].startswith(control_prefix):
            for cytokine, field in cytokines.items():
                controls[row["donor"]][cytokine].append(float(row[field]))
    donors = set(config["expected"]["donors"])
    if set(controls) != donors or any(not values for donor in donors for values in controls[donor].values()):
        raise InputError("every flow donor requires at least one NTC control")

    output = []
    for row in rows:
        if row["perturbation"].startswith(control_prefix):
            continue
        for cytokine, field in cytokines.items():
            baseline = float(np.mean(controls[row["donor"]][cytokine]))
            output.append(
                {
                    "target": row["perturbation"],
                    "donor": row["donor"],
                    "cytokine": cytokine,
                    "percent_positive": float(row[field]),
                    "donor_ntc_mean_percent_positive": baseline,
                    "log2_ratio_to_donor_ntc": float(np.log2(float(row[field]) / baseline)),
                    "n_donor_ntc_measurements": len(controls[row["donor"]][cytokine]),
                }
            )
    return sorted(output, key=lambda row: (row["cytokine"], row["target"], row["donor"]))


def align_profiles(
    screen: dict[str, Any], bulk: dict[str, dict[str, float]], config: dict[str, Any]
) -> dict[str, Any]:
    targets = screen["targets"]
    common = tuple(sorted(set(screen["genes"]) & set(next(iter(bulk.values())))))
    if len(common) != config["expected"]["common_genes"]:
        raise InputError("screen/bulk common-gene count differs from the frozen contract")
    screen_index = {gene: index for index, gene in enumerate(screen["genes"])}
    indices = [screen_index[gene] for gene in common]
    source = screen["profiles"][:, indices]
    arrayed = np.asarray([[bulk[target][gene] for gene in common] for target in targets])
    symbol_to_id = dict(zip(screen["symbols"], screen["genes"]))
    if not config["analysis"]["mask_all_panel_target_genes"]:
        raise InputError("the frozen benchmark requires masking all panel target genes")
    mask = np.ones(source.shape, dtype=bool)
    panel_gene_ids = [symbol_to_id.get(target) for target in targets]
    if any(gene_id not in common for gene_id in panel_gene_ids):
        raise InputError("every panel target gene must occur in the common gene universe")
    for gene_id in panel_gene_ids:
        mask[:, common.index(gene_id)] = False
    scored = np.sum(mask, axis=1)
    if not np.all(scored == config["expected"]["scored_genes_per_profile"]):
        raise InputError("panel-target masking differs from the frozen contract")
    return {"genes": common, "source": source, "arrayed": arrayed, "mask": mask}


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(left @ right / denominator) if denominator else float("nan")


def _correlation(left: np.ndarray, right: np.ndarray) -> float:
    centered_left = left - np.mean(left)
    centered_right = right - np.mean(right)
    return _cosine(centered_left, centered_right)


def _spearman(left: np.ndarray, right: np.ndarray) -> float:
    return _correlation(rankdata(left), rankdata(right))


def profile_metrics(
    aligned: dict[str, Any], targets: tuple[str, ...], scale: str
) -> tuple[list[dict[str, Any]], np.ndarray, np.ndarray]:
    source = np.array(aligned["source"], copy=True)
    arrayed = np.array(aligned["arrayed"], copy=True)
    if scale == "panel_centered":
        source -= np.mean(source, axis=0, keepdims=True)
        arrayed -= np.mean(arrayed, axis=0, keepdims=True)
    elif scale != "raw":
        raise ValueError(scale)
    mask = aligned["mask"]
    common_source = np.mean(source, axis=0)
    rows = []
    similarities = np.empty((len(targets), len(targets)), dtype=float)
    for observed_index, target in enumerate(targets):
        for predicted_index in range(len(targets)):
            pair_mask = mask[observed_index] & mask[predicted_index]
            similarities[observed_index, predicted_index] = _cosine(
                arrayed[observed_index, pair_mask], source[predicted_index, pair_mask]
            )
        values_mask = mask[observed_index]
        observed = arrayed[observed_index, values_mask]
        predicted = source[observed_index, values_mask]
        ordered_similarities = np.sort(similarities[observed_index])
        if np.any(np.diff(ordered_similarities) <= 1e-12):
            raise InputError("retrieval similarities contain a tie")
        rms = float(np.sqrt(np.mean(observed**2)))
        common_cosine = (
            _cosine(observed, common_source[values_mask]) if scale == "raw" else None
        )
        order = np.argsort(-similarities[observed_index], kind="stable")
        rank = int(np.flatnonzero(order == observed_index)[0]) + 1
        rows.append(
            {
                "scale": scale,
                "target": target,
                "genes_scored": int(np.sum(values_mask)),
                "cosine": _cosine(observed, predicted),
                "pearson": _correlation(observed, predicted),
                "spearman": _spearman(observed, predicted),
                "normalized_rmse": float(np.sqrt(np.mean((observed - predicted) ** 2)) / rms),
                "common_source_cosine": common_cosine,
                "cosine_gain_over_common_source": (
                    _cosine(observed, predicted) - common_cosine
                    if common_cosine is not None
                    else None
                ),
                "retrieval_rank": rank,
                "retrieval_reciprocal_rank": 1.0 / rank,
            }
        )
    ranks = np.empty_like(similarities, dtype=int)
    for row in range(len(targets)):
        order = np.argsort(-similarities[row], kind="stable")
        ranks[row, order] = np.arange(1, len(targets) + 1)
    return rows, similarities, ranks


def summarize_flow(
    effects: list[dict[str, Any]], targets: tuple[str, ...]
) -> dict[str, dict[str, Any]]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in effects:
        grouped[(row["cytokine"], row["target"])].append(row["log2_ratio_to_donor_ntc"])
    result = {}
    for cytokine in sorted({row["cytokine"] for row in effects}):
        result[cytokine] = {
            target: {
                "donors": len(grouped[(cytokine, target)]),
                "median_log2_ratio": float(np.median(grouped[(cytokine, target)])),
                "mean_log2_ratio": float(np.mean(grouped[(cytokine, target)])),
                "min_log2_ratio": float(np.min(grouped[(cytokine, target)])),
                "max_log2_ratio": float(np.max(grouped[(cytokine, target)])),
            }
            for target in targets
        }
    return result


def conditional_permutations(
    profile_data: dict[str, tuple[np.ndarray, np.ndarray]],
    rna_cytokine: dict[str, np.ndarray],
    flow_summary: dict[str, dict[str, Any]],
    targets: tuple[str, ...],
) -> dict[str, Any]:
    n = len(targets)
    observed_retrieval = {}
    for scale, (similarities, ranks) in profile_data.items():
        observed_retrieval[scale] = {
            "top1": float(np.mean(ranks[np.arange(n), np.arange(n)] == 1)),
            "mrr": float(np.mean(1.0 / ranks[np.arange(n), np.arange(n)])),
            "top1_extreme": 0,
            "mrr_extreme": 0,
        }

    correlation_vectors: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for assay in ("screen", "bulk"):
        for cytokine in ("IL10", "IL21"):
            left = rankdata(rna_cytokine[f"{assay}_{cytokine}"])
            right = rankdata(
                [flow_summary[cytokine][target]["median_log2_ratio"] for target in targets]
            )
            left_norm = np.linalg.norm(left - np.mean(left))
            right_norm = np.linalg.norm(right - np.mean(right))
            if left_norm == 0 or right_norm == 0:
                raise InputError("RNA/flow rank vectors must have nonzero variation")
            correlation_vectors[f"{assay}_{cytokine}"] = (
                (left - np.mean(left)) / left_norm,
                (right - np.mean(right)) / right_norm,
            )
    observed_correlations = {
        name: float(left @ right) for name, (left, right) in correlation_vectors.items()
    }
    unadjusted_extreme = dict.fromkeys(correlation_vectors, 0)
    max_t_extreme = dict.fromkeys(correlation_vectors, 0)
    total = 0
    for permutation in itertools.permutations(range(n)):
        perm = np.asarray(permutation, dtype=int)
        total += 1
        for scale, (_, ranks) in profile_data.items():
            selected = ranks[np.arange(n), perm]
            top1 = float(np.mean(selected == 1))
            mrr = float(np.mean(1.0 / selected))
            observed = observed_retrieval[scale]
            observed["top1_extreme"] += int(top1 >= observed["top1"] - 1e-15)
            observed["mrr_extreme"] += int(mrr >= observed["mrr"] - 1e-15)
        null_correlations = {
            name: float(left[perm] @ right)
            for name, (left, right) in correlation_vectors.items()
        }
        maximum = max(null_correlations.values())
        for name, observed in observed_correlations.items():
            unadjusted_extreme[name] += int(null_correlations[name] >= observed - 1e-15)
            max_t_extreme[name] += int(maximum >= observed - 1e-15)
    if total != math.factorial(n):
        raise AssertionError("exact target-label permutation enumeration failed")
    retrieval = {
        scale: {
            "top1": values["top1"],
            "mrr": values["mrr"],
            "exhaustive_top1_tail_fraction": values["top1_extreme"] / total,
            "exhaustive_mrr_tail_fraction": values["mrr_extreme"] / total,
        }
        for scale, values in observed_retrieval.items()
    }
    correlations = {
        name: {
            "spearman": observed_correlations[name],
            "exhaustive_one_sided_tail_fraction": unadjusted_extreme[name] / total,
            "exhaustive_one_sided_max_stat_tail_fraction": max_t_extreme[name] / total,
        }
        for name in sorted(observed_correlations)
    }
    return {
        "permutations": total,
        "interpretation": "Exhaustive conditional target-label diagnostics within the source-selected nine-perturbation panel. Unequal donor count/composition and upstream panel selection leave target exchangeability unestablished, so tail fractions are not inferential p-values, multiplicity-adjusted inference, population inference, or held-out-discovery inference.",
        "retrieval": retrieval,
        "cytokine_rank_association": correlations,
    }


def _render_csv(rows: list[dict[str, Any]], fields: list[str]) -> str:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(
            {
                field: (
                    "" if row[field] is None else str(row[field]).lower()
                    if isinstance(row[field], bool)
                    else row[field]
                )
                for field in fields
            }
        )
    return buffer.getvalue()


def run(config: dict[str, Any], *, verify_hash: bool = True) -> tuple[dict[str, Any], str, str]:
    identities = {
        name: verify_input(spec, verify_hash=verify_hash)
        for name, spec in config["inputs"].items()
    }
    screen = load_screen(config)
    bulk = load_bulk(config)
    flow = load_flow(config)
    effects = build_flow_effects(flow, config)
    aligned = align_profiles(screen, bulk, config)
    targets = screen["targets"]

    metric_rows: list[dict[str, Any]] = []
    profile_data = {}
    for scale in config["analysis"]["profile_scales"]:
        rows, similarities, ranks = profile_metrics(aligned, targets, scale)
        metric_rows.extend(rows)
        profile_data[scale] = (similarities, ranks)
    flow_summary = summarize_flow(effects, targets)
    gene_index = {gene: index for index, gene in enumerate(aligned["genes"])}
    symbol_to_id = dict(zip(screen["symbols"], screen["genes"]))
    rna_cytokine = {}
    for cytokine in ("IL10", "IL21"):
        gene_id = symbol_to_id[cytokine]
        index = gene_index[gene_id]
        rna_cytokine[f"screen_{cytokine}"] = aligned["source"][:, index]
        rna_cytokine[f"bulk_{cytokine}"] = aligned["arrayed"][:, index]
    permutation = conditional_permutations(
        profile_data, rna_cytokine, flow_summary, targets
    )
    summaries = {}
    for scale in config["analysis"]["profile_scales"]:
        rows = [row for row in metric_rows if row["scale"] == scale]
        summaries[scale] = {
            field: float(np.median([row[field] for row in rows]))
            for field in ("cosine", "pearson", "spearman", "normalized_rmse")
        }
        if scale == "raw":
            summaries[scale]["median_cosine_gain_over_common_source"] = float(
                np.median([row["cosine_gain_over_common_source"] for row in rows])
            )
    metadata = {
        "schema_version": "1.0.0",
        "status": "PASS",
        "benchmark": config["benchmark"],
        "claim_ceiling": config["claim_ceiling"],
        "provenance": config["provenance"],
        "input_verification": identities,
        "design": {
            "screen_condition": config["screen"]["condition"],
            "perturbations": list(targets),
            "followup_donor_labels": config["expected"]["donors"],
            "common_genes": len(aligned["genes"]),
            "all_panel_target_genes_masked": config["analysis"][
                "mask_all_panel_target_genes"
            ],
            "flow_transform": config["analysis"]["flow_transform"],
        },
        "profile_summary": summaries,
        "flow_summary": flow_summary,
        "conditional_permutation": permutation,
        "limitations": [
            "The nine perturbations were selected upstream using the source screen; this is measured follow-up, not held-out target discovery.",
            "Donor coverage varies by perturbation (three to six), so each target summary gives donors equal weight and remains descriptive.",
            "The conditional label permutations do not undo upstream panel selection and are not donor-population inference.",
            "Unequal donor count/composition and NTC replication leave target exchangeability unestablished; exhaustive tail fractions are not inferential p-values or multiplicity-adjusted inference.",
            "Cytokine-positive percentages do not establish durable cell-state conversion, function, fitness, chromatin remodeling, or intervention efficacy.",
        ],
    }
    profile_fields = [
        "scale", "target", "genes_scored", "cosine", "pearson", "spearman",
        "normalized_rmse", "common_source_cosine", "cosine_gain_over_common_source",
        "retrieval_rank", "retrieval_reciprocal_rank",
    ]
    flow_fields = [
        "target", "donor", "cytokine", "percent_positive",
        "donor_ntc_mean_percent_positive", "log2_ratio_to_donor_ntc",
        "n_donor_ntc_measurements",
    ]
    return metadata, _render_csv(metric_rows, profile_fields), _render_csv(effects, flow_fields)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--skip-hash", action="store_true", help="development only")
    parser.add_argument("--check", action="store_true", help="compare with maintained outputs")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    metadata, profile_csv, flow_csv = run(config, verify_hash=not args.skip_hash)
    metadata["config_sha256"] = sha256_file(args.config)
    payloads = {
        config["outputs"]["metadata"]: json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        config["outputs"]["profile_metrics"]: profile_csv,
        config["outputs"]["flow_effects"]: flow_csv,
    }
    for relative, payload in payloads.items():
        path = ROOT / relative
        if args.check:
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                raise SystemExit(f"maintained output differs: {relative}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload, encoding="utf-8")
            print(f"wrote {relative}")
    if args.check:
        print("arrayed validation outputs match")


if __name__ == "__main__":
    main()
