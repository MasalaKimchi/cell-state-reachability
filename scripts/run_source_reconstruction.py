#!/usr/bin/env python3
"""Reconstruct the frozen Zhu case study from immutable external source files.

This runner is intentionally separate from ``reproduce.sh`` because its primary H5AD
input is 16.8 GB and is not committed. It verifies file bytes, profiles schema/grain,
rebuilds the registered target lineage, freezes a source-bound z-score split protocol,
and runs source-safe Ota-to-Hollbacker and reverse transfer diagnostics on log fold
changes and z-scores. It does not provide donor or guide inference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError, held_out_alignment, project_cone


DEFAULT_CONFIG = ROOT / "configs" / "source_reconstruction.json"
DEFAULT_REPORT = ROOT / "results" / "source_reconstruction.json"


def _progress(message: str) -> None:
    print(f"[source-reconstruction] {message}", file=sys.stderr, flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_input(name: str, spec: dict[str, Any], *, verify_hash: bool) -> dict[str, Any]:
    path = ROOT / spec["path"]
    if not path.is_file():
        raise FileNotFoundError(f"{name}: {path}")
    actual_bytes = path.stat().st_size
    if actual_bytes != spec["bytes"]:
        raise InputError(
            f"{name} byte length differs: expected {spec['bytes']}, found {actual_bytes}"
        )
    actual_hash = sha256(path) if verify_hash else None
    if verify_hash and actual_hash != spec["sha256"]:
        raise InputError(f"{name} SHA-256 differs from the frozen source identity")
    return {
        "path": spec["path"],
        "bytes": actual_bytes,
        "sha256_expected": spec["sha256"],
        "sha256_actual": actual_hash,
        "hash_verified": verify_hash,
    }


def _finite_float(value: str, *, field: str, row_number: int) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(f"row {row_number}: {field} is not numeric") from exc
    if not math.isfinite(result):
        raise InputError(f"row {row_number}: {field} is not finite")
    return result


def load_target_table(
    path: Path, target_config: dict[str, Any]
) -> tuple[int, dict[str, dict[str, dict[str, float]]]]:
    """Load an exact one-row-per-gene/source/value target table."""

    gene_field = target_config["gene_field"]
    contrast_field = target_config["contrast_field"]
    sources = target_config["sources"]
    value_fields = target_config["value_fields"]
    allowed_contrasts = set(sources.values())
    values = {
        field: {source: {} for source in sources}
        for field in value_fields
    }
    source_for_contrast = {contrast: source for source, contrast in sources.items()}
    row_count = 0
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {gene_field, contrast_field, *value_fields}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or ()))
            raise InputError(f"target table is missing required fields: {missing}")
        for row_count, row in enumerate(reader, start=1):
            gene = row[gene_field].strip()
            contrast = row[contrast_field].strip()
            if not gene:
                raise InputError(f"row {row_count}: target gene is empty")
            if contrast not in allowed_contrasts:
                raise InputError(f"row {row_count}: unexpected target contrast {contrast!r}")
            source = source_for_contrast[contrast]
            for field in value_fields:
                if gene in values[field][source]:
                    raise InputError(
                        f"row {row_count}: duplicate target gene/source pair {gene}/{source}"
                    )
                values[field][source][gene] = _finite_float(
                    row[field], field=field, row_number=row_count
                )
    if row_count == 0:
        raise InputError("target table is empty")
    return row_count, values


def build_target_lineage(
    source_values: dict[str, dict[str, float]],
    screen_genes: tuple[str, ...],
    *,
    orientation_multiplier: float,
) -> dict[str, Any]:
    """Build source-safe and sign-selected target universes in deterministic order."""

    sources = tuple(source_values)
    if len(sources) != 2:
        raise InputError("exactly two target sources are required")
    first, second = sources
    first_genes = set(source_values[first])
    second_genes = set(source_values[second])
    union = first_genes | second_genes
    shared = first_genes & second_genes
    concordant = {
        gene
        for gene in shared
        if np.sign(source_values[first][gene]) == np.sign(source_values[second][gene])
    }
    screen = set(screen_genes)
    shared_screen = tuple(sorted(shared & screen))
    final = tuple(sorted(concordant & screen))
    if not shared_screen or not final:
        raise InputError("target alignment produced an empty analysis universe")

    first_direction = orientation_multiplier * np.asarray(
        [source_values[first][gene] for gene in shared_screen], dtype=float
    )
    second_direction = orientation_multiplier * np.asarray(
        [source_values[second][gene] for gene in shared_screen], dtype=float
    )
    merged_direction = orientation_multiplier * np.asarray(
        [
            0.5 * (source_values[first][gene] + source_values[second][gene])
            for gene in final
        ],
        dtype=float,
    )
    shared_order = tuple(sorted(shared))
    first_shared = orientation_multiplier * np.asarray(
        [source_values[first][gene] for gene in shared_order], dtype=float
    )
    second_shared = orientation_multiplier * np.asarray(
        [source_values[second][gene] for gene in shared_order], dtype=float
    )
    return {
        "sources": sources,
        "counts": {
            "union": len(union),
            "shared": len(shared),
            "sign_concordant": len(concordant),
            "screen_overlap": len(union & screen),
            "shared_screen": len(shared_screen),
            "final": len(final),
        },
        "shared_screen_genes": shared_screen,
        "final_genes": final,
        "source_directions": {
            first: first_direction,
            second: second_direction,
        },
        "merged_direction": merged_direction,
        "between_source_cosine": _cosine(first_shared, second_shared),
        "between_source_screen_cosine": _cosine(first_direction, second_direction),
    }


def _decode_strings(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        item.decode("utf-8") if isinstance(item, bytes) else str(item)
        for item in values
    )


def _read_h5_categorical(group: Any) -> tuple[str, ...]:
    categories = _decode_strings(group["categories"][:])
    codes = np.asarray(group["codes"][:], dtype=np.int64)
    if np.any(codes < 0) or np.any(codes >= len(categories)):
        raise InputError("H5AD categorical codes contain missing or invalid values")
    return tuple(categories[index] for index in codes)


def profile_h5ad(path: Path, screen_config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate H5AD grain/schema and return axes needed for reconstruction."""

    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError(
            "source reconstruction requires requirements-external.txt"
        ) from exc

    with h5py.File(path, "r") as handle:
        gene_path = f"var/{screen_config['gene_field']}"
        condition_path = f"obs/{screen_config['condition_field']}"
        perturbation_path = f"obs/{screen_config['perturbation_field']}"
        admission_path = f"obs/{screen_config['admission_field']}"
        for required in (gene_path, condition_path, perturbation_path, admission_path):
            if required not in handle:
                raise InputError(f"H5AD is missing {required}")
        genes = _decode_strings(handle[gene_path][:])
        if any(not gene for gene in genes) or len(set(genes)) != len(genes):
            raise InputError("H5AD screen genes must be non-empty and unique")
        conditions = _read_h5_categorical(handle[condition_path])
        perturbations = _read_h5_categorical(handle[perturbation_path])
        admission = np.asarray(handle[admission_path][:], dtype=bool)
        n_rows = len(conditions)
        if len(perturbations) != n_rows or admission.shape != (n_rows,):
            raise InputError("H5AD observation fields do not share one row grain")
        condition_counts = {
            item: conditions.count(item) for item in sorted(set(conditions))
        }
        selected_rows = np.flatnonzero(
            (np.asarray(conditions) == screen_config["condition"]) & admission
        )
        selected_perturbations = tuple(perturbations[index] for index in selected_rows)
        if len(set(selected_perturbations)) != len(selected_perturbations):
            raise InputError("selected perturbation labels are not unique")
        layer_profile = {}
        for layer in screen_config["layers"]:
            layer_path = f"layers/{layer}"
            if layer_path not in handle:
                raise InputError(f"H5AD is missing {layer_path}")
            dataset = handle[layer_path]
            if dataset.shape != (n_rows, len(genes)):
                raise InputError(f"{layer_path} does not match observation/gene axes")
            layer_profile[layer] = {
                "shape": list(dataset.shape),
                "dtype": str(dataset.dtype),
            }
    profile = {
        "grain": "one perturbation-condition DE profile by measured gene",
        "effect_rows": n_rows,
        "effect_genes": len(genes),
        "condition_counts": condition_counts,
        "selected_condition": screen_config["condition"],
        "selected_admitted_atoms": len(selected_rows),
        "layers": layer_profile,
        "inference_unit_warning": (
            "Rows are donor-collapsed perturbation-condition summaries, not donor replicates."
        ),
    }
    axes = {
        "genes": genes,
        "selected_rows": selected_rows,
        "selected_perturbations": selected_perturbations,
    }
    return profile, axes


def _assert_expected_profile(
    profile: dict[str, Any], lineage: dict[str, Any], expected: dict[str, int]
) -> None:
    observed = {
        "effect_rows": profile["effect_rows"],
        "effect_genes": profile["effect_genes"],
        "rest_admitted_atoms": profile["selected_admitted_atoms"],
        "target_union_genes": lineage["counts"]["union"],
        "target_shared_genes": lineage["counts"]["shared"],
        "target_sign_concordant_genes": lineage["counts"]["sign_concordant"],
        "target_screen_overlap_genes": lineage["counts"]["screen_overlap"],
        "target_shared_screen_genes": lineage["counts"]["shared_screen"],
        "target_final_genes": lineage["counts"]["final"],
    }
    if observed != expected:
        differences = {
            key: {"expected": expected.get(key), "observed": observed.get(key)}
            for key in sorted(set(expected) | set(observed))
            if expected.get(key) != observed.get(key)
        }
        raise InputError(f"source profile differs from frozen contract: {differences}")


def _cosine(first: np.ndarray, second: np.ndarray) -> float:
    first_norm = float(np.linalg.norm(first))
    second_norm = float(np.linalg.norm(second))
    if first_norm == 0 or second_norm == 0:
        return 0.0
    return float((first / first_norm) @ (second / second_norm))


def _metrics(prediction: np.ndarray, target: np.ndarray) -> dict[str, float]:
    target_norm = float(np.linalg.norm(target))
    prediction_norm = float(np.linalg.norm(prediction))
    if target_norm == 0:
        raise InputError("evaluation target has zero norm")
    nonzero = (prediction != 0) & (target != 0)
    sign_agreement = (
        float(np.mean(np.sign(prediction[nonzero]) == np.sign(target[nonzero])))
        if np.any(nonzero)
        else 0.0
    )
    return {
        "cosine": _cosine(prediction, target),
        "normalized_rmse": float(np.linalg.norm(prediction - target) / target_norm),
        "norm_ratio": prediction_norm / target_norm,
        "sign_agreement": sign_agreement,
    }


def _load_effect_subset(
    h5ad_path: Path,
    layer: str,
    rows: np.ndarray,
    all_genes: tuple[str, ...],
    selected_genes: tuple[str, ...],
) -> tuple[np.ndarray, dict[str, Any]]:
    import h5py

    gene_lookup = {gene: index for index, gene in enumerate(all_genes)}
    try:
        columns = np.asarray([gene_lookup[gene] for gene in selected_genes], dtype=np.intp)
    except KeyError as exc:
        raise InputError(f"selected target gene is absent from H5AD: {exc.args[0]}") from exc
    with h5py.File(h5ad_path, "r") as handle:
        full = np.asarray(handle[f"layers/{layer}"][rows, :], dtype=float)
    effects = full[:, columns]
    del full
    if not np.all(np.isfinite(effects)):
        raise InputError(f"selected {layer} effect block contains non-finite values")
    norms = np.linalg.norm(effects, axis=1)
    quality = {
        "shape": list(effects.shape),
        "finite_fraction": 1.0,
        "zero_fraction": float(np.mean(effects == 0)),
        "atom_norm_min": float(np.min(norms)),
        "atom_norm_median": float(np.median(norms)),
        "atom_norm_max": float(np.max(norms)),
        "zero_norm_atoms": int(np.count_nonzero(norms == 0)),
    }
    return effects, quality


def _split(n_genes: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    order = np.random.default_rng(seed).permutation(n_genes)
    middle = n_genes // 2
    return np.sort(order[:middle]), np.sort(order[middle:])


def _partition_hash(genes: tuple[str, ...], indices: np.ndarray) -> str:
    payload = "\n".join(genes[index] for index in indices).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def evaluate_source_bound_splits(
    effects: np.ndarray,
    target: np.ndarray,
    genes: tuple[str, ...],
    split_seeds: list[int],
) -> dict[str, Any]:
    rows = []
    for seed in split_seeds:
        _progress(f"source-bound z-score split seed={seed}")
        fit, score = _split(target.size, seed)
        result = held_out_alignment(effects, target, fit, score)
        rows.append(
            {
                "seed": seed,
                "held_out_cosine": result.held_out_cosine,
                "fit_cosine": result.fit_cosine,
                "support_size": int(np.count_nonzero(result.coefficients > 1e-8)),
                "fit_gene_sha256": _partition_hash(genes, fit),
                "score_gene_sha256": _partition_hash(genes, score),
            }
        )
    held_out = [row["held_out_cosine"] for row in rows]
    return {
        "splits": rows,
        "held_out_cosine_mean": statistics.mean(held_out),
        "held_out_cosine_sd": statistics.stdev(held_out),
        "held_out_cosine_min": min(held_out),
        "held_out_cosine_max": max(held_out),
    }


def _fit_baselines(
    effects: np.ndarray, target: np.ndarray, fit: np.ndarray
) -> dict[str, np.ndarray]:
    common_direction = np.mean(effects, axis=0, keepdims=True)
    common = project_cone(common_direction[:, fit], target[fit])
    common_prediction = common.coefficients @ common_direction

    fit_effects = effects[:, fit]
    target_fit = target[fit]
    target_norm = float(np.linalg.norm(target_fit))
    atom_norms = np.linalg.norm(fit_effects, axis=1)
    scores = np.full(effects.shape[0], -np.inf)
    usable = atom_norms > 0
    scores[usable] = (fit_effects[usable] @ target_fit) / (
        atom_norms[usable] * target_norm
    )
    best_index = int(np.argmax(scores))
    best = project_cone(effects[best_index : best_index + 1, fit], target_fit)
    best_prediction = best.coefficients @ effects[best_index : best_index + 1]
    return {
        "common_response": common_prediction,
        "best_single_atom": best_prediction,
    }


def source_transfer(
    effects: np.ndarray,
    source_directions: dict[str, np.ndarray],
    split_seeds: list[int],
    *,
    genes: tuple[str, ...] | None = None,
    label: str = "source_transfer",
) -> dict[str, Any]:
    if genes is None:
        genes = tuple(str(index) for index in range(effects.shape[1]))
    if len(genes) != effects.shape[1]:
        raise InputError("source-transfer genes do not match the effect coordinates")
    sources = tuple(source_directions)
    directions = ((sources[0], sources[1]), (sources[1], sources[0]))
    output = {}
    for train_source, test_source in directions:
        train_target = source_directions[train_source]
        test_target = source_directions[test_source]
        rows = []
        for seed in split_seeds:
            _progress(f"{label} {train_source}->{test_source} split seed={seed}")
            fit, score = _split(train_target.size, seed)
            cone = project_cone(effects[:, fit], train_target[fit])
            cone_prediction = cone.coefficients @ effects
            predictions = {"cone": cone_prediction}
            predictions.update(_fit_baselines(effects, train_target, fit))
            metrics = {
                name: _metrics(prediction[score], test_target[score])
                for name, prediction in predictions.items()
            }
            best_baseline = max(
                ("common_response", "best_single_atom"),
                key=lambda name: metrics[name]["cosine"],
            )
            rows.append(
                {
                    "seed": seed,
                    "metrics": metrics,
                    "best_baseline": best_baseline,
                    "cosine_improvement_over_best_baseline": (
                        metrics["cone"]["cosine"] - metrics[best_baseline]["cosine"]
                    ),
                    "cone_support_size": int(
                        np.count_nonzero(cone.coefficients > 1e-8)
                    ),
                    "fit_gene_sha256": _partition_hash(genes, fit),
                    "score_gene_sha256": _partition_hash(genes, score),
                }
            )
        improvements = [
            row["cosine_improvement_over_best_baseline"] for row in rows
        ]
        output[f"{train_source}_to_{test_source}"] = {
            "splits": rows,
            "cosine_improvement_mean": statistics.mean(improvements),
            "cosine_improvement_sd": statistics.stdev(improvements),
            "all_improvements_positive": all(value > 0 for value in improvements),
        }
    return output


def build_profile(
    config: dict[str, Any], *, verify_hash: bool
) -> tuple[dict[str, Any], dict[str, Any]]:
    input_status = {
        name: verify_input(name, spec, verify_hash=verify_hash)
        for name, spec in config["inputs"].items()
    }
    h5ad_path = ROOT / config["inputs"]["de_stats"]["path"]
    target_path = ROOT / config["inputs"]["target_table"]["path"]
    h5_profile, axes = profile_h5ad(h5ad_path, config["screen"])
    target_rows, target_values = load_target_table(target_path, config["target"])
    lineages = {
        field: build_target_lineage(
            target_values[field],
            axes["genes"],
            orientation_multiplier=config["target"]["orientation_multiplier"],
        )
        for field in config["target"]["value_fields"]
    }
    _assert_expected_profile(
        h5_profile, lineages["zscore"], config["expected_profile"]
    )
    profile = {
        "input_verification": input_status,
        "h5ad": h5_profile,
        "target_table_rows": target_rows,
        "target_lineage": {
            field: {
                "counts": lineage["counts"],
                "between_source_cosine": lineage["between_source_cosine"],
                "between_source_screen_cosine": lineage[
                    "between_source_screen_cosine"
                ],
            }
            for field, lineage in lineages.items()
        },
    }
    internals = {"axes": axes, "lineages": lineages}
    return profile, internals


def run(config_path: Path, *, verify_hash: bool = True, profile_only: bool = False) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "1.0.0":
        raise InputError("unsupported source-reconstruction config schema")
    _progress("verifying immutable inputs and profiling source grain")
    profile, internals = build_profile(config, verify_hash=verify_hash)
    report = {
        "schema_version": "1.0.0",
        "generated_on": "2026-07-17",
        "config_sha256": sha256(config_path),
        "scope": (
            "source-bound aggregate-screen reconstruction and cross-source diagnostic; "
            "not donor, guide, trajectory, or biological validation"
        ),
        "data_quality": profile,
    }
    if profile_only:
        report["status"] = "PROFILE_PASS" if verify_hash else "PROFILE_UNVERIFIED"
        return report

    h5ad_path = ROOT / config["inputs"]["de_stats"]["path"]
    axes = internals["axes"]
    z_lineage = internals["lineages"]["zscore"]
    source_config = config["source_bound_splits"]
    _progress("loading the source-bound z-score block")
    source_effects, source_quality = _load_effect_subset(
        h5ad_path,
        source_config["layer"],
        axes["selected_rows"],
        axes["genes"],
        z_lineage["final_genes"],
    )
    source_bound = evaluate_source_bound_splits(
        source_effects,
        z_lineage["merged_direction"],
        z_lineage["final_genes"],
        source_config["split_seeds"],
    )
    del source_effects
    expected_splits = source_config["expected_held_out_cosines"]
    observed_splits = [row["held_out_cosine"] for row in source_bound["splits"]]
    if len(observed_splits) != len(expected_splits):
        raise InputError("source-bound split count differs from the frozen contract")
    split_errors = [
        abs(observed - expected)
        for observed, expected in zip(observed_splits, expected_splits)
    ]
    source_bound["max_absolute_error_vs_frozen_splits"] = max(split_errors)
    source_bound["status"] = (
        "PASS"
        if max(split_errors) <= source_config["absolute_tolerance"]
        else "FAIL"
    )
    archived = source_config["retired_archived_comparison"]
    source_bound["effect_quality"] = source_quality
    source_bound["protocol"] = {
        "gene_order": source_config["gene_order"],
        "rng": source_config["rng"],
    }
    source_bound["retired_archived_comparison"] = {
        **archived,
        "mean_delta_source_bound_minus_archived": (
            source_bound["held_out_cosine_mean"] - archived["mean"]
        ),
        "fixed_split_absolute_error": abs(
            source_bound["splits"][0]["held_out_cosine"]
            - archived["fixed_split_cosine"]
        ),
        "status": "RETIRED_UNRECONSTRUCTABLE_PROTOCOL",
    }
    report["source_bound_splits"] = source_bound

    transfer_results = {}
    for layer in config["source_transfer"]["layers"]:
        _progress(f"loading source-transfer layer={layer}")
        lineage = internals["lineages"][layer]
        effects, quality = _load_effect_subset(
            h5ad_path,
            layer,
            axes["selected_rows"],
            axes["genes"],
            lineage["shared_screen_genes"],
        )
        transfer_results[layer] = {
            "effect_quality": quality,
            "between_source_target_cosine": lineage["between_source_cosine"],
            "between_source_screen_target_cosine": lineage[
                "between_source_screen_cosine"
            ],
            "directions": source_transfer(
                effects,
                lineage["source_directions"],
                config["source_transfer"]["split_seeds"],
                genes=lineage["shared_screen_genes"],
                label=layer,
            ),
        }
        del effects
    report["source_transfer"] = transfer_results
    report["status"] = (
        "PASS" if source_bound["status"] == "PASS" and verify_hash else "FAIL"
    )
    report["claim_ceiling"] = (
        "Source transfer uses donor-collapsed effects and shared target-table coordinates. "
        "It tests construction sensitivity, not donor generalization or state conversion."
    )
    return report


def _assert_matches(actual: Any, expected: Any, path: str = "report") -> None:
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise AssertionError(f"{path} keys differ")
        for key in expected:
            _assert_matches(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise AssertionError(f"{path} list length differs")
        for index, value in enumerate(expected):
            _assert_matches(actual[index], value, f"{path}[{index}]")
    elif isinstance(expected, float):
        if not np.isclose(actual, expected, rtol=1e-9, atol=1e-12):
            raise AssertionError(f"{path}: {actual} != {expected}")
    elif actual != expected:
        raise AssertionError(f"{path}: {actual} != {expected}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--profile", action="store_true")
    mode.add_argument("--write", type=Path)
    mode.add_argument("--check", type=Path)
    parser.add_argument(
        "--skip-hash",
        action="store_true",
        help="development-only profile shortcut; a full report cannot pass without hashes",
    )
    args = parser.parse_args()
    report = run(
        args.config,
        verify_hash=not args.skip_hash,
        profile_only=args.profile,
    )
    if args.write:
        args.write.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.write}")
    elif args.check:
        expected = json.loads(args.check.read_text(encoding="utf-8"))
        _assert_matches(report, expected)
        print(f"source reconstruction: {report['status']} (frozen report matches)")
    else:
        print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] == "FAIL":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
