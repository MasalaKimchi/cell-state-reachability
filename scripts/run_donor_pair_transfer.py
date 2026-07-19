#!/usr/bin/env python3
"""Run frozen-weight source-and-donor-pair transfer on the released Zhu H5MU.

The released modalities contain two-donor DE estimates and omit DE-ineligible targets.
This runner therefore measures fixed-cohort robustness under the published eligibility
pipeline. It is not a leakage-free donor holdout or donor-population analysis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError, project_cone  # noqa: E402
from scripts.run_source_reconstruction import (  # noqa: E402
    build_target_lineage,
    load_target_table,
)


DEFAULT_CONFIG = ROOT / "configs" / "donor_pair_transfer.json"


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
    observed = sha256_file(path) if verify_hash else None
    if verify_hash and observed != spec["sha256"]:
        raise InputError(f"{spec['path']} SHA-256 differs from the frozen identity")
    return {
        "path": spec["path"],
        "bytes": size,
        "sha256_expected": spec["sha256"],
        "sha256_actual": observed,
        "hash_verified": verify_hash,
    }


def _decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values],
        dtype=object,
    )


def _categorical(frame: Any, field: str) -> np.ndarray:
    column = frame[field]
    if not hasattr(column, "keys") or "categories" not in column or "codes" not in column:
        raise InputError(f"{field} must use categorical H5 encoding")
    categories = _decode(column["categories"][:])
    codes = np.asarray(column["codes"][:], dtype=np.int64)
    if np.any(codes < 0) or np.any(codes >= len(categories)):
        raise InputError(f"{field} contains missing or invalid categorical codes")
    return categories[codes]


def validate_partitions(config: dict[str, Any]) -> None:
    donor_runs = {values["id"]: values["rest_run"] for values in config["donors"].values()}
    donor_ids = set(donor_runs)
    modalities = set(config["h5mu"]["modalities"])
    seen_ids = set()
    seen_pairings = set()
    for partition in config["partitions"]:
        if partition["id"] in seen_ids:
            raise InputError("donor partition IDs must be unique")
        seen_ids.add(partition["id"])
        if partition["left"] not in modalities or partition["right"] not in modalities:
            raise InputError("donor partition references an undeclared modality")
        left = frozenset(partition["left"].split("_"))
        right = frozenset(partition["right"].split("_"))
        if len(left) != 2 or len(right) != 2 or left & right or left | right != donor_ids:
            raise InputError("donor partition must contain complementary two-donor modalities")
        pairing = frozenset((left, right))
        if pairing in seen_pairings:
            raise InputError("donor partitions must be unique")
        seen_pairings.add(pairing)
        left_runs = {donor_runs[donor] for donor in left}
        right_runs = {donor_runs[donor] for donor in right}
        run_confounded = len(left_runs) == len(right_runs) == 1 and left_runs != right_runs
        if partition["run_confounded"] != run_confounded:
            raise InputError("donor partition run-confounding label differs from donor metadata")
    if len(seen_pairings) != 3:
        raise InputError("exactly three complementary donor partitions are required")


def profile_h5mu(path: Path, config: dict[str, Any]) -> dict[str, Any]:
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError("donor-pair transfer requires requirements-external.txt") from exc

    spec = config["h5mu"]
    expected = config["expected"]
    modalities = tuple(spec["modalities"])
    profiles = {}
    reference_gene_ids = reference_symbols = reference_rest_keys = None
    with h5py.File(path, "r") as handle:
        if handle.attrs.get("encoding-type") != "MuData" or "mod" not in handle:
            raise InputError("donor object is not a MuData H5MU")
        observed_modalities = tuple(sorted(handle["mod"].keys()))
        if observed_modalities != tuple(sorted(modalities)):
            raise InputError("donor H5MU modalities differ from the frozen six pairs")
        for modality in modalities:
            group = handle[f"mod/{modality}"]
            if group.attrs.get("encoding-type") != "anndata":
                raise InputError(f"{modality} is not encoded as AnnData")
            obs = group["obs"]
            var = group["var"]
            index_field = obs.attrs.get("_index")
            if not isinstance(index_field, str) or index_field not in obs:
                raise InputError(f"{modality} observation index is missing")
            keys = _decode(obs[index_field][:])
            conditions = _categorical(obs, spec["condition_field"])
            gene_ids = _decode(var[spec["gene_id_field"]][:])
            symbols = _decode(var[spec["gene_symbol_field"]][:])
            n_rows = len(keys)
            if n_rows != expected["modality_rows"][modality]:
                raise InputError(f"{modality} row count differs")
            if len(gene_ids) != expected["genes"] or len(symbols) != expected["genes"]:
                raise InputError(f"{modality} gene count differs")
            if len(set(keys)) != n_rows:
                raise InputError(f"{modality} has duplicate observation keys")
            if len(set(gene_ids)) != len(gene_ids) or len(set(symbols)) != len(symbols):
                raise InputError(f"{modality} gene IDs/symbols must each be unique")
            for layer in spec["required_layers"]:
                if layer not in group["layers"]:
                    raise InputError(f"{modality} is missing layer {layer}")
                if group[f"layers/{layer}"].shape != (n_rows, len(gene_ids)):
                    raise InputError(f"{modality}/{layer} shape differs")
            rest_rows = np.flatnonzero(conditions == spec["condition"])
            rest_keys = keys[rest_rows]
            if len(rest_keys) != expected["rest_atoms"]:
                raise InputError(f"{modality} Rest atom count differs")
            if len(set(rest_keys)) != len(rest_keys):
                raise InputError(f"{modality} Rest keys are not unique")
            canonical_keys = tuple(sorted(str(item) for item in rest_keys))
            row_for_key = {str(key): int(row) for key, row in zip(rest_keys, rest_rows)}
            aligned_rows = np.asarray([row_for_key[key] for key in canonical_keys], dtype=int)
            for layer in spec["required_layers"]:
                values = np.asarray(group[f"layers/{layer}"][rest_rows, :], dtype=float)
                if not np.all(np.isfinite(values)):
                    raise InputError(f"{modality}/{layer} Rest block contains non-finite values")
            if reference_gene_ids is None:
                reference_gene_ids = gene_ids
                reference_symbols = symbols
                reference_rest_keys = canonical_keys
            elif (
                not np.array_equal(gene_ids, reference_gene_ids)
                or not np.array_equal(symbols, reference_symbols)
                or canonical_keys != reference_rest_keys
            ):
                raise InputError("donor modalities do not share identical genes and Rest atoms")
            profiles[modality] = {
                "rows": n_rows,
                "rest_rows": aligned_rows,
                "condition_counts": {
                    condition: int(np.sum(conditions == condition))
                    for condition in sorted(set(conditions))
                },
            }
    return {
        "modalities": profiles,
        "gene_ids": tuple(str(item) for item in reference_gene_ids),
        "gene_symbols": tuple(str(item) for item in reference_symbols),
        "rest_atom_keys": reference_rest_keys,
    }


def load_modality_matrix(
    path: Path,
    modality: str,
    rows: np.ndarray,
    columns: np.ndarray,
    layer: str,
) -> np.ndarray:
    import h5py

    order = np.argsort(rows)
    sorted_rows = rows[order]
    with h5py.File(path, "r") as handle:
        values = np.asarray(
            handle[f"mod/{modality}/layers/{layer}"][sorted_rows, :], dtype=float
        )
    values = values[np.argsort(order)][:, columns]
    if not np.all(np.isfinite(values)):
        raise InputError(f"{modality} selected matrix contains non-finite values")
    return values


def _vector_hash(values: np.ndarray) -> str:
    array = np.asarray(values, dtype="<f8")
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def _names_hash(names: list[str] | tuple[str, ...]) -> str:
    return hashlib.sha256(("\n".join(names) + "\n").encode("utf-8")).hexdigest()


def fit_training_models(
    train: np.ndarray, target: np.ndarray, atom_keys: tuple[str, ...]
) -> dict[str, Any]:
    train = np.asarray(train, dtype=float)
    target = np.asarray(target, dtype=float)
    if train.shape != (len(atom_keys), len(target)) or not np.all(np.isfinite(train)):
        raise InputError("training matrix, target, and atom keys do not align")
    fit = project_cone(train, target)
    common = np.mean(train, axis=0)
    common_denominator = float(common @ common)
    common_alpha = (
        max(0.0, float(common @ target) / common_denominator)
        if common_denominator > 0
        else 0.0
    )
    atom_norms = np.sum(train * train, axis=1)
    atom_dot = train @ target
    atom_alpha = np.divide(
        np.maximum(atom_dot, 0.0),
        atom_norms,
        out=np.zeros_like(atom_dot),
        where=atom_norms > 0,
    )
    target_norm = float(target @ target)
    atom_loss = target_norm - 2 * atom_alpha * atom_dot + atom_alpha**2 * atom_norms
    best_loss = float(np.min(atom_loss))
    tied = np.flatnonzero(np.isclose(atom_loss, best_loss, rtol=0, atol=1e-12))
    best_index = int(tied[0])
    return {
        "cone_coefficients": fit.coefficients,
        "cone_fit_cosine": fit.cosine,
        "cone_kkt_violation": fit.kkt_violation,
        "common_alpha": common_alpha,
        "best_single_index": best_index,
        "best_single_key": atom_keys[best_index],
        "best_single_alpha": float(atom_alpha[best_index]),
    }


def _metrics(prediction: np.ndarray, target: np.ndarray) -> dict[str, float]:
    prediction = np.asarray(prediction, dtype=float)
    target = np.asarray(target, dtype=float)
    if prediction.shape != target.shape or not np.all(np.isfinite(prediction)):
        raise InputError("prediction and test target do not align")
    target_norm = float(np.linalg.norm(target))
    prediction_norm = float(np.linalg.norm(prediction))
    if target_norm == 0:
        raise InputError("test target has zero norm")
    cosine = (
        float(prediction @ target / (prediction_norm * target_norm))
        if prediction_norm > 0
        else 0.0
    )
    nonzero = target != 0
    return {
        "cosine": cosine,
        "normalized_rmse": float(
            np.sqrt(np.mean((prediction - target) ** 2))
            / np.sqrt(np.mean(target**2))
        ),
        "norm_ratio": prediction_norm / target_norm,
        "sign_agreement": float(
            np.mean(np.sign(prediction[nonzero]) == np.sign(target[nonzero]))
        ),
    }


def score_frozen_models(
    models: dict[str, Any], test: np.ndarray, target: np.ndarray
) -> dict[str, dict[str, float]]:
    test = np.asarray(test, dtype=float)
    if test.ndim != 2 or test.shape[1] != len(target):
        raise InputError("test dictionary and target do not align")
    predictions = {
        "cone": models["cone_coefficients"] @ test,
        "training_common_ray": models["common_alpha"] * np.mean(test, axis=0),
        "training_best_single": models["best_single_alpha"]
        * test[models["best_single_index"]],
        "zero": np.zeros(len(target), dtype=float),
    }
    return {name: _metrics(values, target) for name, values in predictions.items()}


def run_challenge(
    train: np.ndarray,
    test: np.ndarray,
    train_target: np.ndarray,
    test_target: np.ndarray,
    atom_keys: tuple[str, ...],
) -> tuple[dict[str, Any], dict[str, Any]]:
    models = fit_training_models(train, train_target, atom_keys)
    metrics = score_frozen_models(models, test, test_target)
    cone = metrics["cone"]
    comparisons = {}
    for baseline in ("zero", "training_common_ray", "training_best_single"):
        comparisons[baseline] = {
            "cosine_improvement": cone["cosine"] - metrics[baseline]["cosine"],
            "normalized_rmse_improvement": metrics[baseline]["normalized_rmse"]
            - cone["normalized_rmse"],
        }
    serializable = {
        "metrics": metrics,
        "comparisons": comparisons,
        "training": {
            "cone_fit_cosine": models["cone_fit_cosine"],
            "cone_kkt_violation": models["cone_kkt_violation"],
            "cone_support": int(np.count_nonzero(models["cone_coefficients"])),
            "coefficient_sha256": _vector_hash(models["cone_coefficients"]),
            "common_alpha": models["common_alpha"],
            "best_single_key": models["best_single_key"],
            "best_single_alpha": models["best_single_alpha"],
        },
    }
    return serializable, models


def _summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    output = {}
    for scope, selected in (
        ("all_partitions", rows),
        ("run_balanced_only", [row for row in rows if not row["run_confounded"]]),
        ("run_confounded_only", [row for row in rows if row["run_confounded"]]),
    ):
        scope_output = {"challenges": len(selected)}
        for model in ("cone", "training_common_ray", "training_best_single", "zero"):
            for metric in ("cosine", "normalized_rmse"):
                values = [row["metrics"][model][metric] for row in selected]
                scope_output[f"median_{model}_{metric}"] = statistics.median(values)
        for baseline in ("zero", "training_common_ray", "training_best_single"):
            for metric in ("cosine_improvement", "normalized_rmse_improvement"):
                values = [row["comparisons"][baseline][metric] for row in selected]
                scope_output[f"median_{metric}_over_{baseline}"] = statistics.median(values)
                scope_output[f"fraction_{metric}_positive_over_{baseline}"] = statistics.mean(
                    value > 0 for value in values
                )
        output[scope] = scope_output
    return output


def run(config: dict[str, Any], *, verify_hash: bool = True) -> dict[str, Any]:
    validate_partitions(config)
    identities = {
        name: verify_input(spec, verify_hash=verify_hash)
        for name, spec in config["inputs"].items()
    }
    h5_path = ROOT / config["inputs"]["donor_h5mu"]["path"]
    profile = profile_h5mu(h5_path, config)
    _, target_values = load_target_table(
        ROOT / config["inputs"]["target_table"]["path"], config["target"]
    )
    lineage = build_target_lineage(
        target_values["log_fc"],
        profile["gene_symbols"],
        orientation_multiplier=config["target"]["orientation_multiplier"],
    )
    genes = lineage["shared_screen_genes"]
    if len(genes) != config["expected"]["shared_target_genes"]:
        raise InputError("donor/target shared gene count differs")
    column_for_gene = {gene: index for index, gene in enumerate(profile["gene_symbols"])}
    columns = np.asarray([column_for_gene[gene] for gene in genes], dtype=int)
    atom_keys = profile["rest_atom_keys"]
    target_directions = lineage["source_directions"]

    matrices = {}
    for modality in config["h5mu"]["modalities"]:
        matrices[modality] = load_modality_matrix(
            h5_path,
            modality,
            profile["modalities"][modality]["rest_rows"],
            columns,
            config["h5mu"]["layer"],
        )

    rows = []
    for partition in config["partitions"]:
        for train_modality, test_modality in (
            (partition["left"], partition["right"]),
            (partition["right"], partition["left"]),
        ):
            for train_source, test_source in config["analysis"]["directions"]:
                for seed in config["analysis"]["split_seeds"]:
                    rng = np.random.default_rng(seed)
                    permutation = rng.permutation(len(genes))
                    midpoint = len(genes) // 2
                    fit_indices = np.sort(permutation[:midpoint])
                    score_indices = np.sort(permutation[midpoint:])
                    result, _ = run_challenge(
                        matrices[train_modality][:, fit_indices],
                        matrices[test_modality][:, score_indices],
                        target_directions[train_source][fit_indices],
                        target_directions[test_source][score_indices],
                        atom_keys,
                    )
                    result.update(
                        {
                            "partition": partition["id"],
                            "run_confounded": partition["run_confounded"],
                            "train_modality": train_modality,
                            "test_modality": test_modality,
                            "train_target_source": train_source,
                            "test_target_source": test_source,
                            "seed": seed,
                            "fit_genes": len(fit_indices),
                            "score_genes": len(score_indices),
                            "fit_gene_sha256": _names_hash([genes[index] for index in fit_indices]),
                            "score_gene_sha256": _names_hash(
                                [genes[index] for index in score_indices]
                            ),
                        }
                    )
                    rows.append(result)

    return {
        "schema_version": "1.0.0",
        "status": "PASS",
        "benchmark": config["benchmark"],
        "claim_ceiling": config["claim_ceiling"],
        "input_verification": identities,
        "data_quality": {
            "h5mu_encoding": "MuData 0.1.0",
            "modalities": {
                modality: {
                    "rows": values["rows"],
                    "condition_counts": values["condition_counts"],
                }
                for modality, values in profile["modalities"].items()
            },
            "genes": len(profile["gene_symbols"]),
            "rest_atoms_complete_all_modalities": len(atom_keys),
            "shared_target_genes": len(genes),
            "between_source_target_cosine": lineage["between_source_screen_cosine"],
            "eligibility_warning": "Released modalities contain published two-donor DE-eligible targets; complete presence can select on effectiveness.",
        },
        "protocol": {
            "condition": config["h5mu"]["condition"],
            "layer": config["h5mu"]["layer"],
            "target_gene_universe": "both target sources and donor-H5MU gene symbols; no held-out-source sign filter",
            "gene_order": config["analysis"]["gene_order"],
            "split_rng": config["analysis"]["rng"],
            "split_seeds": config["analysis"]["split_seeds"],
            "coefficient_application": "fit on one donor-pair dictionary and source target; apply identical coefficients to the complementary donor-pair dictionary and other source target on held-out genes",
            "baseline_application": "select atom/scalar only on training genes, donor pair, and target source; apply unchanged identity/scalar to test dictionary",
            "inference": config["analysis"]["inference"],
        },
        "challenges": rows,
        "summary": _summary(rows),
        "limitations": [
            "The six modalities are released two-donor summaries, not individual-donor effects.",
            "Published DE eligibility determines modality presence and can select on perturbation effectiveness.",
            "Three donor partitions, both directions, two target-source directions, and gene splits are correlated challenges, not independent replicates.",
            "D1+D2 versus D3+D4 is fully confounded with Rest sequencing run; mixed-run partitions are reported separately.",
            "Four donors cannot support donor-population significance at 0.05 by an exact donor sign-flip test.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--skip-hash", action="store_true", help="development only")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    report = run(config, verify_hash=not args.skip_hash)
    report["config_sha256"] = sha256_file(args.config)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    output = ROOT / config["output"]
    if args.check:
        if not output.is_file() or output.read_text(encoding="utf-8") != payload:
            raise SystemExit(f"maintained output differs: {config['output']}")
        print("donor-pair transfer output matches")
    else:
        output.write_text(payload, encoding="utf-8")
        print(f"wrote {config['output']}")


if __name__ == "__main__":
    main()
