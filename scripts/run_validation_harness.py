#!/usr/bin/env python3
"""Run the deterministic, data-free systemic validation harness."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError, held_out_alignment, project_cone
from validation import (
    LabeledEffects,
    LabeledTarget,
    Provenance,
    active_set_oracle,
    align_labeled_problem,
    binomial_upper_confidence_bound,
    grouped_gene_splits,
    max_t_empirical_p,
)


DEFAULT_CONFIG = ROOT / "configs" / "validation_harness.json"
DEFAULT_REPORT = ROOT / "results" / "validation_harness.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _gate(value: float, threshold: float) -> bool:
    return bool(np.isfinite(value) and value <= threshold)


def _cosine(first: np.ndarray, second: np.ndarray) -> float:
    norm_first = float(np.linalg.norm(first))
    norm_second = float(np.linalg.norm(second))
    if norm_first == 0 or norm_second == 0:
        return 0.0
    return float((first / norm_first) @ (second / norm_second))


def _specificity_scores(
    effects: np.ndarray,
    target: np.ndarray,
    fit: np.ndarray,
    score: np.ndarray,
) -> tuple[float, float, float]:
    cone = project_cone(effects[:, fit], target[fit])
    cone_prediction = cone.coefficients @ effects
    common_direction = np.mean(effects, axis=0, keepdims=True)
    common = project_cone(common_direction[:, fit], target[fit])
    common_prediction = common.coefficients @ common_direction
    cone_cosine = _cosine(cone_prediction[score], target[score])
    common_cosine = _cosine(common_prediction[score], target[score])
    return cone_cosine, common_cosine, cone_cosine - common_cosine


def _oracle_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    fitted_errors = []
    objective_errors = []
    for trial in range(config["oracle_trials"]):
        n_atoms = int(rng.integers(2, 8))
        n_genes = int(rng.integers(n_atoms + 2, 22))
        effects = rng.normal(size=(n_atoms, n_genes))
        if trial % 3 == 0:
            effects *= np.geomspace(1e-12, 1e12, n_atoms)[:, None]
        if trial % 7 == 0 and n_atoms > 2:
            effects[-1] = effects[-2]
        weights = np.exp(rng.uniform(-2, 2, size=n_genes))
        if trial % 2:
            target = rng.uniform(0, 2, size=n_atoms) @ effects
        else:
            target = rng.normal(size=n_genes)
        production = project_cone(effects, target, gene_weights=weights)
        _, fitted, relative_objective = active_set_oracle(
            effects, target, gene_weights=weights, max_atoms=7
        )
        scale = max(float(np.linalg.norm(target)), 1.0)
        fitted_errors.append(float(np.linalg.norm(production.fitted - fitted) / scale))
        objective_errors.append(
            float(abs(production.relative_objective - relative_objective))
        )
    max_fitted = max(fitted_errors)
    max_objective = max(objective_errors)
    thresholds = config["thresholds"]
    passed = _gate(max_fitted, thresholds["max_oracle_fitted_error"]) and _gate(
        max_objective, thresholds["max_oracle_objective_error"]
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "trials": config["oracle_trials"],
        "max_relative_fitted_error": max_fitted,
        "max_relative_objective_error": max_objective,
    }


def _degeneracy_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    effects = rng.normal(size=(6, 30))
    target = rng.normal(size=30)
    base = project_cone(effects, target).fitted
    scales = np.geomspace(1e-8, 1e8, effects.shape[0])
    scaled = project_cone(effects * scales[:, None], target).fitted
    augmented = np.vstack([effects, effects[2], effects[2] * (1 + 1e-12), np.zeros(30)])
    duplicated = project_cone(augmented, target).fitted
    scale = max(float(np.linalg.norm(target)), 1.0)
    errors = [
        float(np.linalg.norm(base - scaled) / scale),
        float(np.linalg.norm(base - duplicated) / scale),
    ]
    maximum = max(errors)
    passed = _gate(maximum, config["thresholds"]["max_degenerate_fitted_error"])
    return {
        "status": "PASS" if passed else "FAIL",
        "max_relative_fitted_error": maximum,
        "challenges": ["atom_rescaling_1e-8_to_1e8", "duplicate_near_duplicate_zero"],
    }


def _label_scenario() -> dict[str, Any]:
    provenance = Provenance(
        dataset_id="synthetic-effects-v1",
        source_sha256="0" * 64,
        gene_namespace="HGNC-symbol",
        units="z_score",
        modality="synthetic-CRISPRi",
        context="synthetic-rest",
        timepoint="synthetic",
        orientation="perturbations_by_genes",
    )
    target_provenance = replace(
        provenance,
        dataset_id="synthetic-target-v1",
        source_sha256="1" * 64,
        modality="RNA-contrast",
        orientation="genes",
    )
    effects = LabeledEffects(
        values=np.eye(4),
        perturbations=("p1", "p2", "p3", "p4"),
        genes=("G1", "G2", "G3", "G4"),
        provenance=provenance,
    )
    target = LabeledTarget(
        values=np.array([1.0, 2.0, 3.0, 4.0]),
        genes=("G1", "G2", "G3", "G4"),
        provenance=target_provenance,
    )
    reordered = replace(target, values=target.values[::-1], genes=target.genes[::-1])
    aligned = align_labeled_problem(effects, reordered)
    reorder_safe = bool(np.array_equal(aligned.target, target.values))

    corruptions: list[Callable[[], object]] = [
        lambda: align_labeled_problem(replace(effects, genes=("G1", "G1", "G3", "G4")), target),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, gene_namespace="Ensembl"))),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, units="log2_fc"))),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, orientation="genes_by_one"))),
        lambda: align_labeled_problem(effects, replace(target, genes=target.genes[:-1], values=target.values[:-1])),
        lambda: align_labeled_problem(effects, replace(target, provenance=replace(target.provenance, source_sha256="unknown"))),
        lambda: align_labeled_problem(replace(effects, values=np.where(np.eye(4) == 1, np.nan, 0)), target),
        lambda: align_labeled_problem(replace(effects, values=np.eye(3)), target),
    ]
    caught = 0
    for corruption in corruptions:
        try:
            corruption()
        except InputError:
            caught += 1
    passed = reorder_safe and caught == len(corruptions)
    return {
        "status": "PASS" if passed else "FAIL",
        "safe_reordering": reorder_safe,
        "faults_caught": caught,
        "faults_injected": len(corruptions),
    }


def _grouped_split_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    n_groups = 14
    genes_per_group = 6
    groups = tuple(f"module_{group}" for group in range(n_groups) for _ in range(genes_per_group))
    module_effects = rng.normal(size=(7, n_groups))
    effects = np.repeat(module_effects, genes_per_group, axis=1) + rng.normal(
        scale=0.05, size=(7, n_groups * genes_per_group)
    )
    target = rng.uniform(0.2, 1.5, size=7) @ effects + rng.normal(
        scale=0.1, size=n_groups * genes_per_group
    )
    scores = []
    leakage = 0
    for fit, score in grouped_gene_splits(
        groups, n_splits=config["grouped_splits"], seed=config["seed"] + 1
    ):
        fit_groups = {groups[index] for index in fit}
        score_groups = {groups[index] for index in score}
        leakage += len(fit_groups & score_groups)
        scores.append(held_out_alignment(effects, target, fit, score).held_out_cosine)
    passed = leakage <= config["thresholds"]["max_group_leakage"] and np.all(
        np.isfinite(scores)
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "splits": config["grouped_splits"],
        "group_leakage_count": leakage,
        "held_out_cosine_mean": float(np.mean(scores)),
        "held_out_cosine_sd": float(np.std(scores, ddof=1)),
        "held_out_cosine_min": float(np.min(scores)),
    }


def _null_scenario(config: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    rejections = 0
    adjusted_not_below_marginal = True
    for _ in range(config["null_trials"]):
        observed = rng.normal(size=config["n_hypotheses"])
        null = rng.normal(
            size=(config["null_resamples"], config["n_hypotheses"])
        )
        adjusted = max_t_empirical_p(observed, null)
        marginal = np.asarray(
            [
                (1 + np.count_nonzero(null[:, index] >= value))
                / (config["null_resamples"] + 1)
                for index, value in enumerate(observed)
            ]
        )
        adjusted_not_below_marginal &= bool(np.all(adjusted >= marginal))
        rejections += int(np.any(adjusted <= 0.05))
    rate = rejections / config["null_trials"]
    upper = binomial_upper_confidence_bound(rejections, config["null_trials"])
    passed = adjusted_not_below_marginal and _gate(
        upper, config["thresholds"]["max_null_familywise_error_upper_95"]
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "trials": config["null_trials"],
        "hypotheses_per_trial": config["n_hypotheses"],
        "resamples_per_trial": config["null_resamples"],
        "familywise_error_rate_at_0.05": rate,
        "familywise_rejections": rejections,
        "familywise_error_upper_95": upper,
        "adjusted_p_never_below_marginal": adjusted_not_below_marginal,
    }


def _structured_specificity_scenario(
    config: dict[str, Any], rng: np.random.Generator
) -> dict[str, Any]:
    """Expose common-response, module-leakage, and sign-selection optimism."""

    trials = config["structured_trials"]
    n_modules = config["structured_modules"]
    genes_per_module = config["structured_genes_per_module"]
    n_atoms = config["structured_atoms"]
    common_strength = config["structured_common_strength"]
    nuisance_strength = config["structured_nuisance_strength"]
    n_genes = n_modules * genes_per_module
    module_ids = np.repeat(np.arange(n_modules), genes_per_module)

    null_raw = []
    null_common = []
    null_random_delta = []
    null_group_delta = []
    alternative_group_delta = []
    selection_inflation = []

    for _ in range(trials):
        common = np.repeat(rng.normal(size=n_modules), genes_per_module)
        common += rng.normal(scale=0.15, size=n_genes)
        common /= np.std(common)

        atom_modules = rng.normal(size=(n_atoms, n_modules))
        atom_modules -= np.mean(atom_modules, axis=0, keepdims=True)
        specific = np.repeat(atom_modules, genes_per_module, axis=1)
        specific += rng.normal(scale=0.05, size=(n_atoms, n_genes))
        effects = common_strength * common + specific

        nuisance = np.repeat(rng.normal(size=n_modules), genes_per_module)
        null_target = common_strength * common + nuisance_strength * nuisance
        null_target += rng.normal(scale=0.10, size=n_genes)

        random_order = rng.permutation(n_genes)
        random_fit = np.sort(random_order[: n_genes // 2])
        random_score = np.sort(random_order[n_genes // 2 :])
        module_order = rng.permutation(n_modules)
        fit_modules = set(module_order[: n_modules // 2])
        group_fit = np.flatnonzero(
            np.fromiter((item in fit_modules for item in module_ids), dtype=bool)
        )
        group_score = np.flatnonzero(
            np.fromiter((item not in fit_modules for item in module_ids), dtype=bool)
        )

        random_metrics = _specificity_scores(
            effects, null_target, random_fit, random_score
        )
        group_metrics = _specificity_scores(
            effects, null_target, group_fit, group_score
        )
        null_raw.append(random_metrics[0])
        null_common.append(random_metrics[1])
        null_random_delta.append(random_metrics[2])
        null_group_delta.append(group_metrics[2])

        support = rng.choice(n_atoms, size=3, replace=False)
        support_weights = rng.dirichlet(np.ones(3))
        alternative_target = support_weights @ effects[support]
        alternative_target += rng.normal(scale=0.10, size=n_genes)
        alternative_group_delta.append(
            _specificity_scores(
                effects, alternative_target, group_fit, group_score
            )[2]
        )

        source_a = rng.normal(size=n_genes)
        source_b = rng.normal(size=n_genes)
        unselected = _cosine(source_a, source_b)
        selected = np.sign(source_a) == np.sign(source_b)
        selected_cosine = _cosine(source_a[selected], source_b[selected])
        selection_inflation.append(selected_cosine - unselected)

    summary = {
        "trials": trials,
        "modules": n_modules,
        "genes_per_module": genes_per_module,
        "atoms": n_atoms,
        "common_strength": common_strength,
        "nuisance_strength": nuisance_strength,
        "null_raw_random_gene_cosine_mean": float(np.mean(null_raw)),
        "null_common_response_cosine_mean": float(np.mean(null_common)),
        "null_random_gene_improvement_mean": float(np.mean(null_random_delta)),
        "null_module_holdout_improvement_mean": float(np.mean(null_group_delta)),
        "random_gene_optimism_gap": float(
            np.mean(null_random_delta) - np.mean(null_group_delta)
        ),
        "alternative_module_holdout_improvement_mean": float(
            np.mean(alternative_group_delta)
        ),
        "sign_selection_cosine_inflation_mean": float(np.mean(selection_inflation)),
    }
    thresholds = config["thresholds"]
    passed = (
        summary["null_raw_random_gene_cosine_mean"]
        >= thresholds["min_nuisance_raw_cosine"]
        and summary["random_gene_optimism_gap"]
        >= thresholds["min_random_gene_optimism_gap"]
        and summary["null_module_holdout_improvement_mean"]
        <= thresholds["max_null_module_improvement"]
        and summary["alternative_module_holdout_improvement_mean"]
        >= thresholds["min_alternative_module_improvement"]
        and summary["sign_selection_cosine_inflation_mean"]
        >= thresholds["min_sign_selection_inflation"]
    )
    return {"status": "PASS" if passed else "FAIL", **summary}


def run(config_path: Path) -> dict[str, Any]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "1.1.0":
        raise ValueError("unsupported harness config schema")
    scenario_seeds = config["scenario_seeds"]
    generators = {
        name: np.random.default_rng(
            np.random.SeedSequence([config["seed"], scenario_seeds[name]])
        )
        for name in scenario_seeds
    }
    scenarios = {
        "active_set_oracle": _oracle_scenario(
            config, generators["active_set_oracle"]
        ),
        "degenerate_cones": _degeneracy_scenario(
            config, generators["degenerate_cones"]
        ),
        "label_and_provenance_faults": _label_scenario(),
        "grouped_gene_holdout": _grouped_split_scenario(
            config, generators["grouped_gene_holdout"]
        ),
        "maxT_null_calibration": _null_scenario(
            config, generators["maxT_null_calibration"]
        ),
        "structured_specificity_calibration": _structured_specificity_scenario(
            config, generators["structured_specificity_calibration"]
        ),
    }
    status = "PASS" if all(item["status"] == "PASS" for item in scenarios.values()) else "FAIL"
    return {
        "schema_version": "1.1.0",
        "generated_on": "2026-07-17",
        "config_sha256": _sha256(config_path),
        "status": status,
        "scope": "deterministic synthetic software/statistical contract; not biological validation",
        "scenarios": scenarios,
    }


def _assert_matches(actual: Any, expected: Any, path: str = "report") -> None:
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise AssertionError(f"{path} keys differ")
        for key in expected:
            _assert_matches(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, float):
        if not np.isclose(actual, expected, rtol=1e-9, atol=1e-12):
            raise AssertionError(f"{path}: {actual} != {expected}")
    elif actual != expected:
        raise AssertionError(f"{path}: {actual} != {expected}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = run(args.config)
    if args.write:
        args.write.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {args.write}")
    elif args.check:
        expected = json.loads(args.check.read_text(encoding="utf-8"))
        _assert_matches(report, expected)
        print(f"validation harness: {report['status']} (frozen report matches)")
    else:
        print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
