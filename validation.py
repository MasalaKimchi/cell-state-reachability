"""Reusable validation contracts for cell-state reachability studies.

This module is deliberately data-agnostic.  It provides strict labelled alignment,
an independent small-problem NNLS oracle, grouped gene splits, and maxT-adjusted
empirical p-values.  None of these functions assigns biological meaning to a fit.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re
from typing import Iterable

import numpy as np
from scipy.stats import beta

from reachability import InputError


_SHA256 = re.compile(r"[0-9a-f]{64}")


def _stable_norm(values: np.ndarray) -> float:
    scale = float(np.max(np.abs(values), initial=0.0))
    return 0.0 if scale == 0.0 else scale * float(np.linalg.norm(values / scale))


@dataclass(frozen=True)
class Provenance:
    """Minimum metadata required to interpret a numerical vector or matrix."""

    dataset_id: str
    source_sha256: str
    gene_namespace: str
    units: str
    modality: str
    context: str
    timepoint: str
    orientation: str


@dataclass(frozen=True)
class LabeledEffects:
    values: np.ndarray
    perturbations: tuple[str, ...]
    genes: tuple[str, ...]
    provenance: Provenance


@dataclass(frozen=True)
class LabeledTarget:
    values: np.ndarray
    genes: tuple[str, ...]
    provenance: Provenance


@dataclass(frozen=True)
class AlignedProblem:
    effects: np.ndarray
    target: np.ndarray
    perturbations: tuple[str, ...]
    genes: tuple[str, ...]
    effects_provenance: Provenance
    target_provenance: Provenance
    effects_gene_fraction: float
    target_gene_fraction: float


@dataclass(frozen=True)
class LabeledNullColumn:
    """One null-score column with explicit synchronized resample identities."""

    values: np.ndarray
    resample_ids: tuple[str, ...]


def _validate_provenance(provenance: Provenance) -> None:
    for name in (
        "dataset_id",
        "gene_namespace",
        "units",
        "modality",
        "context",
        "timepoint",
        "orientation",
    ):
        value = getattr(provenance, name)
        if not isinstance(value, str) or not value.strip():
            raise InputError(f"provenance {name} must be a non-empty string")
    if not isinstance(provenance.source_sha256, str) or not _SHA256.fullmatch(
        provenance.source_sha256
    ):
        raise InputError("provenance source_sha256 must be 64 lowercase hex characters")


def _validate_labels(labels: tuple[str, ...], name: str) -> None:
    if not labels:
        raise InputError(f"{name} must be non-empty")
    if any(not isinstance(label, str) or not label.strip() for label in labels):
        raise InputError(f"{name} must contain non-empty strings")
    if len(set(labels)) != len(labels):
        raise InputError(f"{name} must be unique")


def align_labeled_problem(
    effects: LabeledEffects,
    target: LabeledTarget,
    *,
    allow_intersection: bool = False,
    min_shared_genes: int = 2,
    min_shared_fraction: float = 1.0,
) -> AlignedProblem:
    """Validate and explicitly align labelled effects and target gene axes.

    Reordering is safe because labels, not positions, control alignment. Missing genes
    fail closed unless ``allow_intersection`` is explicit and both coverage gates pass.
    """

    _validate_provenance(effects.provenance)
    _validate_provenance(target.provenance)
    _validate_labels(effects.genes, "effect genes")
    _validate_labels(target.genes, "target genes")
    _validate_labels(effects.perturbations, "perturbations")

    values = np.asarray(effects.values, dtype=float)
    direction = np.asarray(target.values, dtype=float)
    if values.ndim != 2 or values.shape != (
        len(effects.perturbations),
        len(effects.genes),
    ):
        raise InputError("effect values must have perturbations-by-genes shape")
    if direction.ndim != 1 or direction.shape != (len(target.genes),):
        raise InputError("target values must match the labelled target gene axis")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(direction)):
        raise InputError("labelled values must contain only finite numbers")
    if effects.provenance.gene_namespace != target.provenance.gene_namespace:
        raise InputError("gene namespaces differ")
    if effects.provenance.units != target.provenance.units:
        raise InputError("effect and target units differ")
    if effects.provenance.orientation != "perturbations_by_genes":
        raise InputError("effect orientation must be perturbations_by_genes")
    if target.provenance.orientation != "genes":
        raise InputError("target orientation must be genes")
    if not isinstance(min_shared_genes, int) or isinstance(min_shared_genes, bool) or min_shared_genes < 1:
        raise InputError("min_shared_genes must be a positive integer")
    if not np.isfinite(min_shared_fraction) or not 0 < min_shared_fraction <= 1:
        raise InputError("min_shared_fraction must be in (0, 1]")

    effect_lookup = {gene: i for i, gene in enumerate(effects.genes)}
    target_lookup = {gene: i for i, gene in enumerate(target.genes)}
    common = tuple(gene for gene in effects.genes if gene in target_lookup)
    effect_fraction = len(common) / len(effects.genes)
    target_fraction = len(common) / len(target.genes)
    exact_set = len(common) == len(effects.genes) == len(target.genes)
    if not allow_intersection and not exact_set:
        raise InputError("gene sets differ; intersection must be explicitly enabled")
    if (
        len(common) < min_shared_genes
        or effect_fraction < min_shared_fraction
        or target_fraction < min_shared_fraction
    ):
        raise InputError("shared gene coverage is below the declared gate")

    effect_indices = np.fromiter((effect_lookup[g] for g in common), dtype=np.intp)
    target_indices = np.fromiter((target_lookup[g] for g in common), dtype=np.intp)
    return AlignedProblem(
        effects=values[:, effect_indices],
        target=direction[target_indices],
        perturbations=effects.perturbations,
        genes=common,
        effects_provenance=effects.provenance,
        target_provenance=target.provenance,
        effects_gene_fraction=effect_fraction,
        target_gene_fraction=target_fraction,
    )


def active_set_oracle(
    effects: np.ndarray,
    target: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
    nonnegative_tolerance: float = 1e-10,
    max_atoms: int = 16,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Exhaustively solve a small NNLS problem without calling an NNLS solver.

    Every candidate support is fit by ordinary least squares after independent column
    normalization. The returned fitted point and relative objective certify production
    solutions; coefficients may remain non-identifiable.
    """

    matrix = np.asarray(effects, dtype=float)
    direction = np.asarray(target, dtype=float)
    if matrix.ndim != 2 or direction.ndim != 1 or matrix.shape[1] != direction.size:
        raise InputError("oracle inputs must be aligned perturbations-by-genes and target")
    if matrix.shape[0] == 0 or direction.size == 0 or not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(direction)):
        raise InputError("oracle inputs must be non-empty and finite")
    if matrix.shape[0] > max_atoms:
        raise InputError("oracle problem exceeds max_atoms")
    if not np.isfinite(nonnegative_tolerance) or nonnegative_tolerance < 0:
        raise InputError("nonnegative_tolerance must be finite and non-negative")
    weights = np.ones(direction.size) if gene_weights is None else np.asarray(gene_weights, dtype=float)
    if weights.shape != direction.shape or not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise InputError("oracle gene_weights must be finite and positive")
    weights = weights / np.max(weights)
    if np.any(weights == 0):
        raise InputError("oracle gene-weight dynamic range is not representable")

    design = np.sqrt(weights)[:, None] * matrix.T
    response = np.sqrt(weights) * direction
    response_norm = _stable_norm(response)
    if response_norm == 0:
        raise InputError("oracle target is zero or below floating-point resolution")
    response_unit = response / response_norm
    atom_norms = np.asarray(
        [_stable_norm(design[:, index]) for index in range(design.shape[1])]
    )
    active_indices = np.flatnonzero(atom_norms > 0)
    unit_design = design[:, active_indices] / atom_norms[active_indices]
    best_coefficients = np.zeros(matrix.shape[0])
    best_fitted = np.zeros_like(direction)
    best_objective = 0.5
    for size in range(1, active_indices.size + 1):
        for support in combinations(range(active_indices.size), size):
            subset = unit_design[:, support]
            coefficients, _, _, _ = np.linalg.lstsq(subset, response_unit, rcond=None)
            if np.any(coefficients < -nonnegative_tolerance):
                continue
            coefficients = np.maximum(coefficients, 0.0)
            candidate = np.zeros(matrix.shape[0])
            source_indices = active_indices[list(support)]
            candidate[source_indices] = (
                coefficients * response_norm / atom_norms[source_indices]
            )
            fitted = candidate @ matrix
            residual = np.sqrt(weights) * (direction - fitted) / response_norm
            objective = 0.5 * float(np.dot(residual, residual))
            if objective < best_objective:
                best_coefficients = candidate
                best_fitted = fitted
                best_objective = objective
    return best_coefficients, best_fitted, best_objective


def grouped_gene_splits(
    groups: Iterable[str], *, n_splits: int, seed: int
) -> tuple[tuple[np.ndarray, np.ndarray], ...]:
    """Create deterministic half-group fit/score splits with no group leakage."""

    group_array = np.asarray(tuple(groups), dtype=object)
    if group_array.ndim != 1 or group_array.size < 2:
        raise InputError("groups must be a one-dimensional vector")
    if any(not isinstance(group, str) or not group for group in group_array):
        raise InputError("group labels must be non-empty strings")
    unique = np.unique(group_array)
    if unique.size < 2:
        raise InputError("at least two groups are required")
    if not isinstance(n_splits, int) or isinstance(n_splits, bool) or n_splits < 1:
        raise InputError("n_splits must be a positive integer")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise InputError("seed must be an integer")

    rng = np.random.default_rng(seed)
    splits = []
    for _ in range(n_splits):
        shuffled = rng.permutation(unique)
        fit_groups = set(shuffled[: unique.size // 2])
        fit = np.flatnonzero(np.fromiter((g in fit_groups for g in group_array), dtype=bool))
        score = np.flatnonzero(np.fromiter((g not in fit_groups for g in group_array), dtype=bool))
        splits.append((fit, score))
    return tuple(splits)


def max_t_empirical_p(observed: Iterable[float], null_scores: np.ndarray) -> np.ndarray:
    """Single-step maxT familywise p-values with plus-one finite-sample correction."""

    values = np.asarray(tuple(observed), dtype=float)
    null = np.asarray(null_scores, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise InputError("observed must be a non-empty vector")
    if null.ndim != 2 or null.shape[1] != values.size or null.shape[0] == 0:
        raise InputError("null_scores must be resamples-by-hypotheses")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(null)):
        raise InputError("maxT inputs must be finite")
    maxima = np.max(null, axis=1)
    return np.asarray(
        [(1 + np.count_nonzero(maxima >= value)) / (null.shape[0] + 1) for value in values]
    )


def stack_synchronized_nulls(
    columns: Iterable[tuple[str, LabeledNullColumn]],
) -> tuple[tuple[str, ...], tuple[str, ...], np.ndarray]:
    """Stack hypothesis-wise null scores only when resample rows align exactly."""

    materialized = tuple(columns)
    if not materialized:
        raise InputError("at least one labelled null column is required")
    hypotheses = tuple(name for name, _ in materialized)
    if any(not isinstance(name, str) or not name.strip() for name in hypotheses):
        raise InputError("null hypothesis labels must be non-empty strings")
    if len(set(hypotheses)) != len(hypotheses):
        raise InputError("null hypothesis labels must be unique")

    reference_ids: tuple[str, ...] | None = None
    values = []
    for _, column in materialized:
        vector = np.asarray(column.values, dtype=float)
        identifiers = tuple(column.resample_ids)
        if vector.ndim != 1 or vector.size == 0 or not np.all(np.isfinite(vector)):
            raise InputError("each labelled null column must be a non-empty finite vector")
        if len(identifiers) != vector.size:
            raise InputError("null resample IDs must match the score vector")
        if any(not isinstance(item, str) or not item.strip() for item in identifiers):
            raise InputError("null resample IDs must be non-empty strings")
        if len(set(identifiers)) != len(identifiers):
            raise InputError("null resample IDs must be unique")
        if reference_ids is None:
            reference_ids = identifiers
        elif identifiers != reference_ids:
            raise InputError("null resample rows are not synchronized across hypotheses")
        values.append(vector)
    assert reference_ids is not None
    return hypotheses, reference_ids, np.column_stack(values)


def holm_adjusted_p(p_values: Iterable[float]) -> np.ndarray:
    """Holm familywise adjusted p-values for valid marginal p-values."""

    values = np.asarray(tuple(p_values), dtype=float)
    if (
        values.ndim != 1
        or values.size == 0
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
        or np.any(values > 1)
    ):
        raise InputError("p_values must be a non-empty finite vector in [0, 1]")
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted_ranked = np.maximum.accumulate(
        (values.size - np.arange(values.size)) * ranked
    )
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def binomial_upper_confidence_bound(
    successes: int, trials: int, *, confidence: float = 0.95
) -> float:
    """One-sided exact Clopper-Pearson upper confidence bound."""

    if (
        not isinstance(successes, int)
        or isinstance(successes, bool)
        or not isinstance(trials, int)
        or isinstance(trials, bool)
        or trials < 1
        or successes < 0
        or successes > trials
    ):
        raise InputError("successes and trials must be integers with 0 <= successes <= trials")
    if not np.isfinite(confidence) or not 0 < confidence < 1:
        raise InputError("confidence must be in (0, 1)")
    if successes == trials:
        return 1.0
    return float(beta.ppf(confidence, successes + 1, trials - successes))
