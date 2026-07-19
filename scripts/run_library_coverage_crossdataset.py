#!/usr/bin/env python3
"""Run the frozen retrospective effect-dictionary coverage benchmark.

The runner binds three local derived caches by byte length and SHA-256 before
loading any array or computing any claim-bearing metric. Zhu and Replogle are
retrospective compression/self-reconstruction challenges: library, target
catalog, and supplied candidate atoms are rows of the same measured effect
dictionary. Norman is a retrospective single-to-double additivity/coverage
challenge in canonical K562 cells; it is not held-out biological validation.

Strict numerical cone inclusion and soft cosine-threshold coverage are reported
separately. Acquisition candidates are supplied, already-measured effect atoms.
The expensive comparator is realized *mean-cosine gain* after adding each atom;
soft threshold-coverage gain is secondary and explicitly named.

Exact ``--check`` requires the three registered local NPZ byte artifacts. The companion
``build_library_coverage_caches.py`` reconstructs them deterministically from hash-matching
raw inputs. Public URLs remain mutable availability routes, so a hash mismatch fails closed
and durable retrieval of historical source bytes is not implied.
"""

from __future__ import annotations

import argparse
from datetime import date
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from library_coverage import coverage_report, rank_acquisitions  # noqa: E402
from effect_dictionary import load_effect_dictionary  # noqa: E402
from reachability import InputError  # noqa: E402

try:
    from scipy.stats import spearmanr
except ImportError as exc:  # pragma: no cover
    raise SystemExit("scipy is required for the cross-dataset benchmark") from exc


DEFAULT_CONFIG = ROOT / "configs" / "library_coverage_crossdataset.json"
DEFAULT_REPORT = ROOT / "results" / "library_coverage_crossdataset.json"

SCHEMA_VERSION = "3.0.0"
BENCHMARK_ID = "library_coverage_split_first_retrospective_v3"
STANDARD_FEATURE_SELECTION = (
    "stable top-400 variance coordinates from the 50 current-library rows after "
    "the retrospective row split; catalog and candidate rows do not participate"
)
NORMAN_FEATURE_SELECTION = (
    "full additivity description: stable top-400 variance coordinates from all 152 "
    "position-specific single-perturbation rows only; acquisition and split "
    "sensitivity: group by 105 canonical genes, freeze one measured cassette-position "
    "representative per gene, then select a separate stable top-400 set from each "
    "current 40-gene library only; doubles and candidates never select those features"
)
DOMINANCE_TOLERANCE = 1e-8


def _progress(message: str) -> None:
    print(f"[library-coverage-crossdataset] {message}", file=sys.stderr, flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _cache_path(spec: dict[str, Any], cache_dir: Path | None) -> Path:
    configured = Path(spec["path"])
    return ROOT / configured if cache_dir is None else cache_dir / configured.name


def verify_input(
    name: str,
    spec: dict[str, Any],
    *,
    cache_dir: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Verify one cache identity before NumPy is allowed to load it."""
    path = _cache_path(spec, cache_dir)
    if not path.is_file():
        raise FileNotFoundError(f"{name}: missing frozen cache {path}")
    expected_bytes = spec.get("bytes")
    expected_hash = spec.get("sha256")
    if not isinstance(expected_bytes, int) or expected_bytes <= 0:
        raise InputError(f"{name}: config bytes must be a positive integer")
    if (
        not isinstance(expected_hash, str)
        or len(expected_hash) != 64
        or any(ch not in "0123456789abcdef" for ch in expected_hash)
    ):
        raise InputError(f"{name}: config SHA-256 must be 64 lowercase hex characters")
    actual_bytes = path.stat().st_size
    if actual_bytes != expected_bytes:
        raise InputError(
            f"{name}: byte length differs: expected {expected_bytes}, found {actual_bytes}"
        )
    actual_hash = sha256_file(path)
    if actual_hash != expected_hash:
        raise InputError(f"{name}: SHA-256 differs from the frozen cache identity")
    return path, {
        "path": spec["path"],
        "bytes": actual_bytes,
        "sha256_expected": expected_hash,
        "sha256_actual": actual_hash,
        "hash_verified": True,
    }


def _validate_config(config: dict[str, Any]) -> None:
    required = {
        "schema_version",
        "generated_on",
        "benchmark",
        "claim_ceiling",
        "provenance",
        "cache_build",
        "inputs",
        "analysis",
        "datasets",
    }
    if set(config) != required:
        raise InputError(
            "config top-level keys differ: "
            f"missing={sorted(required - set(config))}, "
            f"extra={sorted(set(config) - required)}"
        )
    if config["schema_version"] != SCHEMA_VERSION:
        raise InputError(f"schema_version must be {SCHEMA_VERSION}")
    if config["benchmark"] != BENCHMARK_ID:
        raise InputError(f"benchmark must be {BENCHMARK_ID}")
    try:
        generated_on = date.fromisoformat(config["generated_on"])
    except (TypeError, ValueError) as exc:
        raise InputError("generated_on must be an ISO YYYY-MM-DD date") from exc
    if generated_on.isoformat() != config["generated_on"]:
        raise InputError("generated_on must use canonical YYYY-MM-DD formatting")
    provenance = config["provenance"]
    if not isinstance(provenance, dict) or set(provenance) != {
        "candidate_status",
        "perturbation_modalities",
        "cell_systems",
        "zhu",
        "norman",
        "replogle",
    }:
        raise InputError("provenance must contain the frozen v3 provenance sections")
    cache_build = config["cache_build"]
    if not isinstance(cache_build, dict) or set(cache_build) != {
        "format",
        "writer_contract",
        "zhu",
        "norman",
        "replogle",
    }:
        raise InputError("cache_build must contain the frozen v3 build sections")
    if cache_build["format"] != "portable_effect_dictionary_v1":
        raise InputError("cache_build.format differs from the supported portable format")
    if set(config["inputs"]) != {"zhu", "norman", "replogle"}:
        raise InputError("config inputs must be exactly zhu, norman, and replogle")
    if set(config["datasets"]) != {
        "zhu_crispri_tcell",
        "norman_k562_crispra",
        "replogle_k562_essential_crispri",
    }:
        raise InputError("config dataset IDs differ from the frozen benchmark")
    threshold = config["analysis"].get("soft_cosine_threshold")
    if not isinstance(threshold, (int, float)) or not np.isfinite(threshold):
        raise InputError("analysis.soft_cosine_threshold must be finite")
    if not -1.0 <= float(threshold) <= 1.0:
        raise InputError("analysis.soft_cosine_threshold must lie in [-1, 1]")
    sensitivity = config["analysis"].get("split_sensitivity")
    if not isinstance(sensitivity, dict):
        raise InputError("analysis.split_sensitivity must be configured")
    seeds = sensitivity.get("seeds")
    thresholds = sensitivity.get("soft_cosine_thresholds")
    if (
        not isinstance(seeds, list)
        or len(seeds) != 12
        or any(isinstance(seed, bool) or not isinstance(seed, int) for seed in seeds)
        or len(set(seeds)) != len(seeds)
    ):
        raise InputError("split-sensitivity seeds must be twelve distinct integers")
    if (
        not isinstance(thresholds, list)
        or not thresholds
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not np.isfinite(value)
            or not -1.0 <= float(value) <= 1.0
            for value in thresholds
        )
        or thresholds != sorted(set(thresholds))
    ):
        raise InputError("split-sensitivity thresholds must be sorted unique values in [-1, 1]")

    inputs = config["inputs"]
    expected_input_shapes = {
        "zhu": ([6871, 8950], [6871]),
        "norman": ([283, 5045], [283]),
        "replogle": ([1087, 5000], [1087]),
    }
    for name, (effect_shape, label_shape) in expected_input_shapes.items():
        spec = inputs[name]
        if set(spec) != {
            "path",
            "bytes",
            "sha256",
            "effect_array",
            "label_array",
            "required_arrays",
        }:
            raise InputError(f"{name}: input contract keys differ from portable v1")
        if spec["effect_array"] != "E" or spec["label_array"] != "perts":
            raise InputError(f"{name}: effect/label array names must be E/perts")
        arrays = spec["required_arrays"]
        if arrays["E"]["shape"] != effect_shape or arrays["perts"]["shape"] != label_shape:
            raise InputError(f"{name}: effect/label shapes differ from frozen v3")

    datasets = config["datasets"]
    for name in ("zhu_crispri_tcell", "replogle_k562_essential_crispri"):
        design = datasets[name]
        if design["n_features"] != 400 or design["feature_selection"] != STANDARD_FEATURE_SELECTION:
            raise InputError(f"{name}: split-first feature-selection contract differs")
        if min(design["library_size"], design["catalog_size"], design["candidate_size"]) <= 0:
            raise InputError(f"{name}: reference split sizes must be positive")
        if design["sensitivity_catalog_size"] <= 0 or design["sensitivity_candidate_reserve"] <= 0:
            raise InputError(f"{name}: sensitivity split sizes must be positive")
    norman = datasets["norman_k562_crispra"]
    if norman["n_features"] != 400 or norman["feature_selection"] != NORMAN_FEATURE_SELECTION:
        raise InputError("Norman feature-selection contract differs from frozen v3")
    if norman["acquisition_library_size"] != 40:
        raise InputError("Norman reference and sensitivity library size must be 40")
    if norman["acquisition_candidate_size"] <= 0:
        raise InputError("Norman acquisition candidate size must be positive")
    if not isinstance(norman["sensitivity_catalog_seed"], int):
        raise InputError("Norman sensitivity catalog seed must be an integer")
    if not 0 < norman["sensitivity_catalog_doubles"] <= norman["expected"]["double_effect_rows"]:
        raise InputError("Norman sensitivity catalog size is invalid")


def load_matrix(path: Path, spec: dict[str, Any]) -> dict[str, Any]:
    """Load a hash-verified, pickle-free dictionary and validate every array."""
    effect_key = spec["effect_array"]
    label_key = spec["label_array"]
    expected = spec["required_arrays"]
    try:
        dictionary = load_effect_dictionary(path)
    except (OSError, ValueError) as exc:
        raise InputError(f"{spec['path']}: invalid or unsafe NPZ cache") from exc
    if set(dictionary) != set(expected):
        raise InputError(
            f"{spec['path']}: arrays differ: "
            f"{sorted(dictionary)} != {sorted(expected)}"
        )
    for key, array in dictionary.items():
        array_spec = expected[key]
        if list(array.shape) != array_spec["shape"]:
            raise InputError(
                f"{spec['path']}:{key} shape differs: "
                f"{list(array.shape)} != {array_spec['shape']}"
            )
        if str(array.dtype) != array_spec["dtype"]:
            raise InputError(
                f"{spec['path']}:{key} dtype differs: "
                f"{array.dtype} != {array_spec['dtype']}"
            )
    effects = dictionary[effect_key]
    labels_raw = dictionary[label_key]
    if effects.ndim != 2 or not np.all(np.isfinite(effects)):
        raise InputError(f"{spec['path']} effect matrix must be finite and two-dimensional")
    labels = np.asarray([str(value) for value in labels_raw], dtype=str)
    if labels.ndim != 1 or len(labels) != effects.shape[0]:
        raise InputError(f"{spec['path']} labels do not align with effect rows")
    if np.any(labels == "") or len(set(labels.tolist())) != len(labels):
        raise InputError(f"{spec['path']} labels must be non-empty and unique")
    return {"E": effects, "labels": labels}


def _top_variance_features(effects: np.ndarray, count: int) -> np.ndarray:
    if effects.ndim != 2 or count <= 0 or count > effects.shape[1]:
        raise InputError("feature count must be positive and no larger than the matrix")
    variances = np.var(effects, axis=0, dtype=np.float64)
    if not np.all(np.isfinite(variances)):
        raise InputError("feature variances are not finite")
    return np.argsort(-variances, kind="stable")[:count]


def _library_only_features(
    effects: np.ndarray, library_indices: np.ndarray, count: int
) -> np.ndarray:
    """Rank coordinates using only already-frozen current-library rows."""

    indices = np.asarray(library_indices, dtype=int)
    if indices.ndim != 1 or indices.size == 0:
        raise InputError("feature selection requires non-empty 1D library indices")
    if np.any(indices < 0) or np.any(indices >= effects.shape[0]):
        raise InputError("library indices fall outside the effect matrix")
    return _top_variance_features(np.asarray(effects)[indices], count)


def _index_sha256(indices: np.ndarray) -> str:
    canonical = np.asarray(indices, dtype="<i8")
    return hashlib.sha256(canonical.tobytes(order="C")).hexdigest()


def _row_cosines(targets: np.ndarray, references: np.ndarray) -> np.ndarray:
    target_norms = np.linalg.norm(targets, axis=1)
    reference_norms = np.linalg.norm(references, axis=1)
    denom = target_norms[:, None] * reference_norms[None, :]
    similarities = np.divide(
        targets @ references.T,
        denom,
        out=np.zeros((targets.shape[0], references.shape[0]), dtype=float),
        where=denom > 0,
    )
    return similarities


def _nonnegative_ray_cosines(targets: np.ndarray, rays: np.ndarray) -> np.ndarray:
    """Best cosine after non-negative scaling of each supplied ray."""
    return np.maximum(_row_cosines(targets, rays).max(axis=1), 0.0)


def _alignment_summary(cosines: np.ndarray, threshold: float) -> dict[str, Any]:
    cosines = np.asarray(cosines, dtype=float)
    soft = cosines >= threshold
    return {
        "n_targets": int(cosines.size),
        "mean_cosine": float(cosines.mean()),
        "soft_cosine_threshold": threshold,
        "soft_coverage_count": int(soft.sum()),
        "soft_coverage_fraction": float(soft.mean()),
    }


def _coverage_summary(report: Any) -> dict[str, Any]:
    if report.soft_covered is None or report.soft_coverage_fraction is None:
        raise RuntimeError("frozen benchmark requires explicit soft coverage")
    return {
        "n_targets": int(len(report.cosines)),
        "strict_inside_cone_count": int(report.strict_inside_cone.sum()),
        "strict_inside_cone_fraction": report.strict_inside_cone_fraction,
        "soft_cosine_threshold": report.soft_cosine_threshold,
        "soft_coverage_count": int(report.soft_covered.sum()),
        "soft_coverage_fraction": report.soft_coverage_fraction,
        "mean_cosine": report.mean_cosine,
        "mean_residual_fraction": report.mean_residual_fraction,
    }


def _signed_span_cosines(targets: np.ndarray, library: np.ndarray) -> np.ndarray:
    """Cosines for the orthogonal projection onto the unrestricted row span.

    This is a capacity ceiling, not a biologically admissible comparator:
    coefficients may have either sign.  The rank tolerance is fixed by the
    ordinary LAPACK convention and never selected with target outcomes.
    """

    matrix = np.asarray(library, dtype=float)
    target_matrix = np.asarray(targets, dtype=float)
    if matrix.ndim != 2 or target_matrix.ndim != 2:
        raise InputError("signed-span inputs must be two-dimensional")
    if matrix.shape[1] != target_matrix.shape[1]:
        raise InputError("signed-span inputs must share a feature axis")
    _, singular_values, right_vectors = np.linalg.svd(matrix, full_matrices=False)
    if singular_values.size == 0 or singular_values[0] == 0:
        return np.zeros(target_matrix.shape[0], dtype=float)
    tolerance = (
        max(matrix.shape) * np.finfo(float).eps * float(singular_values[0])
    )
    rank = int(np.count_nonzero(singular_values > tolerance))
    basis = right_vectors[:rank]
    predictions = (target_matrix @ basis.T) @ basis
    return _paired_cosines(predictions, target_matrix)


def _distribution_summary(values: list[float]) -> dict[str, Any]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise InputError("split-sensitivity summaries require finite values")
    return {
        "n_splits": int(array.size),
        "median": float(np.median(array)),
        "q25": float(np.quantile(array, 0.25)),
        "q75": float(np.quantile(array, 0.75)),
        "p05": float(np.quantile(array, 0.05)),
        "p95": float(np.quantile(array, 0.95)),
        "minimum": float(array.min()),
        "maximum": float(array.max()),
        "positive_split_count": int(np.count_nonzero(array > 0)),
    }


def _soft_threshold_fractions(
    cosines: np.ndarray, thresholds: list[float]
) -> dict[str, float]:
    values = np.asarray(cosines, dtype=float)
    return {
        f"cosine_at_least_{threshold:g}": float(np.mean(values >= threshold))
        for threshold in thresholds
    }


def _assert_capacity_dominance(
    cone: np.ndarray,
    best_single: np.ndarray,
    common: np.ndarray,
    signed: np.ndarray,
) -> None:
    """Fail closed if nested model classes violate their mathematical ordering."""

    cone_values = np.asarray(cone, dtype=float)
    best_values = np.asarray(best_single, dtype=float)
    common_values = np.asarray(common, dtype=float)
    signed_values = np.asarray(signed, dtype=float)
    if not (
        cone_values.shape
        == best_values.shape
        == common_values.shape
        == signed_values.shape
    ):
        raise InputError("capacity-comparator cosine axes do not align")
    if np.any(cone_values + DOMINANCE_TOLERANCE < best_values):
        raise InputError("cone fit is worse than its nested best-single comparator")
    if np.any(cone_values + DOMINANCE_TOLERANCE < common_values):
        raise InputError("cone fit is worse than its nested common-response comparator")
    if np.any(signed_values + DOMINANCE_TOLERANCE < cone_values):
        raise InputError("signed-span fit is worse than its nested nonnegative cone")


def standard_split_sensitivity(
    payload: dict[str, Any],
    design: dict[str, Any],
    *,
    seeds: list[int],
    thresholds: list[float],
) -> dict[str, Any]:
    """Coverage-only split sensitivity without repeated acquisition fitting."""

    effects = np.asarray(payload["E"])
    source_indices = np.arange(effects.shape[0], dtype=int)
    feature_count = int(design["n_features"])
    catalog_size = int(design.get("sensitivity_catalog_size", design["catalog_size"]))
    candidate_size = int(
        design.get("sensitivity_candidate_reserve", design["candidate_size"])
    )
    rows = []
    for seed in seeds:
        library_idx, catalog_idx, candidate_idx, unused_idx = _split_design(
            source_indices,
            int(seed),
            int(design["library_size"]),
            catalog_size,
            candidate_size,
        )
        feature_indices = _library_only_features(
            effects, library_idx, feature_count
        )
        selected = np.asarray(effects[:, feature_indices], dtype=float)
        library = selected[library_idx]
        catalog = selected[catalog_idx]
        coverage = coverage_report(library, catalog)
        best_single = _nonnegative_ray_cosines(catalog, library)
        common = _nonnegative_ray_cosines(
            catalog, library.mean(axis=0, keepdims=True)
        )
        signed = _signed_span_cosines(catalog, library)
        _assert_capacity_dominance(coverage.cosines, best_single, common, signed)
        rows.append(
            {
                "seed": int(seed),
                "library_source_row_index_sha256": _index_sha256(library_idx),
                "catalog_source_row_index_sha256": _index_sha256(catalog_idx),
                "candidate_source_row_index_sha256": _index_sha256(candidate_idx),
                "unused_source_row_index_sha256": _index_sha256(unused_idx),
                "feature_index_sha256": _index_sha256(feature_indices),
                "cone_mean_cosine": coverage.mean_cosine,
                "cone_strict_inside_fraction": coverage.strict_inside_cone_fraction,
                "cone_soft_coverage_by_threshold": _soft_threshold_fractions(
                    coverage.cosines, thresholds
                ),
                "best_single_mean_cosine": float(best_single.mean()),
                "common_response_mean_cosine": float(common.mean()),
                "signed_span_mean_cosine": float(signed.mean()),
                "cone_minus_best_single_mean_cosine": float(
                    coverage.mean_cosine - best_single.mean()
                ),
                "cone_minus_common_response_mean_cosine": float(
                    coverage.mean_cosine - common.mean()
                ),
                "signed_span_minus_cone_mean_cosine": float(
                    signed.mean() - coverage.mean_cosine
                ),
            }
        )

    metric_names = (
        "cone_mean_cosine",
        "cone_strict_inside_fraction",
        "best_single_mean_cosine",
        "common_response_mean_cosine",
        "signed_span_mean_cosine",
        "cone_minus_best_single_mean_cosine",
        "cone_minus_common_response_mean_cosine",
        "signed_span_minus_cone_mean_cosine",
    )
    aggregate = {
        name: _distribution_summary([float(row[name]) for row in rows])
        for name in metric_names
    }
    for threshold in thresholds:
        key = f"cosine_at_least_{threshold:g}"
        aggregate[f"cone_soft_coverage_{key}"] = _distribution_summary(
            [float(row["cone_soft_coverage_by_threshold"][key]) for row in rows]
        )
    return {
        "protocol": {
            "purpose": "algorithmic row-split sensitivity; not biological sampling uncertainty",
            "rng": "numpy.default_rng(seed)",
            "seeds": [int(seed) for seed in seeds],
            "library_atoms": int(design["library_size"]),
            "catalog_effect_rows": catalog_size,
            "reserved_candidate_effect_atoms": candidate_size,
            "feature_count": feature_count,
            "feature_selection_rows": "current_library_only",
            "soft_cosine_thresholds": thresholds,
            "acquisition_recomputed": False,
        },
        "splits": rows,
        "aggregate": aggregate,
    }


def _canonical_candidate_label(label: str) -> str:
    constituents = sorted(part for part in str(label).split("+") if part != "ctrl")
    canonical = "+".join(constituents)
    if not canonical:
        raise InputError("candidate effect atom has no non-control perturbation label")
    return canonical


def _candidate_identity(
    position: int,
    source_indices: np.ndarray,
    labels: np.ndarray,
) -> dict[str, Any]:
    label = str(labels[position])
    canonical = _canonical_candidate_label(label)
    return {
        "candidate_pool_position_zero_based": int(position),
        "source_effect_row_index_zero_based": int(source_indices[position]),
        "effect_atom_label": label,
        "canonical_perturbation_label": canonical,
    }


def _acquisition_summary(
    ranking: Any,
    source_indices: np.ndarray,
    labels: np.ndarray,
) -> dict[str, Any]:
    label_values = np.asarray(labels, dtype=str)
    canonical_values = np.asarray(
        [_canonical_candidate_label(label) for label in label_values], dtype=str
    )
    if len(set(label_values.tolist())) != len(label_values):
        raise InputError("acquisition candidate effect-atom labels must be unique")
    if len(set(canonical_values.tolist())) != len(canonical_values):
        raise InputError(
            "acquisition candidates must contain one effect atom per canonical label"
        )
    realized = ranking.realized_mean_cosine_gain
    realized_order = ranking.realized_mean_cosine_order
    soft_gain = ranking.realized_soft_coverage_fraction_gain
    soft_order = ranking.realized_soft_coverage_fraction_order
    strict_gain = ranking.realized_strict_inside_cone_fraction_gain
    if (
        realized is None
        or realized_order is None
        or soft_gain is None
        or soft_order is None
        or strict_gain is None
    ):
        raise RuntimeError("frozen benchmark requires every realized comparator")
    for name, values in (
        ("mean-cosine", realized),
        ("soft-coverage", soft_gain),
        ("strict-inclusion", strict_gain),
    ):
        if np.any(np.asarray(values, dtype=float) < -DOMINANCE_TOLERANCE):
            raise InputError(f"realized {name} gain is negative after adding an atom")
    rho = float(spearmanr(ranking.certificate_score, realized).statistic)
    if not np.isfinite(rho):
        raise InputError("certificate/realized mean-cosine Spearman is undefined")
    certificate_top = int(ranking.certificate_order[0])
    realized_top = int(realized_order[0])
    soft_top = int(soft_order[0])
    certificate_ties = int(
        np.sum(ranking.certificate_score == ranking.certificate_score[certificate_top])
    )
    realized_ties = int(np.sum(realized == realized[realized_top]))
    soft_ties = int(np.sum(soft_gain == soft_gain[soft_top]))
    certificate_identity = _candidate_identity(certificate_top, source_indices, labels)
    realized_identity = _candidate_identity(realized_top, source_indices, labels)
    soft_identity = _candidate_identity(soft_top, source_indices, labels)
    certificate_identity.update(
        {
            "certificate_score": float(ranking.certificate_score[certificate_top]),
            "realized_mean_cosine_gain": float(realized[certificate_top]),
            "realized_soft_coverage_fraction_gain": float(soft_gain[certificate_top]),
            "realized_strict_inside_cone_fraction_gain": float(
                strict_gain[certificate_top]
            ),
        }
    )
    realized_identity.update(
        {
            "certificate_score": float(ranking.certificate_score[realized_top]),
            "realized_mean_cosine_gain": float(realized[realized_top]),
            "realized_soft_coverage_fraction_gain": float(soft_gain[realized_top]),
            "realized_strict_inside_cone_fraction_gain": float(strict_gain[realized_top]),
        }
    )
    soft_identity.update(
        {
            "realized_soft_coverage_fraction_gain": float(soft_gain[soft_top]),
            "realized_mean_cosine_gain": float(realized[soft_top]),
        }
    )
    return {
        "candidate_status": "supplied_already_measured_effect_atoms",
        "certificate_definition": (
            "positive alignment with strict outside-cone dual separators, "
            "weighted by residual fraction"
        ),
        "primary_realized_objective": "mean_catalog_cosine_gain",
        "secondary_realized_objective": "soft_cosine_threshold_coverage_fraction_gain",
        "candidate_pool_effect_atom_labels_unique": True,
        "candidate_pool_canonical_labels_unique": True,
        "top1_agreement_unit": (
            "unique supplied effect atom; exactly one candidate exists per canonical label"
        ),
        "certificate_vs_realized_mean_cosine_gain_spearman": rho,
        "certificate_top1": certificate_identity,
        "certificate_top1_tie_count": certificate_ties,
        "realized_mean_cosine_top1": realized_identity,
        "realized_mean_cosine_top1_tie_count": realized_ties,
        "certificate_matches_realized_mean_cosine_top1": certificate_top == realized_top,
        "realized_soft_coverage_top1": soft_identity,
        "realized_soft_coverage_top1_tie_count": soft_ties,
    }


def _split_design(
    source_indices: np.ndarray,
    seed: int,
    library_size: int,
    catalog_size: int,
    candidate_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    required = library_size + catalog_size + candidate_size
    if min(library_size, catalog_size, candidate_size) <= 0 or required > len(source_indices):
        raise InputError("row split sizes are invalid for the supplied effect pool")
    permutation = np.random.default_rng(seed).permutation(len(source_indices))
    library = permutation[:library_size]
    catalog = permutation[library_size : library_size + catalog_size]
    candidates = permutation[library_size + catalog_size : required]
    unused = permutation[required:]
    return library, catalog, candidates, unused


def standard_split(
    payload: dict[str, Any],
    design: dict[str, Any],
    *,
    source_indices: np.ndarray | None = None,
) -> dict[str, Any]:
    """Retrospective dictionary compression with split-first feature selection.

    Row roles are frozen before any coordinate is ranked.  Only the current
    library rows may determine the feature set; catalog and candidate effects
    remain unseen until that preprocessing choice is fixed.
    """
    effects = payload["E"]
    labels = payload["labels"]
    source_indices = (
        np.arange(effects.shape[0], dtype=int)
        if source_indices is None
        else np.asarray(source_indices, dtype=int)
    )
    if len(source_indices) != effects.shape[0]:
        raise InputError("source row indices do not align with the effect pool")
    seed = int(design["seed"])
    library_idx, catalog_idx, candidate_idx, unused_idx = _split_design(
        source_indices,
        seed,
        int(design["library_size"]),
        int(design["catalog_size"]),
        int(design["candidate_size"]),
    )
    feature_count = int(design["n_features"])
    feature_indices = _library_only_features(effects, library_idx, feature_count)
    selected = np.asarray(effects[:, feature_indices], dtype=float)
    library = selected[library_idx]
    catalog = selected[catalog_idx]
    candidates = selected[candidate_idx]
    threshold = float(design["soft_cosine_threshold"])
    report = coverage_report(
        library, catalog, soft_cosine_threshold=threshold
    )
    ranking = rank_acquisitions(
        library,
        catalog,
        candidates,
        realized=True,
        soft_cosine_threshold=threshold,
    )
    common_ray = library.mean(axis=0, keepdims=True)
    best_single = _nonnegative_ray_cosines(catalog, library)
    common = _nonnegative_ray_cosines(catalog, common_ray)
    candidate_source_indices = source_indices[candidate_idx]
    candidate_labels = labels[candidate_idx]
    return {
        "effect_pool_rows": int(effects.shape[0]),
        "design": {
            "rng": "numpy.default_rng(seed)",
            "seed": seed,
            "library_atoms": int(len(library_idx)),
            "catalog_effect_rows": int(len(catalog_idx)),
            "supplied_candidate_effect_atoms": int(len(candidate_idx)),
            "unused_measured_effect_rows": int(len(unused_idx)),
            "feature_count": feature_count,
            "feature_selection": design["feature_selection"],
            "feature_selection_rows": "current_library_only",
            "feature_index_sha256": _index_sha256(feature_indices),
            "library_source_row_index_sha256": _index_sha256(
                source_indices[library_idx]
            ),
            "catalog_source_row_index_sha256": _index_sha256(
                source_indices[catalog_idx]
            ),
            "candidate_source_row_index_sha256": _index_sha256(
                candidate_source_indices
            ),
        },
        "coverage": _coverage_summary(report),
        "comparators": {
            "best_single_library_atom_nonnegative_ray": _alignment_summary(
                best_single, threshold
            ),
            "common_library_response_nonnegative_ray": _alignment_summary(
                common, threshold
            ),
        },
        "acquisition": _acquisition_summary(
            ranking, candidate_source_indices, candidate_labels
        ),
        "_cosines": report.cosines,
        "_certificate_scores": ranking.certificate_score,
        "_realized_mean_cosine_gains": ranking.realized_mean_cosine_gain,
    }


def build_zhu_from_slice(
    payload: dict[str, Any],
    design: dict[str, Any],
) -> dict[str, Any]:
    """Use every source-admitted Rest effect; do not outcome-filter the pool."""
    effects = payload["E"]
    result = standard_split(payload, design)
    result["source_effect_rows"] = int(effects.shape[0])
    result["preprocessing"] = {
        "effect_definition": design["effect_definition"],
        "eligible_effect_pool": design["eligible_effect_pool"],
        "outcome_based_effect_filtering": False,
    }
    return result


def _norman_membership(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    singles = np.zeros(len(labels), dtype=bool)
    doubles = np.zeros(len(labels), dtype=bool)
    for index, label in enumerate(labels):
        parts = label.split("+")
        if len(parts) != 2 or any(not part for part in parts):
            raise InputError(f"Norman perturbation label is not a two-part pair: {label!r}")
        control_count = parts.count("ctrl")
        if control_count == 1:
            singles[index] = True
        elif control_count == 0:
            doubles[index] = True
        else:
            raise InputError(f"Norman perturbation label has two controls: {label!r}")
    if np.any(singles & doubles) or not np.all(singles | doubles):
        raise InputError("Norman single/double classification is not exhaustive")
    return singles, doubles


def _norman_canonical_representatives(
    effects: np.ndarray,
    labels: np.ndarray,
    source_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Choose one outcome-independent measured row per canonical Norman gene.

    ``GENE+ctrl`` is preferred and ``ctrl+GENE`` is used only when the preferred
    cassette position is absent.  Grouping happens before role assignment, so
    the two measured positions of one gene can never cross library/candidate
    roles or create a tautological family-level top match.
    """

    matrix = np.asarray(effects)
    label_array = np.asarray(labels, dtype=str)
    sources = np.asarray(source_indices, dtype=int)
    if len(matrix) != len(label_array) or len(matrix) != len(sources):
        raise InputError("Norman representative axes do not align")
    label_to_index: dict[str, int] = {}
    genes: set[str] = set()
    for index, label in enumerate(label_array):
        parts = label.split("+")
        noncontrol = [part for part in parts if part != "ctrl"]
        if len(parts) != 2 or len(noncontrol) != 1 or label in label_to_index:
            raise InputError("Norman representative labels are malformed or duplicated")
        label_to_index[label] = index
        genes.add(noncontrol[0])
    selected_indices = []
    fallbacks = 0
    both_positions = 0
    for gene in sorted(genes):
        preferred = f"{gene}+ctrl"
        fallback = f"ctrl+{gene}"
        if preferred in label_to_index and fallback in label_to_index:
            both_positions += 1
        if preferred in label_to_index:
            selected_indices.append(label_to_index[preferred])
        elif fallback in label_to_index:
            selected_indices.append(label_to_index[fallback])
            fallbacks += 1
        else:  # pragma: no cover - genes came from the same labels
            raise InputError(f"Norman gene {gene} has no measured single row")
    chosen = np.asarray(selected_indices, dtype=int)
    return (
        matrix[chosen],
        label_array[chosen],
        sources[chosen],
        {
            "canonical_genes": int(len(chosen)),
            "genes_with_both_cassette_positions": both_positions,
            "opposite_position_fallbacks": fallbacks,
            "representative_rule": (
                "choose measured GENE+ctrl; use ctrl+GENE only when GENE+ctrl is absent; "
                "choice is frozen before row-role assignment and does not use outcomes"
            ),
        },
    )


def norman_split_sensitivity(
    payload: dict[str, Any],
    design: dict[str, Any],
    *,
    seeds: list[int],
    thresholds: list[float],
) -> dict[str, Any]:
    """Vary only the canonical-gene library against one fixed double catalog."""

    effects = np.asarray(payload["E"])
    labels = np.asarray(payload["labels"], dtype=str)
    is_single, is_double = _norman_membership(labels)
    single_source_indices = np.flatnonzero(is_single)
    representatives, representative_labels, representative_sources, representative_meta = (
        _norman_canonical_representatives(
            effects[is_single], labels[is_single], single_source_indices
        )
    )
    doubles = effects[is_double]
    double_sources = np.flatnonzero(is_double)
    library_size = int(design["acquisition_library_size"])
    catalog_size = int(design["sensitivity_catalog_doubles"])
    catalog_seed = int(design["sensitivity_catalog_seed"])
    feature_count = int(design["n_features"])
    if not 0 < library_size <= len(representatives):
        raise InputError("Norman sensitivity library size is invalid")
    if not 0 < catalog_size <= len(doubles):
        raise InputError("Norman sensitivity double-catalog size is invalid")
    catalog_idx = np.random.default_rng(catalog_seed).permutation(len(doubles))[
        :catalog_size
    ]
    catalog_sources = double_sources[catalog_idx]
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(int(seed))
        permutation = rng.permutation(len(representatives))
        library_idx = permutation[:library_size]
        unused_idx = permutation[library_size:]
        feature_indices = _library_only_features(
            representatives, library_idx, feature_count
        )
        library = np.asarray(
            representatives[library_idx][:, feature_indices], dtype=float
        )
        catalog = np.asarray(doubles[catalog_idx][:, feature_indices], dtype=float)
        coverage = coverage_report(library, catalog)
        best_single = _nonnegative_ray_cosines(catalog, library)
        common = _nonnegative_ray_cosines(
            catalog, library.mean(axis=0, keepdims=True)
        )
        signed = _signed_span_cosines(catalog, library)
        _assert_capacity_dominance(coverage.cosines, best_single, common, signed)
        rows.append(
            {
                "seed": int(seed),
                "library_representative_source_row_index_sha256": _index_sha256(
                    representative_sources[library_idx]
                ),
                "unused_representative_source_row_index_sha256": _index_sha256(
                    representative_sources[unused_idx]
                ),
                "library_canonical_labels_sha256": hashlib.sha256(
                    (
                        "\n".join(
                            _canonical_candidate_label(label)
                            for label in representative_labels[library_idx]
                        )
                        + "\n"
                    ).encode("utf-8")
                ).hexdigest(),
                "catalog_double_source_row_index_sha256": _index_sha256(catalog_sources),
                "feature_index_sha256": _index_sha256(feature_indices),
                "cone_mean_cosine": coverage.mean_cosine,
                "cone_strict_inside_fraction": coverage.strict_inside_cone_fraction,
                "cone_soft_coverage_by_threshold": _soft_threshold_fractions(
                    coverage.cosines, thresholds
                ),
                "best_single_mean_cosine": float(best_single.mean()),
                "common_response_mean_cosine": float(common.mean()),
                "signed_span_mean_cosine": float(signed.mean()),
                "cone_minus_best_single_mean_cosine": float(
                    coverage.mean_cosine - best_single.mean()
                ),
                "cone_minus_common_response_mean_cosine": float(
                    coverage.mean_cosine - common.mean()
                ),
                "signed_span_minus_cone_mean_cosine": float(
                    signed.mean() - coverage.mean_cosine
                ),
            }
        )
    metric_names = (
        "cone_mean_cosine",
        "cone_strict_inside_fraction",
        "best_single_mean_cosine",
        "common_response_mean_cosine",
        "signed_span_mean_cosine",
        "cone_minus_best_single_mean_cosine",
        "cone_minus_common_response_mean_cosine",
        "signed_span_minus_cone_mean_cosine",
    )
    aggregate = {
        name: _distribution_summary([float(row[name]) for row in rows])
        for name in metric_names
    }
    for threshold in thresholds:
        key = f"cosine_at_least_{threshold:g}"
        aggregate[f"cone_soft_coverage_{key}"] = _distribution_summary(
            [float(row["cone_soft_coverage_by_threshold"][key]) for row in rows]
        )
    return {
        "protocol": {
            "purpose": (
                "algorithmic current-library partition sensitivity against one fixed "
                "measured double catalog; not biological sampling uncertainty"
            ),
            "rng": "numpy.default_rng(seed)",
            "seeds": [int(seed) for seed in seeds],
            "canonical_single_genes": int(len(representatives)),
            "library_canonical_genes": library_size,
            "canonical_representative_selection": representative_meta,
            "eligible_double_effect_rows": int(doubles.shape[0]),
            "catalog_double_effect_rows": catalog_size,
            "catalog_seed": catalog_seed,
            "catalog_role": "fixed across every library-partition seed",
            "catalog_double_source_row_index_sha256": _index_sha256(catalog_sources),
            "feature_count": feature_count,
            "feature_selection_rows": (
                "current 40-gene measured-representative library only"
            ),
            "soft_cosine_thresholds": thresholds,
            "acquisition_recomputed": False,
        },
        "splits": rows,
        "aggregate": aggregate,
    }


def _paired_cosines(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return np.divide(
        np.sum(left * right, axis=1),
        denom,
        out=np.zeros(left.shape[0], dtype=float),
        where=denom > 0,
    )


def norman_single_to_double(
    payload: dict[str, Any],
    design: dict[str, Any],
) -> dict[str, Any]:
    """Retrospective K562 single-to-double additivity and coverage benchmark."""
    effects = payload["E"]
    labels = payload["labels"]
    is_single, is_double = _norman_membership(labels)
    single_source_indices = np.flatnonzero(is_single)
    double_source_indices = np.flatnonzero(is_double)
    expected = design["expected"]
    if len(single_source_indices) != expected["single_effect_rows"]:
        raise InputError("Norman single-effect count differs from the frozen design")
    if len(double_source_indices) != expected["double_effect_rows"]:
        raise InputError("Norman double-effect count differs from the frozen design")

    # The held-out double catalog never participates in feature selection.
    feature_count = int(design["n_features"])
    raw_singles = np.asarray(effects[is_single])
    raw_doubles = np.asarray(effects[is_double])
    feature_indices = _top_variance_features(raw_singles, feature_count)
    singles = np.asarray(raw_singles[:, feature_indices], dtype=float)
    doubles = np.asarray(raw_doubles[:, feature_indices], dtype=float)
    single_labels = labels[is_single]
    double_labels = labels[is_double]
    threshold = float(design["soft_cosine_threshold"])

    all_single_report = coverage_report(
        singles, doubles, soft_cosine_threshold=threshold
    )
    common = _nonnegative_ray_cosines(
        doubles, singles.mean(axis=0, keepdims=True)
    )
    best_single = _nonnegative_ray_cosines(doubles, singles)

    single_by_label: dict[str, np.ndarray] = {}
    unique_single_genes = set()
    for label, effect in zip(single_labels, singles):
        genes = [part for part in label.split("+") if part != "ctrl"]
        if len(genes) != 1 or label in single_by_label:
            raise InputError("Norman single-effect labels are malformed or duplicated")
        single_by_label[label] = effect
        unique_single_genes.add(genes[0])
    if len(unique_single_genes) != expected["unique_single_genes"]:
        raise InputError("Norman unique single-gene count differs from the frozen design")

    def constituent_effect(gene: str, position: int) -> tuple[np.ndarray, bool]:
        preferred = f"{gene}+ctrl" if position == 0 else f"ctrl+{gene}"
        fallback = f"ctrl+{gene}" if position == 0 else f"{gene}+ctrl"
        if preferred in single_by_label:
            return single_by_label[preferred], False
        if fallback in single_by_label:
            return single_by_label[fallback], True
        raise InputError(f"Norman double lacks a constituent single effect for {gene}")

    constituent_sums = []
    pair_cone_strict = []
    pair_cone_cosines = []
    position_fallbacks = 0
    for label, target in zip(double_labels, doubles):
        genes = label.split("+")
        first, first_fallback = constituent_effect(genes[0], 0)
        second, second_fallback = constituent_effect(genes[1], 1)
        position_fallbacks += int(first_fallback) + int(second_fallback)
        pair = np.vstack([first, second])
        constituent_sums.append(pair.sum(axis=0))
        pair_report = coverage_report(
            pair, target[None, :], soft_cosine_threshold=threshold
        )
        pair_cone_cosines.append(pair_report.mean_cosine)
        pair_cone_strict.append(bool(pair_report.strict_inside_cone[0]))
        if pair_report.soft_covered is None:
            raise RuntimeError("Norman pair-cone soft coverage unexpectedly missing")
    sum_cosines = _paired_cosines(np.vstack(constituent_sums), doubles)
    pair_cone_cosines_array = np.asarray(pair_cone_cosines)
    pair_cone_summary = _alignment_summary(pair_cone_cosines_array, threshold)
    pair_cone_summary.update(
        {
            "strict_inside_cone_count": int(np.sum(pair_cone_strict)),
            "strict_inside_cone_fraction": float(np.mean(pair_cone_strict)),
        }
    )

    seed = int(design["seed"])
    # Collapse the role-assignment axis to one fixed, measured representative per
    # canonical gene before splitting. Cassette-position siblings can therefore
    # never leak across library/candidate roles.
    representatives, representative_labels, representative_source_indices, representative_meta = (
        _norman_canonical_representatives(
            raw_singles, single_labels, single_source_indices
        )
    )
    permutation = np.random.default_rng(seed).permutation(len(representatives))
    library_size = int(design["acquisition_library_size"])
    candidate_size = int(design["acquisition_candidate_size"])
    if library_size + candidate_size > len(representatives):
        raise InputError("Norman acquisition split exceeds canonical single genes")
    library_idx = permutation[:library_size]
    candidate_idx = permutation[library_size : library_size + candidate_size]
    unused_idx = permutation[library_size + candidate_size :]
    library_canonical_labels = {
        _canonical_candidate_label(label)
        for label in representative_labels[library_idx]
    }
    candidate_canonical_labels = {
        _canonical_candidate_label(label)
        for label in representative_labels[candidate_idx]
    }
    role_overlap = library_canonical_labels & candidate_canonical_labels
    if role_overlap:
        raise InputError("Norman canonical genes overlap acquisition roles")
    acquisition_feature_indices = _library_only_features(
        representatives, library_idx, feature_count
    )
    acquisition_singles = np.asarray(
        representatives[:, acquisition_feature_indices], dtype=float
    )
    acquisition_doubles = np.asarray(
        raw_doubles[:, acquisition_feature_indices], dtype=float
    )
    ranking = rank_acquisitions(
        acquisition_singles[library_idx],
        acquisition_doubles,
        acquisition_singles[candidate_idx],
        realized=True,
        soft_cosine_threshold=threshold,
    )
    candidate_source_indices = representative_source_indices[candidate_idx]
    candidate_labels = representative_labels[candidate_idx]
    return {
        "source_effect_rows": int(effects.shape[0]),
        "single_effect_rows": int(len(singles)),
        "unique_single_genes": int(len(unique_single_genes)),
        "double_effect_rows": int(len(doubles)),
        "classification": "retrospective_single_to_double_additivity_and_coverage",
        "preprocessing": {
            "effect_definition": design["effect_definition"],
            "feature_count": feature_count,
            "feature_selection": design["feature_selection"],
            "feature_index_sha256": _index_sha256(feature_indices),
        },
        "coverage": _coverage_summary(all_single_report),
        "comparators": {
            "paired_constituent_sum": {
                **_alignment_summary(sum_cosines, threshold),
                "constituent_selection": (
                    "for A+B use A+ctrl and ctrl+B to match cassette position; "
                    "fall back to the available opposite-position single only when absent"
                ),
                "opposite_position_fallback_count": position_fallbacks,
            },
            "paired_constituent_two_atom_cone": {
                **pair_cone_summary,
                "constituent_selection": "same position-matched rule as paired_constituent_sum",
                "opposite_position_fallback_count": position_fallbacks,
            },
            "best_single_effect_nonnegative_ray": _alignment_summary(
                best_single, threshold
            ),
            "common_single_response_nonnegative_ray": _alignment_summary(
                common, threshold
            ),
        },
        "acquisition_design": {
            "rng": "numpy.default_rng(seed)",
            "seed": seed,
            "role_assignment_unit": "canonical_gene_with_one_fixed_measured_position_representative",
            "canonical_representative_selection": representative_meta,
            "library_candidate_canonical_gene_overlap_count": 0,
            "current_library_canonical_genes": library_size,
            "supplied_candidate_canonical_genes": candidate_size,
            "unused_canonical_genes": int(len(unused_idx)),
            "catalog_double_effect_rows": int(len(doubles)),
            "feature_count": feature_count,
            "feature_selection": (
                "stable top-variance coordinates from current canonical-gene acquisition "
                "library representatives only"
            ),
            "feature_selection_rows": "current_library_only",
            "feature_index_sha256": _index_sha256(acquisition_feature_indices),
            "library_source_row_index_sha256": _index_sha256(
                representative_source_indices[library_idx]
            ),
            "candidate_source_row_index_sha256": _index_sha256(
                candidate_source_indices
            ),
            "catalog_source_row_index_sha256": _index_sha256(double_source_indices),
            "library_canonical_labels_sha256": hashlib.sha256(
                ("\n".join(
                    _canonical_candidate_label(label)
                    for label in representative_labels[library_idx]
                ) + "\n").encode("utf-8")
            ).hexdigest(),
            "candidate_canonical_labels_sha256": hashlib.sha256(
                ("\n".join(
                    _canonical_candidate_label(label)
                    for label in candidate_labels
                ) + "\n").encode("utf-8")
            ).hexdigest(),
        },
        "acquisition": _acquisition_summary(
            ranking, candidate_source_indices, candidate_labels
        ),
        "_cosines": all_single_report.cosines,
        "_certificate_scores": ranking.certificate_score,
        "_realized_mean_cosine_gains": ranking.realized_mean_cosine_gain,
    }


def run(config: dict[str, Any], *, cache_dir: Path | None = None) -> dict[str, Any]:
    _validate_config(config)
    paths: dict[str, Path] = {}
    identities: dict[str, Any] = {}
    # Verify every byte identity before loading any cache or computing metrics.
    for name, spec in config["inputs"].items():
        _progress(f"verifying {name} cache bytes and SHA-256")
        path, identity = verify_input(name, spec, cache_dir=cache_dir)
        paths[name] = path
        identities[name] = identity

    matrices = {
        name: load_matrix(paths[name], config["inputs"][name])
        for name in ("zhu", "norman", "replogle")
    }
    threshold = float(config["analysis"]["soft_cosine_threshold"])
    dataset_config = config["datasets"]

    _progress("Zhu retrospective measured-dictionary compression/self-reconstruction")
    zhu_design = dict(dataset_config["zhu_crispri_tcell"])
    zhu_design["soft_cosine_threshold"] = threshold
    zhu = build_zhu_from_slice(matrices["zhu"], zhu_design)

    _progress("Norman canonical K562 retrospective single-to-double benchmark")
    norman_design = dict(dataset_config["norman_k562_crispra"])
    norman_design["soft_cosine_threshold"] = threshold
    norman = norman_single_to_double(matrices["norman"], norman_design)

    _progress("Replogle K562 retrospective measured-dictionary compression/self-reconstruction")
    replogle_design = dict(dataset_config["replogle_k562_essential_crispri"])
    replogle_design["soft_cosine_threshold"] = threshold
    replogle = standard_split(matrices["replogle"], replogle_design)
    replogle["source_effect_rows"] = int(matrices["replogle"]["E"].shape[0])
    replogle["preprocessing"] = {
        "effect_definition": replogle_design["effect_definition"]
    }

    sensitivity = config["analysis"]["split_sensitivity"]
    seeds = [int(seed) for seed in sensitivity["seeds"]]
    sensitivity_thresholds = [
        float(value) for value in sensitivity["soft_cosine_thresholds"]
    ]
    _progress("Zhu coverage-only sensitivity across 12 row partitions")
    zhu["split_sensitivity"] = standard_split_sensitivity(
        matrices["zhu"],
        zhu_design,
        seeds=seeds,
        thresholds=sensitivity_thresholds,
    )
    _progress("Norman coverage-only sensitivity across 12 canonical-gene libraries")
    norman["split_sensitivity"] = norman_split_sensitivity(
        matrices["norman"],
        norman_design,
        seeds=seeds,
        thresholds=sensitivity_thresholds,
    )
    _progress("Replogle coverage-only sensitivity across 12 row partitions")
    replogle["split_sensitivity"] = standard_split_sensitivity(
        matrices["replogle"],
        replogle_design,
        seeds=seeds,
        thresholds=sensitivity_thresholds,
    )

    return {
        "schema_version": config["schema_version"],
        "generated_on": config["generated_on"],
        "status": "PASS",
        "benchmark": config["benchmark"],
        "claim_ceiling": config["claim_ceiling"],
        "provenance": config["provenance"],
        "input_verification": identities,
        "analysis_contract": config["analysis"],
        "datasets": {
            "zhu_crispri_tcell": zhu,
            "norman_k562_crispra": norman,
            "replogle_k562_essential_crispri": replogle,
        },
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _json_safe(item)
            for key, item in value.items()
            if not key.startswith("_")
        }
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _render_report(results: dict[str, Any]) -> str:
    return json.dumps(
        _json_safe(results), indent=2, sort_keys=True, allow_nan=False
    ) + "\n"


def make_figure(results: dict[str, Any], out_path: Path) -> None:
    """Render an optional descriptive figure from the claim-bearing run."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    datasets = results["datasets"]
    order = [
        "zhu_crispri_tcell",
        "norman_k562_crispra",
        "replogle_k562_essential_crispri",
    ]
    labels = ["Zhu CD4⁺ T", "Norman K562", "Replogle K562"]
    legends = [
        "Zhu — CRISPRi",
        "Norman — CRISPRa",
        "Replogle — essential-gene CRISPRi",
    ]
    colors = ["#B2182B", "#1b7837", "#2166AC"]
    threshold = results["analysis_contract"]["soft_cosine_threshold"]

    fig, (axis_a, axis_b, axis_c) = plt.subplots(1, 3, figsize=(12.6, 4.2))
    for name, legend, color in zip(order, legends, colors):
        sorted_cosines = np.sort(datasets[name]["_cosines"])[::-1]
        axis_a.plot(
            np.linspace(0, 1, len(sorted_cosines)),
            sorted_cosines,
            color=color,
            lw=2.0,
            label=legend,
        )
    axis_a.axhline(threshold, color="#5a5a5a", ls=(0, (4, 2)), lw=1.0)
    axis_a.set_xlabel("Measured target-row rank (normalized)")
    axis_a.set_ylabel("Best cone cosine")
    axis_a.set_title("Retrospective alignment, not strict reachability", loc="left")
    axis_a.legend(loc="lower left", fontsize=6, frameon=False)

    for name, label, color in zip(order, labels, colors):
        certificate = datasets[name]["_certificate_scores"]
        realized = datasets[name]["_realized_mean_cosine_gains"]
        certificate_scaled = (certificate - certificate.min()) / (
            certificate.max() - certificate.min() + 1e-12
        )
        realized_scaled = (realized - realized.min()) / (
            realized.max() - realized.min() + 1e-12
        )
        rho = datasets[name]["acquisition"][
            "certificate_vs_realized_mean_cosine_gain_spearman"
        ]
        axis_b.scatter(
            certificate_scaled,
            realized_scaled,
            s=20,
            color=color,
            alpha=0.65,
            edgecolor="none",
            label=f"{label}: ρ={rho:.2f}",
        )
    axis_b.plot([0, 1], [0, 1], color="#5a5a5a", ls=":", lw=0.8)
    axis_b.set_xlabel("Certificate score (normalized)")
    axis_b.set_ylabel("Realized mean-cosine gain (normalized)")
    axis_b.set_title("Certificate vs explicit retrospective objective", loc="left")
    axis_b.legend(loc="lower right", fontsize=6, frameon=False)

    soft = [datasets[name]["coverage"]["soft_coverage_fraction"] for name in order]
    rhos = [
        datasets[name]["acquisition"][
            "certificate_vs_realized_mean_cosine_gain_spearman"
        ]
        for name in order
    ]
    for index, (fraction, rho, color) in enumerate(zip(soft, rhos, colors)):
        axis_c.bar(index - 0.19, fraction, 0.38, color=color, alpha=0.45)
        axis_c.bar(index + 0.19, rho, 0.38, color=color, alpha=1.0)
    axis_c.set_xticks(range(3))
    axis_c.set_xticklabels(labels, fontsize=6)
    axis_c.set_ylim(0, 1.1)
    axis_c.set_ylabel("value")
    axis_c.set_title("Soft coverage and rank association", loc="left")
    axis_c.legend(
        handles=[
            mpatches.Patch(
                facecolor="0.6", alpha=0.5, label=f"soft coverage (cos ≥ {threshold})"
            ),
            mpatches.Patch(
                facecolor="0.25", label="Spearman ρ vs mean-cosine gain"
            ),
        ],
        loc="lower center",
        fontsize=6,
        frameon=False,
    )
    for axis, letter in zip((axis_a, axis_b, axis_c), "abc"):
        axis.text(
            -0.12,
            1.05,
            letter,
            transform=axis.transAxes,
            fontweight="bold",
            fontsize=11,
            va="bottom",
            ha="right",
        )
    fig.tight_layout(w_pad=2.4)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--slice-dir",
        type=Path,
        default=None,
        help="optional directory override containing the three frozen cache basenames",
    )
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify exact maintained output instead of overwriting it",
    )
    parser.add_argument("--figure", type=Path, default=None)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    observed = run(config, cache_dir=args.slice_dir)
    observed["config_sha256"] = sha256_file(args.config)
    payload = _render_report(observed)
    if args.check:
        if not args.report.is_file():
            raise SystemExit(f"--check requires an existing report at {args.report}")
        if args.report.read_text(encoding="utf-8") != payload:
            raise SystemExit("cross-dataset benchmark output differs from frozen report")
        print("cross-dataset benchmark report matches exactly")
    else:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8")
        print(f"wrote {args.report}")

    if args.figure is not None:
        args.figure.parent.mkdir(parents=True, exist_ok=True)
        make_figure(observed, args.figure)
        print(f"wrote {args.figure}")


if __name__ == "__main__":
    main()
