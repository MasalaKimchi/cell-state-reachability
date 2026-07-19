#!/usr/bin/env python3
"""Run reciprocal transfer across Zhu's released guide-rank modalities.

The official dataset card defines ``guide_1`` and ``guide_2`` as the first and
second alphanumeric sgRNA IDs within each target-condition pair.  The H5MU does
not embed those identifiers, so this benchmark measures robustness across the
released ranked summaries without a verified ID crosswalk.  Upstream
effectiveness selection means it is not leakage-safe guide generalization,
named-sgRNA replication, or evidence of biological state reachability.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError  # noqa: E402
from scripts.run_donor_pair_transfer import (  # noqa: E402
    _decode,
    _names_hash,
    _vector_hash,
    fit_training_models,
    load_modality_matrix,
    score_frozen_models,
    verify_input,
)
from scripts.run_source_reconstruction import (  # noqa: E402
    build_target_lineage,
    load_target_table,
)


DEFAULT_CONFIG = ROOT / "configs" / "guide_pair_transfer.json"
GUIDE_MODALITIES = ("guide_1", "guide_2")
DOCUMENTED_CONDITIONS = ("Rest", "Stim8hr", "Stim48hr")
REQUIRED_LAYERS = (
    "adj_p_value",
    "baseMean",
    "lfcSE",
    "log_fc",
    "p_value",
    "zscore",
)
REQUIRED_CATEGORICALS = (
    "culture_condition",
    "target_contrast",
    "target_contrast_gene_name",
    "chunk",
    "n_total_genes_category",
)
REQUIRED_MODALITY_FIELDS = (
    "X",
    "layers",
    "obs",
    "obsm",
    "obsp",
    "uns",
    "var",
    "varm",
    "varp",
)
REQUIRED_OBS_FIELDS = (
    "chunk",
    "culture_condition",
    "distal_offtarget_flag",
    "donor_correlation_all_mean",
    "donor_correlation_all_min",
    "donor_correlation_hits_mean",
    "donor_correlation_hits_min",
    "guide_correlation_all",
    "guide_correlation_all_pval",
    "guide_correlation_signif",
    "guide_correlation_signif_pval",
    "guide_n_signif_ontarget",
    "low_target_gex",
    "n_cells_target",
    "n_down_genes",
    "n_downstream",
    "n_guides",
    "n_total_de_genes",
    "n_total_genes_category",
    "n_up_genes",
    "neighboring_gene_KD",
    "ontarget_effect_size",
    "ontarget_significant",
    "single_guide_estimate",
    "target_baseMean",
    "target_condition",
    "target_contrast",
    "target_contrast_gene_name",
)
REQUIRED_VAR_FIELDS = ("_index", "gene_ids", "gene_name")
REQUIRED_OBS_COLUMN_ORDER = (
    "target_contrast_gene_name",
    "culture_condition",
    "target_contrast",
    "chunk",
    "n_cells_target",
    "n_up_genes",
    "n_down_genes",
    "n_total_de_genes",
    "ontarget_effect_size",
    "ontarget_significant",
    "target_baseMean",
    "neighboring_gene_KD",
    "n_total_genes_category",
    "distal_offtarget_flag",
    "low_target_gex",
    "n_guides",
    "single_guide_estimate",
    "n_downstream",
    "guide_correlation_signif",
    "guide_correlation_signif_pval",
    "guide_correlation_all",
    "guide_correlation_all_pval",
    "guide_n_signif_ontarget",
    "donor_correlation_all_mean",
    "donor_correlation_all_min",
    "donor_correlation_hits_mean",
    "donor_correlation_hits_min",
)
REQUIRED_ROOT_ATTRIBUTES = {
    "axis": 0,
    "encoding-type": "MuData",
    "encoding-version": "0.1.0",
    "encoder": "mudata",
    "encoder-version": "0.3.1",
}
REQUIRED_MODALITY_ATTRIBUTES = {
    "encoding-type": "anndata",
    "encoding-version": "0.1.0",
    "encoder": "mudata",
    "encoder-version": "0.3.1",
}
REQUIRED_OBS_ATTRIBUTES = {
    "encoding-type": "dataframe",
    "encoding-version": "0.2.0",
    "_index": "target_condition",
    "column-order": list(REQUIRED_OBS_COLUMN_ORDER),
}
REQUIRED_VAR_ATTRIBUTES = {
    "encoding-type": "dataframe",
    "encoding-version": "0.2.0",
    "_index": "_index",
    "column-order": ["gene_ids", "gene_name"],
}
REQUIRED_CATEGORICAL_ATTRIBUTES = {
    "encoding-type": "categorical",
    "encoding-version": "0.2.0",
    "ordered": False,
}
TARGET_SOURCES = {
    "hollbacker": "Th2_vs_Th1 (Hollbacker 2021)",
    "ota": "Th2_vs_Th1 (Ota 2021)",
}
GENE_ORDER = "lexicographic_gene_symbol"
SPLIT_RNG = "numpy.default_rng(seed)"
SPLIT_SEEDS = (0, 1, 2)
SCORE_TARGET_SOURCES = "same source and the opposite source for every frozen fit"
REPRESENTATION_SENSITIVITY = {
    "train_guide_slot": "guide_1",
    "train_target_source": "hollbacker",
    "seed": 0,
    "alternative_atom_order": "reverse_lexicographic",
}
AUTHOR_SELECTION_PROVENANCE = {
    "initial_sample_selection_rule": (
        "keep_min_cells & keep_effective_guides & keep_total_counts"
    ),
    "selection_field": "keep_effective_guides",
    "minimum_passing_replicates_per_guide_condition": 3,
    "minimum_cells_per_guide_per_condition_sample": 5,
    "target_condition_membership_source": (
        "cond_targets from unreleased for_DE_by_guide.csv"
    ),
    "required_testable_guides_per_target_condition": 2,
    "final_sample_selection_rule": "initial keep_for_DE & keep_test_genes",
    "de_formula": "~ log10_n_cells + target",
    "donor_term_in_de_formula": False,
    "unreleased_intermediate": "for_DE_by_guide.csv",
    "official_modality_mapping": (
        "guide_1 is the lowest alphanumeric sgRNA ID and guide_2 the second within "
        "each perturbed-gene/culture-condition pair; targets with one passing guide "
        "appear only in guide_1"
    ),
    "h5mu_identity_field_status": (
        "no guide ID or sequence field is embedded in the guide-level H5MU"
    ),
    "public_identity_sources": [
        "GWCD4i.pseudobulk_merged.h5ad obs/guide_id",
        "sgrna_library_metadata.suppl_table.csv sgRNA",
    ],
    "identity_crosswalk_status": (
        "exact ranked-modality-to-sgRNA IDs not reconstructed or hash-verified in "
        "this benchmark"
    ),
    "guide_identity_status": (
        "OFFICIAL_ALPHANUMERIC_GUIDE_RANKS_IDS_NOT_EMBEDDED_"
        "CROSSWALK_NOT_VERIFIED"
    ),
}
MODELS = ("cone", "training_common_ray", "training_best_single", "zero")
BASELINES = ("zero", "training_common_ray", "training_best_single")
METRICS = ("cosine", "normalized_rmse", "norm_ratio", "sign_agreement")
INFERENCE_CEILING = (
    "descriptive correlated challenges only; no p-values, confidence intervals, "
    "physical-guide generalization, or donor-population inference"
)
MISSING_CATEGORY = "<MISSING_CATEGORY>"


def _resolve(path: str | Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else ROOT / candidate


def _text_attribute(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="strict")
    return str(value)


def _attribute_strings(value: Any) -> tuple[str, ...]:
    array = np.asarray(value, dtype=object).reshape(-1)
    return tuple(_text_attribute(item) for item in array)


def _attribute_matches(observed: Any, expected: Any) -> bool:
    if isinstance(expected, list):
        return _attribute_strings(observed) == tuple(str(value) for value in expected)
    if isinstance(expected, bool):
        return bool(observed) is expected
    if isinstance(expected, int):
        try:
            return int(observed) == expected
        except (TypeError, ValueError):
            return False
    return _text_attribute(observed) == str(expected)


def _encoded_column_length(column: Any) -> int:
    """Return row length without reading unused outcome-bearing values."""

    if hasattr(column, "keys"):
        if "codes" not in column:
            raise InputError("encoded observation group has no codes dataset")
        return int(len(column["codes"]))
    if not hasattr(column, "shape") or len(column.shape) == 0:
        raise InputError("observation field is not a row-aligned encoded column")
    return int(column.shape[0])


def _categorical_allow_missing(frame: Any, field: str) -> tuple[np.ndarray, np.ndarray]:
    """Decode AnnData categoricals while preserving code -1 as explicit missingness."""

    column = frame[field]
    if not hasattr(column, "keys") or set(column.keys()) != {"categories", "codes"}:
        raise InputError(f"{field} must use exact categorical H5 encoding")
    categories = _decode(column["categories"][:])
    if any(not str(value) for value in categories) or len(set(categories)) != len(
        categories
    ):
        raise InputError(f"{field} categories must be nonempty and unique")
    codes = np.asarray(column["codes"][:], dtype=np.int64)
    if codes.ndim != 1 or np.any(codes < -1) or np.any(codes >= len(categories)):
        raise InputError(f"{field} contains an invalid categorical code")
    missing = codes == -1
    values = np.empty(len(codes), dtype=object)
    values[missing] = MISSING_CATEGORY
    values[~missing] = categories[codes[~missing]]
    return values, missing


def _documented_condition_suffix(target_condition: Any) -> str:
    """Read only the documented terminal condition suffix from an opaque key.

    This is a categorical-integrity check, not an attempt to parse or recover an
    sgRNA identity from ``target_condition``.
    """

    key = str(target_condition)
    matches = [
        condition
        for condition in DOCUMENTED_CONDITIONS
        if key.endswith(f"_{condition}")
    ]
    if len(matches) != 1:
        raise InputError(
            "target_condition key must have exactly one documented terminal "
            "condition suffix"
        )
    return matches[0]


def _canonical_hash(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def validate_config(config: dict[str, Any]) -> None:
    """Reject semantic drift before touching either registered input."""

    if config.get("schema_version") != "1.0.0":
        raise InputError("guide benchmark config schema differs")
    for field, expected_value in AUTHOR_SELECTION_PROVENANCE.items():
        if config["author_provenance"].get(field) != expected_value:
            raise InputError(f"author guide-selection provenance {field} differs")
    if config["h5mu"]["required_root_attributes"] != REQUIRED_ROOT_ATTRIBUTES:
        raise InputError("guide H5MU root-attribute contract differs")
    if config["h5mu"]["required_modality_attributes"] != REQUIRED_MODALITY_ATTRIBUTES:
        raise InputError("guide H5MU modality-attribute contract differs")
    if config["h5mu"]["required_obs_attributes"] != REQUIRED_OBS_ATTRIBUTES:
        raise InputError("guide H5MU observation-attribute contract differs")
    if config["h5mu"]["required_var_attributes"] != REQUIRED_VAR_ATTRIBUTES:
        raise InputError("guide H5MU variable-attribute contract differs")
    if (
        config["h5mu"]["required_categorical_attributes"]
        != REQUIRED_CATEGORICAL_ATTRIBUTES
    ):
        raise InputError("guide H5MU categorical-attribute contract differs")
    if tuple(config["h5mu"]["modalities"]) != GUIDE_MODALITIES:
        raise InputError("guide modalities must be exactly guide_1 and guide_2")
    if tuple(config["h5mu"]["required_layers"]) != REQUIRED_LAYERS:
        raise InputError("guide H5MU layer contract differs")
    if tuple(config["h5mu"]["required_categoricals"]) != REQUIRED_CATEGORICALS:
        raise InputError("guide H5MU categorical contract differs")
    if tuple(config["h5mu"]["required_modality_fields"]) != REQUIRED_MODALITY_FIELDS:
        raise InputError("guide H5MU modality-field contract differs")
    if tuple(config["h5mu"]["required_obs_fields"]) != REQUIRED_OBS_FIELDS:
        raise InputError("guide H5MU observation-field contract differs")
    if tuple(config["h5mu"]["required_var_fields"]) != REQUIRED_VAR_FIELDS:
        raise InputError("guide H5MU var-field contract differs")
    if tuple(config["h5mu"]["selected_rest_mapping_fields"]) != (
        "target_contrast",
        "target_contrast_gene_name",
    ):
        raise InputError("selected-Rest mapping contract differs")
    if tuple(config["h5mu"]["required_mod_order"]) != GUIDE_MODALITIES:
        raise InputError("MuData modality order contract differs")
    if config["h5mu"]["observation_index"] != "target_condition":
        raise InputError("observation index must remain target_condition")
    if config["h5mu"]["condition_field"] != "culture_condition":
        raise InputError("condition field must remain culture_condition")
    if config["h5mu"]["layer"] != "log_fc":
        raise InputError("the reciprocal benchmark is frozen to Rest log_fc")
    if config["h5mu"]["condition"] != "Rest":
        raise InputError("the reciprocal benchmark is frozen to Rest")
    if config["target"] != {
        "gene_field": "variable",
        "contrast_field": "contrast",
        "sources": TARGET_SOURCES,
        "value_fields": ["log_fc"],
        "orientation_multiplier": -1.0,
    }:
        raise InputError("guide benchmark target-direction contract differs")
    analysis = config["analysis"]
    if analysis["gene_order"] != GENE_ORDER:
        raise InputError("guide benchmark gene-order contract differs")
    if analysis["rng"] != SPLIT_RNG:
        raise InputError("guide benchmark split-RNG contract differs")
    if analysis["score_target_sources"] != SCORE_TARGET_SOURCES:
        raise InputError("guide benchmark score-source contract differs")
    if analysis["representation_sensitivity"] != REPRESENTATION_SENSITIVITY:
        raise InputError("guide benchmark representation-sensitivity contract differs")
    if tuple(analysis["baselines"]) != BASELINES:
        raise InputError("guide benchmark baseline contract differs")
    if analysis["inference"] != INFERENCE_CEILING:
        raise InputError("guide benchmark inference ceiling differs")
    if tuple(config["analysis"]["train_guide_slots"]) != GUIDE_MODALITIES:
        raise InputError("both guide slots must be trained reciprocally")

    sources = tuple(analysis["train_target_sources"])
    if sources != tuple(TARGET_SOURCES):
        raise InputError("target sources must remain in frozen hollbacker/ota order")
    raw_seeds = analysis["split_seeds"]
    if (
        any(type(value) is not int for value in raw_seeds)
        or tuple(raw_seeds) != SPLIT_SEEDS
    ):
        raise InputError("split seeds must be the exact integer sequence 0, 1, 2")
    seeds = tuple(raw_seeds)

    expected = config["expected"]
    if int(expected["modalities"]) != 2:
        raise InputError("the guide benchmark requires exactly two modalities")
    if set(expected["modality_rows"]) != set(GUIDE_MODALITIES):
        raise InputError("expected modality-row keys differ")
    if set(expected["condition_counts"]) != set(GUIDE_MODALITIES):
        raise InputError("expected condition-count keys differ")
    if set(expected["categorical_missing_rows"]) != set(GUIDE_MODALITIES):
        raise InputError("expected categorical-missing-row keys differ")
    if set(expected["missing_key_condition_suffix_counts"]) != set(
        GUIDE_MODALITIES
    ):
        raise InputError("expected missing-key suffix-count keys differ")
    for modality in GUIDE_MODALITIES:
        counts = expected["condition_counts"][modality]
        if set(counts) != set(DOCUMENTED_CONDITIONS):
            raise InputError(f"{modality} condition contract differs")
        missing_rows = int(expected["categorical_missing_rows"][modality])
        missing_suffix_counts = expected["missing_key_condition_suffix_counts"][
            modality
        ]
        if set(missing_suffix_counts) != set(DOCUMENTED_CONDITIONS) or sum(
            int(value) for value in missing_suffix_counts.values()
        ) != missing_rows:
            raise InputError(f"{modality} missing-key suffix counts differ")
        if missing_rows < 0 or sum(
            int(value) for value in counts.values()
        ) + missing_rows != int(expected["modality_rows"][modality]):
            raise InputError(
                f"{modality} condition and missing counts do not sum to its rows"
            )
    if int(expected["rest_common_atoms"]) != int(
        expected["condition_counts"]["guide_2"]["Rest"]
    ):
        raise InputError("guide_2 Rest rows must equal the common Rest universe")
    if int(expected["condition_counts"]["guide_1"]["Rest"]) != int(
        expected["rest_common_atoms"]
    ) + int(expected["guide_1_only_rest_atoms"]):
        raise InputError("guide_1 Rest attrition arithmetic differs")
    for field in (
        "guide_1_rest_keys_sha256",
        "guide_2_rest_keys_sha256",
        "ordered_gene_ids_sha256",
        "ordered_gene_symbols_sha256",
    ):
        if re.fullmatch(r"[0-9a-f]{64}", str(expected[field])) is None:
            raise InputError(f"{field} must be a lowercase SHA-256")

    fits = len(GUIDE_MODALITIES) * len(sources) * len(seeds)
    if int(expected["fits"]) != fits:
        raise InputError("expected fit count differs from the exhaustive design")
    if int(expected["representation_sensitivity_fits"]) != 1:
        raise InputError("exactly one atom-order sensitivity fit is required")
    if int(expected["challenge_rows"]) != 2 * fits:
        raise InputError("expected challenge-row count differs from the design")
    if config["author_provenance"]["guide_identity_status"] != (
        "OFFICIAL_ALPHANUMERIC_GUIDE_RANKS_IDS_NOT_EMBEDDED_"
        "CROSSWALK_NOT_VERIFIED"
    ):
        raise InputError("guide identity/rank crosswalk status differs")


def verify_registered_input(
    spec: dict[str, Any], *, verify_hash: bool = True
) -> dict[str, Any]:
    """Require a real frozen SHA before an expensive verified file scan."""

    if (
        verify_hash
        and re.fullmatch(r"[0-9a-f]{64}", str(spec.get("sha256", ""))) is None
    ):
        raise InputError(
            f"{spec['path']} has no frozen SHA-256; verified execution is disabled"
        )
    return verify_input(spec, verify_hash=verify_hash)


def profile_h5mu(path: Path, config: dict[str, Any]) -> dict[str, Any]:
    """Validate H5MU structure and align opaque target-condition keys.

    Observation matching uses exact string identity, and the smaller guide_2 key
    set must be a metadata-consistent subset of guide_1.  The only parsed content
    is the documented terminal Rest/Stim8hr/Stim48hr suffix, used to audit
    categorical missingness; no guide identity is parsed or inferred.
    """

    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError(
            "guide-position transfer requires requirements-external.txt"
        ) from exc

    spec = config["h5mu"]
    expected = config["expected"]
    profiles: dict[str, dict[str, Any]] = {}
    reference_gene_ids: np.ndarray | None = None
    reference_symbols: np.ndarray | None = None

    with h5py.File(path, "r") as handle:
        for attribute, expected_value in spec["required_root_attributes"].items():
            if not _attribute_matches(handle.attrs.get(attribute, ""), expected_value):
                raise InputError(f"guide object root attribute {attribute} differs")
        if "mod" not in handle:
            raise InputError("guide object has no modality group")
        if _attribute_strings(handle["mod"].attrs.get("mod-order", ())) != tuple(
            spec["required_mod_order"]
        ):
            raise InputError("guide object mod-order attribute differs")
        if set(handle["mod"].keys()) != set(GUIDE_MODALITIES):
            raise InputError("guide H5MU modalities differ from guide_1/guide_2")

        for modality in GUIDE_MODALITIES:
            group = handle[f"mod/{modality}"]
            for attribute, expected_value in spec[
                "required_modality_attributes"
            ].items():
                if not _attribute_matches(
                    group.attrs.get(attribute, ""), expected_value
                ):
                    raise InputError(f"{modality} attribute {attribute} differs")
            if set(group.keys()) != set(REQUIRED_MODALITY_FIELDS):
                raise InputError(f"{modality} top-level schema differs")
            obs = group["obs"]
            var = group["var"]
            layers = group["layers"]
            for attribute, expected_value in spec["required_obs_attributes"].items():
                if not _attribute_matches(obs.attrs.get(attribute, ""), expected_value):
                    raise InputError(f"{modality} obs attribute {attribute} differs")
            for attribute, expected_value in spec["required_var_attributes"].items():
                if not _attribute_matches(var.attrs.get(attribute, ""), expected_value):
                    raise InputError(f"{modality} var attribute {attribute} differs")
            if set(obs.keys()) != set(REQUIRED_OBS_FIELDS):
                raise InputError(f"{modality} observation schema differs")
            if set(var.keys()) != set(REQUIRED_VAR_FIELDS):
                raise InputError(f"{modality} var schema differs")
            index_field = _text_attribute(obs.attrs.get("_index", ""))
            if index_field != spec["observation_index"] or index_field not in obs:
                raise InputError(f"{modality} observation index differs")

            keys = _decode(obs[index_field][:])
            n_rows = len(keys)
            if n_rows != int(expected["modality_rows"][modality]):
                raise InputError(f"{modality} row count differs")
            if any(not str(key) for key in keys) or len(set(keys)) != n_rows:
                raise InputError(
                    f"{modality} observation keys must be nonempty and unique"
                )
            for field in REQUIRED_OBS_FIELDS:
                if _encoded_column_length(obs[field]) != n_rows:
                    raise InputError(f"{modality}/{field} row count differs")

            categorical_values: dict[str, np.ndarray] = {}
            categorical_missing: dict[str, np.ndarray] = {}
            for field in REQUIRED_CATEGORICALS:
                if field not in obs:
                    raise InputError(f"{modality} is missing categorical {field}")
                for attribute, expected_value in spec[
                    "required_categorical_attributes"
                ].items():
                    if not _attribute_matches(
                        obs[field].attrs.get(attribute, ""), expected_value
                    ):
                        raise InputError(
                            f"{modality}/{field} categorical attribute {attribute} differs"
                        )
                values, missing = _categorical_allow_missing(obs, field)
                if len(values) != n_rows:
                    raise InputError(f"{modality}/{field} length differs")
                categorical_values[field] = values
                categorical_missing[field] = missing
            missing_reference = categorical_missing[spec["condition_field"]]
            for field, missing in categorical_missing.items():
                if not np.array_equal(missing, missing_reference):
                    raise InputError(
                        f"{modality}/{field} categorical missingness mask differs"
                    )
            if int(np.count_nonzero(missing_reference)) != int(
                expected["categorical_missing_rows"][modality]
            ):
                raise InputError(f"{modality} categorical missing-row count differs")
            conditions = categorical_values[spec["condition_field"]]
            condition_counts = {
                condition: int(np.sum(conditions == condition))
                for condition in expected["condition_counts"][modality]
            }
            observed_conditions = {
                str(value) for value in conditions if value != MISSING_CATEGORY
            }
            if observed_conditions != set(condition_counts):
                raise InputError(f"{modality} condition labels differ")
            if condition_counts != {
                key: int(value)
                for key, value in expected["condition_counts"][modality].items()
            }:
                raise InputError(f"{modality} condition counts differ")

            suffix_conditions = np.asarray(
                [_documented_condition_suffix(key) for key in keys], dtype=object
            )
            nonmissing = ~missing_reference
            if not np.array_equal(suffix_conditions[nonmissing], conditions[nonmissing]):
                raise InputError(
                    f"{modality} target_condition suffix disagrees with categorical "
                    "culture_condition"
                )
            missing_key_condition_suffix_counts = {
                condition: int(
                    np.count_nonzero(missing_reference & (suffix_conditions == condition))
                )
                for condition in DOCUMENTED_CONDITIONS
            }
            if missing_key_condition_suffix_counts != {
                key: int(value)
                for key, value in expected["missing_key_condition_suffix_counts"][
                    modality
                ].items()
            }:
                raise InputError(
                    f"{modality} missing-key condition suffix counts differ"
                )

            for field in (spec["gene_id_field"], spec["gene_symbol_field"]):
                if field not in var:
                    raise InputError(f"{modality} is missing var/{field}")
            gene_ids = _decode(var[spec["gene_id_field"]][:])
            symbols = _decode(var[spec["gene_symbol_field"]][:])
            if len(gene_ids) != int(expected["genes"]) or len(symbols) != int(
                expected["genes"]
            ):
                raise InputError(f"{modality} gene count differs")
            if (
                any(not str(value) for value in gene_ids)
                or any(not str(value) for value in symbols)
                or len(set(gene_ids)) != len(gene_ids)
                or len(set(symbols)) != len(symbols)
            ):
                raise InputError(f"{modality} gene IDs and symbols must be unique")
            if (
                _names_hash(tuple(str(value) for value in gene_ids))
                != expected["ordered_gene_ids_sha256"]
            ):
                raise InputError(f"{modality} ordered gene-ID hash differs")
            if (
                _names_hash(tuple(str(value) for value in symbols))
                != expected["ordered_gene_symbols_sha256"]
            ):
                raise InputError(f"{modality} ordered gene-symbol hash differs")
            if reference_gene_ids is None:
                reference_gene_ids = gene_ids
                reference_symbols = symbols
            elif not np.array_equal(gene_ids, reference_gene_ids) or not np.array_equal(
                symbols, reference_symbols
            ):
                raise InputError(
                    "guide modalities do not share identical ordered genes"
                )

            if set(layers.keys()) != set(REQUIRED_LAYERS):
                raise InputError(
                    f"{modality} layer names differ from the exact contract"
                )
            for layer in REQUIRED_LAYERS:
                if layers[layer].shape != (n_rows, len(gene_ids)):
                    raise InputError(f"{modality}/{layer} shape differs")

            key_to_row = {str(key): index for index, key in enumerate(keys)}
            selected_mapping_by_key = {
                str(key): tuple(
                    str(categorical_values[field][index])
                    for field in spec["selected_rest_mapping_fields"]
                )
                for index, key in enumerate(keys)
                if conditions[index] == spec["condition"]
            }
            rest_keys = tuple(
                sorted(
                    str(key)
                    for key, condition in zip(keys, conditions)
                    if condition == spec["condition"]
                )
            )
            profiles[modality] = {
                "rows": n_rows,
                "condition_counts": condition_counts,
                "categorical_missing_rows": int(np.count_nonzero(missing_reference)),
                "missing_key_condition_suffix_counts": (
                    missing_key_condition_suffix_counts
                ),
                "key_to_row": key_to_row,
                "selected_mapping_by_key": selected_mapping_by_key,
                "all_keys": frozenset(str(key) for key in keys),
                "rest_keys": rest_keys,
            }

    if reference_gene_ids is None or reference_symbols is None:
        raise InputError("guide H5MU contains no gene metadata")
    guide_1 = profiles["guide_1"]
    guide_2 = profiles["guide_2"]
    if not guide_2["all_keys"] < guide_1["all_keys"]:
        raise InputError(
            "guide_2 observation keys must be an exact proper subset of guide_1"
        )

    guide_1_rest = guide_1["rest_keys"]
    guide_2_rest = guide_2["rest_keys"]
    if not set(guide_2_rest) < set(guide_1_rest):
        raise InputError(
            "guide_2 Rest keys must be a proper subset of guide_1 Rest keys"
        )
    if len(guide_2_rest) != int(expected["rest_common_atoms"]):
        raise InputError("common Rest key count differs")
    guide_1_only_rest = tuple(sorted(set(guide_1_rest) - set(guide_2_rest)))
    if len(guide_1_only_rest) != int(expected["guide_1_only_rest_atoms"]):
        raise InputError("guide_1-only Rest key count differs")
    if _names_hash(guide_1_rest) != expected["guide_1_rest_keys_sha256"]:
        raise InputError("guide_1 Rest-key hash differs")
    if _names_hash(guide_2_rest) != expected["guide_2_rest_keys_sha256"]:
        raise InputError("guide_2/common Rest-key hash differs")
    for modality, values in (("guide_1", guide_1), ("guide_2", guide_2)):
        selected = [values["selected_mapping_by_key"][key] for key in guide_2_rest]
        for field_index, field in enumerate(spec["selected_rest_mapping_fields"]):
            field_values = [item[field_index] for item in selected]
            if any(
                not value or value == MISSING_CATEGORY for value in field_values
            ) or len(set(field_values)) != len(field_values):
                raise InputError(
                    f"{modality} selected Rest {field} values must be nonempty and unique"
                )
    for key in guide_2_rest:
        if (
            guide_2["selected_mapping_by_key"][key]
            != guide_1["selected_mapping_by_key"][key]
        ):
            raise InputError(f"selected Rest mapping differs for opaque key {key}")

    aligned_rest_rows = {
        modality: np.asarray(
            [profiles[modality]["key_to_row"][key] for key in guide_2_rest],
            dtype=int,
        )
        for modality in GUIDE_MODALITIES
    }
    return {
        "modalities": {
            modality: {
                "rows": profiles[modality]["rows"],
                "condition_counts": profiles[modality]["condition_counts"],
                "categorical_missing_rows": profiles[modality][
                    "categorical_missing_rows"
                ],
                "missing_key_condition_suffix_counts": profiles[modality][
                    "missing_key_condition_suffix_counts"
                ],
                "rest_rows": aligned_rest_rows[modality],
                "rest_row_for_key": {
                    key: int(row)
                    for key, row in zip(guide_2_rest, aligned_rest_rows[modality])
                },
            }
            for modality in GUIDE_MODALITIES
        },
        "gene_ids": tuple(str(value) for value in reference_gene_ids),
        "gene_symbols": tuple(str(value) for value in reference_symbols),
        "rest_atom_keys": guide_2_rest,
        "guide_1_rest_keys": guide_1_rest,
        "guide_1_only_rest_keys": guide_1_only_rest,
        "metadata_hashes": {
            "guide_1_rest_keys_sha256": _names_hash(guide_1_rest),
            "guide_2_common_rest_keys_sha256": _names_hash(guide_2_rest),
            "guide_1_only_rest_keys_sha256": _names_hash(guide_1_only_rest),
            "ordered_gene_ids_sha256": _names_hash(
                tuple(str(value) for value in reference_gene_ids)
            ),
            "ordered_gene_symbols_sha256": _names_hash(
                tuple(str(value) for value in reference_symbols)
            ),
        },
    }


def make_gene_split(
    genes: tuple[str, ...] | list[str], seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return deterministic, disjoint fit/score indices for one seed."""

    if len(genes) < 2 or len(set(genes)) != len(genes):
        raise InputError("gene splits require at least two unique genes")
    rng = np.random.default_rng(int(seed))
    permutation = rng.permutation(len(genes))
    midpoint = len(genes) // 2
    fit_indices = np.sort(permutation[:midpoint])
    score_indices = np.sort(permutation[midpoint:])
    if len(fit_indices) == 0 or len(score_indices) == 0:
        raise InputError("gene split produced an empty side")
    return fit_indices, score_indices


def frozen_predictions(
    models: dict[str, Any], matrix: np.ndarray
) -> dict[str, np.ndarray]:
    """Apply every training-frozen model to one aligned score dictionary."""

    matrix = np.asarray(matrix, dtype=float)
    coefficients = np.asarray(models["cone_coefficients"], dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != len(coefficients):
        raise InputError("frozen model and score dictionary do not align")
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(coefficients)):
        raise InputError("frozen model inputs contain non-finite values")
    best_index = int(models["best_single_index"])
    if not 0 <= best_index < matrix.shape[0]:
        raise InputError("training-selected best-single index is out of range")
    return {
        "cone": coefficients @ matrix,
        "training_common_ray": float(models["common_alpha"]) * np.mean(matrix, axis=0),
        "training_best_single": float(models["best_single_alpha"]) * matrix[best_index],
        "zero": np.zeros(matrix.shape[1], dtype=float),
    }


def _domain_score(
    models: dict[str, Any], matrix: np.ndarray, target: np.ndarray
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    target = np.asarray(target, dtype=float)
    if target.ndim != 1 or not np.all(np.isfinite(target)):
        raise InputError("score target must be a finite vector")
    predictions = frozen_predictions(models, matrix)
    metrics = score_frozen_models(models, matrix, target)
    comparisons = {
        baseline: {
            "cosine_improvement": metrics["cone"]["cosine"]
            - metrics[baseline]["cosine"],
            "normalized_rmse_improvement": metrics[baseline]["normalized_rmse"]
            - metrics["cone"]["normalized_rmse"],
        }
        for baseline in BASELINES
    }
    return {"metrics": metrics, "comparisons": comparisons}, predictions


def _prediction_cosine(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if (
        first.shape != second.shape
        or not np.all(np.isfinite(first))
        or not np.all(np.isfinite(second))
    ):
        raise InputError("prediction-cosine vectors do not align")
    first_norm = float(np.linalg.norm(first))
    second_norm = float(np.linalg.norm(second))
    if first_norm == 0 and second_norm == 0:
        return 1.0
    if first_norm == 0 or second_norm == 0:
        return 0.0
    value = float(first @ second / (first_norm * second_norm))
    return float(np.clip(value, -1.0, 1.0))


def _relative_l2_difference(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if (
        first.shape != second.shape
        or not np.all(np.isfinite(first))
        or not np.all(np.isfinite(second))
    ):
        raise InputError("relative-difference vectors do not align")
    scale = max(float(np.linalg.norm(first)), float(np.linalg.norm(second)))
    if scale == 0:
        return 0.0
    return float(np.linalg.norm(first - second) / scale)


def atom_order_sensitivity(
    canonical_models: dict[str, Any],
    train: np.ndarray,
    within_score: np.ndarray,
    held_score: np.ndarray,
    train_target: np.ndarray,
    atom_keys: tuple[str, ...],
) -> dict[str, Any]:
    """Refit once in reverse atom order and compare mapped predictions."""

    reversed_keys = tuple(reversed(atom_keys))
    reversed_models = fit_training_models(train[::-1], train_target, reversed_keys)
    canonical_coefficients = np.asarray(
        canonical_models["cone_coefficients"], dtype=float
    )
    reversed_coefficients = np.asarray(
        reversed_models["cone_coefficients"], dtype=float
    )[::-1]
    vector_pairs = {
        "coefficients": (canonical_coefficients, reversed_coefficients),
        "fit_prediction": (
            canonical_coefficients @ train,
            reversed_coefficients @ train,
        ),
        "within_score_prediction": (
            canonical_coefficients @ within_score,
            reversed_coefficients @ within_score,
        ),
        "held_score_prediction": (
            canonical_coefficients @ held_score,
            reversed_coefficients @ held_score,
        ),
    }
    return {
        "alternative_atom_order": "reverse_lexicographic",
        "canonical_atom_order_sha256": _names_hash(atom_keys),
        "alternative_atom_order_sha256": _names_hash(reversed_keys),
        "canonical_coefficient_sha256": _vector_hash(canonical_coefficients),
        "alternative_coefficient_mapped_to_canonical_order_sha256": _vector_hash(
            reversed_coefficients
        ),
        "canonical_support": int(np.count_nonzero(canonical_coefficients)),
        "alternative_support": int(np.count_nonzero(reversed_coefficients)),
        "canonical_fit_cosine": float(canonical_models["cone_fit_cosine"]),
        "alternative_fit_cosine": float(reversed_models["cone_fit_cosine"]),
        "canonical_kkt_violation": float(canonical_models["cone_kkt_violation"]),
        "alternative_kkt_violation": float(reversed_models["cone_kkt_violation"]),
        "common_ray_alpha_absolute_difference": abs(
            float(canonical_models["common_alpha"])
            - float(reversed_models["common_alpha"])
        ),
        "best_single_key_match": (
            canonical_models["best_single_key"] == reversed_models["best_single_key"]
        ),
        "best_single_alpha_absolute_difference": abs(
            float(canonical_models["best_single_alpha"])
            - float(reversed_models["best_single_alpha"])
        ),
        "stability": {
            name: {
                "cosine": _prediction_cosine(first, second),
                "relative_l2_difference": _relative_l2_difference(first, second),
            }
            for name, (first, second) in vector_pairs.items()
        },
        "interpretation": (
            "descriptive fixed-solver sensitivity only; one reversed atom order does "
            "not prove coefficient identifiability"
        ),
    }


def score_reciprocal_domains(
    models: dict[str, Any],
    within_test: np.ndarray,
    held_test: np.ndarray,
    score_target: np.ndarray,
) -> dict[str, Any]:
    """Score one target source in both within-slot and reciprocal-slot domains."""

    within, within_predictions = _domain_score(models, within_test, score_target)
    held, held_predictions = _domain_score(models, held_test, score_target)
    deltas = {
        "metrics": {
            model: {
                metric: held["metrics"][model][metric]
                - within["metrics"][model][metric]
                for metric in METRICS
            }
            for model in MODELS
        },
        "comparisons": {
            baseline: {
                metric: held["comparisons"][baseline][metric]
                - within["comparisons"][baseline][metric]
                for metric in ("cosine_improvement", "normalized_rmse_improvement")
            }
            for baseline in BASELINES
        },
    }
    prediction_cosines = {
        model: _prediction_cosine(within_predictions[model], held_predictions[model])
        for model in ("cone", "training_common_ray", "training_best_single")
    }
    return {
        "within_training_guide": within,
        "reciprocal_held_guide": held,
        "reciprocal_minus_within": deltas,
        "prediction_cosine": prediction_cosines,
    }


def run_challenge(
    train: np.ndarray,
    within_test: np.ndarray,
    held_test: np.ndarray,
    train_target: np.ndarray,
    score_target: np.ndarray,
    atom_keys: tuple[str, ...],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Fit once and score both guide-position domains without refitting."""

    models = fit_training_models(train, train_target, atom_keys)
    result = score_reciprocal_domains(models, within_test, held_test, score_target)
    result["training"] = _training_record(models, atom_keys)
    return result, models


def _training_record(
    models: dict[str, Any], atom_keys: tuple[str, ...]
) -> dict[str, Any]:
    best_single_payload = {
        "opaque_key": str(models["best_single_key"]),
        "alpha_hex": float(models["best_single_alpha"]).hex(),
    }
    baseline_hashes = {
        "training_common_ray": _vector_hash(
            np.asarray([models["common_alpha"]], dtype=float)
        ),
        "training_best_single": _canonical_hash(best_single_payload),
        "zero": _vector_hash(np.asarray([0.0], dtype=float)),
    }
    coefficient_hash = _vector_hash(models["cone_coefficients"])
    record = {
        "cone_fit_cosine": float(models["cone_fit_cosine"]),
        "cone_kkt_violation": float(models["cone_kkt_violation"]),
        "cone_support": int(np.count_nonzero(models["cone_coefficients"])),
        "coefficient_sha256": coefficient_hash,
        "baseline_sha256": baseline_hashes,
        "common_alpha": float(models["common_alpha"]),
        "best_single_opaque_key": str(models["best_single_key"]),
        "best_single_alpha": float(models["best_single_alpha"]),
        "atom_key_order_sha256": _names_hash(atom_keys),
    }
    record["training_model_sha256"] = _canonical_hash(
        {
            "coefficient_sha256": coefficient_hash,
            "baseline_sha256": baseline_hashes,
            "atom_key_order_sha256": record["atom_key_order_sha256"],
        }
    )
    return record


def _distribution(values: list[float]) -> dict[str, Any]:
    if not values or not all(np.isfinite(value) for value in values):
        raise InputError("summary values must be nonempty and finite")
    return {
        "median": float(statistics.median(values)),
        "range": [float(min(values)), float(max(values))],
        "fraction_positive": float(statistics.mean(value > 0 for value in values)),
    }


def _domain_summary(rows: list[dict[str, Any]], domain: str) -> dict[str, Any]:
    return {
        "metrics": {
            model: {
                metric: _distribution(
                    [float(row[domain]["metrics"][model][metric]) for row in rows]
                )
                for metric in ("cosine", "normalized_rmse")
            }
            for model in MODELS
        },
        "cone_improvement_over_baselines": {
            baseline: {
                metric: _distribution(
                    [
                        float(row[domain]["comparisons"][baseline][metric])
                        for row in rows
                    ]
                )
                for metric in ("cosine_improvement", "normalized_rmse_improvement")
            }
            for baseline in BASELINES
        },
    }


def _summary_scope(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "challenge_rows": len(rows),
        "unique_fits": len({row["fit_id"] for row in rows}),
        "within_training_guide": _domain_summary(rows, "within_training_guide"),
        "reciprocal_held_guide": _domain_summary(rows, "reciprocal_held_guide"),
        "reciprocal_minus_within": {
            "metrics": {
                model: {
                    metric: _distribution(
                        [
                            float(
                                row["reciprocal_minus_within"]["metrics"][model][metric]
                            )
                            for row in rows
                        ]
                    )
                    for metric in METRICS
                }
                for model in MODELS
            },
            "comparisons": {
                baseline: {
                    metric: _distribution(
                        [
                            float(
                                row["reciprocal_minus_within"]["comparisons"][baseline][
                                    metric
                                ]
                            )
                            for row in rows
                        ]
                    )
                    for metric in (
                        "cosine_improvement",
                        "normalized_rmse_improvement",
                    )
                }
                for baseline in BASELINES
            },
        },
        "prediction_cosine": {
            model: _distribution(
                [float(row["prediction_cosine"][model]) for row in rows]
            )
            for model in ("cone", "training_common_ray", "training_best_single")
        },
    }


def _summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    classes = {
        "same_source_guide_position_transfer": [
            row
            for row in rows
            if row["source_transfer_class"] == "same_source_guide_position_transfer"
        ],
        "joint_guide_position_plus_target_source_transfer": [
            row
            for row in rows
            if row["source_transfer_class"]
            == "joint_guide_position_plus_target_source_transfer"
        ],
    }
    if any(not selected for selected in classes.values()):
        raise InputError("both transfer classes must contain challenge rows")
    return {name: _summary_scope(selected) for name, selected in classes.items()}


def run(config: dict[str, Any], *, verify_hash: bool = True) -> dict[str, Any]:
    validate_config(config)
    # Both complete-file gates intentionally precede HDF5/CSV parsing.
    identities = {
        name: verify_registered_input(spec, verify_hash=verify_hash)
        for name, spec in config["inputs"].items()
    }
    guide_path = _resolve(config["inputs"]["guide_h5mu"]["path"])
    target_path = _resolve(config["inputs"]["target_table"]["path"])
    profile = profile_h5mu(guide_path, config)

    _, target_values = load_target_table(target_path, config["target"])
    lineage = build_target_lineage(
        target_values["log_fc"],
        profile["gene_symbols"],
        orientation_multiplier=config["target"]["orientation_multiplier"],
    )
    genes = tuple(str(value) for value in lineage["shared_screen_genes"])
    expected = config["expected"]
    if len(genes) != int(expected["shared_target_genes"]):
        raise InputError("guide/target shared gene count differs")
    if genes != tuple(sorted(genes)):
        raise InputError("shared target genes must use lexicographic symbol order")
    if not np.isclose(
        float(lineage["between_source_screen_cosine"]),
        float(expected["between_source_target_cosine"]),
        rtol=0,
        atol=1e-12,
    ):
        raise InputError("between-source target cosine differs")

    column_for_gene = {
        gene: index for index, gene in enumerate(profile["gene_symbols"])
    }
    columns = np.asarray([column_for_gene[gene] for gene in genes], dtype=int)
    atom_keys = profile["rest_atom_keys"]
    matrices = {
        modality: load_modality_matrix(
            guide_path,
            modality,
            profile["modalities"][modality]["rest_rows"],
            columns,
            config["h5mu"]["layer"],
        )
        for modality in GUIDE_MODALITIES
    }
    expected_shape = (int(expected["rest_common_atoms"]), len(genes))
    for modality, matrix in matrices.items():
        if matrix.shape != expected_shape or not np.all(np.isfinite(matrix)):
            raise InputError(f"{modality} selected Rest log_fc matrix differs")

    target_directions = lineage["source_directions"]
    rows: list[dict[str, Any]] = []
    representation_sensitivity: dict[str, Any] | None = None
    fit_count = 0
    sources = tuple(config["analysis"]["train_target_sources"])
    sensitivity_spec = config["analysis"]["representation_sensitivity"]
    for train_slot in GUIDE_MODALITIES:
        held_slot = "guide_2" if train_slot == "guide_1" else "guide_1"
        for train_source in sources:
            opposite_source = next(
                source for source in sources if source != train_source
            )
            for seed in config["analysis"]["split_seeds"]:
                fit_indices, score_indices = make_gene_split(genes, int(seed))
                train_matrix = matrices[train_slot][:, fit_indices]
                within_matrix = matrices[train_slot][:, score_indices]
                held_matrix = matrices[held_slot][:, score_indices]
                models = fit_training_models(
                    train_matrix,
                    target_directions[train_source][fit_indices],
                    atom_keys,
                )
                if (
                    train_slot == sensitivity_spec["train_guide_slot"]
                    and train_source == sensitivity_spec["train_target_source"]
                    and int(seed) == sensitivity_spec["seed"]
                ):
                    if representation_sensitivity is not None:
                        raise InputError(
                            "atom-order sensitivity scope matched more than once"
                        )
                    representation_sensitivity = {
                        "train_guide_slot": train_slot,
                        "held_guide_slot": held_slot,
                        "train_target_source": train_source,
                        "seed": int(seed),
                        "fit_gene_sha256": _names_hash(
                            tuple(genes[index] for index in fit_indices)
                        ),
                        "score_gene_sha256": _names_hash(
                            tuple(genes[index] for index in score_indices)
                        ),
                        **atom_order_sensitivity(
                            models,
                            train_matrix,
                            within_matrix,
                            held_matrix,
                            target_directions[train_source][fit_indices],
                            atom_keys,
                        ),
                    }
                fit_count += 1
                training = _training_record(models, atom_keys)
                fit_id = f"{train_slot}|{train_source}|seed={int(seed)}"
                for score_source in (train_source, opposite_source):
                    scored = score_reciprocal_domains(
                        models,
                        within_matrix,
                        held_matrix,
                        target_directions[score_source][score_indices],
                    )
                    transfer_class = (
                        "same_source_guide_position_transfer"
                        if score_source == train_source
                        else "joint_guide_position_plus_target_source_transfer"
                    )
                    rows.append(
                        {
                            "fit_id": fit_id,
                            "train_guide_slot": train_slot,
                            "held_guide_slot": held_slot,
                            "guide_identity_status": config["author_provenance"][
                                "guide_identity_status"
                            ],
                            "train_target_source": train_source,
                            "score_target_source": score_source,
                            "source_transfer_class": transfer_class,
                            "seed": int(seed),
                            "fit_genes": len(fit_indices),
                            "score_genes": len(score_indices),
                            "fit_gene_sha256": _names_hash(
                                tuple(genes[index] for index in fit_indices)
                            ),
                            "score_gene_sha256": _names_hash(
                                tuple(genes[index] for index in score_indices)
                            ),
                            "training": training,
                            "within_training_guide": {
                                "guide_slot": train_slot,
                                **scored["within_training_guide"],
                            },
                            "reciprocal_held_guide": {
                                "guide_slot": held_slot,
                                **scored["reciprocal_held_guide"],
                            },
                            "reciprocal_minus_within": scored[
                                "reciprocal_minus_within"
                            ],
                            "prediction_cosine": scored["prediction_cosine"],
                        }
                    )

    if fit_count != int(expected["fits"]) or len(rows) != int(
        expected["challenge_rows"]
    ):
        raise InputError("observed exhaustive fit/challenge counts differ")
    if representation_sensitivity is None:
        raise InputError("atom-order sensitivity scope was not executed")
    if len({row["fit_id"] for row in rows}) != fit_count:
        raise InputError("fit identities do not map two challenge rows per fit")
    for fit_id in {row["fit_id"] for row in rows}:
        selected = [row for row in rows if row["fit_id"] == fit_id]
        if (
            len(selected) != 2
            or len({row["training"]["training_model_sha256"] for row in selected}) != 1
            or len({_canonical_hash(row["prediction_cosine"]) for row in selected})
            != 1
        ):
            raise InputError(
                "each frozen fit and its target-independent prediction cosines must "
                "be reused for exactly two target scores"
            )

    post_analysis_identities = {
        name: verify_registered_input(spec, verify_hash=verify_hash)
        for name, spec in config["inputs"].items()
    }
    if post_analysis_identities != identities:
        raise InputError("registered input identity changed during analysis")

    report = {
        "schema_version": "1.0.0",
        "generated_on": config["generated_on"],
        "status": "PASS" if verify_hash else "DEVELOPMENT_UNVERIFIED_HASH_SKIPPED",
        "benchmark": config["benchmark"],
        "claim_ceiling": config["claim_ceiling"],
        "source": config["source"],
        "input_verification": identities,
        "source_object_version": {
            key: config["inputs"]["guide_h5mu"][key]
            for key in (
                "url",
                "s3_bucket",
                "s3_key",
                "s3_version_id",
                "etag_header",
                "last_modified_header",
                "last_modified_iso",
            )
        },
        "author_provenance": config["author_provenance"],
        "data_quality": {
            "h5mu_encoding": "MuData with two AnnData modalities",
            "modalities": {
                modality: {
                    "rows": profile["modalities"][modality]["rows"],
                    "condition_counts": profile["modalities"][modality][
                        "condition_counts"
                    ],
                    "categorical_missing_rows": profile["modalities"][modality][
                        "categorical_missing_rows"
                    ],
                    "missing_key_condition_suffix_counts": profile["modalities"][
                        modality
                    ]["missing_key_condition_suffix_counts"],
                }
                for modality in GUIDE_MODALITIES
            },
            "genes": len(profile["gene_symbols"]),
            "shared_source_safe_target_genes": len(genes),
            "between_source_target_cosine": float(
                lineage["between_source_screen_cosine"]
            ),
            "guide_1_rest_atoms": len(profile["guide_1_rest_keys"]),
            "common_rest_atoms": len(atom_keys),
            "guide_1_only_rest_atoms": len(profile["guide_1_only_rest_keys"]),
            "guide_2_keys_exact_subset_of_guide_1": True,
            "selected_rest_target_mapping_exact": True,
            "metadata_hashes": profile["metadata_hashes"],
            "selected_matrix_sha256": {
                modality: _vector_hash(matrix) for modality, matrix in matrices.items()
            },
            "guide_identity_status": config["author_provenance"][
                "guide_identity_status"
            ],
            "eligibility_warning": (
                "The author pipeline uses keep_effective_guides, and the intermediate "
                "for_DE_by_guide.csv is absent. Presence in the released guide-rank "
                "modalities can therefore select on perturbation effectiveness."
            ),
        },
        "protocol": {
            "condition": config["h5mu"]["condition"],
            "layer": config["h5mu"]["layer"],
            "target_gene_universe": (
                "both target sources and guide-H5MU gene symbols; no held-out-source "
                "sign or magnitude filter"
            ),
            "gene_order": config["analysis"]["gene_order"],
            "split_rng": config["analysis"]["rng"],
            "split_seeds": config["analysis"]["split_seeds"],
            "fit_count": fit_count,
            "representation_sensitivity_fits": int(
                expected["representation_sensitivity_fits"]
            ),
            "challenge_rows": len(rows),
            "opaque_key_contract": (
                "target_condition values are matched by exact string identity only; "
                "only the documented terminal Rest/Stim8hr/Stim48hr suffix is parsed "
                "to audit categorical integrity, never to reconstruct an sgRNA ID. "
                "The official dataset card supplies the guide-rank meaning, but exact "
                "IDs are not cross-bound here. Cross-modality target ID/name equality "
                "is required only on the selected common Rest rows; no unsupported "
                "global cross-condition bijection is asserted"
            ),
            "coefficient_application": (
                "fit once on one alphanumeric guide-rank modality, one target source, "
                "and fit genes; apply identical coefficients to both ranked modalities "
                "on disjoint score genes"
            ),
            "baseline_application": (
                "select common-ray scale and best-single opaque key/scale only on the "
                "training modality, training source, and fit genes; freeze them for both "
                "score domains and both target-source scores"
            ),
            "whole_universe_reporting": (
                "all 8,323 common category-labeled Rest atoms are retained; 35 nominal "
                "guide_1 Rest keys with missing categorical metadata are withheld, and "
                "no outcome-ranked top-k subset is selected"
            ),
            "prediction_cosine_zero_policy": (
                "two zero predictions map to 1 and exactly one zero prediction maps to 0"
            ),
            "inference": config["analysis"]["inference"],
        },
        "solver_representation_sensitivity": representation_sensitivity,
        "challenges": rows,
        "summary": _summary(rows),
        "limitations": [
            "The official dataset card defines guide_1 and guide_2 as the first and second alphanumeric sgRNA IDs within each target-condition pair, but those IDs are not embedded or cross-verified in this H5MU benchmark.",
            "guide_2 is a strict subset of guide_1 and the released universe is conditioned on the authors' keep_effective_guides field.",
            "The author pipeline initially required keep_min_cells & keep_effective_guides & keep_total_counts, including >=3 passing replicates per guide-condition and >=5 cells per guide in each condition/sample; target-condition membership came from cond_targets in unreleased for_DE_by_guide.csv, exactly 2 testable guides per target-condition were required, and the final sample mask also required keep_test_genes. The DE model was ~ log10_n_cells + target without a donor term.",
            "The object contains 114 guide_1 and 38 guide_2 rows with code -1 across every required categorical field; documented key suffixes classify these as Rest/Stim8hr/Stim48hr = 35/18/61 and 0/0/38, respectively. They are recorded as missing rather than imputed, and no selected common category-labeled Rest row is missing its target mapping.",
            "Public pseudobulk and guide-library artifacts expose guide IDs, but this report does not hash-bind or reconstruct their exact ranked-modality crosswalk; the absent for_DE_by_guide.csv also prevents auditing the author-selected intermediate from the H5MU alone.",
            "Random gene splits contain correlated coordinates and are descriptive sensitivity challenges, not independent replicates.",
            "The same source study supplies both guide-rank modalities; this is not independent external validation.",
            "Because the fit is underdetermined, transferred coefficients and their held-slot predictions are relative to the frozen SciPy NNLS solver and lexicographic opaque-key order; the fitted cone point need not identify a unique coefficient vector.",
            "No p-values or confidence intervals are emitted; no named-guide, leakage-safe guide-generalization, donor, functional, state-conversion, or intervention claim is supported.",
        ],
    }
    json.dumps(report, allow_nan=False)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--skip-hash", action="store_true", help="development only")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check and args.skip_hash:
        parser.error("--check requires full registered hash verification")

    config_path = _resolve(args.config)
    config_bytes = config_path.read_bytes()
    config = json.loads(config_bytes)
    config_sha256 = hashlib.sha256(config_bytes).hexdigest()
    report = run(config, verify_hash=not args.skip_hash)
    report["config_sha256"] = config_sha256
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    output = _resolve(config["output"])
    if args.check:
        if not output.is_file() or output.read_text(encoding="utf-8") != payload:
            raise SystemExit(f"maintained output differs: {config['output']}")
        print("guide-position transfer output matches")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload, encoding="utf-8")
        print(f"wrote {config['output']}")


if __name__ == "__main__":
    main()
