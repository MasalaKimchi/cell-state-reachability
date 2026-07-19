"""Fail-closed tests for reciprocal released guide-position transfer."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from reachability import InputError  # noqa: E402
from scripts.run_guide_pair_transfer import (  # noqa: E402
    MISSING_CATEGORY,
    _categorical_allow_missing,
    frozen_predictions,
    load_modality_matrix,
    make_gene_split,
    profile_h5mu,
    run,
    run_challenge,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[1]


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
FORBIDDEN_OUTCOME_FIELDS = (
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
    "n_up_genes",
    "neighboring_gene_KD",
    "ontarget_effect_size",
    "ontarget_significant",
    "single_guide_estimate",
    "target_baseMean",
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


def _names_hash(values: tuple[str, ...] | list[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _categorical(
    group,
    name: str,
    values: list[str],
    *,
    invalid_code: bool = False,
    missing_code: bool = False,
):
    string = h5py.string_dtype("utf-8")
    categories = list(dict.fromkeys(values))
    categorical = group.create_group(name)
    categorical.attrs["encoding-type"] = "categorical"
    categorical.attrs["encoding-version"] = "0.2.0"
    categorical.attrs["ordered"] = False
    categorical.create_dataset(
        "categories", data=np.asarray(categories, dtype=object), dtype=string
    )
    codes = np.asarray([categories.index(value) for value in values], dtype=np.int8)
    if invalid_code:
        codes[0] = len(categories)
    if missing_code:
        codes[-1] = -1
    categorical.create_dataset("codes", data=codes)


def _fixture_rows(*, guide2_only: bool = False):
    guide1 = [
        ("opaque:key/b", "Rest", "ENSG_B", "B"),
        ("guide1-only::rest", "Rest", "ENSG_C", "C"),
        ("opaque:key/a", "Rest", "ENSG_A", "A"),
        ("stim-key", "Stim8hr", "ENSG_D", "D"),
        ("stim48-key", "Stim48hr", "ENSG_E", "E"),
    ]
    guide2 = [
        ("opaque:key/a", "Rest", "ENSG_A", "A"),
        ("stim-key", "Stim8hr", "ENSG_D", "D"),
        ("opaque:key/b", "Rest", "ENSG_B", "B"),
        ("stim48-key", "Stim48hr", "ENSG_E", "E"),
    ]
    if guide2_only:
        guide2[2] = ("guide2-only::rest", "Rest", "ENSG_X", "X")
    return {"guide_1": guide1, "guide_2": guide2}


def _write_tiny_h5mu(
    path: Path,
    *,
    modalities: tuple[str, ...] = ("guide_1", "guide_2"),
    root_encoding: str = "MuData",
    modality_encoding: str = "anndata",
    invalid_code: tuple[str, str] | None = None,
    duplicate_key: str | None = None,
    wrong_layer_shape: tuple[str, str] | None = None,
    nonfinite: tuple[str, str] | None = None,
    gene_axis_mismatch: bool = False,
    mapping_mismatch: bool = False,
    guide2_only: bool = False,
    outcome_offset: float = 0.0,
) -> None:
    string = h5py.string_dtype("utf-8")
    rows_by_modality = _fixture_rows(guide2_only=guide2_only)
    vectors = {
        "opaque:key/a": np.asarray([1.0, 0.0, 1.0, 0.0]),
        "opaque:key/b": np.asarray([0.0, 1.0, 0.0, 1.0]),
        "guide1-only::rest": np.asarray([3.0, 3.0, 3.0, 3.0]),
        "guide2-only::rest": np.asarray([7.0, 7.0, 7.0, 7.0]),
        "stim-key": np.asarray([2.0, 2.0, 2.0, 2.0]),
        "stim48-key": np.asarray([4.0, 4.0, 4.0, 4.0]),
    }
    with h5py.File(path, "w") as handle:
        handle.attrs["encoding-type"] = root_encoding
        handle.attrs["encoding-version"] = "0.1.0"
        handle.attrs["axis"] = 0
        handle.attrs["encoder"] = "mudata"
        handle.attrs["encoder-version"] = "0.3.1"
        mod = handle.create_group("mod")
        mod.attrs["mod-order"] = np.asarray([b"guide_1", b"guide_2"])
        for modality in modalities:
            rows = list(rows_by_modality.get(modality, rows_by_modality["guide_2"]))
            if duplicate_key == modality:
                rows[1] = (rows[0][0], rows[1][1], rows[1][2], rows[1][3])
            group = mod.create_group(modality)
            group.attrs["encoding-type"] = modality_encoding
            group.attrs["encoding-version"] = "0.1.0"
            group.attrs["encoder"] = "mudata"
            group.attrs["encoder-version"] = "0.3.1"
            group.create_dataset("X", data=np.zeros((len(rows), 4), dtype=float))
            for empty in ("obsm", "obsp", "uns", "varm", "varp"):
                group.create_group(empty)

            obs = group.create_group("obs")
            obs.attrs["_index"] = "target_condition"
            obs.attrs["encoding-type"] = "dataframe"
            obs.attrs["encoding-version"] = "0.2.0"
            obs.attrs["column-order"] = np.asarray(
                REQUIRED_OBS_COLUMN_ORDER, dtype=string
            )
            obs.create_dataset(
                "target_condition",
                data=np.asarray([row[0] for row in rows], dtype=object),
                dtype=string,
            )
            target_ids = [row[2] for row in rows]
            target_names = [row[3] for row in rows]
            if mapping_mismatch and modality == "guide_2":
                target_names[0] = "WRONG_NAME"
            categorical_values = {
                "culture_condition": [row[1] for row in rows],
                "target_contrast": target_ids,
                "target_contrast_gene_name": target_names,
                "chunk": ["chunk_0"] * len(rows),
                "n_total_genes_category": ["enough"] * len(rows),
            }
            for field, values in categorical_values.items():
                _categorical(
                    obs,
                    field,
                    values,
                    invalid_code=invalid_code == (modality, field),
                )
            for field_index, field in enumerate(FORBIDDEN_OUTCOME_FIELDS):
                values = (
                    outcome_offset + field_index + np.arange(len(rows), dtype=float)
                )
                obs.create_dataset(field, data=values)

            var = group.create_group("var")
            var.attrs["_index"] = "_index"
            var.attrs["encoding-type"] = "dataframe"
            var.attrs["encoding-version"] = "0.2.0"
            var.attrs["column-order"] = np.asarray(
                ["gene_ids", "gene_name"], dtype=string
            )
            gene_ids = ["g1", "g2", "g3", "g4"]
            gene_names = ["G1", "G2", "G3", "G4"]
            if gene_axis_mismatch and modality == "guide_2":
                gene_names[1], gene_names[2] = gene_names[2], gene_names[1]
            var.create_dataset(
                "_index", data=np.asarray(gene_ids, object), dtype=string
            )
            var.create_dataset(
                "gene_ids", data=np.asarray(gene_ids, object), dtype=string
            )
            var.create_dataset(
                "gene_name", data=np.asarray(gene_names, object), dtype=string
            )

            layers = group.create_group("layers")
            base = np.vstack([vectors[row[0]] for row in rows])
            if modality == "guide_2":
                base = 2.0 * base
            for layer in REQUIRED_LAYERS:
                values = base.copy() if layer == "log_fc" else np.ones_like(base)
                if nonfinite == (modality, layer):
                    values[0, 0] = np.nan
                if wrong_layer_shape == (modality, layer):
                    values = values[:, :-1]
                layers.create_dataset(layer, data=values)


def _tiny_config() -> dict:
    common = ("opaque:key/a", "opaque:key/b")
    guide1_rest = ("guide1-only::rest", *common)
    gene_ids = ("g1", "g2", "g3", "g4")
    gene_names = ("G1", "G2", "G3", "G4")
    return {
        "h5mu": {
            "modalities": ["guide_1", "guide_2"],
            "condition": "Rest",
            "condition_field": "culture_condition",
            "observation_index": "target_condition",
            "gene_id_field": "gene_ids",
            "gene_symbol_field": "gene_name",
            "layer": "log_fc",
            "required_root_attributes": {
                "axis": 0,
                "encoding-type": "MuData",
                "encoding-version": "0.1.0",
                "encoder": "mudata",
                "encoder-version": "0.3.1",
            },
            "required_mod_order": ["guide_1", "guide_2"],
            "required_modality_attributes": {
                "encoding-type": "anndata",
                "encoding-version": "0.1.0",
                "encoder": "mudata",
                "encoder-version": "0.3.1",
            },
            "required_modality_fields": list(REQUIRED_MODALITY_FIELDS),
            "required_obs_attributes": {
                "encoding-type": "dataframe",
                "encoding-version": "0.2.0",
                "_index": "target_condition",
                "column-order": list(REQUIRED_OBS_COLUMN_ORDER),
            },
            "required_var_attributes": {
                "encoding-type": "dataframe",
                "encoding-version": "0.2.0",
                "_index": "_index",
                "column-order": ["gene_ids", "gene_name"],
            },
            "required_categorical_attributes": {
                "encoding-type": "categorical",
                "encoding-version": "0.2.0",
                "ordered": False,
            },
            "required_layers": list(REQUIRED_LAYERS),
            "required_categoricals": list(REQUIRED_CATEGORICALS),
            "required_obs_fields": list(REQUIRED_OBS_FIELDS),
            "required_var_fields": ["_index", "gene_ids", "gene_name"],
            "selected_rest_mapping_fields": [
                "target_contrast",
                "target_contrast_gene_name",
            ],
        },
        "expected": {
            "modalities": 2,
            "genes": 4,
            "modality_rows": {"guide_1": 5, "guide_2": 4},
            "condition_counts": {
                "guide_1": {"Rest": 3, "Stim8hr": 1, "Stim48hr": 1},
                "guide_2": {"Rest": 2, "Stim8hr": 1, "Stim48hr": 1},
            },
            "categorical_missing_rows": {"guide_1": 0, "guide_2": 0},
            "missing_key_condition_suffix_counts": {
                "guide_1": {"Rest": 0, "Stim8hr": 0, "Stim48hr": 0},
                "guide_2": {"Rest": 0, "Stim8hr": 0, "Stim48hr": 0},
            },
            "rest_common_atoms": 2,
            "guide_1_only_rest_atoms": 1,
            "guide_1_rest_keys_sha256": _names_hash(list(guide1_rest)),
            "guide_2_rest_keys_sha256": _names_hash(list(common)),
            "ordered_gene_ids_sha256": _names_hash(list(gene_ids)),
            "ordered_gene_symbols_sha256": _names_hash(list(gene_names)),
        },
    }


def _write_target_table(path: Path, *, ota_values: tuple[float, ...]) -> None:
    genes = ("G1", "G2", "G3", "G4")
    hollbacker = (-1.0, -2.0, -3.0, -4.0)
    lines = ["variable,contrast,log_fc"]
    for source, values in (
        ("Th2_vs_Th1 (Hollbacker 2021)", hollbacker),
        ("Th2_vs_Th1 (Ota 2021)", ota_values),
    ):
        lines.extend(
            f'"{gene}","{source}",{value}' for gene, value in zip(genes, values)
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _integration_config(
    tmp_path: Path,
    *,
    ota_values: tuple[float, ...] = (-1.0, -2.0, -3.0, -4.0),
) -> dict:
    tmp_path.mkdir(parents=True, exist_ok=True)
    h5mu_path = tmp_path / "tiny.h5mu"
    target_path = tmp_path / "target.csv"
    _write_tiny_h5mu(h5mu_path)
    _write_target_table(target_path, ota_values=ota_values)

    config = {
        "schema_version": "1.0.0",
        "generated_on": "2026-07-19",
        "benchmark": "tiny_guide_position_transfer",
        "claim_ceiling": "synthetic test only",
        "source": {"citation": "synthetic"},
        "author_provenance": {
            "guide_identity_status": (
                "OFFICIAL_ALPHANUMERIC_GUIDE_RANKS_IDS_NOT_EMBEDDED_"
                "CROSSWALK_NOT_VERIFIED"
            ),
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
                "guide_1 is the lowest alphanumeric sgRNA ID and guide_2 the second "
                "within each perturbed-gene/culture-condition pair; targets with one "
                "passing guide appear only in guide_1"
            ),
            "h5mu_identity_field_status": (
                "no guide ID or sequence field is embedded in the guide-level H5MU"
            ),
            "public_identity_sources": [
                "GWCD4i.pseudobulk_merged.h5ad obs/guide_id",
                "sgrna_library_metadata.suppl_table.csv sgRNA",
            ],
            "identity_crosswalk_status": (
                "exact ranked-modality-to-sgRNA IDs not reconstructed or hash-verified "
                "in this benchmark"
            ),
            "interpretation": "synthetic test",
        },
        "inputs": {
            "guide_h5mu": {
                "path": str(h5mu_path),
                "bytes": h5mu_path.stat().st_size,
                "sha256": _file_sha256(h5mu_path),
                "url": "https://example.invalid/tiny.h5mu",
                "s3_bucket": "synthetic",
                "s3_key": "tiny.h5mu",
                "s3_version_id": "synthetic-v1",
                "etag_header": "synthetic",
                "last_modified_header": "synthetic",
                "last_modified_iso": "2026-07-19T00:00:00Z",
            },
            "target_table": {
                "path": str(target_path),
                "bytes": target_path.stat().st_size,
                "sha256": _file_sha256(target_path),
            },
        },
        "target": {
            "gene_field": "variable",
            "contrast_field": "contrast",
            "sources": {
                "hollbacker": "Th2_vs_Th1 (Hollbacker 2021)",
                "ota": "Th2_vs_Th1 (Ota 2021)",
            },
            "value_fields": ["log_fc"],
            "orientation_multiplier": -1.0,
        },
        "analysis": {
            "gene_order": "lexicographic_gene_symbol",
            "rng": "numpy.default_rng(seed)",
            "split_seeds": [0, 1, 2],
            "train_guide_slots": ["guide_1", "guide_2"],
            "train_target_sources": ["hollbacker", "ota"],
            "representation_sensitivity": {
                "train_guide_slot": "guide_1",
                "train_target_source": "hollbacker",
                "seed": 0,
                "alternative_atom_order": "reverse_lexicographic",
            },
            "score_target_sources": (
                "same source and the opposite source for every frozen fit"
            ),
            "baselines": ["zero", "training_common_ray", "training_best_single"],
            "inference": (
                "descriptive correlated challenges only; no p-values, confidence "
                "intervals, physical-guide generalization, or donor-population "
                "inference"
            ),
        },
        "output": str(tmp_path / "report.json"),
    }
    config.update(copy.deepcopy(_tiny_config()))
    config["expected"].update(
        {
            "shared_target_genes": 4,
            "between_source_target_cosine": float(
                np.dot((-1.0, -2.0, -3.0, -4.0), ota_values)
                / (
                    np.linalg.norm((-1.0, -2.0, -3.0, -4.0))
                    * np.linalg.norm(ota_values)
                )
            ),
            "fits": 12,
            "representation_sensitivity_fits": 1,
            "challenge_rows": 24,
        }
    )
    return config


def _aligned_matrix(path: Path, profile: dict, modality: str) -> np.ndarray:
    rows = profile["modalities"][modality]["rest_rows"]
    return load_modality_matrix(
        path, modality, rows, np.arange(len(profile["gene_symbols"])), "log_fc"
    )


def test_profile_aligns_opaque_keys_independent_of_physical_row_order(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path)
    profile = profile_h5mu(path, _tiny_config())

    assert profile["rest_atom_keys"] == ("opaque:key/a", "opaque:key/b")
    first = _aligned_matrix(path, profile, "guide_1")
    second = _aligned_matrix(path, profile, "guide_2")
    np.testing.assert_array_equal(first, [[1, 0, 1, 0], [0, 1, 0, 1]])
    np.testing.assert_array_equal(second, 2.0 * first)


def test_profile_uses_guide2_rest_subset_and_excludes_guide1_only_rows(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path)
    profile = profile_h5mu(path, _tiny_config())

    assert profile["guide_1_rest_keys"] == (
        "guide1-only::rest",
        "opaque:key/a",
        "opaque:key/b",
    )
    assert profile["guide_1_only_rest_keys"] == ("guide1-only::rest",)
    assert "guide1-only::rest" not in profile["rest_atom_keys"]


@pytest.mark.parametrize("modalities", [("guide_1",), ("guide_1", "guide_2", "extra")])
def test_profile_rejects_missing_or_extra_modality(tmp_path, modalities):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, modalities=modalities)
    with pytest.raises(InputError, match="modalit"):
        profile_h5mu(path, _tiny_config())


@pytest.mark.parametrize(
    ("root_encoding", "modality_encoding", "message"),
    [
        ("not-mudata", "anndata", "MuData|encoding"),
        ("MuData", "not-anndata", "AnnData|encoding"),
    ],
)
def test_profile_rejects_wrong_h5_encoding(
    tmp_path, root_encoding, modality_encoding, message
):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(
        path, root_encoding=root_encoding, modality_encoding=modality_encoding
    )
    with pytest.raises(InputError, match=message):
        profile_h5mu(path, _tiny_config())


def test_profile_rejects_invalid_categorical_code(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, invalid_code=("guide_2", "culture_condition"))
    with pytest.raises(InputError, match="categorical|code"):
        profile_h5mu(path, _tiny_config())


def test_categorical_decoder_preserves_code_minus_one_as_missing(tmp_path):
    path = tmp_path / "missing.h5"
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        _categorical(
            obs,
            "culture_condition",
            ["Rest", "Stim48hr"],
            missing_code=True,
        )
        values, missing = _categorical_allow_missing(obs, "culture_condition")
    assert values.tolist() == ["Rest", MISSING_CATEGORY]
    assert missing.tolist() == [False, True]


def test_profile_rejects_duplicate_observation_keys(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, duplicate_key="guide_1")
    with pytest.raises(InputError, match="unique"):
        profile_h5mu(path, _tiny_config())


def test_profile_rejects_required_layer_shape_drift(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, wrong_layer_shape=("guide_1", "baseMean"))
    with pytest.raises(InputError, match="shape"):
        profile_h5mu(path, _tiny_config())


def test_profile_rejects_nonfinite_selected_log_fc(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, nonfinite=("guide_2", "log_fc"))
    profile = profile_h5mu(path, _tiny_config())
    with pytest.raises(InputError, match="non-finite"):
        _aligned_matrix(path, profile, "guide_2")


def test_profile_does_not_require_unused_layers_to_be_finite(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, nonfinite=("guide_2", "p_value"))
    profile = profile_h5mu(path, _tiny_config())
    assert profile["rest_atom_keys"] == ("opaque:key/a", "opaque:key/b")


def test_profile_rejects_modality_gene_axis_mismatch(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, gene_axis_mismatch=True)
    with pytest.raises(InputError, match="gene.*(axis|differ|match)"):
        profile_h5mu(path, _tiny_config())


def test_profile_rejects_common_key_target_mapping_mismatch(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, mapping_mismatch=True)
    with pytest.raises(InputError, match="mapping|target"):
        profile_h5mu(path, _tiny_config())


def test_profile_rejects_guide2_only_row(tmp_path):
    path = tmp_path / "tiny.h5mu"
    _write_tiny_h5mu(path, guide2_only=True)
    with pytest.raises(InputError, match="subset|guide_2"):
        profile_h5mu(path, _tiny_config())


def test_forbidden_outcome_metadata_cannot_change_universe_or_matrices(tmp_path):
    first_path = tmp_path / "first.h5mu"
    mutated_path = tmp_path / "mutated.h5mu"
    _write_tiny_h5mu(first_path, outcome_offset=0.0)
    _write_tiny_h5mu(mutated_path, outcome_offset=1_000_000.0)

    first = profile_h5mu(first_path, _tiny_config())
    mutated = profile_h5mu(mutated_path, _tiny_config())
    assert first["rest_atom_keys"] == mutated["rest_atom_keys"]
    assert first["gene_ids"] == mutated["gene_ids"]
    assert first["gene_symbols"] == mutated["gene_symbols"]
    for modality in ("guide_1", "guide_2"):
        np.testing.assert_array_equal(
            _aligned_matrix(first_path, first, modality),
            _aligned_matrix(mutated_path, mutated, modality),
        )


def test_training_models_and_baselines_ignore_held_guide_and_score_target():
    train = np.asarray([[2.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    train_target = np.asarray([1.0, 0.25])
    keys = ("A", "B", "C")
    within = np.asarray([[2.0, 1.0], [1.0, 2.0], [3.0, 3.0]])
    held_a = np.asarray([[3.0, 2.0], [7.0, -4.0], [1.0, 9.0]])
    held_b = np.asarray([[-30.0, 20.0], [70.0, -40.0], [10.0, 90.0]])
    first, first_models = run_challenge(
        train,
        within,
        held_a,
        train_target,
        np.asarray([0.5, -0.25]),
        keys,
    )
    second, second_models = run_challenge(
        train,
        within,
        held_b,
        train_target,
        np.asarray([-99.0, 101.0]),
        keys,
    )

    assert first_models["best_single_key"] == second_models["best_single_key"] == "A"
    assert first_models["best_single_index"] == second_models["best_single_index"] == 0
    assert first_models["common_alpha"] == second_models["common_alpha"]
    assert first_models["best_single_alpha"] == second_models["best_single_alpha"]
    np.testing.assert_array_equal(
        first_models["cone_coefficients"], second_models["cone_coefficients"]
    )
    assert first["reciprocal_held_guide"] != second["reciprocal_held_guide"]


def test_run_challenge_reuses_one_frozen_fit_for_within_and_held_domains():
    train = np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    within = np.asarray([[2.0, 0.0], [0.0, 2.0], [2.0, 2.0]])
    held = np.asarray([[3.0, 0.0], [0.0, 3.0], [3.0, 3.0]])
    target = np.asarray([1.0, 0.5])
    result, models = run_challenge(train, within, held, target, target, ("A", "B", "C"))

    assert models["best_single_key"] == "C"
    assert models["common_alpha"] == pytest.approx(1.125)
    within_predictions = frozen_predictions(models, within)
    held_predictions = frozen_predictions(models, held)
    assert set(within_predictions) == set(held_predictions)
    assert set(result) >= {
        "within_training_guide",
        "reciprocal_held_guide",
        "reciprocal_minus_within",
        "prediction_cosine",
    }


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_gene_splits_are_deterministic_disjoint_and_exhaustive(seed):
    genes = tuple(f"G{index:02d}" for index in range(11))
    first_fit, first_score = make_gene_split(genes, seed)
    second_fit, second_score = make_gene_split(genes, seed)

    np.testing.assert_array_equal(first_fit, second_fit)
    np.testing.assert_array_equal(first_score, second_score)
    assert set(first_fit).isdisjoint(set(first_score))
    assert set(first_fit) | set(first_score) == set(range(len(genes)))
    assert len(first_fit) == len(genes) // 2
    assert len(first_score) == len(genes) - len(first_fit)


@pytest.mark.parametrize(
    ("section", "field", "mutated"),
    [
        ("analysis", "gene_order", "input_order"),
        ("analysis", "rng", "legacy RNG"),
        ("analysis", "split_seeds", [0.9, 1, 2]),
        ("analysis", "train_target_sources", ["ota", "hollbacker"]),
        ("analysis", "score_target_sources", "same source only"),
        ("analysis", "baselines", []),
        ("analysis", "representation_sensitivity", {}),
        ("target", "orientation_multiplier", 1.0),
    ],
)
def test_validate_config_rejects_protocol_semantic_drift(
    tmp_path, section, field, mutated
):
    config = _integration_config(tmp_path)
    config[section][field] = mutated
    with pytest.raises(InputError, match="contract|seeds|frozen"):
        validate_config(config)


def test_validate_config_rejects_removed_obs_column_order(tmp_path):
    config = _integration_config(tmp_path)
    del config["h5mu"]["required_obs_attributes"]["column-order"]
    with pytest.raises(InputError, match="observation-attribute"):
        validate_config(config)


def test_validate_config_rejects_config_schema_drift(tmp_path):
    config = _integration_config(tmp_path)
    config["schema_version"] = "999.0.0"
    with pytest.raises(InputError, match="config schema"):
        validate_config(config)


def test_report_has_24_unique_rows_12_fits_and_paired_hash_reuse(tmp_path):
    report = run(_integration_config(tmp_path), verify_hash=True)
    rows = report["challenges"]

    assert len(rows) == 24
    row_keys = {
        (
            row["train_guide_slot"],
            row["held_guide_slot"],
            row["train_target_source"],
            row["score_target_source"],
            row["seed"],
        )
        for row in rows
    }
    assert len(row_keys) == 24
    assert (
        sum(row["train_target_source"] == row["score_target_source"] for row in rows)
        == 12
    )
    assert (
        sum(row["train_target_source"] != row["score_target_source"] for row in rows)
        == 12
    )

    fit_ids = sorted({row["fit_id"] for row in rows})
    assert len(fit_ids) == 12
    for fit_id in fit_ids:
        pair = [row for row in rows if row["fit_id"] == fit_id]
        assert len(pair) == 2
        assert {row["score_target_source"] for row in pair} == {
            row["train_target_source"] for row in pair
        } | ({"hollbacker", "ota"} - {pair[0]["train_target_source"]})
        assert (
            pair[0]["training"]["coefficient_sha256"]
            == pair[1]["training"]["coefficient_sha256"]
        )
        assert (
            pair[0]["training"]["baseline_sha256"]
            == pair[1]["training"]["baseline_sha256"]
        )
        assert (
            pair[0]["training"]["training_model_sha256"]
            == pair[1]["training"]["training_model_sha256"]
        )

    sensitivity = report["solver_representation_sensitivity"]
    assert sensitivity["train_guide_slot"] == "guide_1"
    assert sensitivity["held_guide_slot"] == "guide_2"
    assert sensitivity["train_target_source"] == "hollbacker"
    assert sensitivity["seed"] == 0
    assert sensitivity["canonical_atom_order_sha256"] == _names_hash(
        ["opaque:key/a", "opaque:key/b"]
    )
    assert sensitivity["alternative_atom_order_sha256"] == _names_hash(
        ["opaque:key/b", "opaque:key/a"]
    )
    assert set(sensitivity["stability"]) == {
        "coefficients",
        "fit_prediction",
        "within_score_prediction",
        "held_score_prediction",
    }
    for values in sensitivity["stability"].values():
        assert -1.0 <= values["cosine"] <= 1.0
        assert values["relative_l2_difference"] >= 0.0


def test_held_source_mutation_cannot_change_universe_splits_or_other_source_fit(
    tmp_path,
):
    original = run(_integration_config(tmp_path / "original"), verify_hash=True)
    mutated = run(
        _integration_config(tmp_path / "mutated", ota_values=(9.0, -80.0, -0.7, -60.0)),
        verify_hash=True,
    )

    assert original["data_quality"]["shared_source_safe_target_genes"] == 4
    assert mutated["data_quality"]["shared_source_safe_target_genes"] == 4
    original_splits = {
        (row["seed"], row["fit_gene_sha256"], row["score_gene_sha256"])
        for row in original["challenges"]
    }
    mutated_splits = {
        (row["seed"], row["fit_gene_sha256"], row["score_gene_sha256"])
        for row in mutated["challenges"]
    }
    assert original_splits == mutated_splits

    def unchanged_source_fits(report):
        return {
            row["fit_id"]: row["training"]["training_model_sha256"]
            for row in report["challenges"]
            if row["train_target_source"] == "hollbacker"
        }

    assert unchanged_source_fits(original) == unchanged_source_fits(mutated)


def test_report_is_deterministic_and_records_guide_identity_ceiling(tmp_path):
    config = _integration_config(tmp_path)
    first = run(config, verify_hash=True)
    second = run(config, verify_hash=True)

    def render(value):
        return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"

    assert render(first) == render(second)
    assert (
        first["data_quality"]["guide_identity_status"]
        == "OFFICIAL_ALPHANUMERIC_GUIDE_RANKS_IDS_NOT_EMBEDDED_"
        "CROSSWALK_NOT_VERIFIED"
    )


def test_cli_write_then_exact_check_is_deterministic(tmp_path):
    config = _integration_config(tmp_path)
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    command = [
        sys.executable,
        str(ROOT / "scripts" / "run_guide_pair_transfer.py"),
        "--config",
        str(config_path),
    ]

    written = subprocess.run(command, capture_output=True, text=True, check=False)
    assert written.returncode == 0, written.stderr
    output = Path(config["output"])
    first_bytes = output.read_bytes()
    checked = subprocess.run(
        [*command, "--check"], capture_output=True, text=True, check=False
    )
    assert checked.returncode == 0, checked.stderr
    assert output.read_bytes() == first_bytes


def test_official_guide_rank_mapping_is_scoped_to_unverified_h5mu_crosswalk():
    config = json.loads(
        (ROOT / "configs" / "guide_pair_transfer.json").read_text(encoding="utf-8")
    )
    assert "guide_id_field" not in config["h5mu"]
    assert config["h5mu"]["modalities"] == ["guide_1", "guide_2"]
    assert (
        config["author_provenance"]["guide_identity_status"]
        == "OFFICIAL_ALPHANUMERIC_GUIDE_RANKS_IDS_NOT_EMBEDDED_"
        "CROSSWALK_NOT_VERIFIED"
    )
    assert "lowest alphanumeric sgRNA ID" in config["author_provenance"][
        "official_modality_mapping"
    ]
    assert config["author_provenance"]["public_identity_sources"] == [
        "GWCD4i.pseudobulk_merged.h5ad obs/guide_id",
        "sgrna_library_metadata.suppl_table.csv sgRNA",
    ]
    assert "not reconstructed or hash-verified" in config["author_provenance"][
        "identity_crosswalk_status"
    ]
    assert config["source"]["dataset_card"].startswith(
        "https://virtualcellmodels.cziscience.com/dataset/"
    )
