"""Contracts for clean-source, byte-stable effect-dictionary reconstruction."""

from __future__ import annotations

import hashlib
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from effect_dictionary import save_effect_dictionary
from reachability import InputError
from scripts import build_library_coverage_caches as cache_builder


def _identity(path: Path) -> dict[str, object]:
    return {
        "path": path.name,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def _cell_build_spec(path: Path, **overrides: object) -> dict[str, object]:
    spec: dict[str, object] = {
        "builder": "condition_minus_control_mean",
        "source": _identity(path),
        "condition_field": "condition",
        "control_field": "control",
        "control_label": "ctrl",
        "gene_field": "feature_name",
    }
    spec.update(overrides)
    return spec


def _write_h5ad(tmp_path: Path, *, mismatched_control: bool = False) -> Path:
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    expression = np.asarray(
        [
            [1.0, 2.0],
            [3.0, 0.0],
            [4.0, 1.0],
            [6.0, 3.0],
            [1.0, 5.0],
        ],
        dtype=np.float32,
    )
    conditions = np.asarray(["ctrl", "ctrl", "pert_B", "pert_B", "pert_A"])
    controls = conditions == "ctrl"
    if mismatched_control:
        controls[-1] = True
    # Explicit object-backed strings keep this fixture writable under both pandas 2
    # and pandas 3; AnnData 0.10 cannot serialize pandas 3 ArrowStringArray columns.
    obs_index = pd.Index(
        np.asarray(
            [f"cell_{index}" for index in range(expression.shape[0])], dtype=object
        ),
        dtype=object,
    )
    var_index = pd.Index(
        np.asarray(["id1", "id2"], dtype=object), dtype=object, name="gene_id"
    )
    obs = pd.DataFrame(index=obs_index)
    obs["condition"] = pd.Series(conditions, index=obs_index, dtype=object)
    obs["control"] = pd.Series(controls, index=obs_index, dtype=bool)
    var = pd.DataFrame(index=var_index)
    var["feature_name"] = pd.Series(
        np.asarray(["G1", "G2"], dtype=object), index=var_index, dtype=object
    )
    adata = ad.AnnData(X=expression, obs=obs, var=var)
    path = tmp_path / "processed.h5ad"
    adata.write_h5ad(path)
    return path


def _schema(dictionary: dict[str, np.ndarray]) -> dict[str, dict[str, object]]:
    return {
        key: {"shape": list(value.shape), "dtype": str(value.dtype)}
        for key, value in dictionary.items()
    }


def _full_config(
    raw_path: Path,
    expected_cache: Path,
    dictionary: dict[str, np.ndarray],
) -> dict[str, object]:
    raw_spec = _cell_build_spec(raw_path)
    raw_spec["source"] = {
        **raw_spec["source"],
        "url": "https://example.invalid/processed.h5ad",
        "etag_observed": "synthetic-etag",
        "last_modified_observed": "2026-07-19",
        "version_id_observed": "synthetic-version",
        "retrieved_on": "2026-07-19",
    }
    raw_spec["source_matrix_semantics"] = (
        "synthetic processed expression for a unit test; not raw counts"
    )
    raw_spec["gene_axis_note"] = "synthetic unique gene axis"
    cache_identity = _identity(expected_cache)
    cache_identity["path"] = "slice/norman_effects.npz"
    cache_identity["effect_array"] = "E"
    cache_identity["label_array"] = "perts"
    cache_identity["required_arrays"] = _schema(dictionary)
    return {
        "cache_build": {
            "format": "portable_effect_dictionary_v1",
            "writer_contract": "byte-stable portable NPZ",
            "zhu": {
                "builder": "zhu_de_layer",
                "source_config": "configs/source_reconstruction.json",
                "layer": "log_fc",
                "gene_universe": "shared_screen_genes",
                "effect_rows": "registered synthetic row rule",
            },
            "norman": raw_spec,
            "replogle": {
                **raw_spec,
                "gene_field": "gene_id",
                "feature_name_note": "gene_id is the unique synthetic axis",
            },
        },
        "inputs": {"norman": cache_identity},
    }


def test_registered_source_verification_rejects_wrong_hash_before_loading(tmp_path):
    path = tmp_path / "not-an-h5ad"
    path.write_bytes(b"processed expression placeholder")
    spec = _cell_build_spec(path)
    spec["source"] = {**spec["source"], "sha256": "0" * 64}

    with pytest.raises(InputError, match="SHA-256"):
        cache_builder.build_condition_minus_control_dictionary(
            spec, project_root=tmp_path
        )


def test_control_annotation_requires_boolean_and_exact_ctrl_agreement():
    conditions = ["ctrl", "pert_A"]
    with pytest.raises(InputError, match="boolean or integer"):
        cache_builder.validate_control_annotation(
            conditions,
            ["True", "False"],
            control_label="ctrl",
        )
    accepted_conditions, accepted_controls = cache_builder.validate_control_annotation(
        conditions, np.asarray([1, 0], dtype=np.int64), control_label="ctrl"
    )
    assert accepted_conditions.tolist() == conditions
    assert accepted_controls.dtype == np.bool_
    assert accepted_controls.tolist() == [True, False]
    with pytest.raises(InputError, match="disagrees"):
        cache_builder.validate_control_annotation(
            conditions,
            np.asarray([True, True]),
            control_label="ctrl",
        )


def test_processed_h5ad_build_has_exact_axes_effects_and_counts(tmp_path):
    raw_path = _write_h5ad(tmp_path)
    dictionary = cache_builder.build_condition_minus_control_dictionary(
        _cell_build_spec(raw_path), project_root=tmp_path
    )

    assert dictionary["E"].dtype == np.float32
    assert dictionary["perts"].tolist() == ["pert_A", "pert_B"]
    assert dictionary["genes"].tolist() == ["G1", "G2"]
    assert dictionary["ncells"].tolist() == [1, 2]
    assert np.allclose(
        dictionary["E"],
        np.asarray([[-1.0, 4.0], [3.0, 1.0]], dtype=np.float32),
    )


def test_processed_h5ad_rejects_control_semantic_mismatch(tmp_path):
    raw_path = _write_h5ad(tmp_path, mismatched_control=True)
    with pytest.raises(InputError, match="disagrees"):
        cache_builder.build_condition_minus_control_dictionary(
            _cell_build_spec(raw_path), project_root=tmp_path
        )


def test_missing_or_duplicate_configured_gene_field_fails_closed(tmp_path):
    raw_path = _write_h5ad(tmp_path)
    with pytest.raises(InputError, match="missing configured var field/index"):
        cache_builder.build_condition_minus_control_dictionary(
            _cell_build_spec(raw_path, gene_field="absent"), project_root=tmp_path
        )


def test_configured_named_var_index_is_a_valid_unique_gene_axis(tmp_path):
    raw_path = _write_h5ad(tmp_path)
    dictionary = cache_builder.build_condition_minus_control_dictionary(
        _cell_build_spec(raw_path, gene_field="gene_id"), project_root=tmp_path
    )
    assert dictionary["genes"].tolist() == ["id1", "id2"]


def test_check_rebuilds_byte_identically_without_installing_cache(tmp_path):
    raw_path = _write_h5ad(tmp_path)
    build_spec = _cell_build_spec(raw_path)
    dictionary = cache_builder.build_condition_minus_control_dictionary(
        build_spec, project_root=tmp_path
    )
    expected = save_effect_dictionary(tmp_path / "expected.npz", dictionary)
    config = _full_config(raw_path, expected, dictionary)

    result = cache_builder.run(
        config, dataset="norman", project_root=tmp_path, check=True
    )

    assert result[0]["status"] == "PASS"
    assert result[0]["sha256"] == hashlib.sha256(expected.read_bytes()).hexdigest()
    assert not (tmp_path / "slice" / "norman_effects.npz").exists()


def test_output_dir_build_is_deterministic_and_refuses_overwrite(tmp_path):
    raw_path = _write_h5ad(tmp_path)
    build_spec = _cell_build_spec(raw_path)
    dictionary = cache_builder.build_condition_minus_control_dictionary(
        build_spec, project_root=tmp_path
    )
    expected = save_effect_dictionary(tmp_path / "expected.npz", dictionary)
    config = _full_config(raw_path, expected, dictionary)
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"

    cache_builder.run(
        config, dataset="norman", project_root=tmp_path, output_dir=first_dir
    )
    cache_builder.run(
        config, dataset="norman", project_root=tmp_path, output_dir=second_dir
    )
    first = first_dir / "norman_effects.npz"
    second = second_dir / "norman_effects.npz"
    assert first.read_bytes() == second.read_bytes() == expected.read_bytes()

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        cache_builder.run(
            config, dataset="norman", project_root=tmp_path, output_dir=first_dir
        )


def test_builder_config_rejects_unexpected_keys():
    config = {
        "cache_build": {
            "format": "portable_effect_dictionary_v1",
            "writer_contract": "byte-stable portable NPZ",
            "zhu": {
                "builder": "zhu_de_layer",
                "source_config": "configs/source_reconstruction.json",
                "layer": "log_fc",
                "gene_universe": "shared_screen_genes",
                "effect_rows": "registered synthetic row rule",
                "fallback": "unregistered",
            },
            "norman": {},
            "replogle": {},
        }
    }
    with pytest.raises(InputError, match="unexpected=.*fallback"):
        cache_builder.validate_cache_build_section(config)


def test_output_contract_rejects_misdeclared_axes_and_unexpected_keys(tmp_path):
    dictionary = {
        "E": np.zeros((1, 1), dtype=np.float32),
        "perts": np.array(["p1"]),
        "genes": np.array(["g1"]),
    }
    path = save_effect_dictionary(tmp_path / "portable.npz", dictionary)
    spec = {
        **_identity(path),
        "effect_array": "NOT_E",
        "label_array": "NOT_PERTS",
        "required_arrays": _schema(dictionary),
    }
    with pytest.raises(InputError, match="effect_array='E'"):
        cache_builder.assert_cache_contract(path, spec)
    spec["effect_array"] = "E"
    spec["label_array"] = "perts"
    spec["unexpected"] = True
    with pytest.raises(InputError, match="unexpected=.*unexpected"):
        cache_builder.assert_cache_contract(path, spec)


def test_zhu_builder_respects_custom_project_root_and_reconstructs_axes(tmp_path):
    h5py = pytest.importorskip("h5py")
    data_dir = tmp_path / "data"
    config_dir = tmp_path / "configs"
    data_dir.mkdir()
    config_dir.mkdir()
    h5ad_path = data_dir / "tiny.h5ad"
    target_path = data_dir / "target.csv"
    string = h5py.string_dtype("utf-8")
    log_fc = np.asarray([[1.0, 2.0, 3.0], [-1.0, 0.5, 2.0]], dtype=np.float64)
    with h5py.File(h5ad_path, "w") as handle:
        var = handle.create_group("var")
        var.create_dataset(
            "gene_name",
            data=np.asarray(["G1", "G2", "G3"], dtype=object),
            dtype=string,
        )
        obs = handle.create_group("obs")
        condition = obs.create_group("culture_condition")
        condition.create_dataset(
            "categories", data=np.asarray(["Rest"], dtype=object), dtype=string
        )
        condition.create_dataset("codes", data=np.asarray([0, 0], dtype=np.int8))
        perturbation = obs.create_group("target_contrast_gene_name")
        perturbation.create_dataset(
            "categories", data=np.asarray(["P1", "P2"], dtype=object), dtype=string
        )
        perturbation.create_dataset("codes", data=np.asarray([0, 1], dtype=np.int8))
        obs.create_dataset("ontarget_significant", data=np.asarray([True, True]))
        layers = handle.create_group("layers")
        layers.create_dataset("zscore", data=log_fc + 0.25)
        layers.create_dataset("log_fc", data=log_fc)
    with target_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("gene", "contrast", "zscore", "log_fc")
        )
        writer.writeheader()
        for gene, value in (("G1", 1.0), ("G2", -1.0), ("G3", 2.0)):
            writer.writerow(
                {"gene": gene, "contrast": "A", "zscore": value, "log_fc": value}
            )
            writer.writerow(
                {
                    "gene": gene,
                    "contrast": "B",
                    "zscore": value * 2,
                    "log_fc": value * 2,
                }
            )

    source_config = {
        "schema_version": "1.0.0",
        "inputs": {
            "de_stats": {
                "path": "data/tiny.h5ad",
                "bytes": h5ad_path.stat().st_size,
                "sha256": cache_builder.sha256_file(h5ad_path),
            },
            "target_table": {
                "path": "data/target.csv",
                "bytes": target_path.stat().st_size,
                "sha256": cache_builder.sha256_file(target_path),
            },
        },
        "target": {
            "gene_field": "gene",
            "contrast_field": "contrast",
            "sources": {"first": "A", "second": "B"},
            "value_fields": ["zscore", "log_fc"],
            "orientation_multiplier": -1.0,
        },
        "screen": {
            "condition": "Rest",
            "condition_field": "culture_condition",
            "perturbation_field": "target_contrast_gene_name",
            "admission_field": "ontarget_significant",
            "gene_field": "gene_name",
            "layers": ["zscore", "log_fc"],
        },
        "expected_profile": {
            "effect_rows": 2,
            "effect_genes": 3,
            "rest_admitted_atoms": 2,
            "target_union_genes": 3,
            "target_shared_genes": 3,
            "target_sign_concordant_genes": 3,
            "target_screen_overlap_genes": 3,
            "target_shared_screen_genes": 3,
            "target_final_genes": 3,
        },
    }
    source_config_path = config_dir / "source_reconstruction.json"
    source_config_path.write_text(json.dumps(source_config), encoding="utf-8")
    dictionary = cache_builder.build_zhu_dictionary(
        {
            "builder": "zhu_de_layer",
            "source_config": "configs/source_reconstruction.json",
            "layer": "log_fc",
            "gene_universe": "shared_screen_genes",
            "effect_rows": "tiny registered rows",
        },
        project_root=tmp_path,
    )
    np.testing.assert_array_equal(dictionary["E"], log_fc.astype(np.float32))
    assert dictionary["perts"].tolist() == ["P1", "P2"]
    assert dictionary["genes"].tolist() == ["G1", "G2", "G3"]
