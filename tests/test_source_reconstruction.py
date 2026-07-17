import csv
import hashlib

import numpy as np
import pytest

from reachability import InputError
from scripts.run_source_reconstruction import (
    build_target_lineage,
    load_target_table,
    profile_h5ad,
    source_transfer,
    verify_input,
)


def test_target_lineage_keeps_source_safe_and_sign_selected_universes_separate():
    values = {
        "first": {"G1": 1.0, "G2": -2.0, "G3": 1.0},
        "second": {"G1": 2.0, "G2": -1.0, "G4": 1.0},
    }
    lineage = build_target_lineage(
        values, ("G1", "G2", "G3"), orientation_multiplier=-1.0
    )
    assert lineage["counts"] == {
        "union": 4,
        "shared": 2,
        "sign_concordant": 2,
        "screen_overlap": 3,
        "shared_screen": 2,
        "final": 2,
    }
    assert lineage["shared_screen_genes"] == ("G1", "G2")
    np.testing.assert_allclose(lineage["source_directions"]["first"], [-1.0, 2.0])
    np.testing.assert_allclose(lineage["merged_direction"], [-1.5, 1.5])


def test_target_table_rejects_duplicate_gene_source_rows(tmp_path):
    path = tmp_path / "target.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("gene", "contrast", "zscore", "log_fc")
        )
        writer.writeheader()
        writer.writerow({"gene": "G1", "contrast": "A", "zscore": 1, "log_fc": 1})
        writer.writerow({"gene": "G1", "contrast": "A", "zscore": 2, "log_fc": 2})
    config = {
        "gene_field": "gene",
        "contrast_field": "contrast",
        "sources": {"first": "A", "second": "B"},
        "value_fields": ["zscore", "log_fc"],
    }
    with pytest.raises(InputError, match="duplicate"):
        load_target_table(path, config)


def test_input_verification_binds_declared_hash_to_bytes(tmp_path, monkeypatch):
    path = tmp_path / "source.bin"
    path.write_bytes(b"source bytes")
    monkeypatch.setattr(
        "scripts.run_source_reconstruction.ROOT", tmp_path
    )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    verified = verify_input(
        "source",
        {"path": "source.bin", "bytes": path.stat().st_size, "sha256": digest},
        verify_hash=True,
    )
    assert verified["hash_verified"] is True
    path.write_bytes(b"changed bytes")
    with pytest.raises(InputError):
        verify_input(
            "source",
            {"path": "source.bin", "bytes": path.stat().st_size, "sha256": digest},
            verify_hash=True,
        )


def test_h5ad_profile_enforces_one_row_grain(tmp_path):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "tiny.h5ad"
    string = h5py.string_dtype("utf-8")
    with h5py.File(path, "w") as handle:
        var = handle.create_group("var")
        var.create_dataset("gene_name", data=np.asarray(["G1", "G2"], dtype=object), dtype=string)
        obs = handle.create_group("obs")
        condition = obs.create_group("culture_condition")
        condition.create_dataset("categories", data=np.asarray(["Rest", "Stim"], dtype=object), dtype=string)
        condition.create_dataset("codes", data=np.asarray([0, 1], dtype=np.int8))
        perturbation = obs.create_group("target_contrast_gene_name")
        perturbation.create_dataset("categories", data=np.asarray(["P1", "P2"], dtype=object), dtype=string)
        perturbation.create_dataset("codes", data=np.asarray([0, 1], dtype=np.int8))
        obs.create_dataset("ontarget_significant", data=np.asarray([True, False]))
        layers = handle.create_group("layers")
        layers.create_dataset("zscore", data=np.eye(2))
        layers.create_dataset("log_fc", data=np.eye(2))
    profile, axes = profile_h5ad(
        path,
        {
            "gene_field": "gene_name",
            "condition_field": "culture_condition",
            "perturbation_field": "target_contrast_gene_name",
            "admission_field": "ontarget_significant",
            "condition": "Rest",
            "layers": ["zscore", "log_fc"],
        },
    )
    assert profile["effect_rows"] == 2
    assert profile["selected_admitted_atoms"] == 1
    np.testing.assert_array_equal(axes["selected_rows"], [0])


def test_source_transfer_reports_cone_and_frozen_baselines():
    effects = np.eye(4)
    directions = {
        "first": np.array([1.0, 0.5, -1.0, -0.5]),
        "second": np.array([0.8, 0.4, -0.8, -0.4]),
    }
    result = source_transfer(effects, directions, [0, 1])
    assert set(result) == {"first_to_second", "second_to_first"}
    for direction in result.values():
        assert len(direction["splits"]) == 2
        assert set(direction["splits"][0]["metrics"]) == {
            "cone",
            "common_response",
            "best_single_atom",
        }
