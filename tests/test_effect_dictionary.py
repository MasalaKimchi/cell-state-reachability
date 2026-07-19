"""Contract tests for the reusable effect-dictionary adapter."""

from __future__ import annotations

import hashlib

import numpy as np
import pytest
from scipy import sparse

import effect_dictionary as ed


def _screen() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    expression = np.array(
        [
            [1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0],
            [4.0, 2.0, 3.0],
            [6.0, 4.0, 3.0],
            [2.0, 5.0, 1.0],
            [2.0, 7.0, 3.0],
        ]
    )
    conditions = np.array(["ctrl", "ctrl", "pert_B", "pert_B", "pert_A", "pert_A"])
    genes = np.array(["G1", "G2", "G3"])
    return expression, conditions, genes


def test_build_has_perturbations_by_genes_orientation_and_signed_effects():
    expression, conditions, genes = _screen()
    dictionary = ed.build_effect_dictionary(
        expression, conditions, control_label="ctrl", gene_names=genes
    )

    assert dictionary["E"].shape == (2, 3)
    assert dictionary["perts"].tolist() == ["pert_A", "pert_B"]
    assert dictionary["genes"].tolist() == genes.tolist()
    assert dictionary["ncells"].tolist() == [2, 2]
    control_mean = expression[conditions == "ctrl"].mean(axis=0)
    expected_a = expression[conditions == "pert_A"].mean(axis=0) - control_mean
    np.testing.assert_allclose(dictionary["E"][0], expected_a)
    assert ed.validate_effect_dictionary(dictionary) == []


def test_sparse_cell_matrix_matches_dense_result():
    expression, conditions, genes = _screen()
    dense = ed.build_effect_dictionary(expression, conditions, gene_names=genes)
    sparse_result = ed.build_effect_dictionary(
        sparse.csr_matrix(expression), conditions, gene_names=genes
    )
    np.testing.assert_allclose(sparse_result["E"], dense["E"])
    np.testing.assert_array_equal(sparse_result["perts"], dense["perts"])


def test_sparse_nonfinite_values_fail_precisely():
    expression, conditions, genes = _screen()
    sparse_expression = sparse.csr_matrix(expression)
    sparse_expression.data[0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        ed.build_effect_dictionary(sparse_expression, conditions, gene_names=genes)


@pytest.mark.parametrize(
    "change, message",
    [
        ({"perts": np.array(["same", "same"])}, "duplicate"),
        ({"genes": np.array(["G1", "G1", "G3"])}, "duplicate"),
        ({"genes": np.array(["G1", "", "G3"])}, "non-empty"),
        ({"E": np.array([[1.0, np.nan, 0.0], [0.0, 1.0, 0.0]])}, "non-finite"),
        ({"E": np.zeros((2, 3), dtype=int)}, "floating"),
        ({"ncells": np.array([2, 0])}, "positive"),
    ],
)
def test_validation_rejects_ambiguous_or_nonfinite_contracts(change, message):
    dictionary = {
        "E": np.zeros((2, 3), dtype=float),
        "perts": np.array(["p1", "p2"]),
        "genes": np.array(["G1", "G2", "G3"]),
        "ncells": np.array([2, 2]),
    }
    dictionary.update(change)
    assert any(message in problem for problem in ed.validate_effect_dictionary(dictionary))


def test_validation_rejects_object_label_arrays():
    dictionary = {
        "E": np.zeros((1, 2), dtype=float),
        "perts": np.array(["p1"], dtype=object),
        "genes": np.array(["G1", "G2"]),
    }
    assert any("non-object string" in p for p in ed.validate_effect_dictionary(dictionary))


def test_build_rejects_missing_control_and_duplicate_gene_labels():
    expression, conditions, genes = _screen()
    with pytest.raises(ValueError, match="no control"):
        ed.build_effect_dictionary(
            expression, conditions, control_label="missing", gene_names=genes
        )
    with pytest.raises(ValueError, match="unique"):
        ed.build_effect_dictionary(
            expression, conditions, gene_names=np.array(["G1", "G1", "G3"])
        )


def test_roundtrip_uses_portable_non_object_arrays(tmp_path):
    expression, conditions, genes = _screen()
    dictionary = ed.build_effect_dictionary(expression, conditions, gene_names=genes)
    path = ed.save_effect_dictionary(tmp_path / "effects.npz", dictionary)
    restored = ed.load_effect_dictionary(path)

    np.testing.assert_allclose(restored["E"], dictionary["E"])
    np.testing.assert_array_equal(restored["perts"], dictionary["perts"])
    np.testing.assert_array_equal(restored["genes"], dictionary["genes"])
    with np.load(path, allow_pickle=False) as archive:
        assert archive["perts"].dtype.kind in {"U", "S"}
        assert archive["genes"].dtype.kind in {"U", "S"}


def test_portable_writer_is_byte_stable_and_uses_fixed_zip_metadata(tmp_path):
    expression, conditions, genes = _screen()
    dictionary = ed.build_effect_dictionary(expression, conditions, gene_names=genes)
    first = ed.save_effect_dictionary(tmp_path / "first.npz", dictionary)
    second = ed.save_effect_dictionary(tmp_path / "second.npz", dictionary)

    assert first.read_bytes() == second.read_bytes()
    assert hashlib.sha256(first.read_bytes()).hexdigest() == hashlib.sha256(
        second.read_bytes()
    ).hexdigest()
    assert len(first.read_bytes()) == 1016
    assert (
        hashlib.sha256(first.read_bytes()).hexdigest()
        == "a35dd58395c07240902430cadd5cfd9b586c0bfece3e252f28b1fac4d26ceeff"
    )
    from zipfile import ZIP_STORED, ZipFile

    with ZipFile(first) as archive:
        assert archive.namelist() == ["E.npy", "perts.npy", "genes.npy", "ncells.npy"]
        assert all(info.date_time == (1980, 1, 1, 0, 0, 0) for info in archive.infolist())
        assert all(info.compress_type == ZIP_STORED for info in archive.infolist())


def test_loader_refuses_pickle_backed_labels(tmp_path):
    path = tmp_path / "unsafe.npz"
    np.savez(
        path,
        E=np.zeros((1, 1), dtype=float),
        perts=np.array(["p1"], dtype=object),
        genes=np.array(["G1"]),
    )
    with pytest.raises(ValueError, match="object/pickled"):
        ed.load_effect_dictionary(path)


def test_loader_rejects_unknown_fields(tmp_path):
    path = tmp_path / "extra.npz"
    np.savez(
        path,
        E=np.zeros((1, 1), dtype=float),
        perts=np.array(["p1"]),
        genes=np.array(["G1"]),
        surprise=np.array([1]),
    )
    with pytest.raises(ValueError, match="unexpected key"):
        ed.load_effect_dictionary(path)


def test_conditions_must_be_real_string_labels():
    expression, _, genes = _screen()
    bad = np.array(["ctrl", "ctrl", "p1", "p1", None, "p2"], dtype=object)
    with pytest.raises(ValueError, match="must contain strings"):
        ed.build_effect_dictionary(expression, bad, gene_names=genes)
