"""Tests for the combinatorial-screen ingestion adapter (screen_ingest)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import combicone as cc  # noqa: E402
import screen_ingest as si  # noqa: E402


# --------------------------------------------------------------------------- #
# Label parser
# --------------------------------------------------------------------------- #
def test_parse_single_with_control():
    assert si.parse_condition("AHR+ctrl") == ("AHR",)
    assert si.parse_condition("ctrl+AHR") == ("AHR",)


def test_parse_bare_single():
    assert si.parse_condition("AHR") == ("AHR",)


def test_parse_control_only():
    assert si.parse_condition("ctrl") == ()


def test_parse_double():
    assert si.parse_condition("AHR+KLF1") == ("AHR", "KLF1")


def test_parse_triple():
    assert si.parse_condition("A+B+C") == ("A", "B", "C")


def test_parse_custom_grammar():
    assert si.parse_condition("NT_AHR", control_label="NT", separator="_") == ("AHR",)
    assert si.parse_condition("NT_NT", control_label="NT", separator="_") == ()


def test_parse_rejects_empty():
    with pytest.raises(ValueError):
        si.parse_condition("")
    with pytest.raises(ValueError):
        si.parse_condition("   ")


def test_parse_rejects_empty_arm():
    with pytest.raises(ValueError):
        si.parse_condition("AHR+")


# --------------------------------------------------------------------------- #
# Synthetic screen fixture
# --------------------------------------------------------------------------- #
def _make_screen(seed=0, n_genes=30, n_per=50, emergent_extra=2.0):
    """Build a synthetic cell x gene screen with known singles + one emergent double."""
    rng = np.random.default_rng(seed)
    base = {g: rng.normal(0, 1, n_genes) for g in "ABC"}
    rows, labels = [], []

    def add(label, mu, n=n_per):
        for _ in range(n):
            rows.append(mu + rng.normal(0, 0.3, n_genes))
            labels.append(label)

    add("ctrl", np.zeros(n_genes), 120)
    for g in "ABC":
        add(f"{g}+ctrl", base[g])
    add("A+B", base["A"] + base["B"])  # additive
    add("A+C", base["A"] + base["C"] + emergent_extra * rng.normal(0, 1, n_genes))  # emergent
    X = np.asarray(rows)
    cond = np.asarray(labels)
    genes = np.asarray([f"g{i}" for i in range(n_genes)])
    return X, cond, genes


def test_ingest_arrays_basic():
    X, cond, genes = _make_screen()
    S = si.ingest_screen(expression=X, conditions=cond, gene_names=genes, compute_noise=False)
    assert S.n_atoms == 3
    assert sorted(S.atom_names.tolist()) == ["A", "B", "C"]
    assert S.n_combos == 2
    assert S.coverage() == 1.0
    assert S.n_genes == 30


def test_ingest_combo_resolution():
    X, cond, genes = _make_screen()
    S = si.ingest_screen(expression=X, conditions=cond, gene_names=genes, compute_noise=False)
    by_name = {c.name: c for c in S.combos}
    assert set(by_name) == {"A+B", "A+C"}
    assert by_name["A+B"].constituent_singles == ("A", "B")
    assert by_name["A+B"].has_all_singles


def test_ingest_requires_conditions_with_expression():
    X, cond, genes = _make_screen()
    with pytest.raises(ValueError):
        si.ingest_screen(expression=X)


def test_ingest_no_source_errors():
    with pytest.raises(ValueError):
        si.ingest_screen()


# --------------------------------------------------------------------------- #
# Canonicalization: A+B and B+A collapse; both control arms pool
# --------------------------------------------------------------------------- #
def test_combo_order_collapses():
    X, cond, genes = _make_screen()
    # rename half of A+B cells to B+A; must still be ONE combo
    cond = cond.copy()
    ab = np.where(cond == "A+B")[0]
    cond[ab[: ab.size // 2]] = "B+A"
    S = si.ingest_screen(expression=X, conditions=cond, gene_names=genes, compute_noise=False)
    assert S.n_combos == 2  # A+B (merged from A+B & B+A) and A+C
    names = {c.name for c in S.combos}
    assert "A+B" in names  # canonical sorted form


def test_arm_merge_pools_both_arms():
    X, cond, genes = _make_screen()
    cond = cond.copy()
    # split A's single across both arms
    a = np.where(cond == "A+ctrl")[0]
    cond[a[: a.size // 2]] = "ctrl+A"
    S_merge = si.ingest_screen(
        expression=X, conditions=cond, gene_names=genes,
        arm_handling="merge", compute_noise=False,
    )
    # merge uses all A cells; control_left uses only ctrl+A cells -> different atom
    S_left = si.ingest_screen(
        expression=X, conditions=cond, gene_names=genes,
        arm_handling="control_left", compute_noise=False,
    )
    a_merge = S_merge.atoms[list(S_merge.atom_names).index("A")]
    a_left = S_left.atoms[list(S_left.atom_names).index("A")]
    # both are valid A estimates and highly correlated, but not identical
    assert S_merge.provenance["n_multi_variant_labels"] >= 1
    cos = float(a_merge @ a_left / (np.linalg.norm(a_merge) * np.linalg.norm(a_left)))
    assert cos > 0.8


def test_arm_handling_invalid():
    X, cond, genes = _make_screen()
    with pytest.raises(ValueError):
        si.ingest_screen(expression=X, conditions=cond, gene_names=genes, arm_handling="nonsense")


# --------------------------------------------------------------------------- #
# Split-half noise
# --------------------------------------------------------------------------- #
def test_noise_shape_and_finiteness():
    X, cond, genes = _make_screen(n_per=60)
    S = si.ingest_screen(
        expression=X, conditions=cond, gene_names=genes,
        compute_noise=True, min_cells_per_half=10,
    )
    for c in S.combos:
        assert c.noise_sd is not None
        assert c.noise_sd.shape == (S.n_genes,)
        assert np.all(np.isfinite(c.noise_sd))
        assert np.all(c.noise_sd >= 0)


def test_noise_nan_when_too_few_cells():
    X, cond, genes = _make_screen(n_per=8)  # 8 cells per combo, min_half=10 -> NaN
    S = si.ingest_screen(
        expression=X, conditions=cond, gene_names=genes,
        compute_noise=True, min_cells_per_half=10,
    )
    assert any(np.all(np.isnan(c.noise_sd)) for c in S.combos)


# --------------------------------------------------------------------------- #
# CombiCone handoff
# --------------------------------------------------------------------------- #
def test_triage_ready_runs():
    X, cond, genes = _make_screen()
    S = si.ingest_screen(expression=X, conditions=cond, gene_names=genes, compute_noise=False)
    atoms, names = S.triage_ready()
    tr = cc.triage_combinations(atoms, names, [("A", "B"), ("A", "C")])
    assert list(tr.combos) == [("A", "B"), ("A", "C")]


def test_certify_ready_emergent_vs_additive():
    X, cond, genes = _make_screen(emergent_extra=3.0, n_per=80)
    S = si.ingest_screen(
        expression=X, conditions=cond, gene_names=genes,
        compute_noise=True, min_cells_per_half=15,
    )
    kw_emergent = S.certify_ready("A+C")
    c_em = cc.certify_emergence(
        cone_atoms=kw_emergent["cone_atoms"],
        measured_combo=kw_emergent["measured_combo"],
        noise_sd=kw_emergent["noise_sd"],
        n_boot=100, seed=0,
    )
    # The planted-emergent double should clear a much higher floor than additive A+B
    kw_add = S.certify_ready("A+B")
    c_add = cc.certify_emergence(
        cone_atoms=kw_add["cone_atoms"],
        measured_combo=kw_add["measured_combo"],
        noise_sd=kw_add["noise_sd"],
        n_boot=100, seed=0,
    )
    assert c_em.floor_ratio > c_add.floor_ratio


def test_certify_ready_unknown_combo():
    X, cond, genes = _make_screen()
    S = si.ingest_screen(expression=X, conditions=cond, gene_names=genes, compute_noise=False)
    with pytest.raises(KeyError):
        S.certify_ready("Z+Q")


# --------------------------------------------------------------------------- #
# CSV path
# --------------------------------------------------------------------------- #
def test_ingest_csv_pair(tmp_path):
    X, cond, genes = _make_screen(n_per=40)
    expr_csv = tmp_path / "expr.csv"
    meta_csv = tmp_path / "meta.csv"
    header = ",".join(["cell"] + genes.tolist())
    with expr_csv.open("w") as fh:
        fh.write(header + "\n")
        for i, row in enumerate(X):
            fh.write(f"cell{i}," + ",".join(f"{v:.6g}" for v in row) + "\n")
    with meta_csv.open("w") as fh:
        fh.write("cell,condition\n")
        for i, c in enumerate(cond):
            fh.write(f"cell{i},{c}\n")
    S = si.ingest_screen(str(expr_csv), conditions_csv=str(meta_csv), condition_key="condition", compute_noise=False)
    assert S.n_atoms == 3
    assert S.n_combos == 2
