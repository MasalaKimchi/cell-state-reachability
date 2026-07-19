"""Tests for the CombiCone command-line interface."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import combicone_cli as cli  # noqa: E402


def _write_synthetic_npz(path, seed=0, n_genes=40):
    """A minimal combicone-substrate .npz with atoms + one emergent double."""
    rng = np.random.default_rng(seed)
    A, B, C = (rng.normal(0, 1, n_genes) for _ in range(3))
    atoms = np.vstack([A, B, C])
    single_genes = np.array(["A", "B", "C"])
    genes = np.array([f"g{i}" for i in range(n_genes)])
    ctrl = np.zeros(n_genes)
    ab = A + B
    ac = A + C + 3.0 * rng.normal(0, 1, n_genes)  # emergent
    conditions = np.array(["A+B", "A+C"], dtype=object)
    means = np.vstack([ab, ac])
    # tiny split-half perturbation so noise SE is finite and small
    means1 = means + rng.normal(0, 0.02, means.shape)
    means2 = means - rng.normal(0, 0.02, means.shape)
    doubles = np.array(["A+B", "A+C"])
    np.savez(
        path, atoms=atoms, single_genes=single_genes, genes=genes, ctrl=ctrl,
        conditions=conditions, means=means, means1=means1, means2=means2, doubles=doubles,
    )


@pytest.fixture()
def npz_path(tmp_path):
    p = tmp_path / "sub.npz"
    _write_synthetic_npz(p)
    return str(p)


def test_ingest_reports_structure(npz_path, capsys):
    rc = cli.main(["ingest", npz_path])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["n_atoms"] == 3
    assert out["from_npz"] is True


def test_triage_ranks_and_writes_csv(npz_path, tmp_path, capsys):
    out_csv = tmp_path / "triage.csv"
    rc = cli.main(["triage", npz_path, "--measured-only", "-o", str(out_csv)])
    assert rc == 0
    text = out_csv.read_text().strip().splitlines()
    assert text[0].startswith("rank,combination,score")
    assert len(text) == 3  # header + 2 combos


def test_certify_emergent_vs_additive(npz_path, tmp_path):
    out_csv = tmp_path / "certs.csv"
    rc = cli.main(["certify", npz_path, "--n-boot", "80", "-o", str(out_csv)])
    assert rc == 0
    import csv
    rows = list(csv.DictReader(out_csv.open()))
    by = {r["combination"]: r for r in rows}
    assert set(by) == {"A+B", "A+C"}
    # the planted-emergent A+C should carry a higher floor ratio than additive A+B
    assert float(by["A+C"]["floor_ratio"]) > float(by["A+B"]["floor_ratio"])


def test_certify_named_subset(npz_path):
    rc = cli.main(["certify", npz_path, "--combo", "A+C", "--n-boot", "50"])
    assert rc == 0


def test_certify_unknown_combo_errors(npz_path):
    with pytest.raises(SystemExit):
        cli.main(["certify", npz_path, "--combo", "X+Y"])


def test_npz_missing_keys_errors(tmp_path):
    bad = tmp_path / "bad.npz"
    np.savez(bad, foo=np.zeros(3))
    with pytest.raises(SystemExit):
        cli.main(["ingest", str(bad)])


def test_triage_csv_pair(tmp_path):
    """CLI reads a dense expression CSV + metadata CSV (no AnnData needed)."""
    rng = np.random.default_rng(1)
    n_genes = 30
    base = {g: rng.normal(0, 1, n_genes) for g in "ABC"}
    rows, labels = [], []
    for _ in range(120):
        rows.append(rng.normal(0, 0.3, n_genes)); labels.append("ctrl")
    for g in "ABC":
        for _ in range(50):
            rows.append(base[g] + rng.normal(0, 0.3, n_genes)); labels.append(f"{g}+ctrl")
    for combo, mu in [("A+B", base["A"] + base["B"]), ("A+C", base["A"] + base["C"])]:
        for _ in range(50):
            rows.append(mu + rng.normal(0, 0.3, n_genes)); labels.append(combo)
    X = np.asarray(rows)
    genes = [f"g{i}" for i in range(n_genes)]
    expr = tmp_path / "expr.csv"
    meta = tmp_path / "meta.csv"
    with expr.open("w") as fh:
        fh.write(",".join(["cell"] + genes) + "\n")
        for i, r in enumerate(X):
            fh.write(f"c{i}," + ",".join(f"{v:.6g}" for v in r) + "\n")
    with meta.open("w") as fh:
        fh.write("cell,condition\n")
        for i, lab in enumerate(labels):
            fh.write(f"c{i},{lab}\n")
    rc = cli.main(["triage", str(expr), "--conditions-csv", str(meta),
                   "--condition-key", "condition", "--measured-only"])
    assert rc == 0
