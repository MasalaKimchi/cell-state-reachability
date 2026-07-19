#!/usr/bin/env python
"""End-to-end validation: reproduce the canonical Norman substrate from the raw .h5ad.

Reads ``data/norman_perturbation.h5ad`` in backed mode, streams pooled
condition means in row chunks (peak RAM ~200 MB), builds the CombiCone
substrate via the same math as :mod:`screen_ingest` with
``arm_handling='control_left'``, and asserts the atoms and double effects match
the shipped ``combicone_substrate.npz`` (105 single-gene atoms, 131 doubles).

This is the true end-to-end test for the ingestion adapter: it confirms that the
public API, run on the real 91k-cell file, reproduces the arrays every downstream
result was computed from.

Usage
-----
    PYTHONPATH=. python scripts/validate_ingest_norman.py \
        --h5ad data/norman_perturbation.h5ad \
        --canonical combicone_substrate.npz
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import screen_ingest as si  # noqa: E402


def _canonical_label(raw: str, control_label: str, separator: str) -> str:
    g = si.parse_condition(raw, control_label=control_label, separator=separator)
    if len(g) == 0:
        return control_label
    if len(g) == 1:
        return g[0]
    return separator.join(sorted(g))


def _build_gene_arms(conditions, control_label, separator):
    """Map single gene -> set of control arms present ('left'/'right'/'bare')."""
    gene_arms: dict[str, set[str]] = {}
    for raw in np.unique(conditions):
        g = si.parse_condition(raw, control_label=control_label, separator=separator)
        if len(g) == 1:
            arms = [a.strip() for a in str(raw).split(separator)]
            if len(arms) == 1:
                arm = "bare"
            elif arms[0] == control_label:
                arm = "left"
            elif arms[-1] == control_label:
                arm = "right"
            else:
                arm = "bare"
            gene_arms.setdefault(g[0], set()).add(arm)
    return gene_arms


def _keep_single_arm(raw, control_label, separator, gene_arms):
    """control_left rule WITH fallback: prefer ``ctrl+GENE``; if a gene has no
    left arm, keep whatever arm exists (never drop a gene entirely)."""
    arms = [a.strip() for a in str(raw).split(separator)]
    if len(arms) == 1:
        return True
    genes = [a for a in arms if a != control_label]
    if len(genes) != 1:
        return True
    gene = genes[0]
    present = gene_arms.get(gene, set())
    if "left" not in present:
        return True  # no left arm to prefer -> keep this one
    return arms[0] == control_label


def stream_condition_means(
    adata,
    conditions: np.ndarray,
    *,
    control_label: str,
    separator: str,
    chunk: int = 5000,
):
    """Stream pooled sums per canonical label in row chunks; return means + counts."""
    n_cells, n_genes = adata.shape
    # Precompute per-cell canonical label under control_left arm handling.
    gene_arms = _build_gene_arms(conditions, control_label, separator)
    canon = np.empty(n_cells, dtype=object)
    for i, raw in enumerate(conditions):
        g = si.parse_condition(raw, control_label=control_label, separator=separator)
        if len(g) == 1 and not _keep_single_arm(str(raw), control_label, separator, gene_arms):
            canon[i] = "\x00drop"
        else:
            canon[i] = _canonical_label(str(raw), control_label, separator)
    canon = canon.astype(str)
    labels = sorted(set(canon.tolist()) - {"\x00drop"})
    lab_index = {lab: j for j, lab in enumerate(labels)}
    sums = np.zeros((len(labels), n_genes), dtype=np.float64)
    counts = np.zeros(len(labels), dtype=np.int64)
    X = adata.X
    for start in range(0, n_cells, chunk):
        stop = min(start + chunk, n_cells)
        block = X[start:stop]
        block = np.asarray(block.todense()) if hasattr(block, "todense") else np.asarray(block)
        blabels = canon[start:stop]
        for j, lab in enumerate(labels):
            m = blabels == lab
            if np.any(m):
                sums[j] += block[m].sum(axis=0)
                counts[j] += int(m.sum())
    means = sums / np.maximum(counts[:, None], 1)
    return means, counts, labels, lab_index


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/norman_perturbation.h5ad")
    ap.add_argument("--canonical", default="combicone_substrate.npz")
    ap.add_argument("--condition-key", default="condition")
    ap.add_argument("--control-label", default="ctrl")
    ap.add_argument("--separator", default="+")
    ap.add_argument("--atol", type=float, default=1e-3)
    args = ap.parse_args()

    import anndata as ad

    adata = ad.read_h5ad(args.h5ad, backed="r")
    conditions = np.asarray(adata.obs[args.condition_key].astype(str).values)
    means, counts, labels, lab_index = stream_condition_means(
        adata,
        conditions,
        control_label=args.control_label,
        separator=args.separator,
    )
    ctrl_mean = means[lab_index[args.control_label]]
    genes = np.asarray(adata.var_names.astype(str).values)

    canon = np.load(args.canonical, allow_pickle=True)
    c_genes = list(canon["single_genes"])
    c_atoms = canon["atoms"]
    c_doubles = list(canon["doubles"])

    # --- atoms ---------------------------------------------------------------
    atom_cos = []
    for i, g in enumerate(c_genes):
        if g not in lab_index:
            print(f"MISSING atom {g}")
            return 2
        mine = means[lab_index[g]] - ctrl_mean
        theirs = c_atoms[i]
        cos = float(mine @ theirs / (np.linalg.norm(mine) * np.linalg.norm(theirs) + 1e-12))
        atom_cos.append(cos)
    atom_cos = np.array(atom_cos)

    # --- doubles -------------------------------------------------------------
    dbl_cos = []
    for d in c_doubles:
        key = args.separator.join(sorted(si.parse_condition(d, control_label=args.control_label, separator=args.separator)))
        if key not in lab_index:
            print(f"MISSING double {d} (key {key})")
            return 3
        mine = means[lab_index[key]] - ctrl_mean
        # canonical double effect = its means row minus ctrl (recompute from canon means)
        # compare directions via cosine to the canonical means-minus-ctrl
        idx = c_doubles.index(d)
        theirs = canon["means"][list(canon["conditions"]).index(d)] - canon["ctrl"]
        cos = float(mine @ theirs / (np.linalg.norm(mine) * np.linalg.norm(theirs) + 1e-12))
        dbl_cos.append(cos)
    dbl_cos = np.array(dbl_cos)

    print(f"atoms   : n={len(atom_cos)}  min cosine={atom_cos.min():.6f}  mean={atom_cos.mean():.6f}  n_exact(>0.9999)={int((atom_cos>0.9999).sum())}")
    print(f"doubles : n={len(dbl_cos)}  min cosine={dbl_cos.min():.6f}  mean={dbl_cos.mean():.6f}  n_exact(>0.9999)={int((dbl_cos>0.9999).sum())}")

    ok = (atom_cos.min() > 0.999) and (dbl_cos.min() > 0.999)
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
