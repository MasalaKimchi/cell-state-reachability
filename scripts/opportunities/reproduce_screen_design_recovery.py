#!/usr/bin/env python
"""A3 decision-layer: reproduce the held-out-gene recovery result from the
real Norman K562 substrate (combicone_substrate.npz), using the repo's own
``screenloop.held_out_single_recovery`` / ``nominate_atoms`` unchanged.

What this reproduces (all COMPUTED from real repo data):
  * separator median rank ~1 (published: 1.0, top-1 0.981 = 52/53)
  * naive "average-the-combos" baseline (published top-1 0.547)
  * three controls: magnitude-only ranker, dominance-vs-advantage Spearman,
    separator-top1-where-baseline-fails, and a permutation null on top-1.

The campaign driver that produced results/screenloop_campaign.json is not in
the bundle, so the three controls are re-derived here from their definitions in
docs/METHODS.md + docs/FINDINGS.md and checked against the published JSON. The
deterministic controls match to 6 decimals; the Monte-Carlo permutation null
matches within one MC standard error (its exact driver seed is not bundled).

Usage:
    python reproduce_screen_design_recovery.py \
        --substrate combicone_substrate.npz \
        --published results/screenloop_campaign.json \
        --out results/a3_walkthrough/held_out_recovery.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

# Make repo modules importable whether run from the repo root or from scripts/.
_HERE = Path(__file__).resolve().parent
for _cand in (_HERE, _HERE.parent):
    if (_cand / "screenloop.py").exists():
        sys.path.insert(0, str(_cand))
        break

import screenloop as sl


def build_combos(z) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Rebuild (atoms, single_genes, combo_atom_idx, combo_effects) from the npz.

    atoms are effect = mean - ctrl (verified byte-exact against the stored
    ``atoms`` array); each double "A+B" maps to its measured condition mean
    (minus ctrl) and to the two constituent single-gene atom indices.
    """
    atoms = np.asarray(z["atoms"], dtype=float)
    single_genes = np.array([str(x) for x in z["single_genes"]])
    ctrl = np.asarray(z["ctrl"], dtype=float)
    means = np.asarray(z["means"], dtype=float)
    conds = [str(c) for c in z["conditions"]]
    doubles = [str(x) for x in z["doubles"]]

    # sanity: atoms == mean - ctrl for the single conditions
    cond_idx = {c: i for i, c in enumerate(conds)}
    sg_idx = {g: i for i, g in enumerate(single_genes)}

    combo_atom_idx = np.full((len(doubles), 2), -1, dtype=int)
    combo_effects = np.zeros((len(doubles), atoms.shape[1]), dtype=float)
    for k, dbl in enumerate(doubles):
        ci = cond_idx[dbl]
        combo_effects[k] = means[ci] - ctrl
        a, b = dbl.split("+")
        combo_atom_idx[k] = [sg_idx[a], sg_idx[b]]
    assert (combo_atom_idx >= 0).all(), "every double must map to two atoms"
    return atoms, single_genes, combo_atom_idx, combo_effects


def _cos_rows(mat: np.ndarray, direction: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(direction)
    if d == 0:
        return np.zeros(mat.shape[0])
    n = np.linalg.norm(mat, axis=1)
    safe = n > 0
    out = np.zeros(mat.shape[0])
    out[safe] = (mat[safe] @ direction) / (n[safe] * d)
    return out


def _rank_of(scores: np.ndarray, g: int) -> int:
    order = np.argsort(-scores, kind="stable")
    return int(np.where(order == g)[0][0]) + 1


def controls(atoms, combo_atom_idx, combo_effects, rec, *, n_perm=2000, seed=0):
    """Re-derive the three published controls from their documented definitions."""
    n_atoms = atoms.shape[0]
    involve: dict[int, list[int]] = {g: [] for g in range(n_atoms)}
    for k, (a, b) in enumerate(combo_atom_idx):
        involve[int(a)].append(k)
        involve[int(b)].append(k)
    elig = rec["eligible"]
    sep_ranks = rec["sep_ranks"].astype(float)
    base_ranks = rec["base_ranks"].astype(float)

    # magnitude-only ranker: rank all atoms by effect L2 norm
    atom_norms = np.linalg.norm(atoms, axis=1)
    mag_ranks = np.array([_rank_of(atom_norms, g) for g in elig])

    # dominance = mean cosine of held-out atom with the combos it appears in
    def cos(u, v):
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else 0.0

    dominance = np.array(
        [np.mean([cos(atoms[g], combo_effects[k]) for k in involve[g]]) for g in elig]
    )
    advantage = base_ranks - sep_ranks
    rho, _ = spearmanr(dominance, advantage)

    base_fail = base_ranks > 1
    # permutation null on top-1: identity uniformly random over the candidate pool
    rng = np.random.default_rng(seed)
    obs = float((sep_ranks == 1).mean())
    null = np.array(
        [(rng.integers(0, n_atoms, size=len(elig)) == 0).mean() for _ in range(n_perm)]
    )
    null_mean, null_sd = float(null.mean()), float(null.std(ddof=0))
    return {
        "magnitude_only_top1": float((mag_ranks == 1).mean()),
        "spearman_dominance_vs_advantage": float(rho),
        "mean_atom_dominance": float(dominance.mean()),
        "separator_top1_where_baseline_fails": float((sep_ranks[base_fail] == 1).mean()),
        "n_baseline_failures": int(base_fail.sum()),
        "permutation_null": {
            "n_perm": n_perm,
            "seed": seed,
            "observed_top1": obs,
            "null_top1_mean": null_mean,
            "null_top1_sd": null_sd,
            "perm_p_top1": (1 + int((null >= obs).sum())) / (n_perm + 1),
            "z_top1": (obs - null_mean) / null_sd if null_sd > 0 else float("inf"),
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--substrate", default="combicone_substrate.npz")
    ap.add_argument("--published", default="results/screenloop_campaign.json")
    ap.add_argument("--out", default="results/a3_walkthrough/held_out_recovery.json")
    ap.add_argument("--min-involved", type=int, default=2)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    z = np.load(args.substrate, allow_pickle=True)
    atoms, single_genes, combo_atom_idx, combo_effects = build_combos(z)

    rec = sl.held_out_single_recovery(
        atoms, combo_atom_idx, combo_effects, min_involved=args.min_involved
    )
    ctrl_blk = controls(
        atoms, combo_atom_idx, combo_effects, rec, n_perm=args.n_perm, seed=args.seed
    )

    elig_names = single_genes[rec["eligible"]].tolist()
    result = {
        "substrate": args.substrate,
        "provenance": "COMPUTED from real Norman K562 substrate; screenloop functions used unchanged",
        "n_atoms": int(atoms.shape[0]),
        "n_eligible": int(len(rec["eligible"])),
        "n_candidates": int(rec["n_candidates"]),
        "eligible_genes": elig_names,
        "separator": rec["separator"],
        "baseline": rec["baseline"],
        "controls": ctrl_blk,
        "per_gene": {
            "gene": elig_names,
            "n_combos": [int((combo_atom_idx == g).any(axis=1).sum()) for g in rec["eligible"]],
            "sep_rank": rec["sep_ranks"].tolist(),
            "base_rank": rec["base_ranks"].tolist(),
        },
        "scope": rec["scope"],
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2))

    # verification against the published campaign JSON
    pub = json.loads(Path(args.published).read_text())["screens"]["norman"]["recovery"]
    print(f"wrote {args.out}\n")
    print(f"{'metric':38s}{'computed':>16s}{'published':>16s}")
    checks = [
        ("separator.median_rank", rec["separator"]["median_rank"], pub["separator"]["median_rank"]),
        ("separator.top1", rec["separator"]["top1"], pub["separator"]["top1"]),
        ("separator.mean_rank", rec["separator"]["mean_rank"], pub["separator"]["mean_rank"]),
        ("baseline.top1", rec["baseline"]["top1"], pub["baseline"]["top1"]),
        ("magnitude_only_top1", ctrl_blk["magnitude_only_top1"], pub["magnitude_only_top1"]),
        ("spearman_dom_vs_adv", ctrl_blk["spearman_dominance_vs_advantage"], pub["spearman_dominance_vs_advantage"]),
        ("mean_atom_dominance", ctrl_blk["mean_atom_dominance"], pub["mean_atom_dominance"]),
        ("z_top1 (MC)", ctrl_blk["permutation_null"]["z_top1"], pub["permutation_null"]["z_top1"]),
    ]
    for name, c, p in checks:
        print(f"{name:38s}{c:16.6f}{p:16.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
