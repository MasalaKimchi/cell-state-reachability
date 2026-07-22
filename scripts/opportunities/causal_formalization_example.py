"""Reproduce the worked numeric example in docs/causal_formalization.md.

Design-based causal reading of the CombiCone emergence certificate:
  * single-gene effect atoms  =  average treatment effects (ATEs)
  * additive counterfactual set = non-negative conic hull of the ATEs
  * emergence = combination ATE outside that cone (uplift beyond additive policy)

Everything printed is COMPUTED from the real Norman substrate
(combicone_substrate.npz); no value is hard-coded. Run from the repo root:

    python scripts/causal_formalization_example.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)  # make reachability/combicone importable from any CWD

import reachability as rx  # noqa: E402
import combicone as cc  # noqa: E402

SUBSTRATE = os.path.join(REPO, "combicone_substrate.npz")

FLAGSHIP = "SET+CEBPE"
TOP_CERTIFIED = [
    "SET+CEBPE", "IRF1+SET", "MAPK1+PRTG", "CEBPE+RUNX1T1",
    "TBX3+TBX2", "ETS2+IKZF3", "CEBPE+KLF1",
]


def load_substrate(path: str = SUBSTRATE):
    d = np.load(path, allow_pickle=True)
    return {
        "atoms": d["atoms"],                      # (105, 5045) single-gene ATEs
        "single_genes": d["single_genes"],
        "genes": d["genes"],
        "ctrl": d["ctrl"],                        # control mean = E[Y|do(empty)]
        "means": d["means"],                      # (284, 5045) per-condition means
        "means1": d["means1"], "means2": d["means2"],   # split halves
        "conditions": np.array([str(c) for c in d["conditions"]]),
        "doubles": np.array([str(x) for x in d["doubles"]]),
    }


def _idx(names, name):
    return int(np.where(names == name)[0][0])


def two_agent_uplift(S, sub):
    """Best non-negative 2-agent additive counterfactual and its uplift residual."""
    a, b = S.split("+")
    ia, ib = _idx(sub["single_genes"], a), _idx(sub["single_genes"], b)
    ci = _idx(sub["conditions"], S)
    ate_a, ate_b = sub["atoms"][ia], sub["atoms"][ib]
    double = sub["means"][ci] - sub["ctrl"]          # combination ATE
    proj = rx.project_cone(np.vstack([ate_a, ate_b]), double)
    return proj, ate_a, ate_b, double, ci


def main():
    sub = load_substrate()
    print(f"substrate: {sub['atoms'].shape[0]} single-gene ATEs, "
          f"{len(sub['doubles'])} doubles, {len(sub['genes'])} genes\n")

    # --- atoms ARE ATEs: atoms[g] == means[g+ctrl] - ctrl ---
    g0 = sub["single_genes"][0]
    ci0 = _idx(sub["conditions"], f"{g0}+ctrl")
    d0 = np.max(np.abs(sub["atoms"][0] - (sub["means"][ci0] - sub["ctrl"])))
    print(f"[check] atoms[0] == means['{g0}+ctrl'] - ctrl  (max abs diff = {d0:.1e})\n")

    # --- flagship worked example ---
    proj2, ate_a, ate_b, double, ci = two_agent_uplift(FLAGSHIP, sub)
    uplift = proj2.residual
    print(f"=== {FLAGSHIP} ===")
    print(f"  ||ATE_A||={np.linalg.norm(ate_a):.3f}  ||ATE_B||={np.linalg.norm(ate_b):.3f}  "
          f"||double||={np.linalg.norm(double):.3f}")
    print(f"  cos(ATE_A, ATE_B)              = {cc.single_effect_cosine(ate_a, ate_b):+.3f}")
    print(f"  NNLS doses (c_A, c_B)          = ({proj2.coefficients[0]:.3f}, {proj2.coefficients[1]:.3f})")
    print(f"  cos(double, best mixture)      = {proj2.cosine:.3f}")
    print(f"  2-agent uplift fraction        = {proj2.residual_fraction:.3f}")

    projF = rx.project_cone(sub["atoms"], double)
    noise_sd = np.abs(sub["means1"][ci] - sub["means2"][ci]) / 2.0
    cert = cc.certify_emergence(sub["atoms"], double, noise_sd=noise_sd,
                                method="montecarlo", n_boot=200, seed=0)
    print(f"  full-105-cone unreachable      = {projF.residual_fraction:.3f} "
          f"({int(np.sum(projF.coefficients > 1e-9))} atoms used, KKT={projF.kkt_violation:.1e})")
    print(f"  certificate                    = z={cert.z:.1f} floor={cert.floor_ratio:.2f} "
          f"p={cert.p_value:.3g} -> {cert.verdict}")

    sq = uplift ** 2
    cum = np.cumsum(np.sort(sq)[::-1]) / sq.sum()
    print(f"  uplift concentration           = {int(np.searchsorted(cum, 0.50) + 1)} genes carry 50%, "
          f"{int(np.searchsorted(cum, 0.80) + 1)} carry 80% (direction reported, NOT bench-validated)\n")

    # --- cross-pair table ---
    print("=== top certified pairs: 2-agent uplift decomposition ===")
    print(f"{'pair':>16s} {'cos(A,B)':>9s} {'c_A':>6s} {'c_B':>6s} {'cos(dbl,mix)':>12s} {'uplift_frac':>11s}")
    for S in TOP_CERTIFIED:
        p, *_ = two_agent_uplift(S, sub)
        aa, bb = S.split("+")
        cab = cc.single_effect_cosine(sub["atoms"][_idx(sub["single_genes"], aa)],
                                      sub["atoms"][_idx(sub["single_genes"], bb)])
        print(f"{S:>16s} {cab:+9.3f} {p.coefficients[0]:6.3f} {p.coefficients[1]:6.3f} "
              f"{p.cosine:12.3f} {p.residual_fraction:11.3f}")

    # --- faithfulness: 2-agent uplift vs full-cone certificate across all doubles ---
    uf2, ufF = [], []
    for S in sub["doubles"]:
        a, b = S.split("+")
        if a not in sub["single_genes"] or b not in sub["single_genes"] or S not in sub["conditions"]:
            continue
        p, *_ , dbl, _ = two_agent_uplift(S, sub)
        uf2.append(p.residual_fraction)
        ufF.append(rx.project_cone(sub["atoms"], dbl).residual_fraction)
    uf2, ufF = np.array(uf2), np.array(ufF)
    rho, pv = spearmanr(uf2, ufF)
    print(f"\n=== faithfulness across {len(uf2)} doubles ===")
    print(f"  Spearman(2-agent uplift, full-cone residual) = {rho:.3f} (p={pv:.1e})")
    print(f"  full-cone residual <= 2-agent for all pairs  = {bool(np.all(ufF <= uf2 + 1e-9))}")
    print(f"  median 2-agent {np.median(uf2):.3f} -> full-cone {np.median(ufF):.3f} "
          f"(mean reduction {np.mean(uf2 - ufF):.3f})")


if __name__ == "__main__":
    main()
