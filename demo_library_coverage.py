"""End-to-end demo of the library coverage / gap-map layer — no external data.

Builds a synthetic library that spans part of a feature space, a target catalog
that partly demands directions the library cannot produce, and a candidate pool
in which a few perturbations fill the genuine gaps. Then it shows:

  1. coverage of the catalog by the library,
  2. leave-one-out atom redundancy (a planted duplicate is flagged), and
  3. the cheap certificate acquisition ranking recovering the expensive
     realized ranking.

Run:  python demo_library_coverage.py
"""

from __future__ import annotations

import numpy as np

from library_coverage import coverage_report, atom_redundancy, rank_acquisitions


def build_scenario(seed: int = 7, dim: int = 40):
    rng = np.random.default_rng(seed)
    # Library spans coords 0..4 only; produces nothing in coords >= 5.
    lib = np.zeros((5, dim))
    for i in range(5):
        lib[i, i] = 1.0
    lib[:, :5] = np.clip(lib[:, :5] + 0.15 * rng.standard_normal((5, 5)), 0, None)
    dup = lib[2].copy()                       # planted duplicate -> redundant
    effects = np.vstack([lib, dup])           # 6 atoms

    inside = np.abs(rng.standard_normal((15, dim)))
    inside[:, 5:] = 0.0
    outside = np.abs(rng.standard_normal((15, dim)))
    outside[:, 5:] = 0.0
    outside[:, 32] = 2.0 + rng.random(15)     # demand for a coordinate nothing makes
    targets = np.vstack([inside, outside])

    cand = np.abs(rng.standard_normal((10, dim))) * 0.3
    cand[:, 5:] = 0.0
    cand[0] = 0.0
    cand[0, 32] = 1.0                         # exact gap-filler
    cand[1] = 0.0
    cand[1, [30, 31, 33]] = [0.5, 0.5, 0.7]
    return effects, targets, cand


def main() -> None:
    effects, targets, candidates = build_scenario()

    cov = coverage_report(effects, targets)
    print("== Coverage ==")
    print(f"  reachable fraction (strict): {cov.reachable_fraction:.3f}")
    print(f"  mean best cosine:            {cov.mean_cosine:.3f}")
    print(f"  mean unreachable fraction:   {cov.mean_residual_fraction:.3f}")

    red = atom_redundancy(effects, targets)
    print("\n== Redundancy (leave-one-out marginal cosine loss) ==")
    print(f"  per-atom loss: {np.round(red.marginal_cosine_loss, 4).tolist()}")
    print(f"  most redundant atom index: {int(np.argmin(red.marginal_cosine_loss))} "
          f"(atom 5 is a planted duplicate of atom 2)")

    rank = rank_acquisitions(effects, targets, candidates, realized=True)
    from scipy.stats import spearmanr
    rho, _ = spearmanr(rank.certificate_score, rank.realized_gain)
    print("\n== Acquisition (which perturbation to add next) ==")
    print(f"  certificate order (best first): {rank.order.tolist()}")
    print(f"  Spearman(certificate, realized): {rho:.3f}")
    print(f"  certificate top pick: candidate {int(rank.order[0])} "
          f"(candidate 0 is the exact gap-filler)")
    assert int(rank.order[0]) == 0, "certificate should pick the exact gap-filler first"
    assert int(np.argmin(red.marginal_cosine_loss)) in (2, 5), "duplicate should be flagged"
    print("\nOK: certificate recovers the gap-filler; duplicate flagged as redundant.")


if __name__ == "__main__":
    main()
