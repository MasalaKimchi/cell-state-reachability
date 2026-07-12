"""Repeat the flagship held-out-gene check across fixed train/test splits.

The manuscript's canonical held-out cosine (0.448) is computed from one fixed
half-gene split.  This script measures split sensitivity without changing the
fit, target, or scoring rule.  It deliberately does not estimate a permutation
p-value; the denser 60-shuffle headline null remains the inferential reference.

Output
------
results/headline_heldout_split_stability.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import nnls


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "analysis_cache" / "atlas_work" / "inputs.npz"
OUTPUT = ROOT / "results" / "headline_heldout_split_stability.csv"
N_SPLITS = 12


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / denom) if denom else 0.0


def main() -> None:
    bundle = np.load(INPUT, allow_pickle=True)
    E = bundle["E_Rest"].astype(np.float64)
    d = bundle["t_toward_Th1"].astype(np.float64)
    mask = d != 0
    A = E[:, mask].T
    target = d[mask]

    rows: list[dict[str, float | int]] = []
    for seed in range(N_SPLITS):
        rng = np.random.default_rng(seed)
        order = rng.permutation(target.size)
        split = target.size // 2
        train, test = order[:split], order[split:]
        weights, _ = nnls(A[train], target[train])
        rows.append(
            {
                "seed": seed,
                "n_train_genes": train.size,
                "n_test_genes": test.size,
                "train_cosine": cosine(A[train] @ weights, target[train]),
                "heldout_cosine": cosine(A[test] @ weights, target[test]),
                "support_size": int(np.count_nonzero(weights > 1e-8)),
            }
        )
        print(
            f"seed={seed:2d} heldout={rows[-1]['heldout_cosine']:.4f} "
            f"support={rows[-1]['support_size']}"
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    heldout = np.array([float(row["heldout_cosine"]) for row in rows])
    print(
        "summary: "
        f"mean={heldout.mean():.4f} sd={heldout.std(ddof=1):.4f} "
        f"min={heldout.min():.4f} max={heldout.max():.4f}"
    )
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
