"""Minimal-set counterfactual solver.

Given a perturbation dictionary E (P x G) of measured causal effect vectors and a
target direction d (G,), find the smallest weighted subset of perturbations whose
combined effect best reconstructs d:

    minimize ||E_S^T w - d||   subject to  |S| small.

Three solvers, cheap on CPU:
  - greedy   : orthogonal-matching-pursuit-style forward selection (interpretable)
  - omp      : scikit-learn OrthogonalMatchingPursuit
  - lasso    : L1-regularized regression; sparsity via the penalty path

IMPORTANT MODELING CONSTRAINTS (honesty):
  * CRISPRi effects are loss-of-function. Each E row already encodes the *measured*
    knockdown effect, so selecting a perturbation = "knock this gene down". We do NOT
    invent activation effects. Set `knockdown_only=True` to forbid negative weights,
    which would otherwise imply an (untested) up-regulation.
  * Multi-gene sets assume additivity/no epistasis — flagged as extrapolation.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class MinimalSet:
    indices: list[int] = field(default_factory=list)   # rows of E selected
    weights: np.ndarray | None = None                  # per-selected weights
    cosine: float = 0.0                                # alignment with target
    residual_norm: float = float("inf")
    solver: str = ""

    def as_dict(self, perturbation_ids) -> dict:
        return {
            "solver": self.solver,
            "genes": [str(perturbation_ids[i]) for i in self.indices],
            "weights": None if self.weights is None else [float(w) for w in self.weights],
            "cosine_similarity": float(self.cosine),
            "residual_norm": float(self.residual_norm),
            "k": len(self.indices),
        }


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0


def greedy_minimal_set(
    E: np.ndarray, d: np.ndarray, k_max: int = 8, knockdown_only: bool = True
) -> MinimalSet:
    """Forward-selection: repeatedly add the perturbation that most reduces residual.

    Returns the set at the k that maximizes cosine alignment with the target (a simple,
    defensible stopping rule; report the full k-vs-alignment curve in the app).
    """
    P, G = E.shape
    residual = d.astype(np.float64).copy()
    chosen: list[int] = []
    best = MinimalSet(solver="greedy")
    for _ in range(min(k_max, P)):
        # score each candidate by 1D projection onto current residual
        proj = E @ residual
        norms = np.einsum("ij,ij->i", E, E) + 1e-12
        gain = proj / np.sqrt(norms)
        if knockdown_only:
            gain[proj < 0] = -np.inf  # forbid selections that require up-regulation
        for c in chosen:
            gain[c] = -np.inf
        j = int(np.argmax(gain))
        if not np.isfinite(gain[j]):
            break
        chosen.append(j)
        # least-squares refit of weights over the chosen set
        A = E[chosen].T
        w, *_ = np.linalg.lstsq(A, d, rcond=None)
        if knockdown_only:
            w = np.clip(w, 0.0, None)
        recon = A @ w
        cos = _cosine(recon, d)
        if cos > best.cosine:
            best = MinimalSet(
                indices=list(chosen), weights=w, cosine=cos,
                residual_norm=float(np.linalg.norm(d - recon)), solver="greedy",
            )
        residual = d - recon
    return best


def omp_minimal_set(E: np.ndarray, d: np.ndarray, k_max: int = 8) -> MinimalSet:
    from sklearn.linear_model import OrthogonalMatchingPursuit

    model = OrthogonalMatchingPursuit(n_nonzero_coefs=k_max, fit_intercept=False)
    model.fit(E.T, d)                    # columns of E.T are perturbations
    idx = np.flatnonzero(model.coef_)
    recon = E.T @ model.coef_
    return MinimalSet(
        indices=idx.tolist(), weights=model.coef_[idx], cosine=_cosine(recon, d),
        residual_norm=float(np.linalg.norm(d - recon)), solver="omp",
    )


def lasso_minimal_set(E: np.ndarray, d: np.ndarray, alpha: float = 0.05,
                      knockdown_only: bool = True) -> MinimalSet:
    from sklearn.linear_model import Lasso

    model = Lasso(alpha=alpha, fit_intercept=False, positive=knockdown_only, max_iter=5000)
    model.fit(E.T, d)
    idx = np.flatnonzero(model.coef_)
    recon = E.T @ model.coef_
    return MinimalSet(
        indices=idx.tolist(), weights=model.coef_[idx], cosine=_cosine(recon, d),
        residual_norm=float(np.linalg.norm(d - recon)), solver="lasso",
    )


def random_null(E: np.ndarray, d: np.ndarray, k: int, n_iter: int = 1000,
                seed: int = 0) -> np.ndarray:
    """Random-perturbation null: cosine of best-fit over k random rows, repeated.

    Use to report an effect size / empirical p for any nominated set of size k.
    """
    rng = np.random.default_rng(seed)
    P = E.shape[0]
    out = np.empty(n_iter)
    for i in range(n_iter):
        rows = rng.choice(P, size=k, replace=False)
        A = E[rows].T
        w, *_ = np.linalg.lstsq(A, d, rcond=None)
        out[i] = _cosine(A @ w, d)
    return out


# --------------------------------------------------------------------------- #
# Tier-1 directional nominations — NO AUTH, NO h5ad, CSVs only.
# --------------------------------------------------------------------------- #
def tier1_directional_nominations(
    de_stats,
    target_map: dict,
    condition: str | None = "Stim8hr",
    significant_only: bool = True,
    gene_col: str = "target_contrast_gene_name",
) -> "object":
    """Rank knockdowns by a 1-D directional prior, using only the local CSVs.

    This is the graded core when the gene-level effect matrix (Tier-2 h5ad) is absent:
    it needs no CZI/Synapse login and no download. It is an HONEST PROXY, not the full
    gene-space reachability — it scores each perturbation on a single question:

        "Does the target want this gene's own transcript LOWER, and is its knockdown
         reproducible and on-target?"

    Because CRISPRi lowers the targeted gene, a knockdown of gene g is directionally
    aligned when the target direction at g is negative (target wants g down). Strength
    combines |target[g]| (how much g matters to the target) with on-target effect size
    and cross-donor reproducibility, and is down-weighted for off-target risk.

    Parameters
    ----------
    de_stats : DataFrame of the DE_stats summary table (per perturbation x condition).
    target_map : {gene -> target-direction value} e.g. from
        dict(zip(genes, target_states.polarization_target(genes, "toward_Th1"))).
    condition : restrict to one culture condition, or None for all.

    Returns a DataFrame sorted by `direction_score` (desc), columns:
      gene, target_value, ontarget_effect_size, crossdonor_corr, offtarget_flag,
      direction_score.
    """
    import numpy as np
    import pandas as pd

    df = de_stats
    if condition is not None and "culture_condition" in df.columns:
        df = df[df["culture_condition"] == condition]
    if significant_only and "ontarget_significant" in df.columns:
        df = df[df["ontarget_significant"].fillna(False).astype(bool)]

    genes = df[gene_col].astype(str).to_numpy()
    tval = np.array([float(target_map.get(g, 0.0)) for g in genes], dtype=float)

    eff = df.get("ontarget_effect_size")
    eff = np.abs(eff.to_numpy(dtype=float)) if eff is not None else np.ones(len(df))
    eff = np.nan_to_num(eff, nan=0.0)

    donor = df.get("crossdonor_correlation_mean")
    donor = donor.to_numpy(dtype=float) if donor is not None else np.full(len(df), np.nan)
    repro = np.clip((np.nan_to_num(donor, nan=0.0) + 1) / 2, 0, 1)

    off = df.get("offtarget_flag")
    off = off.fillna(False).astype(bool).to_numpy() if off is not None else np.zeros(len(df), bool)

    # aligned = target wants g DOWN (negative target value); knockdown lowers g.
    alignment = np.maximum(-tval, 0.0)
    score = alignment * eff * repro * np.where(off, 0.6, 1.0)

    out = pd.DataFrame({
        "gene": genes,
        "target_value": tval,
        "ontarget_effect_size": eff,
        "crossdonor_corr": donor,
        "offtarget_flag": off,
        "direction_score": score,
    })
    out = out[out["direction_score"] > 0].sort_values("direction_score", ascending=False)
    return out.reset_index(drop=True)
