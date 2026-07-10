"""
reinforcement_analyses.py
=========================
Helper primitives for the reinforcement analyses that strengthen the
cell-state reachability manuscript (notebook 06). Every function is a thin,
documented wrapper over the *validated* primitives in ``reachability.py`` --
this module adds no new geometry, only the held-out comparison protocols and
the reliability-band bookkeeping that turn the manuscript's limitations
(L1, L2, L4, L5) into runnable, quantitative tests.

Design contract
----------------
* All geometry comes from ``reachability.py`` (reachability, signed_reachability,
  reachability_spectrum, additivity_risk, held_out_gene_validation). We never
  re-derive the NNLS fit, the certificate, or the saturation law here.
* The atlas HVG convention is the *nonzero support of the target vector*
  (``target_mask``) -- this reproduces the published atlas cells exactly
  (verified to 0.00e+00 on every geometry field of the headline card).
* The additivity-risk reference norm ``s`` is the median single-effect norm on
  the FULL readout axis (not the masked subspace); this reproduces the
  published ``atlas_additivity_risk`` values (e.g. Th1/Rest risk 0.082).

Author: reinforcement-analyses notebook build.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Optional, Sequence
from scipy.optimize import nnls
from scipy.linalg import lstsq

import reachability as rx   # the validated core (must be importable on sys.path)

__all__ = [
    "target_mask", "cosine",
    "ablation_cell", "AblationResult",
    "ceiling_cell", "CeilingResult",
    "reliability_curve", "ReliabilityPoint",
    "held_out_modality_test", "ModalityTestResult",
]


# ----------------------------------------------------------------------------- helpers
def target_mask(d_full: np.ndarray) -> np.ndarray:
    """Atlas HVG mask = the nonzero support of the target vector (per run_atlas.py).

    The atlas restricts the readout axis to the genes the target actually moves;
    for ``toward_Th1`` this reduces the 10,282-gene axis to 6,188 genes and
    reproduces the cached headline card exactly.
    """
    m = np.zeros(d_full.shape[0], dtype=bool)
    m[np.where(d_full != 0)[0]] = True
    return m


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na and nb else 0.0


# ----------------------------------------------------------------------------- L4: constraint ablation
@dataclass
class AblationResult:
    """One atlas cell's held-out comparison of NNLS vs unconstrained vs nearest-single."""
    nnls_heldout: float
    nnls_support: int
    unconstrained_heldout: float
    unconstrained_neg_weights: int
    nearest_single_heldout: float
    nearest_single_gene: str
    n_generators: int
    cosine_cost_of_constraint: float   # unconstrained - nnls (>=0 means the cone costs accuracy)


def ablation_cell(E_kd: np.ndarray, gene_names: Sequence[str], d_full: np.ndarray,
                  mask: np.ndarray, *, seed: int = 0) -> AblationResult:
    """L4. Fit three recipes on a gene-half split and score held-out cosine on the other half.

    The three recipes:
      (i)   NNLS -- the non-negative cone (the method).
      (ii)  unconstrained least squares -- signed weights allowed (drops the constraint).
      (iii) nearest single effect -- the best single knockdown, scaled by its LS coefficient.

    The split mirrors ``held_out_gene_validation`` so the NNLS number is comparable to the
    manuscript's honest held-out cosine. Reports the accuracy cost of non-negativity and the
    number of *negative* weights the unconstrained fit uses (biologically unrealizable:
    a knockdown applied in negative amount is an activation the assay never measured).
    """
    d = d_full[mask]
    A = E_kd[:, mask].T                       # (G_mask, P): columns are cone generators
    G = d.shape[0]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(G)
    h1, h2 = perm[: G // 2], perm[G // 2:]
    A1, A2, d1, d2 = A[h1], A[h2], d[h1], d[h2]

    # (i) NNLS
    w_nn, _ = nnls(A1, d1)
    nn_ho = cosine(A2 @ w_nn, d2)
    nn_supp = int((w_nn > 1e-8).sum())

    # (ii) unconstrained least squares
    w_ls, _, _, _ = lstsq(A1, d1, lapack_driver="gelsy")
    ls_ho = cosine(A2 @ w_ls, d2)
    ls_neg = int((w_ls < -1e-8).sum())

    # (iii) nearest single effect
    cn = np.linalg.norm(A1, axis=0)
    cn = np.where(cn > 0, cn, 1.0)
    aln = (A1.T @ d1) / (cn * np.linalg.norm(d1))
    j = int(np.argmax(aln))
    a1 = A1[:, j]
    coef = float(a1 @ d1 / (a1 @ a1)) if a1 @ a1 else 0.0
    ns_ho = cosine(coef * A2[:, j], d2)

    return AblationResult(
        nnls_heldout=nn_ho, nnls_support=nn_supp,
        unconstrained_heldout=ls_ho, unconstrained_neg_weights=ls_neg,
        nearest_single_heldout=ns_ho, nearest_single_gene=str(gene_names[j]),
        n_generators=A.shape[1], cosine_cost_of_constraint=ls_ho - nn_ho,
    )


# ----------------------------------------------------------------------------- L5: reachable-cosine ceiling
@dataclass
class CeilingResult:
    lof_fraction: float
    gof_fraction: float
    neither_fraction: float
    theoretical_lof_ceiling: float   # sqrt(lof_fraction) = best any knockdown-only method could reach
    insample_lof_cosine: float       # reachable_cosine (all generators, seen genes)
    signed_ceiling: float            # signed_cosine (ceiling if activation also allowed)
    achieved_heldout: float          # the honest generalization number
    frac_of_theoretical_ceiling: float


def ceiling_cell(lof_fraction: float, gof_fraction: float, neither_fraction: float,
                 reachable_cosine: float, signed_cosine: float, heldout_cosine: float) -> CeilingResult:
    """L5. Reframe a modest cosine as a fraction of the achievable ceiling.

    The orthogonal signed decomposition gives the theoretical knockdown-only ceiling as
    ``sqrt(lof_fraction)`` -- and this EQUALS the in-sample cone cosine (verified across the
    atlas to <1e-3), because the cone fit achieves the geometric maximum on the LOF-reachable
    subspace. The gap from ``heldout_cosine`` up to this ceiling is the honest generalization
    cost; the gap from the ceiling up to 1.0 is biology (the GOF-locked share the assay can't
    reach by knockdown). Numbers are read from the canonical atlas cells, not recomputed.
    """
    ceil_th = float(np.sqrt(max(lof_fraction, 0.0)))
    return CeilingResult(
        lof_fraction=lof_fraction, gof_fraction=gof_fraction, neither_fraction=neither_fraction,
        theoretical_lof_ceiling=ceil_th, insample_lof_cosine=reachable_cosine,
        signed_ceiling=signed_cosine, achieved_heldout=heldout_cosine,
        frac_of_theoretical_ceiling=(heldout_cosine / ceil_th) if ceil_th else float("nan"),
    )


# ----------------------------------------------------------------------------- L2: magnitude-capped recipes
@dataclass
class ReliabilityPoint:
    k: int
    gene: str
    cosine: float
    additivity_risk: float
    reliability: float               # 1 - risk = expected realized-additive fidelity


def reliability_curve(E_kd: np.ndarray, gene_names: Sequence[str], d_full: np.ndarray,
                      mask: np.ndarray, *, k_max: int = 40) -> list[ReliabilityPoint]:
    """L2. Greedy recipe reliability vs size k, using the validated ``additivity_risk`` primitive.

    At each greedy step we refit NNLS on the active support and score the recipe with
    ``reachability.additivity_risk`` (the magnitude-saturation law calibrated on Norman 2019
    K562 double perturbations, ``M* = 13.9``). The reference norm ``s`` is the median
    single-effect norm on the FULL readout axis -- the convention that reproduces the
    published ``atlas_additivity_risk`` values. Reliability = ``1 - risk`` = the expected
    fraction of the recipe's additively-predicted push that survives saturation.
    """
    s_full = float(np.median([n for n in np.linalg.norm(E_kd, axis=1) if n > 0]))
    spec = rx.reachability_spectrum(E_kd, d_full, k_max=k_max, hvg_mask=mask, refit_full=False)
    order = spec["order"]
    cos_k = spec["cosine"]
    A = E_kd[:, mask].T
    d = d_full[mask]
    out: list[ReliabilityPoint] = []
    for k in range(1, len(order) + 1):
        active = np.array(order[:k])
        w, _ = nnls(A[:, active], d)
        risk = rx.additivity_risk(E_kd[:, mask], w, active=active, median_single_norm=s_full)
        out.append(ReliabilityPoint(k=k, gene=str(gene_names[active[-1]]),
                                    cosine=float(cos_k[k - 1]),
                                    additivity_risk=float(risk), reliability=1.0 - float(risk)))
    return out


def cap_binding_k(curve: Sequence[ReliabilityPoint], threshold: float = 0.90) -> Optional[int]:
    """Smallest recipe size k at which reliability first drops below ``threshold`` (None if never)."""
    for p in curve:
        if p.reliability < threshold:
            return p.k
    return None


# ----------------------------------------------------------------------------- L1: held-out-modality certificate test
@dataclass
class ModalityTestResult:
    n_cert_genes: int
    auroc: float
    precision_at_n: float
    n_true_activation: int
    null_auroc_mean: float
    null_auroc_std: float
    z: float


def held_out_modality_test(E_kd: np.ndarray, d: np.ndarray, activation_hits_mask: np.ndarray, *,
                           hvg_mask: Optional[np.ndarray] = None, cert_top: int = 50,
                           n_null: int = 200, seed: int = 0) -> ModalityTestResult:
    """L1. Does the knockdown-only certificate predict a *held-out* activation arm?

    This is the runnable test for the single most novel -- and currently untested -- claim:
    that the genes the certificate flags as "must be switched ON" (positive entries of the
    NNLS residual, i.e. target demand no non-negative knockdown mix can meet) are genuinely
    the genes an activation (CRISPRa) screen would find effective.

    CONTRACT (requires data not present locally):
      E_kd                 : (P, G) measured KNOCKDOWN effect dictionary (the arm we KEEP).
      d                    : (G,)   target shift, same readout axis.
      activation_hits_mask : (G,)   bool -- genes the HIDDEN activation arm shows to be
                             effective (ground truth). In real use = CRISPRa-responsive genes
                             from a screen with BOTH arms on the same axis.

    Scores the knockdown-only residual (certificate score per gene) against the hidden
    activation labels with AUROC + a shuffled-label null. Demonstrated correct on a synthetic
    dual-modality fixture (AUROC ~1.0, z~9) where the ground truth is known by construction.
    No local screen has both arms (CD4+ is CRISPRi-only; Norman K562 is CRISPRa-only), so this
    is the scaffold that runs the moment such data is supplied.
    """
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, bool)
        E_kd, d, activation_hits_mask = E_kd[:, m], d[m], activation_hits_mask[m]
    res = rx.reachability(E_kd, d)
    cert_score = res.residual                        # d - E^T w*: want-up-but-unreachable direction
    y = activation_hits_mask.astype(int)
    if not (0 < y.sum() < len(y)):
        raise ValueError("activation_hits_mask must have both positive and negative entries")
    auroc = float(roc_auc_score(y, cert_score))
    topn = np.argsort(-cert_score)[:cert_top]
    prec = float(y[topn].mean())
    null = np.array([roc_auc_score(rng.permutation(y), cert_score) for _ in range(n_null)])
    z = float((auroc - null.mean()) / null.std()) if null.std() else float("inf")
    return ModalityTestResult(
        n_cert_genes=int((cert_score > 0).sum()), auroc=auroc, precision_at_n=prec,
        n_true_activation=int(y.sum()), null_auroc_mean=float(null.mean()),
        null_auroc_std=float(null.std()), z=z,
    )
