"""Convex-cone reachability for cell-state engineering.

The headline method for the cell-state-reachability project. Given a dictionary of
*measured* CRISPRi knockdown effect vectors `E` (P perturbations x G genes) and a target
transcriptomic shift `d` (G,), it answers a feasibility question, not a ranking one:

    Is `d` reachable by some NON-NEGATIVE combination of knockdown effects
    (i.e. is `d` inside the convex cone spanned by the rows of E), and if not,
    which genes make it unreachable (the constructive infeasibility certificate)?

Why a cone (the whole method in two facts):
  1. CRISPRi is loss-of-function: you can apply a knockdown or not, you cannot apply a
     "negative knockdown". So combination weights are non-negative: w >= 0.
  2. You can co-perturb: several knockdowns at once => the effect is (assumed) additive.
Together, the reachable set is C = { E^T w : w >= 0 } — a finitely generated convex
cone. Testing d in C is a non-negative least squares (NNLS) problem.

What makes the verdict falsifiable and honest:
  * reachable  -> residual ~ 0; the weights w rank a candidate perturbation mixture.
  * outside    -> residual has a component NO non-negative mix can supply. We return the
                  Farkas separating direction rho = d - E^T w*. Its positive coordinates
                  rank readouts that remain under-delivered at the closest cone point. They
                  motivate, but do not validate, concrete CRISPRa/de-repression hypotheses.
  * "meaningfully outside" is decided against a SHUFFLED-TARGET null, never a hardcoded
    threshold, because with P<<G most random targets are partly outside by construction.

Signed / bidirectional extension (`signed_reachability`):
  The one-sided cone above answers "how far can KNOCKDOWN get?". The signed solver answers
  the modality question: after fitting the measured loss-of-function (knockdown) cone, it
  fits the remaining residual with hypothetical gain-of-function atoms modelled as sign-flips
  of measured knockdown effects. It returns a sequential Pythagorean split of the target norm into
  LOF-reachable / GOF-only / neither (they sum to 1). This turns "31% needs activation" into
  a model-based variance decomposition with a knockdown recipe AND a candidate activation recipe, and is
  the geometric basis for modality triage (inhibit/degrade vs activate/agonize vs undruggable).

Dependencies: numpy, scipy.optimize.nnls. CPU-cheap (seconds at 34k x 2k). No GPU, no torch.

Author: cell-state-reachability. MIT.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import nnls


# Calibrated additivity/epistasis model constant (see the "Additivity / epistasis model"
# section lower in this module for the full derivation). M*: the recipe additive-magnitude
# saturation ceiling, in units of the dictionary's median single-effect norm. Calibrated on
# Norman et al. 2019 measured double perturbations (K562 CRISPRa, GSE133344). Defined here so
# it is available as a default argument of reachability_spectrum().
SATURATION_CEILING = 13.9


# --------------------------------------------------------------------------------------
# NNLS acceleration engine  (output-preserving; opt-in)
# --------------------------------------------------------------------------------------
# Every null / held-out / bootstrap loop refits NNLS on a FIXED dictionary A = E^T while
# only the right-hand side changes. scipy.optimize.nnls re-derives its normal-equation
# quantities from A on every call. Two accelerations below leave the numerical result at
# the SAME optimum (identical support, cosine, KKT certificate) while removing that waste:
#
#   (1) Gram reuse  — precompute G = A^T A once and run a Gram-based Lawson-Hanson active
#       set, so each solve works in the small P-dimensional passive subspace instead of
#       rescanning the (G x P) matrix. Guarded: every Gram solve is verified against the
#       method's own KKT/Farkas certificate (max_p <e_p, rho> <= tol); if it fails the
#       optimality bound (ill-conditioned Gram) or raises, we fall back to exact scipy.nnls.
#
#   (2) Process parallelism — the shuffles/resamples are independent, but scipy's compiled
#       _nnls holds the GIL (threads give no speedup), so parallelism must be process-based.
#       A fork pool COW-inherits the dictionary; one BLAS thread per worker avoids
#       oversubscription. Falls back to serial if fork is unavailable (Windows/spawn) or
#       blocked (sandboxed semaphores) — so reproduction is deterministic on any platform.
#
# The single-solve public API (reachability, signed_reachability) is intentionally NOT
# routed through this engine: its results anchor cross-function equalities asserted to
# <1e-9 in the self-test, so it stays byte-identical to the original scipy path.
_GRAM_CERT_TOL = 1e-5      # the self-test's own KKT optimality bound; fall back if exceeded


def nnls_gram(GtG: np.ndarray, Gtb: np.ndarray, *,
              tol: float = 1e-10, max_iter: Optional[int] = None) -> np.ndarray:
    """Lawson-Hanson NNLS from a PRECOMPUTED Gram matrix GtG = A^T A and Gtb = A^T b.

    Returns x >= 0 minimizing ||A x - b||_2, at the same optimum as scipy.optimize.nnls(A, b).
    Once GtG is cached across a loop, each call costs active-set submatrix solves in P-space,
    not a rescan of the (G x P) matrix A. Raises np.linalg.LinAlgError on a singular passive
    submatrix (caught by _nnls_auto, which then falls back to scipy).
    """
    n = GtG.shape[0]
    if max_iter is None:
        max_iter = 3 * n
    x = np.zeros(n)
    passive = np.zeros(n, dtype=bool)
    grad = Gtb - GtG @ x                       # = A^T (b - A x): the KKT gradient
    it = 0
    while (~passive).any() and grad[~passive].max() > tol:
        it += 1
        if it > max_iter:
            break
        cand = np.where(~passive)[0]
        passive[cand[np.argmax(grad[cand])]] = True
        while True:                            # inner feasibility loop
            idx = np.where(passive)[0]
            s = np.zeros(n)
            s[idx] = np.linalg.solve(GtG[np.ix_(idx, idx)], Gtb[idx])
            if s[idx].min() > 0:
                x = s
                break
            bad = passive & (s <= 0)
            alpha = (x[bad] / (x[bad] - s[bad])).min()
            x = x + alpha * (s - x)
            passive[(np.abs(x) < tol) & passive] = False
        grad = Gtb - GtG @ x
    return x


def _nnls_auto(A: np.ndarray, b: np.ndarray, *,
               GtG: Optional[np.ndarray] = None,
               max_iter: Optional[int] = None,
               cert_tol: float = _GRAM_CERT_TOL) -> np.ndarray:
    """Solve min_{w>=0} ||A w - b|| returning w. Uses the cached-Gram fast path when GtG is
    supplied AND max_iter is unconstrained, but only ACCEPTS that result if it meets the
    KKT/Farkas optimality bound; otherwise (or on a singular Gram) it falls back to the exact
    scipy.optimize.nnls. The returned w is therefore always at the true optimum.
    """
    if GtG is not None and max_iter is None:
        try:
            w = nnls_gram(GtG, A.T @ b)
            rho = b - A @ w
            cert = max(0.0, float(np.max(A.T @ rho))) if A.shape[1] else 0.0
            if cert <= cert_tol:
                return w
        except np.linalg.LinAlgError:
            pass                                # singular passive submatrix -> exact solver
    kw = {} if max_iter is None else {"maxiter": max_iter}
    w, _ = nnls(A, b, **kw)
    return w


def _resolve_n_jobs(n_jobs: int) -> int:
    """Map n_jobs (1=serial, -1=all cores, k=k workers) to a concrete positive worker count."""
    import os
    ncpu = os.cpu_count() or 1
    if n_jobs is None or n_jobs == 0:
        return 1
    if n_jobs < 0:
        return max(1, ncpu + 1 + n_jobs)        # -1 -> ncpu, -2 -> ncpu-1, ...
    return min(n_jobs, ncpu)


# Fork-worker state: the fixed dictionary and its Gram live in a module global so children
# inherit them by copy-on-write instead of re-pickling a 340 MB matrix per task.
_WORKER: dict = {}


def _worker_init(A: np.ndarray, GtG: Optional[np.ndarray]) -> None:
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)                    # one BLAS thread per process
    except Exception:
        pass
    _WORKER["A"] = A
    _WORKER["GtG"] = GtG


def _worker_solve(b: np.ndarray) -> np.ndarray:
    return _nnls_auto(_WORKER["A"], b, GtG=_WORKER["GtG"])


def _nnls_batch_w(A: np.ndarray, rhs_list, *,
                  n_jobs: int = 1, use_gram: bool = True,
                  max_iter: Optional[int] = None):
    """Return [ argmin_{w>=0} ||A w - b|| for b in rhs_list ] for a FIXED matrix A.

    Gram-reused when use_gram and max_iter is None; parallel over the right-hand sides via a
    fork pool when n_jobs != 1. Output is independent of n_jobs and of the Gram/scipy choice
    (both reach the same optimum), so results are deterministic and reproduction-safe. Any
    parallel-backend failure (no fork; sandboxed semaphores) degrades cleanly to serial.
    """
    GtG = (A.T @ A) if (use_gram and max_iter is None) else None
    n_workers = _resolve_n_jobs(n_jobs)
    if n_workers <= 1 or len(rhs_list) <= 1:
        return [_nnls_auto(A, b, GtG=GtG, max_iter=max_iter) for b in rhs_list]
    import multiprocessing as mp
    try:
        ctx = mp.get_context("fork")            # raises ValueError on spawn-only platforms
    except ValueError:
        return [_nnls_auto(A, b, GtG=GtG, max_iter=max_iter) for b in rhs_list]
    try:
        with ctx.Pool(n_workers, initializer=_worker_init, initargs=(A, GtG)) as pool:
            return pool.map(_worker_solve, rhs_list)
    except (OSError, ValueError, ImportError):  # e.g. sandboxed semaphore limits
        return [_nnls_auto(A, b, GtG=GtG, max_iter=max_iter) for b in rhs_list]


# --------------------------------------------------------------------------------------
# Result containers
# --------------------------------------------------------------------------------------
@dataclass
class ReachResult:
    """Geometry of fitting target `d` inside the non-negative cone spanned by rows of E."""
    reachable_cosine: float            # cos(E^T w*, d) at the NNLS optimum, in [-1, 1]
    residual_norm: float               # ||E^T w* - d|| / ||d||  (0 = exactly reachable)
    in_cone_fraction: float            # ||proj onto cone|| / ||d||  (fraction of d explained)
    weights: np.ndarray                # w* >= 0, length P (the ranked knockdown recipe)
    support: np.ndarray                # indices of w* > tol, sorted by weight desc
    fitted: np.ndarray                 # E^T w*  (the closest reachable point), length G
    residual: np.ndarray               # rho = d - E^T w*  (the unmet-demand / certificate direction)
    cert_max_violation: float          # max_p <e_p, rho>  — should be <= tol if w* is optimal
    n_generators: int
    n_genes: int
    # DEG-weighted evaluation (optional; defaults keep every existing construction unchanged).
    # When weighted is True the scalar scores above (reachable_cosine, residual_norm,
    # in_cone_fraction, cert_max_violation) are computed under the gene weights in
    # weight_vector; the arrays (weights, fitted, residual) remain in the ORIGINAL gene space
    # so activation_certificate / directional_decomposition keep working unchanged.
    weighted: bool = False
    weight_vector: Optional[np.ndarray] = None   # (G_fit,) gene weights used, in fit space

    def verdict(self, null: "NullResult | None" = None,
                cos_hi: float = 0.9, resid_lo: float = 0.35) -> str:
        """reachable | partially-reachable | outside.

        If a shuffled-target `null` is supplied the call is data-driven: the label 'outside'
        means no detectable directional reachability above this null, not merely the geometric
        fact that a non-zero residual places a target outside the exact cone. 'reachable'
        means the score clears the null's upper tail AND the residual is small. Without a null
        we fall back to geometry thresholds and you MUST report that the verdict was not
        null-calibrated.
        """
        if null is not None:
            if self.reachable_cosine <= null.p95:
                return "outside"                       # not distinguishable from noise
            if self.reachable_cosine >= null.p99 and self.residual_norm <= resid_lo:
                return "reachable"
            return "partially-reachable"
        # geometry-only fallback (uncalibrated — say so when you report it)
        if self.residual_norm <= resid_lo and self.reachable_cosine >= cos_hi:
            return "reachable"
        if self.reachable_cosine < 0.5:
            return "outside"
        return "partially-reachable"


@dataclass
class NullResult:
    observed: float
    null_cosines: np.ndarray
    p50: float
    p95: float
    p99: float
    percentile_of_observed: float      # where the observed value sits in the null (0-100)
    z: float                           # (observed - null_mean) / null_sd


@dataclass
class ActivationCertificate:
    """Constructive infeasibility certificate for the target direction as a whole."""
    gene_index: np.ndarray             # indices into the gene axis, ranked most-unmet first
    residual_value: np.ndarray         # rho_j at those genes (positive = wants up, unmet)
    target_value: np.ndarray           # d_j at those genes
    explanation: str = field(default="")


@dataclass
class SignedReachResult:
    """Staged signed reachability: how far a target is reachable when BOTH
    loss-of-function (knockdown, measured) AND a *hypothetical* gain-of-function
    (activation, modelled as the sign-flip of a measured effect vector) are allowed.

    The decomposition is an exact *sequential* Pythagorean split of the target norm into
    three parts that sum to 1: the LOF projection is orthogonal to its residual, and the
    GOF-proxy projection is orthogonal to the final residual. The LOF and GOF fitted vectors
    need not be pairwise orthogonal to each other.

        ||d||^2  =  ||fit_lof||^2  +  ||fit_gof||^2  +  ||resid||^2
                     \\_________/     \\_________/       \\______/
                      LOF fraction     GOF fraction      neither

    * lof_fraction : share of the target a KNOCKDOWN screen can reach (measured, real).
    * gof_fraction : ADDITIONAL share unlocked only by activating genes (a CRISPRa/agonist
                     hypothesis -- the modality the knockdown assay structurally cannot test).
    * neither_fraction : share reachable by NO single-gene perturbation in the library, in
                     either direction (the genuinely-outside residual).

    Modelling assumption (state it plainly): activation of gene g is approximated by the
    NEGATED measured knockdown effect vector of g (i.e. -e_g). This is the same additive,
    linear, sign-symmetric assumption the field's linear baseline already licenses; it is a
    hypothesis generator for a CRISPRa arm, not a claim that activation effects are exactly
    the mirror of knockdown effects.
    """
    lof_fraction: float                # ||fit_lof||^2 / ||d||^2  (knockdown-reachable share)
    gof_fraction: float                # ||fit_gof||^2 / ||d||^2  (activation-only share)
    neither_fraction: float            # ||resid||^2   / ||d||^2  (unreachable either way)
    signed_cosine: float               # cos(fit_lof + fit_gof, d) -- best signed alignment
    lof_cosine: float                  # cos(fit_lof, d) -- knockdown-only alignment (== ReachResult.reachable_cosine)
    lof_support: np.ndarray            # knockdown recipe: perturbation indices, weight-desc
    gof_support: np.ndarray            # activation recipe: perturbation indices, weight-desc
    lof_weights: np.ndarray            # w >= 0 over P (knockdown mix)
    gof_weights: np.ndarray            # u >= 0 over P (activation mix, applied to -E rows)
    residual: np.ndarray               # d - fit_lof - fit_gof  (the genuinely-unreachable direction)
    fitted_lof: np.ndarray             # E^T w
    fitted_gof: np.ndarray             # -E_gof^T u
    cert_max_violation: float          # maximum of the two staged NNLS KKT violations
    n_generators: int
    n_genes: int


# --------------------------------------------------------------------------------------
# Core solver
# --------------------------------------------------------------------------------------
def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0.0 or nb == 0.0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def _wcosine(a: np.ndarray, b: np.ndarray, w: np.ndarray) -> float:
    """Weighted cosine: cos in the inner product <x,y>_w = sum_j w_j x_j y_j (w >= 0).

    Equivalent to the plain cosine between sqrt(w)*a and sqrt(w)*b. Reduces to _cosine
    exactly when w is constant. This is the similarity analog of the DEG-weighted MSE
    (WMSE) in Mejia et al. (ICML 2026): genes with larger weight (the perturbation's DEGs)
    dominate the score, so agreement on the quiet, near-zero background can no longer
    inflate it.
    """
    w = np.asarray(w, dtype=float)
    sa, sb = np.sqrt(w) * a, np.sqrt(w) * b
    return _cosine(sa, sb)


def deg_weights(d: np.ndarray, *, scheme: str = "abs", top_k: Optional[int] = None,
                q: float = 0.0, floor: float = 0.0) -> np.ndarray:
    """Build a non-negative DEG-weight vector for a target `d` (length G).

    Turns "score the DEGs, not the background" into an explicit weight per gene, the
    ingredient a DEG-aware metric (WMSE / weighted delta-R^2, Mejia et al. ICML 2026)
    needs. Schemes:

      * "abs"    : w_j = |d_j|              (continuous; a gene's weight is its effect size)
      * "sq"     : w_j = d_j^2             (variance weighting; sharper DEG emphasis)
      * "topk"   : w_j = 1 on the top_k genes by |d_j|, else `floor`  (hard DEG set)
      * "quantile": w_j = 1 where |d_j| >= quantile(|d|, q), else `floor`  (data-driven set)
      * "uniform": w_j = 1 for all j       (recovers the unweighted metric exactly)

    `floor` keeps background genes at a small non-zero weight for the mask schemes so the
    metric never fully ignores them (set 0.0 for a strict DEG-only score). The returned
    vector is normalized to mean 1 so weighted and unweighted magnitudes stay comparable.
    """
    d = np.asarray(d, dtype=float)
    ad = np.abs(d)
    G = d.shape[0]
    if scheme == "uniform":
        return np.ones(G)
    if scheme == "abs":
        w = ad.copy()
    elif scheme == "sq":
        w = d * d
    elif scheme == "topk":
        if top_k is None:
            raise ValueError("scheme='topk' requires top_k")
        w = np.full(G, float(floor))
        idx = np.argsort(-ad)[:min(top_k, G)]
        w[idx] = 1.0
    elif scheme == "quantile":
        thr = np.quantile(ad, q) if q > 0 else 0.0
        w = np.where(ad >= thr, 1.0, float(floor))
    else:
        raise ValueError(f"unknown scheme {scheme!r}")
    m = float(w.mean())
    return w / m if m > 0 else np.ones(G)


def reachability(E: np.ndarray, d: np.ndarray, *,
                 hvg_mask: Optional[np.ndarray] = None,
                 weights: Optional[np.ndarray] = None,
                 weight_tol: float = 1e-8,
                 max_iter: Optional[int] = None) -> ReachResult:
    """Fit target `d` inside the non-negative cone spanned by the rows of `E`.

    Parameters
    ----------
    E : (P, G) array
        Perturbation dictionary; row p is the measured effect vector of knockdown p.
    d : (G,) array
        Target transcriptomic shift (e.g. z-scored `toward_Th1` signature).
    hvg_mask : (G,) bool array, optional
        Restrict the fit to these genes (e.g. highly-variable genes) on BOTH E and d.
    weights : (G,) or (G_masked,) array, optional
        Non-negative per-gene weights for a DEG-aware fit and score (the WMSE analog of
        Mejia et al., ICML 2026). When given: (1) the NNLS objective is reweighted, i.e.
        w* = argmin_{w>=0} sum_j weights_j (E^T w - d)_j^2 (solved by row-scaling A and d
        by sqrt(weights)); (2) reachable_cosine / residual_norm / in_cone_fraction /
        cert_max_violation are all reported in the weighted inner product. `weights` may be
        given in full-gene space (length G) even with hvg_mask set — it is masked to match.
        Default None reproduces the unweighted result BIT-IDENTICALLY (verified in _selftest).
    weight_tol : float
        Weights below this are treated as zero when forming the support.

    Returns
    -------
    ReachResult  (see fields; call .verdict(null) for the classification)
    """
    E = np.asarray(E, dtype=float)
    d = np.asarray(d, dtype=float)
    if E.ndim != 2:
        raise ValueError(f"E must be 2-D (P,G); got {E.shape}")
    if d.shape[0] != E.shape[1]:
        raise ValueError(f"d length {d.shape[0]} != n_genes {E.shape[1]}")

    if hvg_mask is not None:
        hvg_mask = np.asarray(hvg_mask, dtype=bool)
        E = E[:, hvg_mask]
        d = d[hvg_mask]

    wv = None
    if weights is not None:
        wv = np.asarray(weights, dtype=float)
        if hvg_mask is not None and wv.shape[0] == hvg_mask.shape[0]:
            wv = wv[hvg_mask]                  # accept full-gene weights alongside a mask
        if wv.shape[0] != d.shape[0]:
            raise ValueError(f"weights length {wv.shape[0]} != n_genes {d.shape[0]}")
        if np.any(wv < 0):
            raise ValueError("weights must be non-negative")

    A = E.T                                   # (G, P): columns are the cone generators
    kw = {} if max_iter is None else {"maxiter": max_iter}
    if wv is None:
        w, _ = nnls(A, d, **kw)                # min_{w>=0} ||A w - d||_2
    else:
        sw = np.sqrt(wv)                       # weighted NNLS: scale rows of A and d by sqrt(w)
        w, _ = nnls(A * sw[:, None], d * sw, **kw)

    fitted = A @ w                             # closest reachable point E^T w* (ORIGINAL space)
    residual = d - fitted                      # rho: unmet demand / certificate direction

    if wv is None:
        dn = np.linalg.norm(d)
        resid_norm = float(np.linalg.norm(residual) / dn) if dn else 0.0
        in_cone = float(np.linalg.norm(fitted) / dn) if dn else 0.0
        reach_cos = _cosine(fitted, d)
        # KKT / Farkas check: at the optimum, <e_p, rho> <= 0 for every generator p.
        # (equivalently E @ rho <= 0). The max positive violation certifies optimality.
        cert_violation = max(0.0, float(np.max(E @ residual))) if E.shape[0] else 0.0
    else:
        # All scalar scores in the weighted inner product <x,y>_w = sum_j w_j x_j y_j.
        dnw = np.sqrt(float(np.sum(wv * d * d)))
        resid_norm = float(np.sqrt(np.sum(wv * residual * residual)) / dnw) if dnw else 0.0
        in_cone = float(np.sqrt(np.sum(wv * fitted * fitted)) / dnw) if dnw else 0.0
        reach_cos = _wcosine(fitted, d, wv)
        # Weighted KKT: gradient of the weighted objective is E @ (wv * rho); at the
        # weighted optimum its max over generators is <= 0.
        cert_violation = max(0.0, float(np.max(E @ (wv * residual)))) if E.shape[0] else 0.0

    support = np.where(w > weight_tol)[0]
    support = support[np.argsort(-w[support])]

    return ReachResult(
        reachable_cosine=reach_cos,
        residual_norm=resid_norm,
        in_cone_fraction=in_cone,
        weights=w,
        support=support,
        fitted=fitted,
        residual=residual,
        cert_max_violation=cert_violation,
        n_generators=E.shape[0],
        n_genes=E.shape[1],
        weighted=wv is not None,
        weight_vector=wv,
    )


def activation_certificate(res: ReachResult, d: np.ndarray, *,
                           gene_names: Optional[list] = None,
                           top: int = 25,
                           hvg_mask: Optional[np.ndarray] = None) -> ActivationCertificate:
    """Rank positive unmet readouts in the target's separating direction.

    A gene is an 'activation candidate' when the target wants it UP (d_j > 0) yet the
    closest reachable point still under-delivers there (residual rho_j > 0). The full
    residual is a valid Farkas separating direction. An individual coordinate is not, by
    itself, proof that no other cone point can raise that gene; it is a ranked contribution
    to the whole-vector mismatch and therefore a hypothesis for follow-up.
    """
    d = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        d = d[np.asarray(hvg_mask, dtype=bool)]
    rho = res.residual
    # score = unmet upward demand; only where the target itself wants the gene up
    score = np.where((rho > 0) & (d > 0), rho, 0.0)
    order = np.argsort(-score)
    order = order[score[order] > 0][:top]
    expl = (f"{len(order)} genes carry positive unmet demand at the closest cone point. "
            f"Together they contribute to the separating direction; individually they are "
            f"follow-up CRISPRa or de-repression hypotheses, not validated interventions.")
    return ActivationCertificate(
        gene_index=order,
        residual_value=rho[order],
        target_value=d[order],
        explanation=expl,
    )


def directional_decomposition(res: ReachResult, d: np.ndarray, *,
                              hvg_mask: Optional[np.ndarray] = None) -> dict:
    """Split the target norm into knockdown-reachable vs activation-requiring fractions.

    Returns fractions that sum to ~1 (Pythagorean split of ||d||^2 into the fitted
    component inside the cone and the orthogonal-ish residual outside it).
    """
    d = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        d = d[np.asarray(hvg_mask, dtype=bool)]
    dn2 = float(np.dot(d, d))
    if dn2 == 0.0:
        return {"reachable_fraction": 0.0, "residual_fraction": 0.0}
    reachable = float(np.dot(res.fitted, res.fitted)) / dn2
    residual = float(np.dot(res.residual, res.residual)) / dn2
    # unmet UPWARD demand specifically (the activation story)
    up = res.residual[(res.residual > 0) & (d > 0)]
    up_frac = float(np.dot(up, up)) / dn2 if up.size else 0.0
    return {
        "reachable_fraction": reachable,
        "residual_fraction": residual,
        "activation_required_fraction": up_frac,
    }


def signed_reachability(E: np.ndarray, d: np.ndarray, *,
                        hvg_mask: Optional[np.ndarray] = None,
                        gof_mask: Optional[np.ndarray] = None,
                        weight_tol: float = 1e-8,
                        max_iter: Optional[int] = None) -> SignedReachResult:
    """Fit `d` in the JOINT cone of loss-of-function (rows of E) AND gain-of-function
    (rows of -E) generators, and return an exact three-way decomposition of the target
    norm into LOF-reachable / GOF-only / neither.

    The naive way to do this -- stack [E; -E] and run one NNLS -- is degenerate: for any
    gene the solver can trivially use either +e or -e, the combined generator set spans a
    much larger space, and the "how much needs activation" question dissolves because
    everything looks reachable. We instead answer the decision question directly:

        1. LOF fit:  w* = argmin_{w>=0} ||E^T w - d||         (the real knockdown reach)
        2. residual after LOF:  r = d - E^T w*
        3. GOF fit on the residual:  u* = argmin_{u>=0} ||(-E_gof)^T u - r||
           -- i.e. how much of what knockdown COULD NOT reach is reachable by activation.
        4. final residual: rho = r - (-E_gof)^T u*  (reachable by neither direction)

    Because step 3 fits the LOF residual, the three squared-norm shares are (near-)exactly
    orthogonal and sum to ~1. This staged fit encodes the biological prior that knockdown is
    the assay we actually have, and activation is the additional modality we are QUANTIFYING
    the need for -- not a free generator to be used interchangeably.

    Parameters
    ----------
    E : (P, G)   measured knockdown effect dictionary.
    d : (G,)     target shift.
    hvg_mask : (G,) bool, optional   restrict fit to these genes (applied to E and d).
    gof_mask : (P,) bool, optional   which perturbations are allowed as ACTIVATION
        generators (e.g. only reproducible, strongly-on-target knockdowns whose sign-flip is
        a credible CRISPRa hypothesis). Defaults to all rows.
    """
    E = np.asarray(E, dtype=float)
    d = np.asarray(d, dtype=float)
    if E.ndim != 2:
        raise ValueError(f"E must be 2-D (P,G); got {E.shape}")
    if d.shape[0] != E.shape[1]:
        raise ValueError(f"d length {d.shape[0]} != n_genes {E.shape[1]}")
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, dtype=bool)
        E = E[:, m]; d = d[m]
    P = E.shape[0]
    if gof_mask is not None:
        gof_mask = np.asarray(gof_mask, dtype=bool)
        if gof_mask.shape[0] != P:
            raise ValueError(f"gof_mask length {gof_mask.shape[0]} != n_generators {P}")
    else:
        gof_mask = np.ones(P, dtype=bool)

    kw = {} if max_iter is None else {"maxiter": max_iter}
    dn2 = float(np.dot(d, d))
    if dn2 == 0.0:
        z = np.zeros(P)
        return SignedReachResult(0.0, 0.0, 0.0, 0.0, 0.0,
                                 np.array([], int), np.array([], int), z, z,
                                 d.copy(), np.zeros_like(d), np.zeros_like(d), 0.0, P, E.shape[1])

    # 1. LOF fit
    A = E.T
    w, _ = nnls(A, d, **kw)
    fit_lof = A @ w
    r = d - fit_lof

    # 2. GOF fit on the LOF residual, using the sign-flipped allowed generators
    Egof = -E[gof_mask]              # activation generators (rows)
    Agof = Egof.T
    u_sub, _ = nnls(Agof, r, **kw)
    fit_gof = Agof @ u_sub
    rho = r - fit_gof                # reachable by neither

    # scatter GOF weights back to full-P indexing
    u = np.zeros(P)
    u[np.where(gof_mask)[0]] = u_sub

    lof_frac = float(np.dot(fit_lof, fit_lof) / dn2)
    gof_frac = float(np.dot(fit_gof, fit_gof) / dn2)
    neither = float(np.dot(rho, rho) / dn2)

    fitted_signed = fit_lof + fit_gof
    # Per-stage optimality certificate. The fit is STAGED (LOF on d, then GOF on the LOF
    # residual r), so the certificate verifies each stage reached its own optimum -- NOT a
    # single joint-KKT condition on the final residual (which would not hold, because rho is
    # post-GOF and the LOF condition is defined on the pre-GOF residual r):
    #   * LOF optimal on d:  max_p <e_p, r>  <= tol   (Farkas certificate for the knockdown fit)
    #   * GOF optimal on r:  min gradient of NNLS(-E_gof, r)  >= -tol
    viol_lof = max(0.0, float(np.max(E @ r))) if P else 0.0
    grad_gof = Agof.T @ (fit_gof - r)                  # = Agof^T (Agof u* - r), >= 0 at optimum
    viol_gof = max(0.0, float(-grad_gof.min())) if grad_gof.size else 0.0
    cert_viol = max(viol_lof, viol_gof)

    lof_support = np.where(w > weight_tol)[0]
    lof_support = lof_support[np.argsort(-w[lof_support])]
    gof_support = np.where(u > weight_tol)[0]
    gof_support = gof_support[np.argsort(-u[gof_support])]

    return SignedReachResult(
        lof_fraction=lof_frac,
        gof_fraction=gof_frac,
        neither_fraction=neither,
        signed_cosine=_cosine(fitted_signed, d),
        lof_cosine=_cosine(fit_lof, d),
        lof_support=lof_support,
        gof_support=gof_support,
        lof_weights=w,
        gof_weights=u,
        residual=rho,
        fitted_lof=fit_lof,
        fitted_gof=fit_gof,
        cert_max_violation=cert_viol,
        n_generators=P,
        n_genes=E.shape[1],
    )


# --------------------------------------------------------------------------------------
# Honest nulls
# --------------------------------------------------------------------------------------
def shuffled_target_null(E: np.ndarray, d: np.ndarray, *,
                         n_iter: int = 1000, seed: int = 0,
                         hvg_mask: Optional[np.ndarray] = None,
                         weights: Optional[np.ndarray] = None,
                         n_jobs: int = 1) -> NullResult:
    """Permute the gene labels of `d` and refit; the null band is the achievable-by-chance
    reachability cosine. The observed cosine is meaningful only above this band.

    weights : (G,) or (G_masked,), optional
        Non-negative DEG weights. When given, both the observed cosine and every shuffled
        cosine are computed under the weighted objective/metric, with the weights permuted
        together with the target values, so the null is the correct negative control for a
        DEG-weighted statistic. Default None reproduces the unweighted null exactly.
    """
    rng = np.random.default_rng(seed)
    dd = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, dtype=bool)
        E_use, d_use = E[:, m], dd[m]
    else:
        E_use, d_use = E, dd
    wv = None
    if weights is not None:
        wv = np.asarray(weights, dtype=float)
        if hvg_mask is not None and wv.shape[0] == np.asarray(hvg_mask).shape[0]:
            wv = wv[np.asarray(hvg_mask, dtype=bool)]
        if wv.shape[0] != d_use.size:
            raise ValueError(f"weights length {wv.shape[0]} != n_genes {d_use.size}")
    observed = reachability(E_use, d_use, weights=wv).reachable_cosine
    A = E_use.T
    G = d_use.size
    cos = np.empty(n_iter)
    # Draw every shuffle up front so the rng stream (and thus the null) is IDENTICAL to the
    # original serial loop regardless of n_jobs — permutation is the only per-iter rng draw.
    perms = [rng.permutation(G) for _ in range(n_iter)]
    if wv is None:
        # Fixed dictionary A, varying RHS -> Gram-reused, optionally process-parallel.
        ws = _nnls_batch_w(A, [d_use[pp] for pp in perms], n_jobs=n_jobs)
        for i, (pp, w) in enumerate(zip(perms, ws)):
            dp = d_use[pp]
            cos[i] = _cosine(A @ w, dp)
    else:
        # Weighted: A is rescaled by sqrt(weights) each iter (not fixed) -> unchanged path.
        for i, pp in enumerate(perms):
            dp = d_use[pp]
            wvp = wv[pp]
            s = np.sqrt(wvp)
            w, _ = nnls(A * s[:, None], dp * s)
            cos[i] = _wcosine(A @ w, dp, wvp)
    p50, p95, p99 = np.percentile(cos, [50, 95, 99])
    pct = float((cos < observed).mean() * 100.0)
    z = float((observed - cos.mean()) / cos.std()) if cos.std() else np.inf
    return NullResult(observed=observed, null_cosines=cos,
                      p50=float(p50), p95=float(p95), p99=float(p99),
                      percentile_of_observed=pct, z=z)


# --------------------------------------------------------------------------------------
# Metric calibration: positive control (ceiling) + dynamic-range fraction
# (Mejia et al., "Needles in the Haystack", ICML 2026; bioRxiv 10.1101/2025.10.20.683304)
# --------------------------------------------------------------------------------------
@dataclass
class CalibrationResult:
    """Control-anchored calibration of a reachability score.

    A raw reachable cosine is uninterpretable on its own: different targets have different
    sparsity, so the SAME cosine can mean different things. This anchors it between two
    controls and reports a scale-free placement, following the metric-calibration framework
    of Mejia et al. (ICML 2026):

      * floor   (negative control): the shuffled-target null mean — what the metric awards
                to a target with no real structure. A well-calibrated metric puts the mean
                baseline here.
      * ceiling (positive control): the score achieved by a KNOWN-reachable target (the
                'interpolated duplicate'). A well-calibrated metric awards near-max here.
      * dynamic_range_fraction = (observed - floor) / (ceiling - floor), clipped to [0, 1]:
                0 = indistinguishable from noise, 1 = as good as a known-reachable target.
    """
    observed: float
    floor: float
    ceiling: float
    dynamic_range_fraction: float
    weighted: bool
    ceiling_kind: str
    n_ceiling: int


def positive_control_ceiling(E: np.ndarray, *,
                             hvg_mask: Optional[np.ndarray] = None,
                             weights_scheme: Optional[str] = None,
                             kind: str = "interpolated_duplicate",
                             n_targets: int = 40, support_size: int = 5,
                             noise: float = 0.25, seed: int = 0) -> dict:
    """Positive-control ceiling: the score a KNOWN-reachable target earns under this metric.

    The 'interpolated duplicate' of Mejia et al. (ICML 2026): construct targets that ARE, by
    construction, reachable, then measure how high the metric scores them. A trustworthy
    metric must award near-max here — this is the upper anchor the raw cosine is placed
    against (the lower anchor being the shuffled-target null).

    kind :
      * "interpolated_duplicate" : d = E^T w + noise for a random non-negative sparse w
        (a genuine cone member, jittered) — reachable by the exact generators that built it.
      * "single"                 : d = a single random generator's effect vector + noise
        (the sparsest reachable target).
    weights_scheme : if given (e.g. "abs"), each target is scored with its own DEG weights,
        so the ceiling is measured under the SAME weighting as the observed statistic.
    Returns dict with the per-target cosines and summary ceiling (median) + spread.
    """
    E = np.asarray(E, dtype=float)
    if hvg_mask is not None:
        E = E[:, np.asarray(hvg_mask, dtype=bool)]
    P, G = E.shape
    A = E.T
    rng = np.random.default_rng(seed)
    cosines = np.empty(n_targets)
    for i in range(n_targets):
        if kind == "single":
            j = rng.integers(P)
            base = A[:, j].copy()
        else:  # interpolated_duplicate
            idx = rng.choice(P, size=min(support_size, P), replace=False)
            w_true = rng.uniform(0.5, 2.0, size=idx.size)
            base = A[:, idx] @ w_true
        bn = np.linalg.norm(base)
        d_t = base + noise * bn / np.sqrt(G) * rng.standard_normal(G) if bn > 0 else base
        wv = deg_weights(d_t, scheme=weights_scheme) if weights_scheme else None
        cosines[i] = reachability(E, d_t, weights=wv).reachable_cosine
    return {
        "cosines": cosines,
        "ceiling": float(np.median(cosines)),
        "ceiling_p25": float(np.percentile(cosines, 25)),
        "ceiling_p75": float(np.percentile(cosines, 75)),
        "kind": kind, "n_targets": int(n_targets),
        "weighted": weights_scheme is not None,
    }


def calibrate_reachability(observed: float, floor: float, ceiling: float, *,
                           weighted: bool = False, ceiling_kind: str = "interpolated_duplicate",
                           n_ceiling: int = 0) -> CalibrationResult:
    """Place an observed reachable cosine between a negative floor and a positive ceiling.

    dynamic_range_fraction = (observed - floor) / (ceiling - floor), clipped to [0, 1].
    Scale-free and comparable across targets of differing sparsity — the calibrated number
    to report instead of a raw cosine (Mejia et al., ICML 2026).
    """
    denom = ceiling - floor
    drf = (observed - floor) / denom if denom > 1e-12 else 0.0
    drf = float(min(1.0, max(0.0, drf)))
    return CalibrationResult(observed=float(observed), floor=float(floor),
                             ceiling=float(ceiling), dynamic_range_fraction=drf,
                             weighted=weighted, ceiling_kind=ceiling_kind,
                             n_ceiling=int(n_ceiling))


@dataclass
class AnalyticNullResult:
    """Closed-form anisotropy-corrected null for the reachable cosine.

    The shuffled-target null permutes the gene labels of `d`, which PRESERVES the target's
    value multiset -- hence its mean, its norm, and therefore the fraction of its energy that
    lies along the uniform ("all-genes") direction. That surviving component is what lets a
    non-negative cone score a high cosine by chance whenever the dictionary has a shared,
    same-sign transcriptional signature (e.g. a common activation/stress axis). Decomposing a
    shuffled target into (uniform component) + (mean-zero residual), which the cone reaches
    through two ORTHOGONAL channels, gives

        E[null cosine]  ~  sqrt( (a . rho1)^2  +  (1 - a^2) . kappa )

    where
        a^2   = G * mean(d)^2 / ||d||^2      the target's DC / anisotropy fraction (shuffle-invariant),
        rho1  = reachable cosine of the uniform direction 1 (a fixed property of the dictionary),
        kappa = chance reachable cosine^2 for a random MEAN-ZERO target (the isotropic AC floor).

    Validated against shuffled-target nulls on synthetic dictionaries (Pearson(pred, empirical)
    ~ 0.999 held-out / 0.995 in-sample; RMSE ~ 0.03-0.046) and on the primary CD4 Rest
    dictionary across four target axes (max abs error ~0.05 in-sample, ~0.07 held-out): it
    reproduces both the ~0 null of a sign-balanced target (Th1/Th2) and the elevated ~0.26-0.34
    null of an up-dominated target (CD4 aging) that a naive z-against-zero would mistake for
    signal. `rho1` and `kappa` are estimated with a handful of NNLS solves, so the analytic
    null replaces ~1000 shuffles by ~10-20 fits.
    """
    observed: float
    anisotropy: float                   # a^2, the DC/uniform energy fraction of the target
    rho_uniform: float                  # rho1, reachable cosine of the uniform direction
    kappa: float                        # AC chance reachable cosine^2
    null_mean: float                    # analytic E[null cosine]
    null_sd: float                      # analytic SD of the null cosine
    z: float                            # anisotropy-corrected z of the observed cosine
    held_out: bool                      # invariants measured held-out (True) or in-sample (False)


def analytic_anisotropy_null(E: np.ndarray, d: np.ndarray, *,
                             hvg_mask: Optional[np.ndarray] = None,
                             held_out: bool = False, n_probe: int = 16,
                             seed: int = 0,
                             observed: Optional[float] = None) -> AnalyticNullResult:
    """Closed-form anisotropy-corrected null (see AnalyticNullResult for the derivation).

    Parameters
    ----------
    E, d, hvg_mask   as in `reachability` / `shuffled_target_null`.
    held_out : bool  if True, measure the dictionary invariants (rho1, kappa) with a held-out
                     gene split -- the honest regime that matches held_out_gene_validation, and
                     the one to use when comparing an observed HELD-OUT cosine. If False, they
                     are measured in-sample (matches the in-sample reachable_cosine).
    n_probe : int    number of random mean-zero probes used to estimate kappa (the AC floor).
    observed : float optional; the observed reachable cosine to z-score. Defaults to the
                     in-sample reachable cosine of `d` (or its held-out analogue when held_out).
    """
    dd = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, dtype=bool)
        E_use, d_use = E[:, m], dd[m]
    else:
        E_use, d_use = E, dd
    A = E_use.T
    G = d_use.size
    rng = np.random.default_rng(seed)

    # target anisotropy a^2 = G*mean^2 / ||d||^2  (== cos^2(d, 1); shuffle-invariant)
    dnorm2 = float(d_use @ d_use)
    a2 = float(G * d_use.mean() ** 2 / dnorm2) if dnorm2 > 0 else 0.0
    a = np.sqrt(a2)

    one = np.ones(G) / np.sqrt(G)

    def _fit_cos_insample(target):
        w, _ = nnls(A, target)
        return _cosine(A @ w, target)

    def _fit_cos_heldout(target, s):
        r = np.random.default_rng(s)
        idx = r.permutation(G); h = G // 2
        tr, te = idx[:h], idx[h:]
        w, _ = nnls(A[tr], target[tr])
        return _cosine(A[te] @ w, target[te])

    # rho1: reachable cosine of the uniform direction
    if held_out:
        rho1 = float(np.mean([_fit_cos_heldout(one, 7000 + s) for s in range(max(4, n_probe // 2))]))
    else:
        rho1 = float(_fit_cos_insample(one))

    # kappa: chance reachable cosine^2 for random mean-zero targets (AC floor)
    kv = []
    for s in range(n_probe):
        v = rng.standard_normal(G); v -= v.mean()
        if held_out:
            c = _fit_cos_heldout(v, 8500 + s)
        else:
            c = _fit_cos_insample(v)
        kv.append(max(0.0, c) ** 2)
    kappa = float(np.mean(kv))
    kappa_sd = float(np.std(kv))

    null_mean = float(np.sqrt(max(0.0, (a * rho1) ** 2 + (1.0 - a2) * kappa)))
    # SD propagated from the AC-floor channel (the dominant source of shuffle-to-shuffle
    # variance); the anisotropy channel is fixed under shuffling, so it does not fluctuate.
    null_sd = float(np.sqrt((1.0 - a2) * kappa_sd / 2.0)) if null_mean > 0 else 0.0

    if observed is None:
        observed = _fit_cos_insample(d_use) if not held_out else float(
            np.mean([_fit_cos_heldout(d_use, 9000 + s) for s in range(max(4, n_probe // 2))]))
    z = float((observed - null_mean) / null_sd) if null_sd > 0 else np.inf

    return AnalyticNullResult(observed=float(observed), anisotropy=a2, rho_uniform=rho1,
                              kappa=kappa, null_mean=null_mean, null_sd=null_sd, z=z,
                              held_out=held_out)


def reachability_spectrum(E: np.ndarray, d: np.ndarray, *, k_max: int = 12,
                          hvg_mask: Optional[np.ndarray] = None,
                          refit_full: bool = False,
                          epistasis_penalty: float = 0.0,
                          median_single_norm: Optional[float] = None,
                          ceiling: float = SATURATION_CEILING) -> dict:
    """Greedy forward selection under non-negativity: reachable cosine vs. sparsity k.
    The knee is a compact candidate panel, not a proof of globally minimum cardinality.
    Returns arrays k, cosine, residual, order.

    Two selection strategies (identical results on realistic P<G, mixed-sign dictionaries —
    verified in `_selftest`):

    * refit_full=False (default, FAST): orthogonal-matching-pursuit style. At each step,
      score every candidate generator by its correlation with the current residual, add the
      best, then refit NNLS once on the small active set. Cost ~ O(k_max) NNLS solves plus
      k_max cheap (P x G)@(G,) scorings — seconds at P~7000.
    * refit_full=True  (EXHAUSTIVE ONE-STEP GREEDY, SLOW): at each step, try adding each
      remaining generator and
      pick the one that maximises the refit cosine. Cost ~ O(k_max x P) NNLS solves —
      minutes-to-hours at P~7000. It is exact for each greedy step, not for the global
      best subset of size k. Use only as a correctness reference on small problems.

    The fast path is what makes the spectrum usable on the full Tier-2 dictionary; the exact
    path is retained so the equivalence can be re-checked whenever the data shape changes.

    Epistasis / additivity penalty
    ------------------------------
    epistasis_penalty : float, default 0.0
        0.0 reproduces the additive greedy selection EXACTLY (bit-identical order). When > 0,
        each candidate's directional benefit is discounted by the saturation it would incur
        under the magnitude model calibrated on measured double perturbations (see
        `additivity_risk` / SATURATION_CEILING): keep = 1 - penalty * risk(recipe magnitude).
        This steers the recipe toward reaching the target with less magnitude inflation --
        i.e. toward a set whose *additive* prediction is more likely to hold in real
        co-perturbation. 1.0 is a full saturation discount; ~0.5 is a moderate nudge. The
        reported `cosine` is always the TRUE (unpenalised) reachable cosine of the selected
        set, so curves remain comparable across penalty settings.
    median_single_norm : float, optional
        `s` in the risk model; defaults to the median row norm of E in the fit space.
    ceiling : float
        `M*`, the calibrated dimensionless saturation ceiling.
    """
    E = np.asarray(E, dtype=float)
    d = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, dtype=bool)
        E, d = E[:, m], d[m]
    A = E.T                                        # (G, P)
    P = E.shape[0]
    dn = np.linalg.norm(d)
    k_max = min(k_max, P)
    chosen: list[int] = []
    ks, coss, resids = [], [], []

    # Saturation-penalty scaffolding (inert when epistasis_penalty == 0.0).
    # `s` is the dictionary's median single-effect norm in the FIT space, so the scale-free
    # ceiling M* transfers; risk(a) = a / (a + M*.s) is the expected magnitude de-rating of a
    # recipe whose additive magnitude is `a`. A candidate's directional score is discounted by
    # the saturation it would incur, so a generator that reaches alignment with less magnitude
    # inflation is preferred over one that merely piles onto an already-saturated recipe.
    if epistasis_penalty:
        _rn = np.linalg.norm(E, axis=1)
        _s = median_single_norm if median_single_norm is not None else (
            float(np.median(_rn[_rn > 0])) if np.any(_rn > 0) else 1.0)
        def _risk(a):                              # magnitude a -> expected fraction lost
            return a / (a + ceiling * _s) if a > 0 else 0.0

    if refit_full:
        remaining = set(range(P))
        for _ in range(k_max):
            best_j, best_cos, best_fit = None, -np.inf, None
            for j in remaining:
                idx = chosen + [j]
                w, _ = nnls(A[:, idx], d)
                fit = A[:, idx] @ w
                c = _cosine(fit, d)
                if epistasis_penalty:              # discount by the refit recipe's saturation
                    a = float(np.linalg.norm((E[idx].T * w).sum(axis=1))) if len(idx) else 0.0
                    c = c * (1.0 - epistasis_penalty * _risk(a))
                if c > best_cos:
                    best_j, best_cos, best_fit = j, c, fit
            if best_j is None:
                break
            chosen.append(best_j)
            remaining.discard(best_j)
            ks.append(len(chosen)); coss.append(_cosine(best_fit, d))   # report the TRUE cosine
            resids.append(float(np.linalg.norm(best_fit - d) / dn) if dn else 0.0)
    else:
        col_norms = np.linalg.norm(A, axis=0) + 1e-12
        col_sq = col_norms ** 2
        residual = d.copy()
        cur_add = np.zeros_like(d)                 # running additive vector Sum_i w_i E_i (fit space)
        for _ in range(k_max):
            score = (A.T @ residual) / col_norms   # correlation of each generator with residual
            if epistasis_penalty:
                # candidate's 1-D non-negative step wj, tentative additive magnitude & risk
                wj = np.maximum(0.0, (A.T @ residual) / col_sq)          # (P,) 1-D NNLS coefficient
                a_tent = np.linalg.norm(cur_add[:, None] + wj[None, :] * A, axis=0)  # (P,)
                keep = 1.0 - epistasis_penalty * (a_tent / (a_tent + ceiling * _s))
                score = score * keep               # discount directional benefit by saturation
            if chosen:
                score[chosen] = -np.inf
            j = int(np.argmax(score))
            if score[j] <= 0:                      # no generator can further reduce the residual
                break
            chosen.append(j)
            w, _ = nnls(A[:, chosen], d)
            fit = A[:, chosen] @ w
            residual = d - fit
            cur_add = (E[chosen].T * w).sum(axis=1) if epistasis_penalty else cur_add
            ks.append(len(chosen)); coss.append(_cosine(fit, d))
            resids.append(float(np.linalg.norm(residual) / dn) if dn else 0.0)

    return {"k": np.array(ks), "cosine": np.array(coss),
            "residual": np.array(resids), "order": np.array(chosen)}


# --------------------------------------------------------------------------------------
# Additivity / epistasis model  (calibrated on measured double perturbations)
# --------------------------------------------------------------------------------------
# The reachable-cone method composes single-knockdown effect vectors ADDITIVELY. That
# assumption was calibrated against the one public screen that measures both singles and
# their doubles: Norman et al. 2019 (K562 CRISPRa, GSE133344), 126 doubles with both
# singles present. Two candidate epistasis mechanisms were tested against the measured
# non-additivity of each double:
#
#   * effect-vector COLLINEARITY (the intuitive "shared-program => sub-additive" prior):
#     Spearman rho ~ -0.16 with directional non-additivity (n.s.). REFUTED as a predictor
#     -- collinear pairs are, if anything, slightly MORE additive. Do not penalise on it.
#   * combined MAGNITUDE (saturation ceiling): Spearman rho ~ +0.58 (p<0.01) with the
#     magnitude deficit 1 - ||measured|| / ||additive||. ROBUST and biologically sensible:
#     you cannot push a gene's expression arbitrarily far, so stacking large effects yields
#     less than their sum. This is the calibrated penalty.
#
# The saturation law   achieved = a / (1 + a / (M* . s))   fits the Norman doubles with
# R^2 ~ 0.57, where `a` is the additive-predicted magnitude ||Sum_i w_i E_i||, `s` is the
# dictionary's median single-effect (row) norm, and M* is a dimensionless ceiling. Dividing
# by `s` makes the law SCALE-FREE, so the ceiling calibrated on Norman log-fold-change
# transfers to a z-score dictionary. Directional (angular) non-additivity, by contrast,
# improves with effect size (median cos(measured, additive) 0.64 -> 0.81 from the weakest to
# the strongest magnitude quartile), i.e. it is largely measurement noise at low SNR, so the
# calibrated model corrects only the magnitude channel and leaves the fit DIRECTION intact.
# (SATURATION_CEILING is defined near the top of the module so it can be a default argument
#  of reachability_spectrum.)


def additivity_risk(E: np.ndarray, weights: np.ndarray, *,
                    median_single_norm: Optional[float] = None,
                    ceiling: float = SATURATION_CEILING,
                    active: Optional[np.ndarray] = None) -> float:
    """Expected fraction of a recipe's additively-predicted magnitude lost to saturation.

    A per-recipe *additivity-risk* score in [0, 1): 0 = the recipe is small enough that
    additive composition is safe; ->1 = the recipe stacks so much effect that most of the
    predicted push is expected to be lost to the saturation ceiling measured in real
    double-perturbation data (Norman 2019). It scores the MAGNITUDE channel only -- the
    direction (hence the reachable cosine and the verdict) is unchanged; what it flags is
    how much to trust that the recipe's *realised* effect will match its additive prediction.

    risk = a / (a + M* . s),   a = ||Sum_i w_i E_i||,   s = median single-effect norm.

    Parameters
    ----------
    E : (P, G) array          the perturbation dictionary (rows = single-knockdown effects).
    weights : (P,) array      non-negative recipe weights (e.g. ReachResult.weights, or the
                              refit weights of a greedy support). Only nonzero entries matter.
    median_single_norm : float, optional
                              `s`; defaults to the median L2 norm of the rows of E. Pass an
                              explicit value to score against a reference dictionary.
    ceiling : float           `M*`; the calibrated dimensionless ceiling.
    active : (k,) int, optional
                              restrict to these row indices (a sparse support); weights are
                              then taken as `weights[active]` aligned to `active`.
    """
    E = np.asarray(E, dtype=float)
    w = np.asarray(weights, dtype=float)
    if active is not None:
        active = np.asarray(active, dtype=int)
        rows = E[active]
        wv = w[active] if w.shape[0] == E.shape[0] else w
    else:
        rows = E
        wv = w
    if median_single_norm is None:
        row_norms = np.linalg.norm(E, axis=1)
        median_single_norm = float(np.median(row_norms[row_norms > 0])) if np.any(row_norms > 0) else 1.0
    additive = wv @ rows                        # Sum_i w_i E_i   (G,)
    a = float(np.linalg.norm(additive))
    if a == 0.0:
        return 0.0
    return a / (a + ceiling * median_single_norm)


@dataclass
class HeldOutResult:
    """Held-out-GENE validation: does the recipe generalise to genes it was not fit on?"""
    in_sample_cosine: float
    held_out_cosine: float              # THE honest number — fit on half the genes, scored on the other half
    null_mean: float                    # shuffled-target held-out cosine (chance level for THIS target)
    null_std: float
    null_max: float
    z: float                            # (held_out - null_mean) / null_std
    n_shuffles: int


def held_out_gene_validation(E: np.ndarray, d: np.ndarray, *,
                             hvg_mask: Optional[np.ndarray] = None,
                             weights: Optional[np.ndarray] = None,
                             n_shuffles: int = 60, seed: int = 0,
                             n_jobs: int = 1) -> HeldOutResult:
    """Split the gene axis in half: fit NNLS weights on genes H1, score cosine on genes H2.

    This is the load-bearing honesty test. A dense non-negative fit over thousands of
    correlated generators can reach a target in-sample by overfitting; only structure that
    *generalises to held-out genes* is real. The paired shuffled-target null gives the
    chance level for THIS specific target — essential because a direction-biased target
    (e.g. mostly-up) is partly reachable by chance and has a null well above zero.

    weights : (G,) or (G_masked,) array, optional
        Non-negative per-gene DEG weights (WMSE analog). When given, BOTH the fit on H1 and
        the score on H2 use the weighted objective / weighted cosine, restricted to the
        genes in each half. The shuffled-target null permutes d AND carries its weights with
        the permuted values, so the null answers "chance level under the SAME weighting" —
        the correct control for a DEG-weighted statistic. Default None = unchanged behaviour.
    """
    E = np.asarray(E, dtype=float)
    d = np.asarray(d, dtype=float)
    if hvg_mask is not None:
        m = np.asarray(hvg_mask, dtype=bool)
        E, d = E[:, m], d[m]
    G = d.shape[0]
    wv = None
    if weights is not None:
        wv = np.asarray(weights, dtype=float)
        if hvg_mask is not None and wv.shape[0] == np.asarray(hvg_mask).shape[0]:
            wv = wv[np.asarray(hvg_mask, dtype=bool)]
        if wv.shape[0] != G:
            raise ValueError(f"weights length {wv.shape[0]} != n_genes {G}")
    rng = np.random.default_rng(seed)
    perm = rng.permutation(G)
    h1, h2 = perm[:G // 2], perm[G // 2:]
    A1, A2 = E[:, h1].T, E[:, h2].T

    def _fit_score(dfit, dscore):
        if wv is None:
            w, _ = nnls(A1, dfit[h1])
            return _cosine(A1 @ w, dfit[h1]), _cosine(A2 @ w, dscore[h2])
        s1 = np.sqrt(wv[h1])
        w, _ = nnls(A1 * s1[:, None], dfit[h1] * s1)
        return (_wcosine(A1 @ w, dfit[h1], wv[h1]),
                _wcosine(A2 @ w, dscore[h2], wv[h2]))

    in_cos, ho_cos = _fit_score(d, d)

    null = np.empty(n_shuffles)
    # Draw every shuffle up front so the null is IDENTICAL to the original serial loop
    # regardless of n_jobs (permutation is the only per-iter rng draw).
    perms = [rng.permutation(G) for _ in range(n_shuffles)]
    if wv is None:
        # Fixed fit matrix A1, varying RHS -> Gram-reused, optionally process-parallel.
        ws = _nnls_batch_w(A1, [d[pp][h1] for pp in perms], n_jobs=n_jobs)
        for i, (pp, wj) in enumerate(zip(perms, ws)):
            dp = d[pp]
            null[i] = _cosine(A2 @ wj, dp[h2])
    else:
        for i, pp in enumerate(perms):
            dp = d[pp]
            # weights travel with the permuted target values (weighting is a property of d)
            wvp = wv[pp]
            s1 = np.sqrt(wvp[h1])
            wj, _ = nnls(A1 * s1[:, None], dp[h1] * s1)
            null[i] = _wcosine(A2 @ wj, dp[h2], wvp[h2])
    z = float((ho_cos - null.mean()) / null.std()) if null.std() else np.inf
    return HeldOutResult(
        in_sample_cosine=in_cos, held_out_cosine=ho_cos,
        null_mean=float(null.mean()), null_std=float(null.std()),
        null_max=float(null.max()), z=z, n_shuffles=n_shuffles,
    )


def run_reachability(E: np.ndarray, d: np.ndarray, *,
                     gene_names: Optional[list] = None,
                     hvg_mask: Optional[np.ndarray] = None,
                     n_shuffles: int = 60, k_max: int = 15, seed: int = 0) -> dict:
    """One-call reachability report that REFUSES to give a verdict without a target-specific
    null. Returns the cone fit, held-out-gene validation, directional decomposition, greedy
    spectrum, activation certificate, and a null-calibrated verdict — everything the write-up
    needs, computed for THIS target so the raw cosine is never reported unqualified.
    """
    res = reachability(E, d, hvg_mask=hvg_mask)
    ho = held_out_gene_validation(E, d, hvg_mask=hvg_mask, n_shuffles=n_shuffles, seed=seed)
    dec = directional_decomposition(res, d, hvg_mask=hvg_mask)
    spec = reachability_spectrum(E, d, k_max=k_max, hvg_mask=hvg_mask)
    cert = activation_certificate(res, d, gene_names=gene_names, hvg_mask=hvg_mask)
    # verdict is driven by the held-out z-score, not the in-sample cosine
    if ho.z < 3:
        verdict = "outside"           # not distinguishable from a shuffled target
    elif dec["reachable_fraction"] >= 0.6 and res.residual_norm <= 0.35:
        verdict = "reachable"
    else:
        verdict = "partially-reachable"
    return {"result": res, "held_out": ho, "decomposition": dec,
            "spectrum": spec, "certificate": cert, "verdict": verdict,
            "cert_max_violation": res.cert_max_violation}


# --------------------------------------------------------------------------------------
# Experiment-design convenience API
# --------------------------------------------------------------------------------------
def _knee(k: np.ndarray, y: np.ndarray) -> int:
    """Return the `k` at the knee of a concave-increasing curve (Kneedle-lite).

    The optimal library size is the point of diminishing returns on the greedy
    reachability spectrum: cosine rises steeply, then flattens. For a concave-increasing
    curve the knee is the point of maximum vertical distance ABOVE the straight chord
    joining the first and last points. Dependency-free (no `kneed` package needed).
    Returns the *k value* (not the array index); falls back to the last k for degenerate
    (flat or <3-point) curves.
    """
    k = np.asarray(k, dtype=float); y = np.asarray(y, dtype=float)
    if k.size < 3:
        return int(k[-1]) if k.size else 0
    x0, x1 = k[0], k[-1]
    y0, y1 = y[0], y[-1]
    if x1 == x0 or y1 == y0:
        return int(k[-1])
    xn = (k - x0) / (x1 - x0)
    yn = (y - y0) / (y1 - y0)
    dist = yn - xn                       # >0 where the curve bows above the chord
    return int(k[int(np.argmax(dist))])


@dataclass
class DesignResult:
    """A researcher-facing experiment-design card for one current->target transition.

    Bundles everything needed to decide (a) whether the state is reachable at all,
    (b) which perturbations to apply and in which direction, and (c) how big a screen
    to run — with every number qualified by a target-specific null so the reachable
    cosine is never reported uncalibrated.
    """
    # --- verdict & headline geometry -------------------------------------------------
    verdict: str                       # outside | weakly reachable | partially reachable | reachable
    null_calibrated: bool              # False if held-out significance was skipped (n_shuffles=0)
    reachable_cosine: float            # LOF-only cosine (knockdown reach; == signed.lof_cosine)
    signed_cosine: float               # best LOF+GOF alignment
    residual_norm: float               # ||d - fitted_lof|| / ||d||
    lof_fraction: float                # knockdown-reachable share of the target norm
    gof_fraction: float                # activation-only share
    neither_fraction: float            # reachable by neither direction
    # --- significance ----------------------------------------------------------------
    held_out_cosine: float             # honest number: fit on half the genes, scored on the other half
    held_out_z: float                  # (held_out - shuffled_null_mean) / null_sd
    # --- the recipe (ranked, top-N) --------------------------------------------------
    knockdown_recipe: list             # [{gene, weight, rank}] — CRISPRi / inhibit / degrade moves
    activation_recipe: list            # [{gene, weight, rank}] — CRISPRa / agonize moves (from -E fit)
    n_knockdown_support: int           # full LOF support size (recipe is the top-N of this)
    n_activation_support: int          # full GOF support size
    # --- the library (optimal next-screen design) ------------------------------------
    optimal_k: int                     # knee of the greedy spectrum: how many perturbations to screen
    library: list                      # [{gene, k, cosine_at_k, marginal_gain}] up to k_max
    spectrum: dict                     # raw {k, cosine, residual, order} from reachability_spectrum
    # --- the constructive certificate ------------------------------------------------
    activation_certificate: list       # [{gene, residual, target}] readout genes knockdown cannot reach
    # --- bookkeeping -----------------------------------------------------------------
    n_generators: int
    n_readout: int
    cert_max_violation: float          # signed-cone KKT residual (should be ~0)

    def summary(self) -> str:
        """One-line human-readable headline for a design card."""
        return (f"[{self.verdict.upper()}] reach_cos={self.reachable_cosine:.3f} "
                f"(held-out {self.held_out_cosine:.3f}, z={self.held_out_z:.1f})  "
                f"LOF={self.lof_fraction:.2f}/GOF={self.gof_fraction:.2f}/"
                f"neither={self.neither_fraction:.2f}  optimal_k={self.optimal_k}")


def _grade_verdict(reachable_cosine: float, residual_norm: float,
                   held_out_z: float, null_calibrated: bool) -> str:
    """Four-level null-calibrated verdict, matching the atlas labelling convention.

    Significance gate first: a target whose held-out reach is indistinguishable from a
    shuffled target (z<3) is 'outside' regardless of its in-sample cosine. Above the gate,
    grade by knockdown reach — thresholds chosen so the real atlas transitions reproduce
    (reach_cos>=0.60 & z>=3 -> partially reachable; else weakly reachable), with a top
    'reachable' bin reserved for genuinely in-cone targets (small residual, high cosine).
    """
    if null_calibrated and held_out_z < 3:
        return "outside"
    if residual_norm <= 0.35 and reachable_cosine >= 0.85:
        return "reachable"
    if reachable_cosine >= 0.60:
        return "partially reachable"
    return "weakly reachable"


def design_experiment(E: np.ndarray, d: np.ndarray, *,
                      perturbation_names: Optional[list] = None,
                      readout_names: Optional[list] = None,
                      hvg_mask: Optional[np.ndarray] = None,
                      gof_mask: Optional[np.ndarray] = None,
                      k_max: int = 12,
                      top: int = 15,
                      n_shuffles: int = 20,
                      cert_top: int = 25,
                      seed: int = 0) -> DesignResult:
    """Design a cell-state-engineering experiment for one current->target transition.

    This is the one-call researcher-facing API. Given a measured knockdown effect
    dictionary `E` (P perturbations x G readout genes) and a target transcriptomic shift
    `d` (G,), it returns a complete, self-qualifying *design card*:

      1. VERDICT  — is the target reachable? (null-calibrated: 'outside' when held-out
         reach is indistinguishable from a shuffled target, never a hardcoded cosine cut).
      2. RECIPE   — the ranked knockdown recipe AND the ranked activation recipe, split by
         the signed LOF/GOF/neither decomposition (so 'what to do' comes with 'in which
         direction', i.e. CRISPRi vs CRISPRa).
      3. LIBRARY  — the optimal k-perturbation set for the NEXT screen, from the greedy
         reachability spectrum, with the knee (`optimal_k`) detected automatically.
      4. CERTIFICATE — the readout genes the target wants UP that no non-negative knockdown
         mix can deliver (the Farkas separating direction; concrete CRISPRa hypotheses).

    Parameters
    ----------
    E : (P, G) array          measured knockdown effect dictionary; row p is perturbation p.
    d : (G,)  array           target shift (current->target), on the same G gene axis as E.
    perturbation_names : length-P list, optional
        Name of the perturbed gene for each generator row -> names the recipe/library.
    readout_names : length-G list, optional
        Name of each readout gene -> names the activation certificate. If `hvg_mask` is
        given these are indexed on the FULL G axis (masking is applied internally).
    hvg_mask : (G,) bool, optional      restrict the fit to these readout genes.
    gof_mask : (P,) bool, optional      which perturbations may act as activation generators.
    k_max : int               how far to extend the greedy library spectrum.
    top : int                 how many recipe entries to return per direction.
    n_shuffles : int          held-out-gene shuffles for the significance gate; 0 skips it
                              (verdict then geometry-only, null_calibrated=False).
    cert_top : int            how many activation-certificate genes to return.
    seed : int                RNG seed for the held-out split & null.

    Returns
    -------
    DesignResult   (call .summary() for a one-liner; see the dataclass fields).

    Notes
    -----
    Modelling assumption (stated plainly, inherited from `signed_reachability`): activation
    of a gene is approximated by the sign-flip of its measured knockdown effect vector. The
    activation recipe and certificate are therefore falsifiable CRISPRa *hypotheses*, not a
    claim that activation exactly mirrors knockdown. The held-out cosine is the honest
    generalisation number; the in-sample `reachable_cosine` typically overstates it, which is
    exactly what the calibration/reliability analysis quantifies.
    """
    E = np.asarray(E, dtype=float)
    d = np.asarray(d, dtype=float)
    P, G = E.shape
    if d.shape[0] != G:
        raise ValueError(f"d length {d.shape[0]} != n_genes {G}")
    if perturbation_names is not None and len(perturbation_names) != P:
        raise ValueError(f"perturbation_names length {len(perturbation_names)} != P {P}")
    if readout_names is not None and len(readout_names) != G:
        raise ValueError(f"readout_names length {len(readout_names)} != G {G}")

    def pname(i: int):
        return str(perturbation_names[i]) if perturbation_names is not None else int(i)

    # 1. signed decomposition — the modality-aware fit (LOF + GOF + neither)
    s = signed_reachability(E, d, hvg_mask=hvg_mask, gof_mask=gof_mask)
    dm = d[np.asarray(hvg_mask, dtype=bool)] if hvg_mask is not None else d
    dmn = float(np.linalg.norm(dm))
    residual_norm = float(np.linalg.norm(dm - s.fitted_lof) / dmn) if dmn else 0.0

    # 2. significance gate — held-out-gene validation against a shuffled-target null
    if n_shuffles and n_shuffles > 0:
        ho = held_out_gene_validation(E, d, hvg_mask=hvg_mask, n_shuffles=n_shuffles, seed=seed)
        held_out_cosine, held_out_z, null_cal = ho.held_out_cosine, ho.z, True
    else:
        held_out_cosine, held_out_z, null_cal = float("nan"), float("nan"), False

    verdict = _grade_verdict(s.lof_cosine, residual_norm, held_out_z, null_cal)

    # 3. ranked recipes (top-N of each direction's support)
    knockdown_recipe = [{"gene": pname(int(i)), "weight": float(s.lof_weights[int(i)]),
                         "rank": r + 1} for r, i in enumerate(s.lof_support[:top])]
    activation_recipe = [{"gene": pname(int(i)), "weight": float(s.gof_weights[int(i)]),
                          "rank": r + 1} for r, i in enumerate(s.gof_support[:top])]

    # 4. optimal-k library from the greedy spectrum (fast OMP path) + knee
    spec = reachability_spectrum(E, d, k_max=k_max, hvg_mask=hvg_mask, refit_full=False)
    ks, coss, order = spec["k"], spec["cosine"], spec["order"]
    optimal_k = _knee(ks, coss) if ks.size else 0
    library = []
    prev = 0.0
    for j, kk in enumerate(ks):
        gi = int(order[j])
        library.append({"gene": pname(gi), "k": int(kk),
                        "cosine_at_k": float(coss[j]),
                        "marginal_gain": float(coss[j] - prev)})
        prev = float(coss[j])

    # 5. activation certificate — readout genes knockdown cannot push up
    base = reachability(E, d, hvg_mask=hvg_mask)
    cert = activation_certificate(base, d, hvg_mask=hvg_mask, top=cert_top)
    # certificate gene indices are into the (masked) readout axis; map back to names
    if readout_names is not None:
        if hvg_mask is not None:
            masked_names = list(np.asarray(readout_names)[np.asarray(hvg_mask, dtype=bool)])
        else:
            masked_names = list(readout_names)
        cert_names = [str(masked_names[int(g)]) for g in cert.gene_index]
    else:
        cert_names = [int(g) for g in cert.gene_index]
    activation_cert = [{"gene": cert_names[i], "residual": float(cert.residual_value[i]),
                        "target": float(cert.target_value[i])} for i in range(len(cert_names))]

    return DesignResult(
        verdict=verdict, null_calibrated=null_cal,
        reachable_cosine=float(s.lof_cosine), signed_cosine=float(s.signed_cosine),
        residual_norm=residual_norm,
        lof_fraction=float(s.lof_fraction), gof_fraction=float(s.gof_fraction),
        neither_fraction=float(s.neither_fraction),
        held_out_cosine=float(held_out_cosine), held_out_z=float(held_out_z),
        knockdown_recipe=knockdown_recipe, activation_recipe=activation_recipe,
        n_knockdown_support=int(s.lof_support.size), n_activation_support=int(s.gof_support.size),
        optimal_k=optimal_k, library=library, spectrum=spec,
        activation_certificate=activation_cert,
        n_generators=P, n_readout=int(dm.shape[0]),
        cert_max_violation=float(s.cert_max_violation),
    )


# --------------------------------------------------------------------------------------
# Self-test on synthetic data (runnable now, no Tier-2 download needed)
# --------------------------------------------------------------------------------------
def _selftest(seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    P, G = 300, 800
    E = rng.standard_normal((P, G))                       # mixed-sign effect vectors

    # (1) an IN-cone target: a known non-negative mix of a few generators
    true_idx = rng.choice(P, size=4, replace=False)
    true_w = rng.uniform(0.5, 2.0, size=4)
    d_in = (true_w[:, None] * E[true_idx]).sum(0)
    r_in = reachability(E, d_in)
    null_in = shuffled_target_null(E, d_in, n_iter=200, seed=1)

    # (2) an OUT-of-cone target: add a direction orthogonal to the row space is hard in
    #     P<G, so instead demand the NEGATIVE of an in-cone vector (needs "un-mixing").
    d_out = -d_in + 0.3 * rng.standard_normal(G)
    r_out = reachability(E, d_out)
    null_out = shuffled_target_null(E, d_out, n_iter=200, seed=2)
    cert = activation_certificate(r_out, d_out, top=10)
    dec_in = directional_decomposition(r_in, d_in)
    dec_out = directional_decomposition(r_out, d_out)

    print("=== IN-cone target (should be reachable) ===")
    print(f"  cosine={r_in.reachable_cosine:.4f}  residual={r_in.residual_norm:.4f}  "
          f"|support|={r_in.support.size}  recovered_true={set(true_idx) <= set(r_in.support.tolist())}")
    print(f"  KKT/Farkas max violation (should be <= ~1e-6): {r_in.cert_max_violation:.2e}")
    print(f"  null p95={null_in.p95:.3f} p99={null_in.p99:.3f}  observed pct={null_in.percentile_of_observed:.1f}")
    print(f"  verdict(null)={r_in.verdict(null_in)}")
    print(f"  decomposition reachable={dec_in['reachable_fraction']:.3f} "
          f"residual={dec_in['residual_fraction']:.3f}")

    print("=== OUT-of-cone target (should be less reachable) ===")
    print(f"  cosine={r_out.reachable_cosine:.4f}  residual={r_out.residual_norm:.4f}")
    print(f"  KKT/Farkas max violation (should be <= ~1e-6): {r_out.cert_max_violation:.2e}")
    print(f"  verdict(null)={r_out.verdict(null_out)}")
    print(f"  decomposition reachable={dec_out['reachable_fraction']:.3f} "
          f"activation_required={dec_out['activation_required_fraction']:.3f}")
    print(f"  activation certificate: {cert.gene_index.size} genes; top rho="
          f"{np.round(cert.residual_value[:5], 2).tolist()}")

    # (3) fast (OMP-style) spectrum must select the SAME generators as the exact scan
    spec_fast = reachability_spectrum(E, d_out, k_max=8, refit_full=False)
    spec_exact = reachability_spectrum(E, d_out, k_max=8, refit_full=True)
    same_order = spec_fast["order"].tolist() == spec_exact["order"].tolist()
    print("=== spectrum: fast (OMP) vs exact greedy ===")
    print(f"  same selection order: {same_order}")
    print(f"  fast cosine[-1]={spec_fast['cosine'][-1]:.4f}  exact cosine[-1]={spec_exact['cosine'][-1]:.4f}")

    # (4) held-out-gene validation + one-call report
    ho = held_out_gene_validation(E, d_in, n_shuffles=50, seed=3)
    print("=== held-out-gene validation (in-cone target) ===")
    print(f"  in-sample={ho.in_sample_cosine:.3f}  held-out={ho.held_out_cosine:.3f}  "
          f"null_mean={ho.null_mean:+.3f}  z={ho.z:.1f}")
    report = run_reachability(E, d_out, n_shuffles=50, k_max=8, seed=4)
    print(f"=== run_reachability(out-of-cone) verdict: {report['verdict']} "
          f"(held-out z={report['held_out'].z:.1f}) ===")

    # (5) SIGNED / bidirectional reachability — three-way orthogonal decomposition
    print("=== signed_reachability (LOF + GOF joint cone) ===")
    s_in = signed_reachability(E, d_in)                 # pure-LOF target
    s_gof = signed_reachability(E, -d_in)               # purely-negated => activation target
    tot_in = s_in.lof_fraction + s_in.gof_fraction + s_in.neither_fraction
    tot_gof = s_gof.lof_fraction + s_gof.gof_fraction + s_gof.neither_fraction
    print(f"  pure-LOF target : LOF={s_in.lof_fraction:.3f} GOF={s_in.gof_fraction:.3f} "
          f"neither={s_in.neither_fraction:.3f} sum={tot_in:.4f}")
    print(f"  negated target  : LOF={s_gof.lof_fraction:.3f} (chance) GOF={s_gof.gof_fraction:.3f} "
          f"neither={s_gof.neither_fraction:.3f} sum={tot_gof:.4f}")
    print(f"  LOF-half matches one-sided solver: "
          f"{abs(s_in.lof_cosine - r_in.reachable_cosine) < 1e-9}")

    # (6) DESIGN_EXPERIMENT convenience API — end-to-end design card
    print("=== design_experiment (one-call design card) ===")
    names = [f"P{i}" for i in range(P)]
    gnames = [f"g{j}" for j in range(G)]
    dz_in = design_experiment(E, d_in, perturbation_names=names, readout_names=gnames,
                              k_max=8, top=6, n_shuffles=30, seed=5)
    dz_out = design_experiment(E, d_out, perturbation_names=names, readout_names=gnames,
                               k_max=8, top=6, n_shuffles=30, seed=6)
    dz_geo = design_experiment(E, d_in, k_max=6, n_shuffles=0, seed=7)   # geometry-only path
    print(f"  in-cone  : {dz_in.summary()}")
    print(f"    recipe head={[r['gene'] for r in dz_in.knockdown_recipe[:4]]} "
          f"optimal_k={dz_in.optimal_k} library_k={[l['k'] for l in dz_in.library]}")
    print(f"  out-cone : {dz_out.summary()}")
    print(f"    activation recipe head={[r['gene'] for r in dz_out.activation_recipe[:4]]} "
          f"cert genes={len(dz_out.activation_certificate)}")
    print(f"  geometry-only null_calibrated={dz_geo.null_calibrated} verdict={dz_geo.verdict}")
    # (7) ANALYTIC ANISOTROPY-CORRECTED NULL must match the empirical shuffled-target null,
    #     and must correctly ELEVATE for an anisotropic (DC-biased) target. Elevation is only
    #     POSSIBLE when the dictionary reaches the uniform direction better than a random one
    #     (rho1^2 > kappa) -- i.e. when it carries a shared same-sign signature, as real
    #     effect-vector dictionaries do. The iid mixed-sign `E` above has no such signature
    #     (rho1^2 ~ kappa), so section (7) builds a dictionary WITH one to exercise both
    #     regimes honestly.
    print("=== analytic anisotropy-corrected null ===")
    E_sig = E + 1.2 * np.abs(rng.standard_normal((P, 1))) * (np.ones(G) / np.sqrt(G))[None, :] * np.sqrt(G)
    ac = rng.standard_normal(G); ac -= ac.mean()          # a fixed mean-zero direction
    ac = ac / np.linalg.norm(ac)
    uni = np.ones(G) / np.sqrt(G)
    # balanced (mean-zero) target: anisotropy ~0, analytic null ~ empirical shuffled null
    d_bal = ac.copy()
    emp_bal = shuffled_target_null(E_sig, d_bal, n_iter=200, seed=11)
    an_bal = analytic_anisotropy_null(E_sig, d_bal, n_probe=24, seed=11,
                                      observed=emp_bal.observed)
    # anisotropic target: same AC part, but a large uniform component => large DC fraction
    a_dc = 0.8
    d_ani = a_dc * uni + np.sqrt(1 - a_dc ** 2) * ac
    emp_ani = shuffled_target_null(E_sig, d_ani, n_iter=200, seed=12)
    an_ani = analytic_anisotropy_null(E_sig, d_ani, n_probe=24, seed=12,
                                      observed=emp_ani.observed)
    print(f"  balanced   a^2={an_bal.anisotropy:.3f}  analytic null={an_bal.null_mean:.3f}  "
          f"empirical null={emp_bal.null_cosines.mean():.3f}  (rho1={an_bal.rho_uniform:.3f}, kappa={an_bal.kappa:.3f})")
    print(f"  anisotropic a^2={an_ani.anisotropy:.3f}  analytic null={an_ani.null_mean:.3f}  "
          f"empirical null={emp_ani.null_cosines.mean():.3f}  (rho1={an_ani.rho_uniform:.3f})")

    # (8) DEG-WEIGHTED EVALUATION — the WMSE analog of Mejia et al. (ICML 2026).
    #     Non-breaking guarantees + a demonstration that weighting removes background dilution.
    print("=== DEG-weighted evaluation (Needles-in-the-Haystack calibration) ===")
    # (8a) weights=None and weights=uniform must reproduce the unweighted result BIT-IDENTICALLY
    r_uwt = reachability(E, d_out)
    r_none = reachability(E, d_out, weights=None)
    r_uni = reachability(E, d_out, weights=deg_weights(d_out, scheme="uniform"))
    identical_none = (r_none.reachable_cosine == r_uwt.reachable_cosine
                      and r_none.residual_norm == r_uwt.residual_norm
                      and np.array_equal(r_none.weights, r_uwt.weights))
    uniform_matches = (abs(r_uni.reachable_cosine - r_uwt.reachable_cosine) < 1e-9
                       and abs(r_uni.residual_norm - r_uwt.residual_norm) < 1e-9)
    print(f"  weights=None bit-identical to unweighted: {identical_none}")
    print(f"  weights=uniform matches unweighted (<1e-9): {uniform_matches}")

    # (8b) SIGNAL DILUTION at the METRIC level (the core Mejia et al. point, isolated from
    #      the NNLS fit). Take a target d and a prediction p that errs by the SAME amount in
    #      two scenarios: (i) error on the high-|d| DEG genes, (ii) error on the quiet
    #      background genes. An unweighted relative error cannot tell these apart well — the
    #      many background genes dilute it — whereas a DEG-weighted error AMPLIFIES the DEG
    #      error and DE-EMPHASISES the background error. This is exactly why an unweighted
    #      metric rewards a mean/null prediction when the real signal is sparse.
    rng2 = np.random.default_rng(seed + 99)
    dt = d_in.copy()
    ad = np.abs(dt)
    deg = np.argsort(-ad)[:20]                                 # the 20 DEG genes (high |d|)
    bg = np.argsort(ad)[:20]                                   # 20 quiet background genes
    emag = 3.0 * dt.std()
    def _rel(p, w=None):
        r = p - dt
        if w is None:
            return float(np.linalg.norm(r) / np.linalg.norm(dt))
        return float(np.sqrt(np.sum(w * r * r)) / np.sqrt(np.sum(w * dt * dt)))
    w_abs = deg_weights(dt, scheme="abs")
    p_deg = dt.copy(); p_deg[deg] += emag * rng2.choice([-1.0, 1.0], 20)   # error ON the DEGs
    p_bg = dt.copy(); p_bg[bg] += emag * rng2.choice([-1.0, 1.0], 20)      # error on background
    err_deg_unw, err_deg_wt = _rel(p_deg), _rel(p_deg, w_abs)
    err_bg_unw, err_bg_wt = _rel(p_bg), _rel(p_bg, w_abs)
    print(f"  metric dilution: error ON DEGs    unweighted={err_deg_unw:.3f} DEG-weighted={err_deg_wt:.3f}")
    print(f"                   error on backgrnd unweighted={err_bg_unw:.3f} DEG-weighted={err_bg_wt:.3f}")

    # (8c) positive-control ceiling + dynamic-range fraction
    pc = positive_control_ceiling(E, n_targets=30, support_size=4, noise=0.2, seed=1)
    null_out2 = shuffled_target_null(E, d_out, n_iter=200, seed=2)
    cal = calibrate_reachability(r_out.reachable_cosine, null_out2.null_cosines.mean(),
                                 pc["ceiling"], n_ceiling=pc["n_targets"])
    cal_in = calibrate_reachability(r_in.reachable_cosine, null_in.null_cosines.mean(),
                                    pc["ceiling"], n_ceiling=pc["n_targets"])
    print(f"  positive-control ceiling (interpolated duplicate): {pc['ceiling']:.3f} "
          f"[{pc['ceiling_p25']:.3f},{pc['ceiling_p75']:.3f}]")
    print(f"  dynamic-range fraction: in-cone={cal_in.dynamic_range_fraction:.3f} "
          f"out-cone={cal.dynamic_range_fraction:.3f}")

    # (8d) weighted held-out-gene validation still runs and stays honest
    ho_w = held_out_gene_validation(E, d_in, weights=deg_weights(d_in, scheme="abs"),
                                    n_shuffles=50, seed=3)
    print(f"  weighted held-out: held-out={ho_w.held_out_cosine:.3f} z={ho_w.z:.1f}")

    # knee helper sanity: a saturating curve knees before its last point
    kk = np.arange(1, 11); yy = 1 - np.exp(-kk / 1.5)
    knee_k = _knee(kk, yy)

    # assertions that must hold if the math is right
    assert r_in.reachable_cosine > 0.999, "in-cone target must be essentially fully reachable"
    assert r_in.residual_norm < 1e-3, "in-cone residual must be ~0"
    assert r_in.cert_max_violation < 1e-5, "KKT violation must be ~0 at optimum"
    assert r_out.cert_max_violation < 1e-5, "KKT violation must be ~0 at optimum"
    assert r_out.reachable_cosine < r_in.reachable_cosine, "out-of-cone must fit worse"
    assert same_order, "fast and exact spectrum must select the same generators"
    assert ho.held_out_cosine > 0.5, "in-cone target must generalise to held-out genes"
    # signed-reachability invariants
    assert s_in.lof_fraction > 0.99, "pure-LOF target must be ~fully knockdown-reachable"
    assert s_in.gof_fraction < 0.01, "pure-LOF target must need ~no activation"
    assert abs(tot_in - 1.0) < 0.02, "three-way decomposition must sum to ~1 (orthogonality)"
    assert abs(tot_gof - 1.0) < 0.02, "three-way decomposition must sum to ~1 (orthogonality)"
    assert s_gof.neither_fraction < 0.01, "negated target is fully reachable using BOTH directions"
    assert s_gof.gof_fraction > s_gof.lof_fraction, "activation must dominate for a negated target"
    assert abs(s_in.lof_cosine - r_in.reachable_cosine) < 1e-9, "LOF half must match one-sided solver"
    assert s_in.cert_max_violation < 1e-5, "signed KKT violation must be ~0 at optimum"
    # design_experiment invariants
    assert abs(dz_in.reachable_cosine - s_in.lof_cosine) < 1e-9, "design reach_cos must match signed LOF cosine"
    assert dz_in.reachable_cosine > dz_out.reachable_cosine, "in-cone target must design a stronger recipe"
    assert dz_in.null_calibrated and not dz_geo.null_calibrated, "n_shuffles flag must toggle calibration"
    assert 1 <= dz_in.optimal_k <= 8, "optimal_k must be a real knee within k_max"
    assert len(dz_in.library) == dz_in.spectrum["k"].size, "library must mirror the spectrum length"
    assert all(dz_in.library[i]["cosine_at_k"] >= dz_in.library[i-1]["cosine_at_k"] - 1e-9
               for i in range(1, len(dz_in.library))), "spectrum cosine must be non-decreasing"
    assert dz_out.n_activation_support > 0 and len(dz_out.activation_recipe) > 0, "negated target must yield an activation recipe"
    assert 1 <= knee_k <= 5, f"knee of a saturating curve must be early, got {knee_k}"
    assert dz_in.cert_max_violation < 1e-5, "design KKT violation must be ~0"
    # analytic-null invariants
    assert an_bal.anisotropy < 0.02, "mean-zero target must have ~0 DC anisotropy"
    assert an_ani.anisotropy > 0.5, "strongly DC-shifted target must be highly anisotropic"
    assert abs(an_bal.null_mean - emp_bal.null_cosines.mean()) < 0.06, \
        "analytic null must match empirical shuffled null for a balanced target"
    assert abs(an_ani.null_mean - emp_ani.null_cosines.mean()) < 0.08, \
        "analytic null must match empirical shuffled null for an anisotropic target"
    assert an_ani.null_mean > an_bal.null_mean + 0.1, \
        "anisotropy channel must ELEVATE the null for a DC-biased target"
    # DEG-weighted evaluation invariants (Needles-in-the-Haystack calibration)
    assert identical_none, "weights=None must reproduce the unweighted result BIT-IDENTICALLY"
    assert uniform_matches, "weights=uniform must match the unweighted result to <1e-9"
    # Signal dilution (metric level): the SAME-magnitude error is scored far higher when it
    # lands on the DEGs than on the background under DEG weighting, but the unweighted metric
    # barely distinguishes them — the artifact Mejia et al. (ICML 2026) identify.
    assert err_deg_wt > err_deg_unw, "DEG-weighting must AMPLIFY error that lands on the DEGs"
    assert err_bg_wt < err_bg_unw, "DEG-weighting must DE-EMPHASISE error on the background"
    assert err_deg_wt > err_bg_wt, \
        "under DEG weighting, DEG error must score worse than identical background error"
    assert 0.0 <= cal.dynamic_range_fraction <= 1.0 and 0.0 <= cal_in.dynamic_range_fraction <= 1.0, \
        "dynamic-range fraction must be in [0,1]"
    assert cal_in.dynamic_range_fraction > cal.dynamic_range_fraction, \
        "an in-cone target must calibrate ABOVE an out-of-cone one"
    assert pc["ceiling"] > 0.9, \
        "interpolated-duplicate positive control must score near-max (metric can reward truth)"
    assert np.isfinite(ho_w.z), "weighted held-out z must be finite"
    print("\nALL SELF-TESTS PASSED")


if __name__ == "__main__":
    _selftest()
