"""conformal_certificate: distribution-free false-certification control for CombiCone.

A **selective-prediction ("reject option")** layer that wraps the CombiCone cone
residual / two-bar emergence verdict in a split-conformal calibration. It does NOT
modify :mod:`combicone`; it imports it. The single thing it buys you:

    Turn the machine-precision separator + hand-picked ``floor_ratio >= 1.9`` bar
    into a rule that either **certifies emergent** or **abstains**, and carries a
    *distribution-free* bound on the false-certification rate at a target level
    ``alpha`` — under exchangeability with a calibration set of known
    non-emergent (additive) combinations.

Why this matters
----------------
The shipped certificate answers "is this residual above measurement noise?" with a
noise-injection p-value and an effect-size floor ratio, and reports a
false-positive RATE measured post-hoc on additive negative controls (0/100 in the
dossier). That is an *observed* specificity, not a *guaranteed* one, and the
``1.9`` floor threshold is a hand-picked operating point. Split-conformal
prediction upgrades this: pick a target error ``alpha`` first, calibrate a
threshold on the additive negative controls, and the realized
false-certification rate is provably ``<= alpha`` in finite samples for any
exchangeable non-emergent input, with **no distributional assumption** on the
score beyond exchangeability.

The nonconformity score
-----------------------
Two scores are supported (higher = more emergent = more nonconforming w.r.t. the
"non-emergent" null); both are read straight out of the existing geometry, nothing
new is fitted:

* ``"floor_ratio"`` (default, recommended): the noise-aware effect-size ratio
  ``residual_fraction / noise_floor`` from :func:`combicone.certify_emergence`
  (analytic null). This is *exactly the quantity the two-bar verdict thresholds at
  1.9*. It is decorrelated from raw effect magnitude, and — the property that makes
  it an ideal conformal score — a truly additive combination has residual ~ noise
  floor, so its floor ratio concentrates at ~1.0 by construction. Requires a
  per-combination noise model (the split-half ``|t1 - t2| / 2``).

* ``"residual_fraction"``: the raw cone residual fraction from
  :func:`reachability.project_cone` — pure geometry, deterministic, and needs **no
  noise model at all**. The purest distribution-free story, at the cost of leaving
  the score confounded with per-combination noise magnitude (a noisier combination
  reads higher for reasons unrelated to emergence), which weakens exchangeability
  between synthetic controls and real measurements.

The calibration set
-------------------
Additive negative controls: ``add = atom_i + atom_j`` for real single-gene atoms,
plus real per-gene measurement noise (a split-half ``|t1 - t2| / 2`` sampled from
the measured doubles). The deterministic part is a non-negative mix of two atoms,
hence *inside the cone by construction* — non-emergent by design. This is the same
negative-control construction the repo's ``certificate_dossier`` already uses to
report specificity; we reuse it as a calibration set so the certificate inherits a
*guaranteed* rate rather than an observed one.

The guarantee (honest statement)
--------------------------------
Let the calibration scores be ``s_1..s_n`` on ``n`` additive negative controls and
``s*`` the score of a test combination. The (one-class / "conformal anomaly
detection") conformal p-value is

    p(s*) = (1 + #{i : s_i >= s*}) / (n + 1).

We **certify emergent iff p(s*) <= alpha**, else **abstain**. If the test
combination is truly non-emergent and *exchangeable* with the calibration
controls, then ``p(s*)`` is super-uniform, so

    P(falsely certify) = P(p(s*) <= alpha) <= alpha,

with the exact finite-sample expectation ``floor(alpha * (n + 1)) / (n + 1)``. No
assumption on the shape of the score distribution is used — only exchangeability.

What exchangeability does and does NOT cover (read this before quoting a number):
  * It holds *by construction* among the additive negative controls (they are
    i.i.d. draws of the same generative process), so the bound is exact for the
    synthetic-additive null. That is what the empirical coverage curve verifies.
  * A *real* measured double that happens to be additive is only *approximately*
    exchangeable with the synthetic controls (its idiosyncratic noise structure and
    the exactness of the synthetic additive part differ). The FCR bound is anchored
    to the additive-NC null; on real doubles we can only *report* the
    certify/abstain split, not certify ground-truth emergence (which is unknown).
  * ``"floor_ratio"`` additionally assumes the split-half ``|t1 - t2| / 2`` is the
    true per-gene measurement SE. ``"residual_fraction"`` makes no noise assumption.

All results are **model-relative** in exactly the sense the rest of CombiCone means
it: "emergent" = "outside the non-negative cone of THESE atoms under THIS metric",
never a claim of biological impossibility.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Callable, Literal, Sequence

import numpy as np

import combicone as cc
import reachability as rx

__all__ = [
    "NonconformityScore",
    "ConformalCalibration",
    "SelectiveVerdict",
    "residual_nonconformity",
    "floor_ratio_nonconformity",
    "build_additive_controls",
    "calibrate_from_scores",
    "calibrate_from_additive_controls",
    "selective_certify",
]

_EPS = 1e-12
_SCOPE = (
    "model-relative: 'emergent' = outside the non-negative cone of the supplied "
    "effect atoms under the chosen metric; false-certification rate is guaranteed "
    "only against exchangeable non-emergent (additive) inputs, not biological truth"
)

NonconformityScore = Literal["floor_ratio", "residual_fraction"]


# --------------------------------------------------------------------------- #
# Nonconformity scores (higher = more emergent). Both read out of the existing
# CombiCone geometry; nothing is fitted here.
# --------------------------------------------------------------------------- #
def residual_nonconformity(
    cone_atoms: np.ndarray,
    measured_combo: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
) -> float:
    """Raw cone residual fraction — pure geometry, no noise model.

    ``residual_fraction`` from :func:`reachability.project_cone`: the fraction of
    the (metric-whitened) measured-combo norm that lies outside the non-negative
    cone of ``cone_atoms``. Deterministic; needs no per-combination noise estimate.
    Higher = further outside the cone = more nonconforming with the additive null.
    """
    res = rx.project_cone(
        np.asarray(cone_atoms, dtype=float),
        np.asarray(measured_combo, dtype=float),
        gene_weights=gene_weights,
    )
    return float(res.residual_fraction)


def floor_ratio_nonconformity(
    cone_atoms: np.ndarray,
    measured_combo: np.ndarray,
    noise_sd: np.ndarray | float,
    *,
    gene_weights: np.ndarray | None = None,
    method: str = "analytic",
    n_boot: int = 200,
    seed: int = 0,
) -> float:
    """Noise-aware effect-size ratio — the two-bar verdict's floor quantity.

    ``floor_ratio = residual_fraction / noise_floor`` from
    :func:`combicone.certify_emergence`. This is the exact quantity the shipped
    two-bar verdict thresholds at ``1.9``; conformal calibration replaces that
    hand-picked threshold with a distribution-free one. Defaults to the analytic
    (deterministic, conservative) null so a large calibration set is cheap. A truly
    additive combination has residual ~ noise floor, so this score concentrates at
    ~1.0 on the additive negative controls by construction.
    """
    cert = cc.certify_emergence(
        cone_atoms=np.asarray(cone_atoms, dtype=float),
        measured_combo=np.asarray(measured_combo, dtype=float),
        noise_sd=noise_sd,
        gene_weights=gene_weights,
        method=method,
        n_boot=n_boot,
        seed=seed,
    )
    return float(cert.floor_ratio)


# --------------------------------------------------------------------------- #
# Additive negative-control construction (the calibration set)
# --------------------------------------------------------------------------- #
def build_additive_controls(
    atoms: np.ndarray,
    *,
    n_controls: int,
    noise_pool: np.ndarray,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Synthesize additive negative controls (known non-emergent) + their noise.

    Each control is ``atoms[i] + atoms[j] + N(0, se)`` where ``(i, j)`` is a random
    distinct atom pair and ``se`` is a per-gene noise vector drawn (with
    replacement) from ``noise_pool`` (the real split-half ``|t1 - t2| / 2`` of the
    measured doubles). The deterministic part is a non-negative combination of two
    atoms, so it is inside the cone by construction — additive, non-emergent.

    Returns ``(controls, noises)`` with shapes ``(n_controls, n_genes)`` each;
    ``noises[k]`` is the ground-truth per-gene SE of ``controls[k]`` (needed only by
    the ``"floor_ratio"`` score; ignored by ``"residual_fraction"``).
    """
    atoms = np.asarray(atoms, dtype=float)
    noise_pool = np.asarray(noise_pool, dtype=float)
    n_atoms, n_genes = atoms.shape
    if n_atoms < 2:
        raise rx.InputError("need at least 2 atoms to form additive controls")
    if noise_pool.ndim != 2 or noise_pool.shape[1] != n_genes:
        raise rx.InputError("noise_pool must be (n_pool, n_genes) matching atoms")
    if n_controls < 1:
        raise rx.InputError("n_controls must be >= 1")
    rng = np.random.default_rng(seed)
    controls = np.empty((n_controls, n_genes), dtype=float)
    noises = np.empty((n_controls, n_genes), dtype=float)
    n_pool = noise_pool.shape[0]
    for k in range(n_controls):
        i, j = rng.choice(n_atoms, size=2, replace=False)
        se = noise_pool[rng.integers(n_pool)]
        controls[k] = atoms[i] + atoms[j] + rng.normal(0.0, se)
        noises[k] = se
    return controls, noises


# --------------------------------------------------------------------------- #
# The calibration object + selective rule
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SelectiveVerdict:
    """The selective-prediction outcome for one measured combination.

    Attributes
    ----------
    decision : str
        ``"certify_emergent"`` or ``"abstain"``.
    conformal_p : float
        One-class conformal p-value ``(1 + #{cal >= s*}) / (n_cal + 1)``. Small =
        the score is extreme relative to the additive null.
    score : float
        The nonconformity score of this combination.
    alpha : float
        Target false-certification level the decision was taken at.
    threshold : float
        Calibrated score threshold at ``alpha`` (certify iff ``score >= threshold``).
    score_name : str
        Which nonconformity score was used.
    scope : str
        Model-relative scope disclaimer.
    """

    decision: str
    conformal_p: float
    score: float
    alpha: float
    threshold: float
    score_name: str
    scope: str = _SCOPE


@dataclass(frozen=True)
class ConformalCalibration:
    """A split-conformal calibration on additive negative controls.

    Holds the sorted calibration scores and turns any test score into a conformal
    p-value / selective decision with a distribution-free false-certification bound.
    Construct via :func:`calibrate_from_additive_controls` (does the geometry) or
    :func:`calibrate_from_scores` (if you already have negative-control scores).

    Attributes
    ----------
    cal_scores : np.ndarray
        Sorted (ascending) nonconformity scores of the ``n`` additive negative
        controls.
    score_name : str
        ``"floor_ratio"`` or ``"residual_fraction"``.
    alpha : float
        Default target false-certification level.
    """

    cal_scores: np.ndarray
    score_name: NonconformityScore
    alpha: float = 0.05
    scope: str = _SCOPE
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.cal_scores.size)

    # ---- core conformal machinery ---------------------------------------- #
    def conformal_p(self, score: float | np.ndarray) -> float | np.ndarray:
        """One-class conformal p-value(s): ``(1 + #{cal >= s*}) / (n + 1)``.

        Super-uniform under exchangeability of a non-emergent test point with the
        calibration controls, so ``P(p <= alpha) <= alpha``.
        """
        s = np.asarray(score, dtype=float)
        # #{cal >= s*} via the sorted calibration scores.
        ge = self.n - np.searchsorted(self.cal_scores, s, side="left")
        p = (1.0 + ge) / (self.n + 1.0)
        return float(p) if np.ndim(score) == 0 else p

    def threshold(self, alpha: float | None = None) -> float:
        """Smallest score that gets certified at ``alpha`` (certify iff score >= it).

        The conformal rule ``p(s*) <= alpha`` is equivalent to ``s* >= tau`` where
        ``tau`` is the ``ceil((1 - alpha) * (n + 1))``-th smallest calibration
        score (the standard split-conformal quantile). Returns ``+inf`` when
        ``alpha`` is too small for the calibration size to ever certify (i.e.
        ``ceil((1 - alpha)(n + 1)) > n``), which is the honest "cannot certify at
        this alpha with this many controls" outcome.
        """
        a = self.alpha if alpha is None else float(alpha)
        if not (0.0 < a < 1.0):
            raise rx.InputError("alpha must be in (0, 1)")
        rank = int(np.ceil((1.0 - a) * (self.n + 1)))
        if rank > self.n:
            return float("inf")
        # rank-th smallest (1-based) -> index rank-1
        return float(self.cal_scores[rank - 1])

    def expected_fcr_bound(self, alpha: float | None = None) -> float:
        """Exact finite-sample expected FCR ``floor(alpha (n+1)) / (n+1)`` (<= alpha)."""
        a = self.alpha if alpha is None else float(alpha)
        return float(np.floor(a * (self.n + 1)) / (self.n + 1))

    def certify(
        self, score: float, alpha: float | None = None
    ) -> SelectiveVerdict:
        """Selective decision for a precomputed nonconformity ``score``."""
        a = self.alpha if alpha is None else float(alpha)
        p = float(self.conformal_p(score))
        tau = self.threshold(a)
        decision = "certify_emergent" if p <= a else "abstain"
        return SelectiveVerdict(
            decision=decision,
            conformal_p=p,
            score=float(score),
            alpha=a,
            threshold=tau,
            score_name=self.score_name,
        )


def calibrate_from_scores(
    cal_scores: Sequence[float],
    score_name: NonconformityScore = "floor_ratio",
    *,
    alpha: float = 0.05,
    meta: dict | None = None,
) -> ConformalCalibration:
    """Build a calibration directly from precomputed negative-control scores."""
    arr = np.sort(np.asarray(cal_scores, dtype=float))
    if arr.ndim != 1 or arr.size < 1 or not np.all(np.isfinite(arr)):
        raise rx.InputError("cal_scores must be a non-empty finite 1D vector")
    return ConformalCalibration(
        cal_scores=arr, score_name=score_name, alpha=alpha, meta=meta or {}
    )


def calibrate_from_additive_controls(
    atoms: np.ndarray,
    *,
    noise_pool: np.ndarray,
    score: NonconformityScore = "floor_ratio",
    n_controls: int = 500,
    alpha: float = 0.05,
    gene_weights: np.ndarray | None = None,
    method: str = "analytic",
    n_boot: int = 200,
    seed: int = 0,
) -> ConformalCalibration:
    """Calibrate on freshly-synthesized additive negative controls.

    Synthesizes ``n_controls`` additive negative controls from ``atoms`` (+ real
    noise sampled from ``noise_pool``), scores them with the chosen nonconformity
    score against the ``atoms`` cone, and returns a :class:`ConformalCalibration`.
    ``"floor_ratio"`` defaults to the analytic null so this stays cheap even for
    large ``n_controls``.
    """
    controls, noises = build_additive_controls(
        atoms, n_controls=n_controls, noise_pool=noise_pool, seed=seed
    )
    if score == "residual_fraction":
        cal = np.array(
            [
                residual_nonconformity(atoms, controls[k], gene_weights=gene_weights)
                for k in range(n_controls)
            ]
        )
    elif score == "floor_ratio":
        cal = np.array(
            [
                floor_ratio_nonconformity(
                    atoms,
                    controls[k],
                    noises[k],
                    gene_weights=gene_weights,
                    method=method,
                    n_boot=n_boot,
                    seed=seed + k,
                )
                for k in range(n_controls)
            ]
        )
    else:
        raise rx.InputError(
            f"score must be 'floor_ratio' or 'residual_fraction', got {score!r}"
        )
    meta = {
        "n_controls": int(n_controls),
        "score": score,
        "method": method if score == "floor_ratio" else "geometry",
        "seed": int(seed),
        "cal_mean": float(cal.mean()),
        "cal_p95": float(np.percentile(cal, 95)),
    }
    return calibrate_from_scores(cal, score, alpha=alpha, meta=meta)


def selective_certify(
    calibration: ConformalCalibration,
    cone_atoms: np.ndarray,
    measured_combo: np.ndarray,
    *,
    noise_sd: np.ndarray | float | None = None,
    alpha: float | None = None,
    gene_weights: np.ndarray | None = None,
    method: str = "analytic",
    n_boot: int = 200,
    seed: int = 0,
) -> SelectiveVerdict:
    """End-to-end selective certificate for a measured combination.

    Scores ``measured_combo`` against the ``cone_atoms`` cone with the calibration's
    nonconformity score, then applies the conformal certify/abstain rule at
    ``alpha``. ``noise_sd`` is required for the ``"floor_ratio"`` score and ignored
    for ``"residual_fraction"``.
    """
    if calibration.score_name == "residual_fraction":
        s = residual_nonconformity(cone_atoms, measured_combo, gene_weights=gene_weights)
    else:
        if noise_sd is None:
            raise rx.InputError("noise_sd is required for the 'floor_ratio' score")
        s = floor_ratio_nonconformity(
            cone_atoms,
            measured_combo,
            noise_sd,
            gene_weights=gene_weights,
            method=method,
            n_boot=n_boot,
            seed=seed,
        )
    return calibration.certify(s, alpha=alpha)


def _demo() -> None:
    """Deterministic, data-free smoke run of the conformal layer."""
    rng = np.random.default_rng(0)
    # A small 4-atom library; additive controls are pairwise sums + noise.
    atoms = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0, 0.0],
        ]
    )
    noise_pool = np.full((3, 4), 0.03)
    calib = calibrate_from_additive_controls(
        atoms, noise_pool=noise_pool, score="residual_fraction", n_controls=50, seed=0
    )
    print("n_cal:", calib.n, "score:", calib.score_name)
    print("threshold@0.1:", round(calib.threshold(0.1), 4))
    emergent = np.array([1.0, 1.0, 0.0, 3.0])  # big off-cone 4th axis
    additive = np.array([1.0, 1.0, 0.0, 0.0])  # inside the cone
    for label, combo in [("emergent", emergent), ("additive", additive)]:
        v = selective_certify(calib, atoms, combo)
        print(f"[{label}] p={v.conformal_p:.3g} -> {v.decision}")


if __name__ == "__main__":
    _demo()
