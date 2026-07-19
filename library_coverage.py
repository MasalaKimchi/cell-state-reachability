"""Library-level coverage, redundancy, and acquisition geometry.

A thin layer over ``reachability.project_cone``. Where the core answers one
question — *is this one target reachable from this library?* — this layer
aggregates the same projection over a **catalog** of targets to answer three
portfolio questions about a whole perturbation **library**:

1. **Coverage.**   What fraction of the target catalog is strictly inside the
   measured cone, and what fraction clears an optional soft cosine bar?
   (`coverage_report`)
2. **Redundancy.** Which atoms carry no weight — i.e. removing them does not
   shrink the reachable cone over the catalog? (`atom_redundancy`)
3. **Acquisition.** Which supplied, already-measured candidate effect atom is
   predicted to improve mean catalog cosine the most? (`rank_acquisitions`)

Nothing here is new geometry. It emits no biological verdict, recipe, dose, or
candidate ranking of *genes*; it ranks supplied *perturbation effect atoms*
using the certificates ``project_cone`` already emits. The module has no notion
of whether an atom has or has not been measured. Callers must not describe a
ranking over supplied effect vectors as a recommendation for an unmeasured
experiment.

The acquisition ranking has two paths that this module deliberately keeps
distinct so one can validate the other:

* **certificate (cheap).** For every catalog target that is *outside* the
  current cone, ``project_cone`` returns a ``residual`` — the component of the
  target no non-negative combination of the library can produce — and certifies
  it with a ``dual_separator``: every current atom projects ``<= 0`` onto that
  residual direction while the target projects ``> 0``. A candidate atom that
  aligns with the residual is therefore one that crosses the separating
  hyperplane and can begin to close the gap. Scoring a candidate is one dot
  product per outside-target. Cost: ``O(|catalog|)`` projections total, then
  ``O(|candidates| * |catalog|)`` dot products.

* **realized (expensive, retrospective comparator).** Actually add each supplied
  candidate effect atom to the library and recompute the whole catalog. The
  primary realized objective is explicitly **mean-cosine gain**. Strict
  inside-cone fraction gain and, when requested, soft cosine-threshold coverage
  gain are returned separately. Cost: ``O(|candidates| * |catalog|)``
  *projections*. This compares two computations on already-observed vectors; it
  is not prospective experimental ground truth.

Matrices use ``(atoms, features)`` orientation throughout, matching
``reachability``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from reachability import InputError, ProjectionResult, project_cone


# --------------------------------------------------------------------------- #
# 1. Coverage
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CoverageReport:
    """Strict cone inclusion and optional soft alignment for one catalog.

    ``strict_inside_cone`` is true only when the core reports
    ``inside_tolerance``. ``soft_covered`` is a distinct application-chosen
    cosine rule and is ``None`` unless a threshold was supplied.

    The compatibility properties ``reachable_fraction``, ``reachable``, and
    ``mean_reachable_cosine`` preserve the former active-rule behavior: they
    refer to soft coverage when a soft threshold is present and strict cone
    inclusion otherwise. New code should use the explicit fields.
    """

    strict_inside_cone_fraction: float
    soft_coverage_fraction: float | None
    soft_cosine_threshold: float | None
    mean_cosine: float
    mean_strict_inside_cone_cosine: float
    mean_soft_covered_cosine: float | None
    mean_residual_fraction: float
    cosines: np.ndarray
    strict_inside_cone: np.ndarray
    soft_covered: np.ndarray | None
    residual_fractions: np.ndarray
    statuses: list[str]
    results: list[ProjectionResult] = field(repr=False)

    @property
    def reachable_fraction(self) -> float:
        """Compatibility alias for the formerly active reachability rule."""
        if self.soft_coverage_fraction is not None:
            return self.soft_coverage_fraction
        return self.strict_inside_cone_fraction

    @property
    def reachable(self) -> np.ndarray:
        """Compatibility alias for the formerly active reachability mask."""
        if self.soft_covered is not None:
            return self.soft_covered
        return self.strict_inside_cone

    @property
    def mean_reachable_cosine(self) -> float:
        """Compatibility alias for the formerly active reached-target mean."""
        if self.mean_soft_covered_cosine is not None:
            return self.mean_soft_covered_cosine
        return self.mean_strict_inside_cone_cosine


def _resolve_soft_cosine_threshold(
    soft_cosine_threshold: float | None,
    reach_cosine: float | None,
) -> float | None:
    """Resolve the explicit threshold and its backwards-compatible alias."""
    if soft_cosine_threshold is not None and reach_cosine is not None:
        raise InputError(
            "use soft_cosine_threshold or the compatibility alias reach_cosine, not both"
        )
    value = soft_cosine_threshold if soft_cosine_threshold is not None else reach_cosine
    if value is None:
        return None
    try:
        threshold = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError("soft cosine threshold must be a finite scalar in [-1, 1]") from exc
    if not np.isfinite(threshold) or not -1.0 <= threshold <= 1.0:
        raise InputError("soft cosine threshold must be a finite scalar in [-1, 1]")
    return threshold


def _project_all(
    effects: np.ndarray,
    targets: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
) -> list[ProjectionResult]:
    """Project every row of ``targets`` onto the cone of ``effects``.

    This helper computes geometry only. Soft coverage is derived separately so
    it cannot change or be confused with the core's strict status.
    """
    targets = np.asarray(targets, dtype=float)
    if targets.ndim != 2:
        raise InputError("targets must be 2D (catalog_size, features)")
    if targets.shape[0] == 0 or targets.shape[1] == 0:
        raise InputError("target catalog and feature axis must be non-empty")
    return [
        project_cone(effects, t, gene_weights=gene_weights)
        for t in targets
    ]


def coverage_report(
    effects: np.ndarray,
    targets: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
    soft_cosine_threshold: float | None = None,
    reach_cosine: float | None = None,
) -> CoverageReport:
    """Summarize how well a library covers a target catalog.

    Parameters
    ----------
    effects : (atoms, features) array   the perturbation library.
    targets : (catalog, features) array  the target directions to cover.
    soft_cosine_threshold : optional float
        If set, separately report whether each best cosine clears this soft,
        application-chosen bar. This never changes strict cone inclusion.
    reach_cosine : optional float
        Backwards-compatible alias for ``soft_cosine_threshold``. Supplying
        both names fails closed.
    """
    threshold = _resolve_soft_cosine_threshold(soft_cosine_threshold, reach_cosine)
    results = _project_all(effects, targets, gene_weights=gene_weights)
    cosines = np.array([r.cosine for r in results])
    residual_fractions = np.array([r.residual_fraction for r in results])
    statuses = [r.geometry_status for r in results]
    strict_inside = np.array([s == "inside_tolerance" for s in statuses])
    soft_covered = None if threshold is None else cosines >= threshold
    strict_cosines = cosines[strict_inside]
    soft_cosines = None if soft_covered is None else cosines[soft_covered]
    return CoverageReport(
        strict_inside_cone_fraction=float(strict_inside.mean()),
        soft_coverage_fraction=(
            None if soft_covered is None else float(soft_covered.mean())
        ),
        soft_cosine_threshold=threshold,
        mean_cosine=float(cosines.mean()),
        mean_strict_inside_cone_cosine=(
            float(strict_cosines.mean()) if strict_cosines.size else 0.0
        ),
        mean_soft_covered_cosine=(
            None
            if soft_cosines is None
            else float(soft_cosines.mean()) if soft_cosines.size else 0.0
        ),
        mean_residual_fraction=float(residual_fractions.mean()),
        cosines=cosines,
        strict_inside_cone=strict_inside,
        soft_covered=soft_covered,
        residual_fractions=residual_fractions,
        statuses=statuses,
        results=results,
    )


# --------------------------------------------------------------------------- #
# 2. Redundancy
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class RedundancyReport:
    """Per-atom leave-one-out changes under explicitly named objectives."""

    marginal_mean_cosine_loss: np.ndarray
    marginal_strict_inside_cone_fraction_loss: np.ndarray
    marginal_soft_coverage_fraction_loss: np.ndarray | None
    full_mean_cosine: float
    full_strict_inside_cone_fraction: float
    full_soft_coverage_fraction: float | None
    soft_cosine_threshold: float | None

    @property
    def marginal_cosine_loss(self) -> np.ndarray:
        """Compatibility alias for marginal mean-cosine loss."""
        return self.marginal_mean_cosine_loss

    @property
    def marginal_reach_loss(self) -> np.ndarray:
        """Compatibility alias for the formerly active coverage rule."""
        if self.marginal_soft_coverage_fraction_loss is not None:
            return self.marginal_soft_coverage_fraction_loss
        return self.marginal_strict_inside_cone_fraction_loss

    @property
    def full_reachable_fraction(self) -> float:
        """Compatibility alias for the formerly active coverage rule."""
        if self.full_soft_coverage_fraction is not None:
            return self.full_soft_coverage_fraction
        return self.full_strict_inside_cone_fraction


def atom_redundancy(
    effects: np.ndarray,
    targets: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
    soft_cosine_threshold: float | None = None,
    reach_cosine: float | None = None,
) -> RedundancyReport:
    """Leave-one-atom-out marginal contribution to coverage.

    An atom whose removal changes none of these descriptive metrics adds no
    catalog-level distinction under the declared rules. This is not an
    experimental drop recommendation. It recomputes the catalog for every
    held-out, already-measured atom: ``O(atoms * catalog)`` projections.
    """
    effects = np.asarray(effects, dtype=float)
    if effects.ndim != 2 or effects.shape[0] < 2:
        raise InputError("atom redundancy requires at least two effect atoms")
    threshold = _resolve_soft_cosine_threshold(soft_cosine_threshold, reach_cosine)
    full = coverage_report(
        effects,
        targets,
        gene_weights=gene_weights,
        soft_cosine_threshold=threshold,
    )
    n_atoms = effects.shape[0]
    cos_loss = np.zeros(n_atoms)
    strict_loss = np.zeros(n_atoms)
    soft_loss = None if threshold is None else np.zeros(n_atoms)
    for i in range(n_atoms):
        keep = np.arange(n_atoms) != i
        sub = coverage_report(
            effects[keep],
            targets,
            gene_weights=gene_weights,
            soft_cosine_threshold=threshold,
        )
        cos_loss[i] = full.mean_cosine - sub.mean_cosine
        strict_loss[i] = (
            full.strict_inside_cone_fraction - sub.strict_inside_cone_fraction
        )
        if soft_loss is not None:
            if full.soft_coverage_fraction is None or sub.soft_coverage_fraction is None:
                raise RuntimeError("soft coverage unexpectedly missing from redundancy pass")
            soft_loss[i] = full.soft_coverage_fraction - sub.soft_coverage_fraction
    return RedundancyReport(
        marginal_mean_cosine_loss=cos_loss,
        marginal_strict_inside_cone_fraction_loss=strict_loss,
        marginal_soft_coverage_fraction_loss=soft_loss,
        full_mean_cosine=full.mean_cosine,
        full_strict_inside_cone_fraction=full.strict_inside_cone_fraction,
        full_soft_coverage_fraction=full.soft_coverage_fraction,
        soft_cosine_threshold=threshold,
    )


# --------------------------------------------------------------------------- #
# 3. Acquisition
# --------------------------------------------------------------------------- #
def gap_directions(
    effects: np.ndarray,
    targets: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
    report: CoverageReport | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the unreachable residual directions and their target weights.

    For every target strictly outside the cone, the core emits a certified dual
    separator. Its unit normal is a direction that every current atom scores
    non-positively against while the target scores positively. These normals,
    weighted by how much of each target is unreachable (``residual_fraction``),
    are what a candidate atom is scored against. Using the separator itself is
    essential when ``gene_weights`` is supplied: the raw residual is not the
    separating normal in a weighted metric.

    Returns
    -------
    directions : (n_outside, features)  unit dual-separator normals.
    weights    : (n_outside,)           residual_fraction per outside target.
    """
    rep = report or coverage_report(effects, targets, gene_weights=gene_weights)
    dirs = []
    wts = []
    for r in rep.results:
        if r.geometry_status != "outside_model_cone":
            continue
        if r.dual_separator is None:  # fail closed if the core did not certify one
            raise RuntimeError("outside_model_cone result is missing a dual separator")
        separator = np.asarray(r.dual_separator, dtype=float)
        n = np.linalg.norm(separator)
        if n == 0.0:
            raise RuntimeError("outside_model_cone result has a zero dual separator")
        dirs.append(separator / n)
        wts.append(r.residual_fraction)
    if not dirs:
        features = np.asarray(targets, dtype=float).shape[1]
        return np.empty((0, features)), np.empty((0,))
    return np.vstack(dirs), np.array(wts)


@dataclass(frozen=True)
class AcquisitionRanking:
    """Supplied effect atoms ranked against explicit retrospective objectives."""

    certificate_order: np.ndarray
    certificate_score: np.ndarray
    realized_mean_cosine_gain: np.ndarray | None
    realized_mean_cosine_order: np.ndarray | None
    realized_strict_inside_cone_fraction_gain: np.ndarray | None
    realized_strict_inside_cone_fraction_order: np.ndarray | None
    realized_soft_coverage_fraction_gain: np.ndarray | None
    realized_soft_coverage_fraction_order: np.ndarray | None
    soft_cosine_threshold: float | None

    @property
    def order(self) -> np.ndarray:
        """Compatibility alias; always the deployable certificate order."""
        return self.certificate_order

    @property
    def realized_gain(self) -> np.ndarray | None:
        """Compatibility alias for realized mean-cosine gain."""
        return self.realized_mean_cosine_gain

    @property
    def realized_order(self) -> np.ndarray | None:
        """Compatibility alias for realized mean-cosine order."""
        return self.realized_mean_cosine_order


def rank_acquisitions(
    effects: np.ndarray,
    targets: np.ndarray,
    candidates: np.ndarray,
    *,
    gene_weights: np.ndarray | None = None,
    realized: bool = False,
    soft_cosine_threshold: float | None = None,
    reach_cosine: float | None = None,
) -> AcquisitionRanking:
    """Rank supplied, already-measured effect atoms by a gap certificate.

    certificate score (always computed, cheap): for each candidate c, sum over
    outside targets of ``max(0, <c_unit, gap_dir>) * residual_fraction``. A
    candidate scores high when it points into the directions the current library
    cannot reach, across many under-covered targets. Cost beyond the one
    coverage pass is ``O(|candidates| * |outside targets|)`` dot products.

    realized mean-cosine gain (optional, expensive retrospective comparator):
    for each supplied candidate effect vector, append it and recompute the mean
    catalog cosine. Strict inside-cone fraction gain and optional soft coverage
    fraction gain are returned as separate objectives. This validates a
    computational approximation on observed vectors; it is not experimental
    ground truth and does not imply that candidates are unmeasured.
    """
    effects = np.asarray(effects, dtype=float)
    candidates = np.asarray(candidates, dtype=float)
    if candidates.ndim != 2:
        raise InputError("candidates must be 2D (n_candidates, features)")
    if candidates.shape[0] == 0 or candidates.shape[1] == 0:
        raise InputError("candidate pool and feature axis must be non-empty")
    if effects.ndim != 2 or effects.shape[1] != candidates.shape[1]:
        raise InputError("effects and candidates must have aligned feature axes")
    if not np.all(np.isfinite(candidates)):
        raise InputError("candidates must contain only finite values")

    threshold = _resolve_soft_cosine_threshold(soft_cosine_threshold, reach_cosine)

    base = coverage_report(
        effects,
        targets,
        gene_weights=gene_weights,
        soft_cosine_threshold=threshold,
    )
    dirs, wts = gap_directions(effects, targets, gene_weights=gene_weights, report=base)

    n_cand = candidates.shape[0]
    cert = np.zeros(n_cand)
    if dirs.shape[0] > 0:
        cand_norms = np.linalg.norm(candidates, axis=1)
        safe = cand_norms > 0
        unit_cand = np.zeros_like(candidates)
        unit_cand[safe] = candidates[safe] / cand_norms[safe, None]
        # (n_cand, n_outside) alignment of each candidate with each gap direction
        align = unit_cand @ dirs.T
        cert = (np.maximum(align, 0.0) * wts[None, :]).sum(axis=1)

    realized_mean_cosine_gain = None
    realized_strict_gain = None
    realized_soft_gain = None
    if realized:
        realized_mean_cosine_gain = np.zeros(n_cand)
        realized_strict_gain = np.zeros(n_cand)
        realized_soft_gain = None if threshold is None else np.zeros(n_cand)
        base_mean = base.mean_cosine
        for j in range(n_cand):
            augmented = np.vstack([effects, candidates[j][None, :]])
            aug = coverage_report(
                augmented,
                targets,
                gene_weights=gene_weights,
                soft_cosine_threshold=threshold,
            )
            realized_mean_cosine_gain[j] = aug.mean_cosine - base_mean
            realized_strict_gain[j] = (
                aug.strict_inside_cone_fraction - base.strict_inside_cone_fraction
            )
            if realized_soft_gain is not None:
                if aug.soft_coverage_fraction is None or base.soft_coverage_fraction is None:
                    raise RuntimeError("soft coverage unexpectedly missing from realized pass")
                realized_soft_gain[j] = (
                    aug.soft_coverage_fraction - base.soft_coverage_fraction
                )

    # The compatibility ``order`` is always the certificate ranking. Keeping the
    # realized comparator separate prevents evaluation code from accidentally
    # reporting a tautological top-pick match when ``realized=True``.
    certificate_order = np.argsort(-cert, kind="stable")
    realized_mean_cosine_order = (
        None
        if realized_mean_cosine_gain is None
        else np.argsort(-realized_mean_cosine_gain, kind="stable")
    )
    realized_strict_order = (
        None
        if realized_strict_gain is None
        else np.argsort(-realized_strict_gain, kind="stable")
    )
    realized_soft_order = (
        None
        if realized_soft_gain is None
        else np.argsort(-realized_soft_gain, kind="stable")
    )
    return AcquisitionRanking(
        certificate_order=certificate_order,
        certificate_score=cert,
        realized_mean_cosine_gain=realized_mean_cosine_gain,
        realized_mean_cosine_order=realized_mean_cosine_order,
        realized_strict_inside_cone_fraction_gain=realized_strict_gain,
        realized_strict_inside_cone_fraction_order=realized_strict_order,
        realized_soft_coverage_fraction_gain=realized_soft_gain,
        realized_soft_coverage_fraction_order=realized_soft_order,
        soft_cosine_threshold=threshold,
    )
