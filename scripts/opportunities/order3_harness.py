"""order3_harness: order-aware (k=3) certification benchmark for CombiCone.

A thin, reusable harness over :mod:`combicone` / :mod:`reachability` that
answers the two order-3 questions the repo raises but validates only on
synthetic data (the Norman screen measures no triples):

1. **Certified-set monotonicity (I3).** The two-bar certified-emergent set
   *shrinks monotonically* as the reference cone is enriched
   (singles -> singles + doubles). The richer cone's certified set is a strict
   subset of the sparser one, and every triple's residual fraction is
   non-increasing. This generalizes the repo's reported order-2 shrink
   (real Norman doubles: 105 -> 235 atoms, certified 40 -> 16) to order 3.

2. **Order-3 triage recovery (E1).** Running the *unchanged* two-bar verdict
   (bar-a: noise-injection p < alpha; bar-b: floor_ratio >= 1.9) with per-triple
   ``noise_sd`` discriminates planted 3-way-emergent triples from additive /
   2-way-reducible ones. Reported as AUROC and top-k precision vs a random
   baseline, for the raw residual, the noise-aware z, and the floor ratio.

SCOPE / HONESTY
---------------
Everything computed here is on ``synth_triple_screen.npz`` -- a SYNTHETIC
planted-epistasis substrate. It is a validation of the *code path and its
discrimination*, NEVER evidence about real biological 3-way epistasis. See
:func:`load_synth_triple_screen` and the module ``__main__`` banner. The
geometry is model-relative in exactly the sense :mod:`reachability` defines
("outside the non-negative cone of THESE atoms under THIS metric").

The monotonicity is not merely empirical: :func:`reachability.project_cone`
enforces the Pythagorean projection identity ``||t||^2 = ||fit||^2 + ||r||^2``
(residual orthogonal to the fitted point) as a KKT-certified invariant, so the
denominator ``||t||`` is fixed and the residual fraction ``||r|| / ||t||`` can
only decrease when atoms are added to the cone. The certified-set nesting is a
consequence, up to the effect-size bar. :func:`monotonicity_report` verifies
both the per-triple residual monotonicity and the set nesting empirically.

Dependencies: numpy + scipy only (via combicone/reachability). No sklearn; the
AUROC and rank utilities are dependency-free (fail-closed, deterministic).

Usage
-----
    python order3_harness.py --npz synth_triple_screen.npz --method analytic

or programmatically::

    import order3_harness as h
    data = h.load_synth_triple_screen("synth_triple_screen.npz")
    mono = h.monotonicity_report(data, method="analytic")
    tri  = h.triage_report(data, method="analytic")
    print(mono.n_certified_sparse, "->", mono.n_certified_rich)
    print(tri.auroc["z"])
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field, asdict
from typing import Sequence

import numpy as np

import combicone as cc
import reachability as rx

__all__ = [
    "TripleScreen",
    "CertificationTable",
    "MonotonicityReport",
    "TriageReport",
    "load_synth_triple_screen",
    "certify_triples",
    "two_bar_mask",
    "monotonicity_report",
    "triage_report",
    "auroc",
    "top_k_precision",
    "random_top_k_precision",
]

_FLOOR_THRESHOLD = 1.9
_ALPHA = 0.05


# --------------------------------------------------------------------------- #
# Substrate
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TripleScreen:
    """A parsed synthetic order-3 screen.

    Attributes mirror ``synth_triple_screen.npz``. Everything is SYNTHETIC;
    ``triple_label_emergent`` / ``triple_class`` are the planted ground truth,
    available here precisely because no real screen supplies measured triples.
    """

    singles: np.ndarray            # (n_singles, n_genes) single-gene atoms
    single_names: np.ndarray       # (n_singles,)
    double_atoms: np.ndarray       # (n_doubles, n_genes) measured double atoms
    double_pairs: np.ndarray       # (n_doubles, 2) member single indices
    triples_measured: np.ndarray   # (n_triples, n_genes) measured triple effects
    triple_members: np.ndarray     # (n_triples, 3) member single indices
    triple_label_emergent: np.ndarray  # (n_triples,) 1 = planted 3-way emergent
    triple_class: np.ndarray       # (n_triples,) additive|emergent|reducible
    noise_sd: np.ndarray           # (n_triples, n_genes) per-gene split-half SE
    tt_norm: np.ndarray            # (n_triples,) planted extra-term magnitude
    synthetic: bool = True         # ALWAYS true; this substrate is synthetic

    @property
    def cone_singles(self) -> np.ndarray:
        """The sparse reference cone: single-gene atoms only."""
        return self.singles

    @property
    def cone_singles_doubles(self) -> np.ndarray:
        """The enriched reference cone: singles + all measured doubles."""
        return np.vstack([self.singles, self.double_atoms])


def load_synth_triple_screen(path: str) -> TripleScreen:
    """Load ``synth_triple_screen.npz`` into a :class:`TripleScreen`.

    Fail-closed on the expected keys and on shape agreement across arrays.
    """
    d = np.load(path, allow_pickle=True)
    required = (
        "singles", "single_names", "double_atoms", "double_pairs",
        "triples_measured", "triple_members", "triple_label_emergent",
        "triple_class", "noise_sd",
    )
    missing = [k for k in required if k not in d.files]
    if missing:
        raise rx.InputError(f"synth screen missing keys: {missing}")

    singles = np.asarray(d["singles"], dtype=float)
    doubles = np.asarray(d["double_atoms"], dtype=float)
    triples = np.asarray(d["triples_measured"], dtype=float)
    noise_sd = np.asarray(d["noise_sd"], dtype=float)
    n_genes = singles.shape[1]
    if not (doubles.shape[1] == triples.shape[1] == noise_sd.shape[1] == n_genes):
        raise rx.InputError("gene axis mismatch among singles/doubles/triples/noise")
    if triples.shape[0] != noise_sd.shape[0]:
        raise rx.InputError("triples_measured and noise_sd disagree on n_triples")

    tt_norm = (
        np.asarray(d["tt_norm"], dtype=float)
        if "tt_norm" in d.files
        else np.full(triples.shape[0], np.nan)
    )
    return TripleScreen(
        singles=singles,
        single_names=np.asarray(d["single_names"]),
        double_atoms=doubles,
        double_pairs=np.asarray(d["double_pairs"], dtype=int),
        triples_measured=triples,
        triple_members=np.asarray(d["triple_members"], dtype=int),
        triple_label_emergent=np.asarray(d["triple_label_emergent"], dtype=int),
        triple_class=np.asarray(d["triple_class"]),
        noise_sd=noise_sd,
        tt_norm=tt_norm,
    )


# --------------------------------------------------------------------------- #
# Certification
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CertificationTable:
    """Per-triple certificate fields for one cone (order-agnostic)."""

    residual: np.ndarray       # unreachable_fraction per triple
    z: np.ndarray              # noise-aware emergence z
    p_value: np.ndarray        # bar-a noise-injection p
    floor_ratio: np.ndarray    # bar-b effect-size ratio
    geometry_status: np.ndarray  # 'outside_model_cone' / 'inside_tolerance'
    cone_size: int             # n atoms in the reference cone
    method: str


def certify_triples(
    cone_atoms: np.ndarray,
    triples: np.ndarray,
    noise_sd: np.ndarray | float,
    *,
    method: str = "analytic",
    gene_weights: np.ndarray | None = None,
    seed: int = 0,
    n_boot: int = 200,
) -> CertificationTable:
    """Certify every measured triple against a single cone of effect atoms.

    Thin loop over :func:`combicone.certify_emergence` (order-agnostic; the cone
    may hold atoms of any order). ``noise_sd`` is either a per-triple ``(n_triples,
    n_genes)`` array, a per-gene ``(n_genes,)`` vector, or a scalar.
    """
    triples = np.asarray(triples, dtype=float)
    n = triples.shape[0]
    noise_sd = np.asarray(noise_sd, dtype=float) if not np.isscalar(noise_sd) else noise_sd

    def _row_noise(i: int):
        if np.isscalar(noise_sd):
            return float(noise_sd)
        if noise_sd.ndim == 2:
            return noise_sd[i]
        return noise_sd  # per-gene vector shared across triples

    resid = np.empty(n)
    z = np.empty(n)
    p = np.empty(n)
    floor = np.empty(n)
    geom = np.empty(n, dtype=object)
    for i in range(n):
        cert = cc.certify_emergence(
            cone_atoms=cone_atoms,
            measured_combo=triples[i],
            noise_sd=_row_noise(i),
            gene_weights=gene_weights,
            method=method,
            seed=seed,
            n_boot=n_boot,
            floor_threshold=_FLOOR_THRESHOLD,
            alpha=_ALPHA,
        )
        resid[i] = cert.unreachable_fraction
        z[i] = cert.z
        p[i] = cert.p_value
        floor[i] = cert.floor_ratio
        geom[i] = cert.geometry_status
    return CertificationTable(
        residual=resid, z=z, p_value=p, floor_ratio=floor,
        geometry_status=geom, cone_size=int(np.asarray(cone_atoms).shape[0]),
        method=method,
    )


def two_bar_mask(
    table: CertificationTable, *, alpha: float = _ALPHA, floor_threshold: float = _FLOOR_THRESHOLD
) -> np.ndarray:
    """Boolean certified-emergent mask: bar-a (p < alpha) AND bar-b (floor >= thr)."""
    return (table.p_value < alpha) & (table.floor_ratio >= floor_threshold)


# --------------------------------------------------------------------------- #
# (I3) Monotonicity
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MonotonicityReport:
    """Certified-set monotonicity under cone enrichment (singles -> singles+doubles)."""

    n_triples: int
    cone_sparse_size: int
    cone_rich_size: int
    n_certified_sparse: int
    n_certified_rich: int
    is_subset: bool             # certified(rich) subset-of certified(sparse)
    n_newly_certified: int      # in rich, not in sparse (must be 0 for monotone)
    n_flipped_to_reachable: int  # certified sparse, not rich
    shrink_fraction: float
    residual_non_increasing: bool
    max_residual_increase: float  # <= tol confirms geometric monotonicity
    certified_by_class_sparse: dict
    certified_by_class_rich: dict
    flipped_by_class: dict
    method: str
    tol: float = 1e-9

    def to_dict(self) -> dict:
        return asdict(self)


def monotonicity_report(
    data: TripleScreen, *, method: str = "analytic", seed: int = 0, n_boot: int = 200
) -> MonotonicityReport:
    """Show the certified-emergent set shrinks monotonically as the cone is enriched.

    Certifies the triples against cone(singles) and cone(singles+doubles), then
    verifies: (a) certified(rich) is a subset of certified(sparse); (b) no triple
    is newly certified under the richer cone; (c) every residual fraction is
    non-increasing (the geometric guarantee).
    """
    sparse = certify_triples(
        data.cone_singles, data.triples_measured, data.noise_sd,
        method=method, seed=seed, n_boot=n_boot,
    )
    rich = certify_triples(
        data.cone_singles_doubles, data.triples_measured, data.noise_sd,
        method=method, seed=seed, n_boot=n_boot,
    )
    cs = two_bar_mask(sparse)
    cr = two_bar_mask(rich)

    delta = rich.residual - sparse.residual  # <= 0 up to fp for monotone
    classes = list(np.unique(data.triple_class))

    def _by_class(mask):
        return {c: int(np.sum(mask & (data.triple_class == c))) for c in classes}

    flipped = cs & ~cr
    return MonotonicityReport(
        n_triples=int(data.triples_measured.shape[0]),
        cone_sparse_size=sparse.cone_size,
        cone_rich_size=rich.cone_size,
        n_certified_sparse=int(cs.sum()),
        n_certified_rich=int(cr.sum()),
        is_subset=bool(np.all(cr <= cs)),
        n_newly_certified=int(np.sum(cr & ~cs)),
        n_flipped_to_reachable=int(flipped.sum()),
        shrink_fraction=float(1.0 - cr.sum() / max(cs.sum(), 1)),
        residual_non_increasing=bool(np.all(delta <= 1e-9)),
        max_residual_increase=float(delta.max()),
        certified_by_class_sparse=_by_class(cs),
        certified_by_class_rich=_by_class(cr),
        flipped_by_class=_by_class(flipped),
        method=method,
    )


# --------------------------------------------------------------------------- #
# (E1) Triage recovery
# --------------------------------------------------------------------------- #
def auroc(scores: Sequence[float], labels: Sequence[int]) -> float:
    """Area under ROC via the rank (Mann-Whitney U) identity. Ties averaged.

    Dependency-free; returns NaN if a class is absent.
    """
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)
    n1 = int((labels == 1).sum())
    n0 = int((labels == 0).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(scores, kind="stable")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1, dtype=float)
    # average tied ranks so ties contribute 0.5
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    ranks = (sums / counts)[inv]
    r_pos = ranks[labels == 1].sum()
    return float((r_pos - n1 * (n1 + 1) / 2) / (n1 * n0))


def top_k_precision(scores: Sequence[float], labels: Sequence[int], k: int) -> float:
    """Fraction of the top-k highest-scoring items that are positive."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)
    k = min(k, len(scores))
    idx = np.argsort(-scores, kind="stable")[:k]
    return float(labels[idx].mean())


def random_top_k_precision(
    labels: Sequence[int], k: int, *, reps: int = 5000, seed: int = 0
) -> float:
    """Expected top-k precision under uniformly random ranking (== base rate)."""
    labels = np.asarray(labels, dtype=int)
    rng = np.random.default_rng(seed)
    n = len(labels)
    k = min(k, n)
    return float(np.mean([labels[rng.permutation(n)[:k]].mean() for _ in range(reps)]))


@dataclass(frozen=True)
class TriageReport:
    """Order-3 triage recovery of planted emergence, on the enriched cone."""

    n_triples: int
    n_emergent: int
    base_rate: float
    cone_size: int
    auroc: dict                 # {'residual','z','floor'} -> AUROC
    top_k_precision: dict       # {k: {'z':..,'random':..,'enrichment':..}}
    two_bar_confusion: dict     # TP/FP/TN/FN + sens/spec/prec on planted label
    certified_by_class: dict
    method: str

    def to_dict(self) -> dict:
        return asdict(self)


def triage_report(
    data: TripleScreen,
    *,
    method: str = "analytic",
    cone: str = "rich",
    ks: Sequence[int] = (20, 40),
    seed: int = 0,
    n_boot: int = 200,
) -> TriageReport:
    """Measure recovery of ``triple_label_emergent`` by order-3 certification.

    ``cone`` selects the reference: ``"rich"`` (singles + doubles, the order-aware
    cone that isolates genuine 3-way structure) or ``"sparse"`` (singles only).
    Reports AUROC for raw residual / noise-aware z / floor ratio, top-k precision
    vs random, and the two-bar confusion matrix against the planted label.
    """
    cone_atoms = data.cone_singles_doubles if cone == "rich" else data.cone_singles
    table = certify_triples(
        cone_atoms, data.triples_measured, data.noise_sd,
        method=method, seed=seed, n_boot=n_boot,
    )
    labels = data.triple_label_emergent
    base = float(labels.mean())

    au = {
        "residual": auroc(table.residual, labels),
        "z": auroc(table.z, labels),
        "floor": auroc(table.floor_ratio, labels),
    }
    tk = {}
    for k in ks:
        prec_z = top_k_precision(table.z, labels, k)
        rnd = random_top_k_precision(labels, k, seed=seed)
        tk[int(k)] = {
            "z": prec_z,
            "random": rnd,
            "enrichment": float(prec_z / base) if base > 0 else float("nan"),
        }

    cert = two_bar_mask(table)
    tp = int(np.sum(cert & (labels == 1)))
    fp = int(np.sum(cert & (labels == 0)))
    fn = int(np.sum(~cert & (labels == 1)))
    tn = int(np.sum(~cert & (labels == 0)))
    confusion = {
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else float("nan"),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else float("nan"),
        "precision": float(tp / (tp + fp)) if (tp + fp) else float("nan"),
    }
    classes = list(np.unique(data.triple_class))
    by_class = {c: int(np.sum(cert & (data.triple_class == c))) for c in classes}

    return TriageReport(
        n_triples=int(len(labels)),
        n_emergent=int(labels.sum()),
        base_rate=base,
        cone_size=table.cone_size,
        auroc=au,
        top_k_precision=tk,
        two_bar_confusion=confusion,
        certified_by_class=by_class,
        method=method,
    )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _banner() -> str:
    return (
        "order3_harness: SYNTHETIC order-3 certification benchmark\n"
        "  substrate = synth_triple_screen.npz (planted epistasis; NOT real triples)\n"
        "  this validates the code path + discrimination, never real 3-way biology\n"
    )


def main(argv: Sequence[str] | None = None) -> dict:
    ap = argparse.ArgumentParser(description="Order-3 CombiCone certification harness.")
    ap.add_argument("--npz", default="synth_triple_screen.npz", help="path to synthetic triple screen")
    ap.add_argument("--method", default="analytic", choices=["analytic", "montecarlo"])
    ap.add_argument("--cone", default="rich", choices=["rich", "sparse"], help="triage reference cone")
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--json", default=None, help="optional path to write the full report JSON")
    args = ap.parse_args(argv)

    print(_banner())
    data = load_synth_triple_screen(args.npz)
    mono = monotonicity_report(data, method=args.method, seed=args.seed, n_boot=args.n_boot)
    tri = triage_report(data, method=args.method, cone=args.cone, seed=args.seed, n_boot=args.n_boot)

    print(f"[I3] MONOTONICITY ({args.method}):")
    print(f"     cone(singles, {mono.cone_sparse_size}) certified = {mono.n_certified_sparse}")
    print(f"     cone(singles+doubles, {mono.cone_rich_size}) certified = {mono.n_certified_rich}")
    print(f"     strict nesting (rich subset-of sparse) : {mono.is_subset}")
    print(f"     newly certified under rich cone        : {mono.n_newly_certified}")
    print(f"     flipped emergent->reachable            : {mono.n_flipped_to_reachable}")
    print(f"     shrink fraction                        : {mono.shrink_fraction:.1%}")
    print(f"     residual non-increasing per triple     : {mono.residual_non_increasing} "
          f"(max increase {mono.max_residual_increase:.2e})")
    print(f"     certified by class sparse -> rich      : {mono.certified_by_class_sparse} -> "
          f"{mono.certified_by_class_rich}")

    print(f"\n[E1] TRIAGE RECOVERY (cone={args.cone}, {args.method}):")
    print(f"     n_triples={tri.n_triples} emergent={tri.n_emergent} base_rate={tri.base_rate:.3f}")
    print(f"     AUROC  residual={tri.auroc['residual']:.4f}  z={tri.auroc['z']:.4f}  "
          f"floor={tri.auroc['floor']:.4f}")
    for k, v in tri.top_k_precision.items():
        print(f"     top-{k}: z-prec={v['z']:.3f}  random={v['random']:.3f}  "
              f"enrichment={v['enrichment']:.2f}x")
    c = tri.two_bar_confusion
    print(f"     two-bar vs planted: TP={c['TP']} FP={c['FP']} TN={c['TN']} FN={c['FN']} "
          f"(sens={c['sensitivity']:.3f} spec={c['specificity']:.3f} prec={c['precision']:.3f})")
    print(f"     certified by class : {tri.certified_by_class}")

    report = {"synthetic": True, "npz": args.npz, "monotonicity": mono.to_dict(),
              "triage": tri.to_dict()}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2, default=lambda o: o.tolist()
                      if isinstance(o, np.ndarray) else o)
        print(f"\nwrote {args.json}")
    return report


if __name__ == "__main__":
    main()
