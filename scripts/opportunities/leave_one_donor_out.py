"""Leakage-free leave-one-donor-out (LODO) certification for CombiCone.

WHAT THIS IS
------------
A donor-holdout harness for multi-donor perturbation screens. For each donor
``d`` it builds the CombiCone reference **only from the OTHER donors**, then
certifies / tests donor ``d`` against that reference. Donor ``d``'s cells never
enter the reference it is scored against -- that disjointness is what makes the
estimate *leakage-free*.

It answers, per fold:

  * **Cone-transfer (singles screen)** -- for each single-gene effect measured in
    the held-out donor, is it reachable from the non-negative conic hull of the
    OTHER donors' single-gene effects (:func:`reachability.project_cone`)? This is
    the leakage-free donor-holdout version of *state reachability* -- the precise
    quantity the shipped ``results/donor_pair_transfer.json`` explicitly lists in
    its ``claim_ceiling`` as **NOT** covered (see ``CONTRAST_WITH_DONOR_PAIR``).

  * **Combination certification (combinatorial screen)** -- if the held-out donor
    has *measured combinations*, certify each against the leave-one-donor
    single-gene cone with :func:`combicone.certify_emergence` (the genuine
    CombiCone emergence certificate, run leakage-free across donors).

HOW IT DIFFERS FROM ``results/donor_pair_transfer.json``
--------------------------------------------------------
``donor_pair_transfer.json`` is a robustness sensitivity over *published,
two-donor-group summaries* with a fixed eligibility filter. Its own
``claim_ceiling`` reads:

    "fixed-four-donor, published-eligibility robustness sensitivity; not
     leakage-free donor holdout, donor-population inference, predictive utility,
     or state reachability"

i.e. it is explicitly **not** a leakage-free donor holdout, because the released
modalities are two-donor *group* summaries, so a single held-out donor cannot be
isolated. THIS harness closes exactly that gap: it consumes **individual
per-donor** effect vectors and holds one donor fully out of the reference. When
you only have the group summaries, you cannot run this -- and that is the honest
boundary, not a thing to paper over.

INPUT CONTRACT
--------------
The canonical input is a mapping ``per_donor_atoms: {donor_key -> (n_atoms,
n_genes) array}`` on a shared, aligned gene axis, with a shared ``atom_names``
(one perturbed-gene name per atom row). Build it however your screen dictates
(per-donor pseudobulk minus that donor's own control baseline is the standard
recipe; :func:`build_per_donor_atoms_from_pseudobulk` does this from an AnnData).
Optionally pass ``per_donor_combos: {donor_key -> list[ComboRecord]}`` to run the
combination-certification path.

Nothing is faked: if a donor lacks an atom, it is marked absent and skipped in
the folds where it would be needed; nothing is imputed.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field, asdict
from typing import Any, Mapping, Sequence

import numpy as np

import reachability as rc

try:
    import combicone as cc
    _HAVE_CC = True
except Exception:  # combicone optional for the pure cone-transfer path
    _HAVE_CC = False


CONTRAST_WITH_DONOR_PAIR = (
    "results/donor_pair_transfer.json is a fixed-four-donor, published-eligibility "
    "ROBUSTNESS SENSITIVITY over two-donor-GROUP summaries; its own claim_ceiling "
    "states it is NOT a leakage-free donor holdout and NOT state reachability. "
    "This harness IS a leakage-free donor holdout: the reference cone for fold d is "
    "built only from donors != d, using individual per-donor effect vectors, so "
    "donor d never leaks into the reference it is scored against."
)


# --------------------------------------------------------------------------- #
# Result containers
# --------------------------------------------------------------------------- #
@dataclass
class AtomTransferResult:
    """One held-out single-gene atom projected onto the leave-one-donor cone."""

    atom_name: str
    held_out_donor: str
    reference_donors: tuple[str, ...]
    n_cone_atoms: int
    self_gene_in_cone: bool
    geometry_status: str          # "inside_tolerance" | "outside_model_cone"
    residual_fraction: float      # unreachable fraction (0 = fully reachable)
    cosine: float                 # weighted cosine of fitted cone point to target
    kkt_violation: float


@dataclass
class ComboCertifyResult:
    """One held-out measured combination certified against the leave-one-donor cone."""

    combo_name: str
    held_out_donor: str
    reference_donors: tuple[str, ...]
    n_cone_atoms: int
    emergent: bool
    residual_fraction: float
    floor_ratio: float | None
    p_value: float | None
    verdict: str


@dataclass
class FoldResult:
    held_out_donor: str
    reference_donors: tuple[str, ...]
    mode: str                                     # "cone_transfer" | "combo_certify"
    n_scored: int
    atom_results: list[AtomTransferResult] = field(default_factory=list)
    combo_results: list[ComboCertifyResult] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "held_out_donor": self.held_out_donor,
            "reference_donors": list(self.reference_donors),
            "mode": self.mode,
            "n_scored": self.n_scored,
        }
        if self.atom_results:
            rf = np.array([a.residual_fraction for a in self.atom_results])
            cos = np.array([a.cosine for a in self.atom_results])
            outside = np.array(
                [a.geometry_status == "outside_model_cone" for a in self.atom_results]
            )
            d.update(
                median_residual_fraction=float(np.median(rf)),
                mean_residual_fraction=float(np.mean(rf)),
                median_cosine=float(np.median(cos)),
                frac_outside_cone=float(np.mean(outside)),
                max_kkt=float(max(a.kkt_violation for a in self.atom_results)),
            )
        if self.combo_results:
            em = np.array([c.emergent for c in self.combo_results])
            d.update(
                n_emergent=int(em.sum()),
                frac_emergent=float(np.mean(em)),
            )
        return d


# --------------------------------------------------------------------------- #
# Core: cone-transfer LODO (singles screen)
# --------------------------------------------------------------------------- #
def _align_reference_cone(
    per_donor_atoms: Mapping[str, np.ndarray],
    per_donor_present: Mapping[str, np.ndarray],
    atom_names: Sequence[str],
    reference_donors: Sequence[str],
    *,
    include_atom_index: int | None,
    exclude_self_gene: bool,
) -> tuple[np.ndarray, list[str]]:
    """Stack the reference donors' present atoms into one cone.

    If ``exclude_self_gene`` and ``include_atom_index`` is given, the atom for the
    *same gene* is dropped from the cone (leave-one-gene-and-donor-out), so a high
    reachability reflects cone geometry rather than trivial same-gene cross-donor
    reproducibility.
    """
    rows: list[np.ndarray] = []
    src: list[str] = []
    self_name = atom_names[include_atom_index] if include_atom_index is not None else None
    for d in reference_donors:
        atoms = per_donor_atoms[d]
        present = per_donor_present[d]
        for i, name in enumerate(atom_names):
            if not present[i]:
                continue
            if exclude_self_gene and self_name is not None and name == self_name:
                continue
            rows.append(atoms[i])
            src.append(f"{d}:{name}")
    if not rows:
        raise ValueError("reference cone is empty after alignment/exclusion")
    return np.asarray(rows, dtype=float), src


def leave_one_donor_out_cone_transfer(
    per_donor_atoms: Mapping[str, np.ndarray],
    atom_names: Sequence[str],
    held_out_donor: str,
    *,
    per_donor_present: Mapping[str, np.ndarray] | None = None,
    exclude_self_gene: bool = True,
    gene_weights: np.ndarray | None = None,
    max_atoms: int | None = None,
) -> FoldResult:
    """One LODO fold: project the held-out donor's atoms onto the OTHER donors' cone.

    Parameters
    ----------
    per_donor_atoms : {donor -> (n_atoms, n_genes)}
        Per-donor single-gene effect matrices on a shared gene axis.
    atom_names : (n_atoms,)
        Perturbed-gene name per atom row (shared across donors).
    held_out_donor : str
        The donor to score; its cells never enter the reference cone.
    per_donor_present : {donor -> (n_atoms,) bool}, optional
        Coverage mask; atoms flagged False are treated as not measured in that
        donor and skipped (never imputed). Defaults to all-present.
    exclude_self_gene : bool
        Drop the same-gene atom from the reference cone (leave-one-gene-out within
        the donor holdout). Recommended True for a genuine geometry test.
    """
    donors = list(per_donor_atoms)
    if held_out_donor not in donors:
        raise KeyError(f"{held_out_donor!r} not among donors {donors}")
    ref = [d for d in donors if d != held_out_donor]
    if not ref:
        raise ValueError("need >= 2 donors for a leave-one-donor-out fold")
    n_atoms = len(atom_names)
    if per_donor_present is None:
        per_donor_present = {d: np.ones(n_atoms, dtype=bool) for d in donors}

    held_atoms = per_donor_atoms[held_out_donor]
    held_present = per_donor_present[held_out_donor]
    score_idx = [i for i in range(n_atoms) if held_present[i]]
    if max_atoms is not None:
        score_idx = score_idx[:max_atoms]

    results: list[AtomTransferResult] = []
    for i in score_idx:
        cone, src = _align_reference_cone(
            per_donor_atoms, per_donor_present, atom_names, ref,
            include_atom_index=i, exclude_self_gene=exclude_self_gene,
        )
        target = held_atoms[i]
        self_in_cone = any(s.split(":", 1)[1] == atom_names[i] for s in src)
        try:
            pr = rc.project_cone(cone, target, gene_weights=gene_weights)
        except rc.InputError as exc:
            # honest failure: record it, do not fabricate a verdict
            results.append(AtomTransferResult(
                atom_name=str(atom_names[i]), held_out_donor=held_out_donor,
                reference_donors=tuple(ref), n_cone_atoms=cone.shape[0],
                self_gene_in_cone=self_in_cone, geometry_status=f"error:{exc}",
                residual_fraction=float("nan"), cosine=float("nan"),
                kkt_violation=float("nan"),
            ))
            continue
        results.append(AtomTransferResult(
            atom_name=str(atom_names[i]), held_out_donor=held_out_donor,
            reference_donors=tuple(ref), n_cone_atoms=cone.shape[0],
            self_gene_in_cone=self_in_cone, geometry_status=pr.geometry_status,
            residual_fraction=float(pr.residual_fraction),
            cosine=float(pr.cosine), kkt_violation=float(pr.kkt_violation),
        ))
    return FoldResult(
        held_out_donor=held_out_donor, reference_donors=tuple(ref),
        mode="cone_transfer", n_scored=len(results), atom_results=results,
    )


# --------------------------------------------------------------------------- #
# Core: combination-certification LODO (combinatorial screen)
# --------------------------------------------------------------------------- #
def leave_one_donor_out_combos(
    per_donor_atoms: Mapping[str, np.ndarray],
    atom_names: Sequence[str],
    per_donor_combos: Mapping[str, Sequence[Any]],
    held_out_donor: str,
    *,
    per_donor_present: Mapping[str, np.ndarray] | None = None,
    gene_weights: np.ndarray | None = None,
    floor_threshold: float = 1.9,
    alpha: float = 0.05,
    method: str = "analytic",
) -> FoldResult:
    """One LODO fold for a COMBINATORIAL screen: certify the held-out donor's
    measured combinations against the single-gene cone built from the OTHER donors.

    ``per_donor_combos[d]`` is a list of records each exposing ``.name``,
    ``.effect`` (n_genes,), and optionally ``.noise_sd`` -- e.g.
    :class:`screen_ingest.ComboRecord`. The two-bar verdict (significance AND
    floor-ratio >= ``floor_threshold``) is CombiCone's own; nothing new is claimed.
    """
    if not _HAVE_CC:
        raise ImportError("combicone is required for the combination-certification path")
    donors = list(per_donor_atoms)
    ref = [d for d in donors if d != held_out_donor]
    if not ref:
        raise ValueError("need >= 2 donors for a leave-one-donor-out fold")
    n_atoms = len(atom_names)
    if per_donor_present is None:
        per_donor_present = {d: np.ones(n_atoms, dtype=bool) for d in donors}

    cone, _src = _align_reference_cone(
        per_donor_atoms, per_donor_present, atom_names, ref,
        include_atom_index=None, exclude_self_gene=False,
    )
    combos = per_donor_combos.get(held_out_donor, [])
    results: list[ComboCertifyResult] = []
    for rec in combos:
        cert = cc.certify_emergence(
            cone_atoms=cone, measured_combo=np.asarray(rec.effect, dtype=float),
            noise_sd=getattr(rec, "noise_sd", None), gene_weights=gene_weights,
            method=method, floor_threshold=floor_threshold, alpha=alpha,
        )
        # EmergenceCertificate fields: unreachable_fraction, floor_ratio, p_value,
        # geometry_status, verdict (there is no boolean; derive it from the two bars).
        rf = float(getattr(cert, "unreachable_fraction", float("nan")))
        floor = getattr(cert, "floor_ratio", None)
        pval = getattr(cert, "p_value", None)
        emergent = bool(
            getattr(cert, "geometry_status", "") == "outside_model_cone"
            and pval is not None and pval < alpha
            and floor is not None and floor >= floor_threshold
        )
        results.append(ComboCertifyResult(
            combo_name=str(rec.name), held_out_donor=held_out_donor,
            reference_donors=tuple(ref), n_cone_atoms=cone.shape[0],
            emergent=emergent, residual_fraction=rf,
            floor_ratio=(float(floor) if floor is not None else None),
            p_value=(float(pval) if pval is not None else None),
            verdict=str(getattr(cert, "verdict", "emergent" if emergent else "reachable")),
        ))
    return FoldResult(
        held_out_donor=held_out_donor, reference_donors=tuple(ref),
        mode="combo_certify", n_scored=len(results), combo_results=results,
    )


# --------------------------------------------------------------------------- #
# Driver over all folds
# --------------------------------------------------------------------------- #
def run_all_folds(
    per_donor_atoms: Mapping[str, np.ndarray],
    atom_names: Sequence[str],
    *,
    per_donor_present: Mapping[str, np.ndarray] | None = None,
    per_donor_combos: Mapping[str, Sequence[Any]] | None = None,
    exclude_self_gene: bool = True,
    gene_weights: np.ndarray | None = None,
    max_atoms: int | None = None,
    floor_threshold: float = 1.9,
) -> dict[str, Any]:
    """Run every donor as the held-out fold; return a JSON-serializable report."""
    donors = list(per_donor_atoms)
    folds: list[FoldResult] = []
    for d in donors:
        if per_donor_combos is not None:
            folds.append(leave_one_donor_out_combos(
                per_donor_atoms, atom_names, per_donor_combos, d,
                per_donor_present=per_donor_present, gene_weights=gene_weights,
                floor_threshold=floor_threshold,
            ))
        else:
            folds.append(leave_one_donor_out_cone_transfer(
                per_donor_atoms, atom_names, d,
                per_donor_present=per_donor_present,
                exclude_self_gene=exclude_self_gene,
                gene_weights=gene_weights, max_atoms=max_atoms,
            ))
    report = {
        "harness": "leave_one_donor_out",
        "leakage_free": True,
        "mode": folds[0].mode if folds else None,
        "donors": donors,
        "exclude_self_gene": exclude_self_gene,
        "contrast_with_donor_pair_transfer": CONTRAST_WITH_DONOR_PAIR,
        "per_fold": [f.summary() for f in folds],
        "per_fold_detail": [
            {
                "held_out_donor": f.held_out_donor,
                "atoms": [asdict(a) for a in f.atom_results],
                "combos": [asdict(c) for c in f.combo_results],
            }
            for f in folds
        ],
    }
    return report


# --------------------------------------------------------------------------- #
# Builder: per-donor atoms from a pseudobulk AnnData
# --------------------------------------------------------------------------- #
def build_per_donor_atoms_from_pseudobulk(
    adata: Any,
    *,
    donor_key: str = "donor_id",
    condition_key: str | None = "culture_condition",
    condition: str | None = None,
    guide_type_key: str | None = "guide_type",
    control_value: str = "non-targeting",
    gene_name_key: str = "perturbed_gene_name",
    keep_key: str | None = "keep_for_DE",
    normalize: str = "log1p_cp10k",
    min_donors: int = 2,
) -> dict[str, Any]:
    """Construct ``per_donor_atoms`` / ``atom_names`` / ``per_donor_present`` from a
    pseudobulk AnnData whose ``.obs`` carries donor, condition, guide-type and
    perturbed-gene columns (the schema of ``GWCD4i.pseudobulk_merged.h5ad``).

    For each donor: normalize every pseudobulk, average within (donor, gene) over
    guides, and subtract that donor's own control (NTC) mean -> the donor's
    single-gene effect atoms. No cross-donor pooling anywhere. Returns a dict ready
    to splat into :func:`run_all_folds`.
    """
    import scipy.sparse as sp

    obs = adata.obs
    X = adata.X
    def _dense(rows):
        m = X[rows]
        return np.asarray(m.todense()) if sp.issparse(m) else np.asarray(m)

    def _norm(mat):
        if normalize == "log1p_cp10k":
            tot = mat.sum(axis=1, keepdims=True)
            tot[tot == 0] = 1.0
            return np.log1p(mat / tot * 1e4)
        return mat  # already normalized

    mask = np.ones(adata.n_obs, dtype=bool)
    if condition_key is not None and condition is not None:
        mask &= (obs[condition_key].astype(str).values == condition)
    if keep_key is not None and keep_key in obs:
        mask &= obs[keep_key].astype(bool).values

    donors = sorted(set(obs[donor_key].astype(str).values[mask]))
    gene_names_all = obs[gene_name_key].astype(str).values
    gtype = obs[guide_type_key].astype(str).values if guide_type_key else None

    # union of perturbed genes (targeting) present in >= min_donors
    from collections import defaultdict
    gene_donor = defaultdict(set)
    for i in np.where(mask)[0]:
        if gtype is not None and gtype[i] == control_value:
            continue
        gene_donor[gene_names_all[i]].add(str(obs[donor_key].values[i]))
    atom_names = sorted(g for g, ds in gene_donor.items() if len(ds) >= min_donors)
    name_idx = {g: i for i, g in enumerate(atom_names)}
    NG = adata.n_vars

    per_donor_atoms = {d: np.zeros((len(atom_names), NG)) for d in donors}
    per_donor_present = {d: np.zeros(len(atom_names), dtype=bool) for d in donors}
    for d in donors:
        dmask = mask & (obs[donor_key].astype(str).values == d)
        # control baseline
        if gtype is not None:
            ctrl_rows = np.where(dmask & (gtype == control_value))[0]
        else:
            ctrl_rows = np.array([], dtype=int)
        if ctrl_rows.size == 0:
            raise ValueError(f"donor {d} has no control ({control_value}) pseudobulks")
        base = _norm(_dense(ctrl_rows)).mean(axis=0)
        for g, gi in name_idx.items():
            grows = np.where(dmask & (gene_names_all == g) &
                             ((gtype != control_value) if gtype is not None else True))[0]
            if grows.size == 0:
                continue
            eff = _norm(_dense(grows)).mean(axis=0) - base
            per_donor_atoms[d][gi] = eff
            per_donor_present[d][gi] = True
    return {
        "per_donor_atoms": per_donor_atoms,
        "atom_names": atom_names,
        "per_donor_present": per_donor_present,
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _load_substrate_npz(path: str) -> dict[str, Any]:
    """Load a per-donor substrate .npz (keys atoms_<D>, coverage_<D>, panel_genes)."""
    z = np.load(path, allow_pickle=True)
    donors = sorted({k.split("_", 1)[1] for k in z.files if k.startswith("atoms_")})
    per_donor_atoms = {d: z[f"atoms_{d}"] for d in donors}
    per_donor_present = (
        {d: z[f"coverage_{d}"] for d in donors}
        if all(f"coverage_{d}" in z.files for d in donors) else None
    )
    atom_names = list(z["panel_genes"]) if "panel_genes" in z.files else \
        [f"atom_{i}" for i in range(next(iter(per_donor_atoms.values())).shape[0])]
    return {"per_donor_atoms": per_donor_atoms, "atom_names": atom_names,
            "per_donor_present": per_donor_present}


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--substrate", required=True,
                    help="per-donor substrate .npz (atoms_<D>, coverage_<D>, panel_genes)")
    ap.add_argument("--held-out-donor", default=None,
                    help="run a single fold for this donor key (default: all folds)")
    ap.add_argument("--no-exclude-self-gene", action="store_true",
                    help="keep same-gene atoms in the reference cone (measures "
                         "same-gene cross-donor reproducibility instead of geometry)")
    ap.add_argument("--max-atoms", type=int, default=None)
    ap.add_argument("--out", default="results/tcell_lodo_result.json")
    args = ap.parse_args(argv)

    sub = _load_substrate_npz(args.substrate)
    kwargs = dict(
        per_donor_present=sub["per_donor_present"],
        exclude_self_gene=not args.no_exclude_self_gene,
        max_atoms=args.max_atoms,
    )
    if args.held_out_donor:
        fold = leave_one_donor_out_cone_transfer(
            sub["per_donor_atoms"], sub["atom_names"], args.held_out_donor, **kwargs
        )
        report = {"harness": "leave_one_donor_out", "leakage_free": True,
                  "mode": "cone_transfer",
                  "contrast_with_donor_pair_transfer": CONTRAST_WITH_DONOR_PAIR,
                  "per_fold": [fold.summary()],
                  "per_fold_detail": [{"held_out_donor": fold.held_out_donor,
                                       "atoms": [asdict(a) for a in fold.atom_results],
                                       "combos": []}]}
    else:
        report = run_all_folds(sub["per_donor_atoms"], sub["atom_names"], **kwargs)

    with open(args.out, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"wrote {args.out}")
    for fs in report["per_fold"]:
        print(json.dumps(fs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
