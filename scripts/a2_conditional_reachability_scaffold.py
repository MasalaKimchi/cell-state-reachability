"""A2 — Conditional (subtype-stratified) reachability.  SCAFFOLD ONLY.

STATUS: blocked in-sandbox. This is the code path and expected-output schema for the
conditional-reachability analysis; it does NOT run here because the raw single-cell count
matrix is not available in the repo (only the DE-statistics summary GWCD4i.DE_stats.h5ad
is present, whose rows are already <gene>_<condition> perturbations with the cells
aggregated away — see CAUSAL.md agenda item A2). To execute, re-pull the raw counts
from the CZI Virtual Cell Platform source and point RAW_COUNTS at them.

WHAT IT DOES (once data is available):
The published effect matrix E[g,:] averages over a heterogeneous CD4+ population (naive
vs memory, cell-cycle phase, incipient Treg). A population-average counterfactual can
hold for no individual subtype, and the minimal knockdown recipe may differ by subtype.
A2 stratifies effect ESTIMATION by a baseline covariate, re-derives a stratified E per
stratum, re-solves the oracle, and reports (a) which verdicts are stratum-stable and
(b) whether the recipe reorders across strata.

This is CATE (conditional average treatment effect) estimation fused with the cone:
instead of one ATE per perturbation, a per-stratum ATE, hence a per-stratum reachability
verdict.
"""
import os, sys, json
import numpy as np

# ---- data requirement (unmet in-sandbox) -----------------------------------
RAW_COUNTS = os.environ.get("A2_RAW_COUNTS", "")   # AnnData .h5ad: cells x genes, RAW
GUIDE_CALLS = os.environ.get("A2_GUIDE_CALLS", "") # per-cell guide assignment
if not (RAW_COUNTS and os.path.exists(RAW_COUNTS)):
    raise SystemExit(
        "A2 requires the raw single-cell count matrix (cells x genes) + per-cell guide "
        "calls, which are not in the repo. Set A2_RAW_COUNTS / A2_GUIDE_CALLS to the CZI "
        "VCP source and re-run. This scaffold intentionally stops here.")

import os as _os
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_os.chdir(_REPO_ROOT); sys.path.insert(0, _REPO_ROOT)
import reachability as R
import scanpy as sc

# ---- 1. baseline covariate -> strata ----------------------------------------
def assign_strata(adata):
    """Return a categorical stratum label per cell from BASELINE (pre-perturbation) state.
    Candidate covariates, in order of preference:
      * cell-cycle phase   (sc.tl.score_genes_cell_cycle on control cells)
      * naive/memory score (CCR7/SELL/IL7R hi = naive; effector/memory markers = memory)
      * baseline polarization score (project control cells onto the Th1/Th2 axis)
    Coarse 2-3 strata are sufficient; the point is heterogeneity, not resolution.
    """
    raise NotImplementedError("choose covariate; see docstring")

# ---- 2. per-stratum effect matrix -------------------------------------------
def stratified_effects(adata, guide_calls, stratum):
    """E_s[g,:] = mean(X | guide g, stratum s) - mean(X | control, stratum s), on the SAME
    log-normalized readout as the pooled E (build_effect_matrices.py convention). Returns
    {stratum: (E_s, gene_axis)}."""
    raise NotImplementedError

# ---- 3. per-stratum verdict + recipe ----------------------------------------
def conditional_reachability(strat_effects, targets, hvg_mask):
    rows = []
    for s, (E_s, gaxis) in strat_effects.items():
        for tname, d in targets.items():
            res = R.reachability(E_s, d, hvg_mask=hvg_mask)
            spec = R.reachability_spectrum(E_s, d, k_max=12, hvg_mask=hvg_mask, refit_full=True)
            rows.append(dict(stratum=s, target=tname,
                             reachable_cosine=float(res.reachable_cosine),
                             residual_norm=float(res.residual_norm),
                             recipe=[str(gaxis[i]) for i in spec["order"][:12]]))
    return rows

# ---- expected output schema (documented so downstream can be built now) ------
EXPECTED_OUTPUT_SCHEMA = {
    "file": "a2_conditional_reachability.csv",
    "columns": ["stratum", "target", "reachable_cosine", "residual_norm",
                "verdict_grade", "recipe_top12", "recipe_jaccard_vs_pooled",
                "cosine_delta_vs_pooled"],
    "headline_metrics": {
        "stratum_stability": "fraction of (target) verdicts whose grade is invariant across strata",
        "recipe_reordering": "mean 1 - Jaccard(recipe_s, recipe_pooled) over strata",
    },
    "interpretation": "A stratum-stable verdict with a reordering recipe => the target is "
                      "reachable in every subtype but via subtype-specific levers (actionable "
                      "for delivery). A stratum-UNstable verdict => the pooled ATE hides a "
                      "subtype where the target is provably outside (a hidden-heterogeneity "
                      "caveat the pooled oracle cannot see).",
}

if __name__ == "__main__":
    json.dump(EXPECTED_OUTPUT_SCHEMA, open("results/a2_expected_output_schema.json", "w"), indent=1)
    print("A2 scaffold — data unmet; wrote results/a2_expected_output_schema.json")
