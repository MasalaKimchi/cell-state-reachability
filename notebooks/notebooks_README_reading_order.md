# Notebooks — reading order

**Start here.** The nine notebooks are numbered in *build* order (how the project was assembled),
but a first-time reader should meet the **payoff before the plumbing**. This overlay reorders them by
the question you're probably asking. Work top to bottom: see *what the method tells you*, then whether
it's *real*, whether you can *trust* it, whether it *generalizes*, and finally *how to run one yourself*.
The physical file names are unchanged — this is a reading route, not a rename.

| # | Reader question | Notebook | What you get |
|---|---|---|---|
| 1 | **What can it tell me?** | [`05_target_id_showcase.ipynb`](05_target_id_showcase.ipynb) | The teaser. A disease-relevant transition walked end to end: verdict → *reach vs. must-activate* split → ranked recipe → druggability + human-genetics triage → a one-page target dossier. Read this first to see the whole service in one pass. |
| 2 | **Is it real?** | [`02_reachability_on_tier2.ipynb`](02_reachability_on_tier2.ipynb) | The headline result plus the honesty tests. Builds the effect dictionary `E` and target `d`, runs the cone fit (the verdict), then held-out-gene validation and a shuffled-target null so the reachable cosine is never reported uncalibrated. Includes the minimal recipe and the activation certificate. |
| 3 | **Can I trust it?** | [`09_causal_validation_dossier.ipynb`](09_causal_validation_dossier.ipynb) | The reviewer's notebook. Organized by a six-assumption identifying stack (unbiased ATEs, additivity, SUTVA, homogeneous effect, exclusion restriction, transportability): each assumption made explicit, the testable ones tested, and — where a test is impossible in silico — a quantified *how far the assumption must fail to flip the verdict*. |
| 4 | **Does it generalize?** | [`03_generalizability_and_impact.ipynb`](03_generalizability_and_impact.ipynb) + [`07_cross_celltype_transfer.ipynb`](07_cross_celltype_transfer.ipynb) | Two transfer tests on unchanged code. **03**: the identical method on a second dataset (Norman K562 CRISPRa) — orientation gate, held-out target, recipe, spectrum, certificate, and the additivity test only combinatorial data enables. **07**: cross-cell-type (K562 vs RPE1) — verdict is *direction transfers, recipe does not*. |
| 5 | **How do I run one?** | [`04_experimental_design_toolkit.ipynb`](04_experimental_design_toolkit.ipynb) + [`bring_your_own_target.ipynb`](bring_your_own_target.ipynb) | **04**: turns a verdict into a runnable screen — `design_experiment()` across all 12 atlas transitions with signed decomposition, calibration, ranked recipe, optimal-k library, and modality triage. **`bring_your_own_target`**: the same as a one-call service on *your* target — paste a gene signature, get a live verdict, recipe, library, and certificate back (a worked example ships in-notebook). |
| — | **Appendix / provenance** | [`01_exploratory_data_analysis.ipynb`](01_exploratory_data_analysis.ipynb), [`06_reinforcement_analyses.ipynb`](06_reinforcement_analyses.ipynb), [`08_deg_weighted_evaluation.ipynb`](08_deg_weighted_evaluation.ipynb) | The supporting record. **01**: the two-tier data inventory and a first look at `E`. **06**: the reinforcement/limitations battery (L1/L2/L4/L5) that hardens the manuscript against its own stated caveats. **08**: re-scores every verdict under a DEG-weighted readout to check the geometry survives emphasis on the most-perturbed genes. |

### Notes

- **The method lives in [`../reachability.py`](../reachability.py).** Every notebook calls into it; none
  re-implements the geometry. The batch drivers (`../run_atlas.py`, `../run_nulls.py`,
  `../run_bootstrap.py`, `../run_a1_sensitivity.py`, `../run_iv_compliance.py`) reproduce the atlas, nulls,
  CIs, and the causal-layer result files headlessly.
- **Reading order ≠ dependency order.** A few notebooks read frozen result files written by their batch
  drivers (notably 09, which depends on the `a{1,3,4,5,6}_*.csv` and `iv_compliance_verdicts.csv` files).
  Notebooks 06 and 07 are independent reinforcement passes and each runs standalone.
- **Heavy inputs are gitignored.** The 809 MB atlas matrix and the large `.npz` bundles are not in the
  repo; notebooks work from small cached tables and bundles (e.g. `cache/norman_effect_bundle.npz`, 3.3 MB).
- Figures are written to `figures/` (and `../analysis_cache/czi_fig/` for 07); result tables to `../results/` and `../analysis_cache/nb_out/`.

*This file is a reading-order overlay. For the per-notebook build-order description, see
[`README.md`](README.md).*
