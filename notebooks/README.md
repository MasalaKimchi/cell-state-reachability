# Notebooks

The nine notebooks are numbered in *build* order (how the project was assembled). A first-time
reader should meet the **payoff before the plumbing**, so start with the reading route below;
the per-notebook detail follows in file order. File names are unchanged — this is a route, not
a rename.

The method itself lives in [`../reachability.py`](../reachability.py). Every notebook calls into
it; none re-implements the geometry.

## Reading route

| Reader question | Notebook | What you get |
|---|---|---|
| **What can it tell me?** | [`05_target_id_showcase`](05_target_id_showcase.ipynb) | The teaser. A disease-relevant transition walked end to end: verdict → *reach vs. must-activate* split → ranked recipe → druggability + human-genetics triage → a one-page target dossier. Read this first to see the whole service in one pass. |
| **Is it real?** | [`02_reachability_on_tier2`](02_reachability_on_tier2.ipynb) | The headline result plus the honesty tests. Builds the effect dictionary `E` and target `d`, runs the cone fit (the verdict), then held-out-gene validation and a shuffled-target null so the reachable cosine is never reported uncalibrated. Includes the minimal recipe and the activation certificate. |
| **Can I trust it?** | [`09_causal_validation_dossier`](09_causal_validation_dossier.ipynb) | The reviewer's notebook. Organized by a six-assumption identifying stack (unbiased ATEs, additivity, SUTVA, homogeneous effect, exclusion restriction, transportability): each assumption made explicit, the testable ones tested, and — where a test is impossible in silico — a quantified *how far the assumption must fail to flip the verdict*. |
| **Does it generalize?** | [`03_generalizability_and_impact`](03_generalizability_and_impact.ipynb) + [`07_cross_celltype_transfer`](07_cross_celltype_transfer.ipynb) | Two transfer tests on unchanged code. **03**: the identical method on a second dataset (Norman K562 CRISPRa) — orientation gate, held-out target, recipe, spectrum, certificate, and the additivity test only combinatorial data enables. **07**: cross-cell-type (K562 vs RPE1) — *direction transfers, recipe does not*. |
| **How do I run one?** | [`04_experimental_design_toolkit`](04_experimental_design_toolkit.ipynb) + [`bring_your_own_target`](bring_your_own_target.ipynb) | **04**: turns a verdict into a runnable screen — `design_experiment()` across all 12 atlas transitions with signed decomposition, calibration, ranked recipe, optimal-k library, and modality triage. **`bring_your_own_target`**: the same as a one-call service on *your* target — paste a gene signature, get a live verdict, recipe, library, and certificate back. |
| Appendix / provenance | [`01_exploratory_data_analysis`](01_exploratory_data_analysis.ipynb), [`06_reinforcement_analyses`](06_reinforcement_analyses.ipynb), [`08_deg_weighted_evaluation`](08_deg_weighted_evaluation.ipynb) | The supporting record. **01**: the two-tier data inventory and a first look at `E`. **06**: the reinforcement/limitations battery (L1/L2/L4/L5). **08**: re-scores every verdict under a DEG-weighted readout. |

## In file order

1. **`01_exploratory_data_analysis.ipynb`** — the two-tier data inventory, the Tier-1
   perturbation-level DE summary, a first look at the Tier-2 effect matrix `E`, the Th1/Th2 and
   aging target states, a reachability preview (why knockdown-only changes the geometry), and the
   autoimmune-enrichment disease linkage.
2. **`02_reachability_on_tier2.ipynb`** — the headline. Loads the Tier-2 matrix, builds the
   dictionary `E` and target `d`, runs the cone fit (the verdict), the honesty tests (held-out-gene
   validation + shuffled-target null), the reachability spectrum / minimal recipe, the
   GATA3↓/TBX21↓ positive control, the activation certificate, and the condition-resolved
   sensitivity sweep. Writes `../results/table*.csv` and the `fig1–fig5` figures.
3. **`03_generalizability_and_impact.ipynb`** — transfers the identical code to a second dataset
   (Norman K562 CRISPRa): orientation gate, held-out CEBPA target, minimal recipe, spectrum,
   activation certificate, and the additivity test that only combinatorial data enables. Writes the
   `norman_table*.csv` tables and the `nb03_fig*` figures.
4. **`04_experimental_design_toolkit.ipynb`** — turns a reachability verdict into a runnable
   screen. Runs `design_experiment()` across all 12 atlas transitions: signed decomposition,
   held-out calibration (in-sample vs held-out reach), ranked knockdown/activation recipe,
   optimal-k library, and modality triage. Writes `nb04_fig*` and the design cards. See the
   experimental-design toolkit section of [`../docs/RESULTS.md`](../docs/RESULTS.md).
5. **`05_target_id_showcase.ipynb`** — "From screen to shortlist," the pharma target-ID
   walkthrough: a disease-relevant transition → verdict → reach vs. must-activate split → ranked
   recipe → druggability + human-genetics triage → one-page target dossier.
6. **`06_reinforcement_analyses.ipynb`** — the reinforcement/limitations battery. Four analyses
   that harden the manuscript against its own stated limitations (**L4** constraint ablation,
   **L5** achievable ceiling, **L2** per-recipe additivity reliability, **L1** synthetic certificate
   test). Writes `../analysis_cache/nb_out/L*.{csv,json}` and `figR1_reinforcement_analyses.png`.
   Summary in [`../docs/REINFORCEMENT_RESULTS.md`](../docs/REINFORCEMENT_RESULTS.md).
7. **`07_cross_celltype_transfer.ipynb`** — cross-cell-type robustness. Uses the CZI Virtual Cell
   Models re-host of the Replogle 2022 essential-gene screens in **K562** and **RPE1** to test
   whether the cone geometry is a property of the biology or of one cell type's basis. Builds
   row/column-aligned `E_K562` / `E_RPE1` over 843 shared perturbations, then runs three tests
   against the unchanged `reachability.py`: single-perturbation effect correspondence (with
   shuffled-gene and mismatched-perturbation nulls), within- vs cross-cell-type reachability, and
   minimal-recipe overlap. Verdict: direction transfers, recipe does not. Writes
   `../analysis_cache/czi_fig/nb07_fig*` and `../analysis_cache/czi_data/per_perturbation_transfer.csv`.
   Summary in [`../docs/CROSS_CELLTYPE_TRANSFER.md`](../docs/CROSS_CELLTYPE_TRANSFER.md).
8. **`08_deg_weighted_evaluation.ipynb`** — DEG-weighted verdict evaluation: re-scores the
   reachability verdicts under a differentially-expressed-gene weighting of the readout to test
   whether the cone geometry survives emphasis on the most-perturbed genes. Writes the
   `deg_weighted_*` result tables.
9. **`09_causal_validation_dossier.ipynb`** — the causal-validation dossier. Organized by a
   **six-assumption identifying stack**, it makes each assumption behind a reachability verdict
   explicit, tests the testable ones, and — where a test is impossible in silico — quantifies **how
   far the assumption must fail** to flip the verdict. Ports the IV/compliance layer and the Part-A
   analyses: **A1** verdict sensitivity radius (the headline), **A5** negative-control outcomes,
   **A6** signed construct validity, **A4** weak-instrument-robust recipe intervals, and the **A2**
   conditional-reachability scaffold. Reads the frozen `a{1,3,4,5,6}_*.csv` result files and
   `iv_compliance_verdicts.csv`; renders `fig_a1_sensitivity.png`. Companion prose:
   [`../docs/CAUSAL.md`](../docs/CAUSAL.md).

Plus **`bring_your_own_target.ipynb`** — the one-call service on an arbitrary target signature.

## Notes

- **Reading order ≠ dependency order.** Notebook 09 reads frozen result files written by its batch
  drivers (`../scripts/run_a1_sensitivity.py`, `../scripts/run_iv_compliance.py`). Notebooks 06 and
  07 are independent reinforcement passes; each runs standalone.
- **Batch drivers.** `../scripts/run_atlas.py`, `../scripts/run_nulls.py`,
  `../scripts/run_bootstrap.py`, `../scripts/run_a1_sensitivity.py` and
  `../scripts/run_iv_compliance.py` reproduce the atlas, nulls, CIs and causal-layer result files
  headlessly.
- **Heavy inputs are gitignored** — the Tier-2 h5ad (16.8 GB), `analysis_cache/**/inputs.npz`
  (772 MB), and the three per-condition `cache/E_*.npz` matrices (150–160 MB each). Notebooks work
  from small cached tables and bundles (e.g. `cache/norman_effect_bundle.npz`, 3.2 MB).
- **Figures** are written to `figures/` (and `../analysis_cache/czi_fig/` for 07); result tables to
  `../results/` and `../analysis_cache/nb_out/`. `notebooks/figures/` is gitignored; the figures the
  narrative docs embed are copied into `../docs/figures/`, which is tracked.
