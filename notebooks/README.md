# Notebooks

The list below is in **build order**. A first-time reader should instead follow the
**reading order** — payoff before plumbing. The full version with one-line "what you get"
summaries is in [`notebooks_README_reading_order.md`](notebooks_README_reading_order.md);
the short route:

| Reader question | Notebook(s) |
|---|---|
| **What can it tell me?** | `05_target_id_showcase` |
| **Is it real?** | `02_reachability_on_tier2` |
| **Can I trust it?** | `09_causal_validation_dossier` |
| **Does it generalize?** | `03_generalizability_and_impact` + `07_cross_celltype_transfer` |
| **How do I run one?** | `04_experimental_design_toolkit` + `bring_your_own_target` |
| Appendix / provenance | `01_exploratory_data_analysis`, `06_reinforcement_analyses`, `08_deg_weighted_evaluation` |

Physical file names are unchanged — this is a reading route, not a rename.

---

The method itself lives in `../reachability.py`; these
notebooks call it, and the batch drivers (`../scripts/run_atlas.py`, `../scripts/run_nulls.py`,
`../scripts/run_bootstrap.py`) reproduce the atlas and CIs headlessly. The per-notebook
build-order descriptions follow.

1. **`01_exploratory_data_analysis.ipynb`** — the two-tier data inventory, the Tier-1
   perturbation-level DE summary, a first look at the Tier-2 effect matrix `E`, the
   Th1/Th2 and aging target states, a reachability preview (why knockdown-only changes
   the geometry), and the autoimmune-enrichment disease linkage.
2. **`02_reachability_on_tier2.ipynb`** — the headline. Loads the Tier-2 matrix, builds
   the dictionary `E` and target `d`, runs the cone fit (the verdict), the honesty tests
   (held-out-gene validation + shuffled-target null), the reachability spectrum / minimal
   recipe, the GATA3↓/TBX21↓ positive control, the activation certificate, and the
   condition-resolved sensitivity sweep. Writes `results/table*.csv` and the `fig1–fig5`
   figures.
3. **`03_generalizability_and_impact.ipynb`** — transfers the identical code to a second
   dataset (Norman K562 CRISPRa): orientation gate, held-out CEBPA target, minimal recipe,
   spectrum, activation certificate, and the additivity test that only combinatorial data
   enables. Writes the `norman_table*.csv` tables and the `nb03_fig*` figures.
4. **`04_experimental_design_toolkit.ipynb`** — turns a reachability verdict into a runnable
   screen. Runs `design_experiment()` across all 12 atlas transitions: signed decomposition,
   held-out calibration (in-sample vs held-out reach), ranked knockdown/activation recipe,
   optimal-k library, and modality triage. Writes `nb04_fig*` and the design cards. See
   the experimental-design toolkit section of `../docs/RESULTS.md`.
5. **`05_target_id_showcase.ipynb`** — "From screen to shortlist," the pharma target-ID
   walkthrough: a disease-relevant transition → verdict → reach vs. must-activate split →
   ranked recipe → druggability + human-genetics triage → one-page target dossier.
6. **`06_reinforcement_analyses.ipynb`** — the reinforcement/limitations battery. Four analyses
   that harden the manuscript against its own stated limitations (**L4** constraint ablation,
   **L5** achievable ceiling, **L2** per-recipe additivity reliability, **L1** synthetic
   certificate test). Writes `../analysis_cache/nb_out/L*.{csv,json}` and `figR1_reinforcement_analyses.png`.
   Summary in `../analysis_cache/nb_out/REINFORCEMENT_RESULTS.md`.
7. **`07_cross_celltype_transfer.ipynb`** — cross-cell-type robustness. Uses the CZI
   Virtual Cell Models re-host of the Replogle 2022 essential-gene screens in **K562** and
   **RPE1** to test whether the cone geometry is a property of the biology or of one cell
   type's basis. Builds row/column-aligned `E_K562` / `E_RPE1` over 843 shared perturbations,
   then runs three tests against the unchanged `reachability.py`: single-perturbation effect
   correspondence (with shuffled-gene and mismatched-perturbation nulls), within- vs
   cross-cell-type reachability, and minimal-recipe overlap. Verdict: direction transfers,
   recipe does not. Writes `../analysis_cache/czi_fig/nb07_fig*` and `../analysis_cache/czi_data/per_perturbation_transfer.csv`.
   Summary in `../analysis_cache/czi_data/CROSS_CELLTYPE_TRANSFER.md`.
8. **`08_deg_weighted_evaluation.ipynb`** — DEG-weighted verdict evaluation: re-scores the
   reachability verdicts under a differentially-expressed-gene weighting of the readout to test
   whether the cone geometry survives emphasis on the most-perturbed genes. Writes the
   `deg_weighted_*` result tables.
9. **`09_causal_validation_dossier.ipynb`** — **"Can I trust it?"**, the reviewer's notebook and
   the causal-validation dossier. Organized by a **six-assumption identifying stack** (unbiased
   ATEs, additivity, SUTVA, homogeneous effect, exclusion restriction, transportability), it
   makes each assumption behind a reachability verdict explicit, tests the testable ones, and —
   where a test is impossible in silico — quantifies **how far the assumption must fail** to flip
   the verdict. Ports the IV/compliance layer and the Part-A analyses: **A1** verdict sensitivity
   radius (in measured-SE units, the headline), **A5** negative-control outcomes, **A6** signed
   construct validity, **A4** weak-instrument-robust recipe intervals, and the **A2** conditional-
   reachability scaffold. Reads the frozen `a{1,3,4,5,6}_*.csv` result files and
   `iv_compliance_verdicts.csv`; renders `fig_a1_sensitivity.png`. Companion prose:
   `../docs/CAUSAL.md` (the causal reframe, the A1–A6/B1–B4 agenda, and the validation ledger). This is the methods/validation
   dossier — distinct from the separately-planned certificate deep-dive (NB-C).

(Notebooks 06 and 07 are independent reinforcement passes; either runs standalone. Notebook 09 is
the causal-validation dossier and depends on the Part-A result files produced by its own batch
drivers `../scripts/run_a1_sensitivity.py` and `../scripts/run_iv_compliance.py`.)

Figures are written to `figures/` (and `../analysis_cache/czi_fig/` for nb07); result tables to `../results/`
and `../analysis_cache/nb_out/`. The heavy inputs (the h5ad matrix, cached `.npz` bundles) are gitignored.
