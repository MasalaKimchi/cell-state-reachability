# Generality on a second public dataset: the verdict transfers, the recipe does not

> **Provenance / scope note.** The Replogle-2022 K562/RPE1 cross-cell-type effect
> matrices and the headline transfer numbers below (verdict agreement 99.3% / 100%,
> recipe Jaccard median 0.11 vs null 0.053) were **first established in earlier work in
> this project** (notebook `07_cross_celltype_transfer`, cached in
> `czi_data/transfer_summary.json` and `per_perturbation_transfer.csv`). What this
> five-day-push step adds is **independent reproduction and consolidation**, not a new
> download: the full held-out reachability sweep was re-run from the raw effect vectors
> and reproduces the cached per-perturbation table to `max|Δ| = 0.0000` on every
> reach-cosine column, and the result is re-expressed here as an explicit
> verdict-transfers / recipes-don't generality statement with a publication figure and a
> tidy summary table. A genuinely new genome-wide Perturb-seq download was not feasible
> under this session's compute budget (no GPU, ~0.8 GB free RAM); the plan's feasibility
> note anticipated this fallback. The value here is the reproducibility check plus the
> packaged, figure-backed generality claim — treat it as consolidation of an existing
> in-repo result on a second, independent dataset, not as a first observation.


**Dataset.** Replogle et al. 2022 genome-wide Perturb-seq (CRISPRi), essential-gene
screens in two cell types — K562 (chronic myeloid leukemia) and RPE1 (retinal pigment
epithelium). We use the 843 perturbations measured in **both** cell types, reduced to
pseudobulk log-fold-change effect vectors over 2,832 shared genes (mean over perturbed
cells minus mean over control cells, on the source files' pre-log-normalized counts; no
re-normalization). Source: `cz-benchmarks-data` (Replogle 2022 K562/RPE1 essential).
This is a genuinely independent test: a different laboratory, different perturbation
technology mix, and two different cell types from the primary Tier-2 CD4 T-cell CRISPRi
atlas and the Norman 2019 K562 CRISPRa doubles used elsewhere in this work.

**Operator.** For every one of the 843 perturbations we hold it out as the target `d`
and fit the reachability cone over the remaining 842 perturbation vectors, both
**within** cell type (K562 target on K562 basis; RPE1 target on RPE1 basis) and
**cross-basis** (K562 target on the RPE1 cone, and vice versa). All fits are exact NNLS;
the reproduction from raw effect vectors reproduces the cached table to `max|Δ| = 0.0000`
on every reach-cosine column.

## Result 1 — the reachability verdict is portable across cell types

The per-perturbation reach cosine in K562 and in RPE1 are strongly rank-correlated
(**Spearman ρ = 0.57, p = 1.8 × 10⁻⁷³, n = 843**): perturbations whose target direction
is well inside the achievable cone in one cell type tend to be well inside it in the
other. Cross-basis reach — fitting a K562 target with only RPE1 perturbation vectors, and
vice versa — stays far above the shuffled-gene null (95th percentile **0.058**): median
cross-basis reach cosine **0.50** (K562 target / RPE1 basis) and **0.73** (RPE1 target /
K562 basis). Reduced to a binary reachable / not-reachable call, the verdict agrees across
bases for **99.3%** of K562 targets and **100%** of RPE1 targets. The *direction* of an
achievable state shift is a portable property of the biology.

## Result 2 — the minimal recipe is basis-specific

The same-target recipe overlap (Jaccard of the greedy minimal knockdown sets in K562 vs
RPE1) has median **0.11**, only marginally above the shuffled-recipe null (95th percentile
**0.053**); just 65% of perturbations clear that null at all, and cross-basis recipe
overlap sits **at** the null (0.053). So while the feasibility verdict transfers, the
*specific set of perturbations* that realizes it does not — it is a property of the
available basis in each cell type, not a portable prescription.

## Why this is the honest scope statement

"Reachability transfers, recipes don't" is the correct generalization claim for the
method. It says a practitioner can trust a **feasibility verdict** computed in one cellular
context as a guide to another, but must **re-derive the recipe** in the target context
against its own measured effects. This is exactly the division of labor the method is built
for: the verdict is the portable scientific claim; the recipe is the context-specific
engineering answer.

## Artifacts

- `fig_generality_second_dataset.png` — 3 panels: (A) reach cosine K562 vs RPE1 (verdict
  transfers, ρ = 0.57); (B) within- and cross-basis reach cosine vs the shuffled-gene null;
  (C) same-target recipe-overlap distribution vs the shuffled-recipe null (recipes don't
  transfer).
- `generality_second_dataset_summary.csv` — the headline metrics above.
- `generality_second_dataset_per_perturbation.csv` — all 843 perturbations, every within-
  and cross-basis reach cosine, verdict, and recipe-overlap value.
- `cross_celltype_effects.npz` — the effect matrices (K562 & RPE1, 843 × 2,832).
