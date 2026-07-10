# Cross-cell-type reachability transfer: what generalizes, what doesn't

*Extension analysis for the cell-state-reachability project. Companion to notebook
`07_cross_celltype_transfer.ipynb`.*

## The question

The project's headline result — that a Th2→Th1 polarization target state is *reachable* by a
short recipe of knockdowns — is a statement about the geometry of a **convex cone** spanned by
perturbation-effect vectors. That cone was measured in CD4⁺ T cells. A natural and load-bearing
question for any reviewer is:

> Is reachability a property of the **biology** of the target state, or of the **specific cell
> type's measured perturbation basis**?

If the former, the cone geometry should transfer to other cell types. If the latter, it is a
per-dataset curiosity. We now have the data to test this directly.

## The data (and why it's newly available)

The CZI Virtual Cell Models platform — the same platform hosting this project's CD4⁺ T anchor —
publicly re-hosts the **Replogle et al. 2022** genome-scale essential-gene CRISPRi Perturb-seq
screens in two human cell types: **K562** (chronic myelogenous leukemia) and **RPE1** (retinal
pigment epithelium). These were previously unreachable in this environment (figshare-blocked);
CZI's `cz-benchmarks-data` bucket resolves that.

We restricted to the **843 single-gene perturbations measured in both cell types** and the
**2,832 shared readout genes**, and built two effect matrices, `E_K562` and `E_RPE1`, that are
row- and column-aligned so every comparison is like-for-like. Each perturbation's effect is a
pseudobulk log fold-change versus the non-targeting control pool, computed identically in both
cell types. On-target self-effects are negative for **100 %** of perturbations in both cell
types — CRISPRi knocks down its own target, confirming the effect vectors are correctly built
and oriented.

We then reused the project's `reachability.py` **unchanged**.

## Three tests, three answers

### 1. Single-perturbation effect direction — *transfers (moderately)*
Matched knockdowns produce correlated effect vectors across cell types (median cosine
**+0.35**), far above a shuffled-gene null (≈0). There is a real shared component — a common
essential-stress / growth-arrest direction that any essential-gene knockdown induces (mean-effect
cross-type cosine **+0.64**) — but gene-specific structure survives removing it: matched cosine
stays positive for **96 %** of perturbations after deflation, and **69 %** are individually
gene-specific at p<0.05. *The biology of a single knockdown partly transfers.*

### 2. Reachability verdict — *transfers (nearly completely, but for a subtle reason)*
Held-out targets are more reachable **within** their own cell type than **across**:
cone residuals are K562 0.61 (within) vs 0.87 (cross), RPE1 0.37 vs 0.68. So the cone geometry
is measurably cell-type-specific. Yet cross-type reachable cosines (0.50–0.73) stay well above
the shuffled-target null, so the binary reachable/outside **verdict agrees ~100 %** within vs
cross. This near-perfect agreement is partly a consequence of an expressive over-complete cone
(843 non-negative generators over 2,832 genes) — the honest, discriminating signal is the
**graded residual**, which does show a cell-type penalty.

### 3. Minimal recipe — *does NOT transfer*
The sharpest test. Even for targets reachable in both cell types, the NNLS minimal recipe
recruits **substantially different generators**. Same-target recipes overlap at median Jaccard
**0.11** — ≈20× above a random-subset null (0.006), so *not* random, but far from identity.
And when the *same* target is reached from the *other* cell type's basis, recipe overlap
collapses to the null (**0.05**). *Knowing which knockdowns reach a state requires that cell
type's own perturbation responses.*

## What it means for the Th2→Th1 headline

This is a **robustness result with a sharp boundary condition**:

- **It defends the method.** The convex-cone framework produces coherent, above-null, partly
  cell-type-invariant structure in three independent human cell types (CD4⁺ T, K562, RPE1). The
  cone is not an artifact of one dataset. A reviewer asking "does this only work in T cells?" is
  answered.
- **It bounds the prescription.** The *direction* toward a target state transfers; the *specific
  minimal recipe* does not. The GATA3↓/TBX21↓-style recipe the paper derives for Th2→Th1 is
  validated **for CD4⁺ T cells** and should not be assumed to be the recipe in another cell type
  without re-fitting on that cell type's basis. This is exactly what the project's philosophy
  (trust null-calibrated, held-out metrics; never the raw in-sample cosine) would predict.

## Caveats

1. K562/RPE1 are essential-gene loss-of-function screens; the shared basis is enriched for core
   cellular machinery, which raises the common-stress baseline and inflates raw cosine — hence
   the gene-specificity controls throughout.
2. Near-100 % binary verdict agreement reflects an expressive cone, not perfect biological
   transfer; the graded residual and recipe overlap are the honest discriminators.
3. Effect magnitudes differ ~2× between cell types (RPE1 larger); all cross-type comparisons are
   cosine/rank-based and therefore scale-invariant.
4. This tests transfer *between two non-T-cell lines*; it is evidence about the method's
   cell-type-invariance in general, not a direct re-measurement of the Th2→Th1 recipe in a
   second T-cell system (no second genome-scale CD4⁺ T screen is available).

## Reproducibility

- `build_effect_matrices.py` — builds the aligned effect matrices from the two CZI h5ad files.
- `czi_data/cross_celltype_effects.npz` — the checkpoint (E_K562, E_RPE1, shared perturbations
  and genes, per-perturbation cell counts).
- `notebooks/07_cross_celltype_transfer.ipynb` — runs all three tests against the unchanged
  `reachability.py`. Its code cells were re-run as a standalone script (`_nb07_verify.py`,
  extracted verbatim from the notebook) in the `cellreach` environment to confirm the numbers
  in this writeup regenerate; `jupyter nbconvert --execute` itself is unavailable here because
  the sandbox blocks the kernel's TCP-socket bind, so the committed `.ipynb` carries the code
  and markdown but not embedded cell outputs.
- `czi_data/per_perturbation_transfer.csv` — every metric, per perturbation (843 rows).
