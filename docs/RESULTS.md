# RESULTS — the reachability method on real Tier-2 data

*This is the consolidated results record for the cell-state reachability project. It absorbs
five former root-level files, folded into the sections noted: the original **RESULTS.md**
(spine — what was run and what it found, §§1–4, §8.2, §9); **MODALITY.md** (signed
loss-/gain-of-function decomposition, the 12-cell atlas, and the druggability×genetics triage,
§5); **GENERALIZABILITY.md** (cross-dataset and cross-cell-type transfer, the input-contract
framing, and the application map, §6); **DESIGN.md** (the `design_experiment()` toolkit and its
calibration caveat, §7); and **METHOD_IMPROVEMENTS.md** (the "novel algorithm vs. novel
formulation" answer, the nine post-hackathon in-silico results, and the methodological +
experimental agendas, §8.1 and §8.3–§8.7). Cross-references have been repointed to the merged
section numbers; every number, table pointer, and figure from the five sources is preserved.*

*What was run, what it found, what succeeded, what failed, why, and what to improve.
Every number here is computed from the real `GWCD4i.DE_stats.h5ad` effect matrix
(33,983 CRISPRi knockdowns × 10,282 genes, 16.8 GB) this session — reproduced by
`notebooks/02_reachability_on_tier2.ipynb`, which runs end-to-end and was validated
cell-by-cell. The method module is `reachability.py`.*

---

## TL;DR

The convex-cone reachability method **works on real data and returns a scientifically
meaningful, honest verdict.** The Th2→Th1 polarization shift is **partially reachable** by
CRISPRi knockdown, and the method's distinctive contribution — the part no differential-
expression ranking produces — is that it says *which part* and *why*:

> **Knockdown can remove the Th2 program but cannot install the Th1 program.**
> GATA3 (the Th2 master driver) is reached; TBX21/IFN-γ (the Th1 drivers) land in the
> activation certificate because you cannot raise a gene by knocking things down.

The headline number — a held-out reachable cosine of **0.45**, which clears its shuffled-target
null at **z ≈ 24** (60 shuffles) — is the honest one: it survives held-out-gene validation, so it
is not an artifact of fitting 700+ non-negative weights. (The same cell's held-out-gene z is
reported as 45 in the atlas below; that run used a tighter 8-shuffle null band, so z ≈ 24 is the
conservative figure and the two are the same 0.45 cosine — see §5.7.)

**Post-hackathon (§8):** three methodological deepenings were added and validated against data —
(1) the positive control broadened from 2 genes to a regulator panel (master-TF drivers give
AUROC = 1.00, p = 0.014; induced receptors correctly do *not*); (2) an epistasis penalty
calibrated on measured Norman doubles — collinearity is refuted, **magnitude saturation** is the
real mechanism (scale-free ceiling M\* = 13.9), and the atlas recipes are shown to sit safely in
the additive regime (risk 0.04–0.08 ≪ 0.5); (3) a **closed-form anisotropy-corrected null**,
`E[null cos] ≈ √((a·ρ₁)² + (1−a²)·κ)`, validated against the empirical shuffle null (Pearson
0.995) and replacing ~1000 refits with ~10–20. All three are wired into `reachability.py` behind
its self-test.

---

## 1. What was run

| Stage | What | Status |
|---|---|---|
| Data integrity | Structure/dtype/NaN audit + self-knockdown orientation gate | ✅ passed |
| QC reproduction | Cross-source concordance, per-condition significant counts | ✅ reproduced exactly |
| Target construction | `toward_Th1` and `aging` signatures mapped into E's gene space | ✅ validated by markers |
| Cone fit | NNLS projection onto the knockdown cone + KKT/Farkas certificate | ✅ certified optimal |
| Honesty | Held-out-gene validation + shuffled-target null (60 shuffles) | ✅ signal real (z≈24) |
| Spectrum | Greedy minimal recipe + null band | ✅ above p99 at every k |
| Positive control | GATA3↓ / TBX21↓ placement | ✅ both correct |
| Activation certificate | Ranked unreachable-upward genes (CRISPRa hypotheses) | ✅ immunologically credible |
| Conditions + sensitivity | Rest/Stim8hr/Stim48hr × {z-score, log-fc} × {sig, all} | ✅ verdict stable |
| Second axis | CD4 aging signature | ⚠️ ran; null correctly discounts it |
| IV / compliance | Guide non-compliance as instrumental variables (ITT vs LATE; valid-instrument subset) | ✅ verdict invariant (max \|Δcos\|=2.2e-16 rescale, 4e-4 drop) |
| A1 · sensitivity radius | Verdict robustness in measured-SE units (measurement-error bootstrap + coordinated-bias radius) | ✅ Th1/Th2 robust to noise, flip at ≈0.03 SE-units; aging axes fragile |
| A3 · k-way additivity | Directional additivity retention vs recipe size (Norman-calibrated) | ✅ 0.71 (k=2) → 0.29 (k=12) |
| A4 · weak instruments | Anderson–Rubin-style 1/π recipe-weight intervals | ✅ headline recipe clean; catches SNX5 (π≈0.07) |
| A5 · negative controls | Housekeeping/proteasome negative-control-outcome projection | ✅ 4.1–5.5× positive-over-negative enrichment |
| A6 · construct validity | Signed cosine of master-regulator KD vectors vs Th1 axis | ✅ TBX21↓/GATA3↓ correct-signed |
| A2 · conditional reachability | Subtype-stratified re-solve | ⏳ scaffold — needs raw single-cell counts |

---

## 2. What succeeded, with the numbers

### 2.1 The data is real and correctly oriented
The decisive gate: in 300 sampled significant perturbations, the knocked-down gene's own
z-score in its own row has **median −9.93, is negative in 100% of rows, and below z = −2 in
99.3%.** A CRISPRi knockdown lowers its own target — the matrix means what it claims. All six
layers (`log_fc`, `zscore`, `p_value`, `adj_p_value`, `baseMean`, `lfcSE`) are dense and
sane. Documented cross-source Th2/Th1 concordance reproduced to the digit: **11,616 shared
genes, 68.5% sign-concordant, Spearman ρ = 0.562 (z-score) / 0.533 (log-fc).**

### 2.2 The verdict is meaningful and certified
For Rest / `toward_Th1` over 6,188 signature genes with 6,871 significant knockdowns:

| Quantity | Value | Meaning |
|---|---|---|
| reachable cosine (in-sample) | 0.627 | best cone alignment |
| **held-out cosine** | **0.448** | **honest — generalizes to unseen genes** |
| residual norm (relative) | 0.779 | most of the shift is *not* reached |
| reachable fraction | 0.393 | ~39% of target norm is knockdown-reachable |
| activation-required fraction | 0.309 | ~31% needs genes to go *up* |
| KKT/Farkas violation | 1.1 × 10⁻¹¹ | the "outside the cone" claim is certified |

The KKT/Farkas violation near machine-zero matters: it means "partially outside the cone" is
**proved** (the solver found the true cone projection), not asserted. The proof is inherited
from convex-optimization theory — we verify the certificate numerically rather than proving
theorems.

### 2.3 The signal is real — it survives the honesty tests
This is the result that separates a method from a coincidence. A 6,871-generator non-negative
fit in 6,188-dimensional space *could* reach a lot by chance. It doesn't:

- **Held-out-gene validation:** fit weights on half the genes, score on the other half →
  cosine **0.448**. In-sample was 0.686, so there is overfitting, but the held-out value is
  where the real structure lives.
- **Shuffled-target null (60 permutations):** null mean **−0.005**, SD 0.019, **max 0.029**.
  The observed 0.448 is **z ≈ 24** and ~15× the largest random draw.
- **Greedy spectrum vs null band:** the observed reachable cosine exceeds the 99th-percentile
  null at **every** sparsity k (k=1: 0.226 vs 0.059; k=15: 0.457 vs 0.161).

### 2.4 The geometry respects known biology (positive control)
Both ends of the Th1/Th2 axis land where immunology predicts:

- **GATA3↓** (remove the Th2 master TF): cosine **+0.101, rank 155 / 6,871 (97.7th pctile)** —
  top ~2%, pointing toward Th1. ✅
- **TBX21↓** (remove the Th1 master TF, T-bet): cosine **−0.099, rank 6,775 / 6,871
  (1.4th pctile)** — bottom ~1%, correctly pointing *away* from Th1. ✅

GATA3 was **not** hand-picked into this result; it emerges from the geometry. Note it is a
positive *control*, not a discovery — the source screen already reports it.

### 2.5 The activation certificate is the constructive payoff
The residual names **2,481 genes with positive unmet upward demand** — genes `toward_Th1`
wants raised that no non-negative knockdown mix delivers. The top of that list is
**LYAR, IKZF3 (Aiolos), CRTAM, CBLB, GBP5, LAG3, IRF8** — and it must be read carefully.
A Jul-10 literature/genetics cross-reference (20 genes; 4 SUPPORTS / 7 PLAUSIBLE / 9 UNKNOWN,
with 8/20 carrying a strong ≥0.5 Open Targets link to a Th1-driven autoimmune disease) shows
that the highest-demand genes are **state-markers** (LYAR, CRTAM, GBP5) and **negative
regulators** (IKZF3/Aiolos, CBLB, LAG3) — genes whose *deletion*, not activation, is the known
way to boost T-cell function. So for those entries an "activate this gene" reading may be
**directionally backwards**. The canonical Th1 master TFs are present but rank deep
(**IFN-γ #35, TBX21 #144, STAT4 #364, STAT1 #473**), and GATA3, which must go *down*, is
correctly **absent** (it is reached). The certificate is therefore best read as a **reproducible
ranking of where activation is required** (stable under a gene-axis split: cross-half Spearman
ρ ≈ 0.65, half-vs-full ≈ 0.84, z ≈ 35), *not* as a validated list of activation targets. It is
the "what to look at instead" that a bare infeasibility verdict lacks — each entry a falsifiable
CRISPRa hypothesis whose sign must still be checked at the bench. See manuscript Fig. S9 and
§8.3.6.

### 2.6 The verdict is robust
| Variant | held-out cosine |
|---|---|
| z-score + significant (headline) | 0.448 |
| log-fc + significant | 0.334 |
| z-score + all perturbations | 0.424 |
| Rest / Stim8hr / Stim48hr | 0.448 / 0.283 / 0.268 |

The verdict does not hinge on a single knob. z-score is the strongest metric (as the source
documentation states); resting cells are the most reachable, and stimulation raises the
activation-required fraction (0.31 → 0.38) — biologically sensible, since stimulation adds an
activation-driven component that loss-of-function cannot supply.

It also survives a harder stress test. Under **DEG-magnitude weighting** (w = |d|, the
signal-dilution correction of Mejia et al. 2026), the Rest verdict *strengthens* rather than
collapses (reachable cosine **0.627 → 0.803**; held-out z stays far above the z=3 floor,
28.3 → 14.1), and the calibrated dynamic-range placement is stable because the shuffled floor
rises in step (0.346 → 0.619). The unweighted metric remains the default and reproduces these
published numbers bit-for-bit. See manuscript Figs. S5–S6.

---

## 3. What failed, or fell short — and why

### 3.1 The dense fit is not the recipe — it spreads weight across hundreds of generators
The full NNLS puts weight on **715 knockdowns**, and GATA3 lands at **weight-rank 116** in that
dense solution — not because GATA3 is unimportant, but because 6,871 effect vectors are highly
correlated, so the solver distributes weight across many near-equivalent directions.
**Why it matters:** the dense weight vector is *not* a usable experimental recipe, and reading
biology off it is misleading. **The fix is already in place:** the *sparse greedy spectrum* is
the object to read for the minimal set and the positive control — there, the recipe is
interpretable (LAT2, ICOS, RARA … at the top) and the knee is at k≈7. Report the spectrum, not
the dense support, as "the recipe."

### 3.2 The greedy spectrum's top genes are not the canonical regulators
The k=1–7 greedy picks are LAT2, APPBP2, RARA, ICOS, SNAP23, SNX4, VTI1A — **not** GATA3.
**Why:** reachability of the *whole multi-gene state shift* is a different question from "which
single knockdown is the best-known regulator." The greedy set is chosen to cover the largest
share of the target vector, which rewards knockdowns with broad, well-aligned effects over the
textbook master regulator. This is a real limitation for interpretability, and an honest one to
state: **the minimal reachable set is an engineering answer, not a named-regulator answer.**
GATA3's role is captured by the positive-control ranking (§2.4), not the greedy recipe.

### 3.3 The absolute cosine is modest (0.45), not near 1
Th2→Th1 is **not** cleanly inside the cone. Only ~39% of the target norm is reachable.
**Why:** by construction, ~55% of the target's norm lives in genes that must go *up*, and
CRISPRi is loss-of-function. The method is not failing here — it is correctly reporting that a
knockdown-only screen structurally cannot achieve full polarization. That *is* the headline
finding, but it means the project should be framed as "reachable-fraction + certificate," never
as "we can drive Th2→Th1 by knockdown."

### 3.4 The aging axis exposed a real statistical subtlety
Run on the CD4 aging signature (86% of its norm in up-genes), the held-out cosine is 0.41 with
z ≈ 11.5 — but the **shuffled-target null mean is 0.256, not ≈ 0.** 
**Why:** when a target is dominated by a single shared direction (here, "most genes up"), a
high-dimensional non-negative cone can reach *any* such biased target partway by chance, so the
chance ceiling rises. **This is the method working, not breaking:** the null band automatically
discounts the aging result to its true above-chance margin. The lesson — worth stating in any
write-up — is that **the raw cosine is meaningless without the target-specific null**, and a
mostly-one-direction target needs the null even more than a balanced one.

### 3.5 Documentation bugs found in the data
- The h5ad is **16.8 GB**, not the documented "~1.4 GB" — a 12× underestimate that forces
  streaming/selective reads on an 18 GB-RAM laptop (never a full load).
- The sgRNA library table has **26,504 rows**, not the documented 31,110 (confirms a prior
  Day-0 finding; the schema gate should be corrected to 26,504).
- The aging signature keys genes by symbol in `gene_name` (the `variable` column is Ensembl
  IDs) — the naïve join returns 0 genes; use `gene_name`.

---

## 4. What to improve (ranked)

1. ~~**Fold the fast greedy into `reachability.py`.**~~ **✅ DONE.** `reachability_spectrum()`
   now defaults to the OMP-style fast path (`refit_full=False`): score candidates by residual
   correlation, refit NNLS on the small active set. `refit_full=True` keeps the exact reference.
   Verified same selection order, ~100× faster.
2. **Make the target-specific null a required output, not an option.** §3.4 shows the raw cosine
   can mislead. `reachability()` should refuse to print a verdict without the paired shuffled-
   target null for *that* target, and report the z-score as the primary statistic. *(Still open;
   the atlas runner enforces this at the pipeline level — every cell carries a held-out z.)*
3. **Report the greedy spectrum as "the recipe," retire the dense support for interpretation.**
   Keep the dense fit for the reachable-fraction; use the sparse path for anything a biologist
   reads. *(Adopted in the atlas: nominations come from the greedy path only.)*
4. ~~**Add a signed-reachability mode.**~~ **✅ DONE.** `signed_reachability()` fits knockdown
   (LOF) and sign-flipped activation (GOF) generators in two stages, reports the LOF/GOF/neither
   split, and certifies each stage. See §5. This became the centerpiece of the
   expansion.
5. ~~**Cache the reduced matrices as a shipped artifact.**~~ **✅ DONE.** `analysis_cache/atlas_work/inputs.npz`
   bundles all three per-condition effect matrices + the four transition targets; the atlas
   reruns from it in minutes without touching the 16.8 GB file.
6. **Broaden the positive control.** GATA3/TBX21 is a two-point check; a fuller validation would
   correlate single-knockdown alignment against an external Th1/Th2-regulator gene set and
   report an enrichment statistic. *(Still open.)*

---

## 5. The expansion — signed reachability, the atlas, and modality triage

The single Th2→Th1 verdict above is one cell of a much larger, cross-verified
result. This round added the signed cone, ran a 12-cell atlas, and crossed every
nomination against druggability and human genetics. The modality framing and full triage are detailed in the subsections below
(§5.1–§5.8, absorbing the former `MODALITY.md`).

### 5.1 The modality question — can the state change be produced at all, and by adding or removing function?

Target discovery for immune and inflammatory disease is dominated by two
questions asked in isolation: **"is this gene associated with disease?"**
(genetics) and **"can this gene be drugged?"** (tractability). Neither asks the
prior question that a cell-engineering program actually faces: **can the desired
cell-state change be produced at all, and if so, by adding function or by
removing it?** A knockdown-reachable state is a small-molecule / degrader /
siRNA program; a state that requires *gaining* function is an agonist,
cytokine, or cell-engineering program; a state reachable by *neither* is not a
drug target at any modality.

This project adds that missing axis. The signed reachability cone decides, for a
measured cell-state transition, how much of the required change is:

- **LOF-reachable** — producible by knocking genes *down* (the perturbation
  library is a knockdown screen, so this is the directly-measured cone);
- **GOF-required** — only producible by *activating* genes (reachable after
  allowing sign-flipped generators);
- **neither** — orthogonal to everything the library can do in either direction
  (irreducible; not a drug target at any modality).

Crossing that per-node modality requirement with Open Targets tractability and
human genetics turns a ranked gene list into a **triage**: green-light
(reachable + genetic + druggable), and — the finding that motivates the whole
contribution — **required-but-undruggable** (reachable + genetic, but no handle
at any conventional modality).

### 5.2 Signed reachability decomposes every verdict into LOF / GOF / neither

`signed_reachability()` in `reachability.py` performs a staged non-negative fit:

1. **LOF stage** — `w* = argmin_{w≥0} ‖Eᵀw − d‖`, the knockdown cone (identical
   to the original one-sided reachability; the LOF cosine reproduces the
   published one-sided reachable cosine exactly).
2. **GOF stage** — refit the LOF residual `r = d − Eᵀw*` using sign-flipped
   generators `−E`, i.e. `u* = argmin_{u≥0} ‖(−E)ᵀu − r‖`. Only perturbations in
   an optional `gof_mask` may act as activation generators.
3. The three squared-norm shares (LOF / GOF / neither) are orthogonal and sum
   to 1.

Optimality is certified **per stage** (not by a joint KKT check on the final
residual, which is not the right condition): the LOF fit satisfies the Farkas
certificate `max(E @ r) ≤ tol`, and the GOF fit satisfies its own NNLS gradient
condition. Both hold to machine precision (‖violation‖ ~1e-11) on every atlas
cell.

**Validation that the sign is real, not an artifact.** On Th2→Th1 (Rest), the
master regulators land where biology demands: GATA3 (Th2 driver, target
component d = −15.2) is assigned to **LOF/remove** (w_LOF = 0.061, w_GOF = 0);
STAT4 (Th1 driver, d = +11.6) is assigned to **GOF/activate** (w_GOF = 0.062,
w_LOF = 0). No gene is mis-signed: GATA3 never receives activation weight, STAT4
never receives knockdown weight.

On Th2→Th1 Rest the signed split is **LOF 0.39 / GOF
0.25 / neither 0.35**, and the LOF cosine (0.627) reproduces the published
one-sided reachable cosine *exactly* — the signed fit is a strict superset of the
original method.

*Reconciliation with the one-sided numbers in §2.2/§3.3.* The one-sided report
splits the target into reachable **0.39** / activation-required **0.31** /
orthogonal **0.30**; the signed report splits it into LOF **0.39** / GOF **0.25**
/ neither **0.35**. The reachable share is identical (0.393 either way). The
difference is only in the *not-reached* remainder (0.61 in both): the one-sided
"activation-required 0.31" is a heuristic upper bound (the LOF residual projected
onto genes that must go up), whereas the signed GOF **0.25** is what a *non-negative
activation cone actually reaches*, leaving the true irreducible remainder as
neither **0.35**. The signed split is the more honest of the two — it does not
credit activation with more than a real activation basis could supply.

### 5.3 A 12-cell reachability atlas, every cell validated

Four measured cell-state transitions — `toward_Th1`, `toward_Th2`,
`toward_younger`, `toward_older` — across three culture conditions (Rest,
Stim8hr, Stim48hr). Each cell is a full-readout signed decomposition with a
held-out-gene validation (n = 8 shuffles). Full table:
`results/atlas_reachability.csv`.

| Finding | Value |
|---|---|
| Cells statistically significant (held-out z ≥ 3) | **12 / 12** (z range 4.7–45.0) |
| Verdict = PARTIALLY REACHABLE | 3 (toward_Th1, toward_Th2, toward_younger — all **Rest**) |
| Verdict = WEAKLY REACHABLE | 9 |
| Mean LOF fraction (knockdown-reachable) | **0.34** |
| Mean GOF fraction (activation-required) | 0.24 |
| Mean *neither* fraction (irreducible) | **0.42** |

**Two robust patterns.**

1. **Resting cells are the most reachable for every axis.** Rest gives the
   highest reachable cosine in all four transitions; the Th axes fall most
   steeply under stimulation (~0.63 → ~0.53), the aging axes fall less
   (~0.61 → ~0.58). An activated T cell's transcriptional program is harder to
   redirect than a resting one — for polarization more than for aging.
2. **Knockdown is never the majority modality.** Across all 12 cells, LOF
   explains ~34% of the target demand on average; activation (~24%) plus the
   irreducible component (~42%) always account for the rest. Even the best
   verdict (Th2 Rest, reachable cosine 0.64) is only 41% LOF. **A pure
   knockdown/degrader program cannot deliver any of these transitions on its
   own** — this is the quantitative core of the modality argument.

Table: `results/atlas_reachability.csv`. Figure:
`notebooks/figures/fig_atlas_decomposition.png`.

### 5.4 Modality triage — the required knockdown nodes are often undruggable

The greedy knockdown recipe for every atlas cell yields a union of **102 unique
LOF-node nominations**. Each was profiled through Open Targets (tractability +
human genetics). Full table: `results/modality_intervention_map.csv`.

**Druggability of the required knockdown nodes:**

| Modality tier | Count |
|---|---|
| Clinical-grade drug/candidate | 10 |
| Small-molecule tractable | 30 |
| Antibody-tractable (surface) | 17 |
| Degrader-only (*predicted* ubiquitination, no real handle) | 33 |
| Conventionally undruggable | 12 |

**45 of 102 (44%) of the required knockdown nodes are hard-to-drug**, and only
10 have a clinical-grade drug. The reachability cone repeatedly nominates genes
the druggable genome cannot yet reach.

**The collision — genetically-supported but undruggable.** Crossing reachability
priority with immune-disease human genetics exposes nominations that are
well-supported yet have no conventional handle:

| Gene | Immune-genetic assoc. | Top disease | Modality tier | Atlas axis |
|---|---|---|---|---|
| **IRF1** | 17 | asthma (0.83) | degrader-only (predicted) | toward_Th2 (rank 1) |
| ATXN2 | 7 | Hashimoto thyroiditis | degrader-only (predicted) | toward_younger (rank 2) |
| FAM98A | 5 | asthma | degrader-only (predicted) | toward_younger (rank 1) |
| C1D | 4 | seasonal allergic rhinitis | degrader-only (predicted) | toward_older (rank 1) |

IRF1 is the **strongest immune-genetics nomination among the hard-to-drug
subset** (17 associations) and a top-ranked knockdown node for the Th2 axis, yet
its only Open Targets tractability is a *predicted* ubiquitination annotation —
no ligand, no antibody, no drug. It is exactly the kind of target this triage is
built to flag: do not start a small-molecule campaign; it needs a degrader
discovery effort or a genetic/cell-engineering approach.

The genetically-supported nominations split cleanly into three actionable
buckets by how far along the modality is:

**Green-light — reachable + genetic + a clinical-grade drug already exists.**
Nine nominations have both immune-disease genetics and ≥1 approved drug or
clinical candidate: **JAK2** (14 immune-genetic assoc., ulcerative colitis; 31
drugs/candidates), **ICOS** (13; common variable immunodeficiency; 4),
**MAPK14** (8; ulcerative colitis; 28), **CD3D** (3; immunodeficiency 19; 16),
plus CD5, IFNGR1, CASP7, KEAP1, VKORC1. These are where a knockdown /
small-molecule program is both justified and feasible today.

**Tractable handle but not yet drugged — reachable + genetic + a plausible
modality, zero drugs so far.** These are the highest-value *new* nominations,
because the genetics is strong and a modality exists but no one has drugged them:
**IL7R** (the single highest immune-genetic count at **22** associations;
antibody-tractable surface receptor; 0 drugs), **ZAP70** (12; SM-tractable; top
knockdown node for toward_older; 0 drugs), and **TET2** (12; SM-tractable; 0
drugs). IL7R is therefore *not* a collision (it has an antibody handle) and *not*
green-light (nothing is drugged yet) — it is precisely the reachable,
genetically-anchored, druggable-in-principle-but-untried target a discovery team
most wants surfaced.

Tables: `results/modality_intervention_map.csv`,
`results/genetics_crossverification.csv`. Figure:
`notebooks/figures/fig_modality_triage.png`.

### 5.5 Negative result — a disease GWAS set is not a reachability target

Before pivoting to transition axes, we tested whether a disease's own GWAS gene
set defines a reachable target. **Four constructions, all failed:** unit-
suppression on the module subspace is *degenerate* (reach cosine 1.000 trivially,
because a small subspace is always fittable when P ≫ subspace); the mean of the
disease genes' measured knockdown vectors is *circular* (LOF 1.000 by
construction — the target is a non-negative mix of E's own rows); a flat
suppression over a variable readout with disease perturbations held out is *worse
than random* (IBD z = −2.6, asthma z = −4.4 vs. a matched random gene set); and
polarization values restricted to disease genes are *modest and inconsistent*
(IBD z = 1.3, SLE z = 1.8, asthma z = 3.0). **Root cause:** reachability needs a
*coherent cell-state-transition signature* (Th2→Th1 gets z ≈ 45 because it *is*
one); a GWAS list is a bag of heterogeneous risk loci with no coordinated
transcriptional direction. This fixed a validity rule now enforced throughout:
**the target `d` must be specified independently of the effect matrix `E`.**
Disease relevance therefore enters the atlas *per nomination*, via genetics —
never by building a disease-shaped target. (The four-construction test was a
one-off diagnostic run this round; its per-construction z-scores are quoted
above.)

### 5.6 What the modality triage does *not* claim

- The transition signatures (Th1/Th2 regulator identities, the aging axis) are
  taken from the source papers; the reachability verdict is the contribution,
  not the biology of the signatures.
- Tractability tiers and drug counts are Open Targets facts. Transcription-factor
  identity (GATA3/TBX21/STAT4/IRF1 as TFs) is domain knowledge — Open Targets
  returned an empty target class for these genes this round, so the "TFs are not
  a small-molecule modality" argument rests on their biology, not on an OT
  classification.
- **A disease's GWAS gene set is *not* a validly-reachable target.** We tested
  four ways of turning a disease's risk genes into a reachability target; all
  failed (degenerate, circular, or non-significant vs. a matched random gene
  set — see §5.5). Reachability requires
  a *coherent transition signature*; a bag of heterogeneous risk loci does not
  have one. Disease relevance therefore enters this atlas **honestly, per
  nomination**, through Open Targets genetics — never by constructing a
  disease-shaped target.

### 5.7 Robustness

- **Held-out-gene validation on all 12 atlas cells (z 4.7–45.0)** — the primary
  significance test, and the one that carries the "signal is real" claim: fit
  weights on half the signature genes, score on the held-out half, above a
  paired shuffled-target null. Every cell clears it.
- **Random-perturbation null** (headline pert z = −18.4; atlas −23 to −55). Note
  the **sign is negative**: this null shuffles gene labels *within each
  perturbation column*, which decorrelates the dictionary and inflates the cone,
  so the shuffled generators fit the target *better* (null cosine ≈0.83) than the
  real ones (≈0.63). It therefore does **not** corroborate the verdict — it is a
  diagnostic that the measured cone is structured (correlated columns) rather than
  a random spanning set, and a caution that a fatter, decorrelated basis would
  over-reach. The load-bearing nulls are the shuffled-*target* and held-out-*gene*
  tests above, both strongly positive.
- Gene-panel subsampling stability on the headline verdict (toward_Th1 Rest):
  point estimate reachable cosine 0.627 / LOF fraction 0.393; resampling 85% of
  the signature genes (×15) gives subsampling mean 0.637 (95% interval
  [0.629, 0.648]) and LOF mean 0.406 ([0.396, 0.420]). The intervals are tight
  (width ~0.02) and sit just above the point estimate — the verdict is stable to
  the gene panel (the small upward shift is expected: subsampling drops the
  least-fittable genes, so each 85% draw fits marginally better than the full set).
  (`analysis_cache/atlas_work/bootstrap_ci.json`.)
- **Leave-one-donor-out is not possible** — the effect matrix is donor-collapsed
  (no per-donor effect vectors are local). Stated as a limitation; the gene-panel
  bootstrap is the honest substitute rather than a fabricated donor split.
- **Additivity risk is calibrated, and the atlas recipes are safe.** The cone
  composes single-knockdown effects additively. That assumption was calibrated
  against the 126 measured Norman double perturbations: the intuitive
  collinearity → sub-additivity mechanism is **refuted** (Spearman ρ ≈ −0.16,
  n.s.), and the real mechanism is **magnitude saturation** (deficit grows with
  combined magnitude, ρ ≈ +0.58; scale-free ceiling M\* = 13.9). Scoring each of
  the 12 atlas recipes with the resulting `additivity_risk()` gives risk
  **0.04–0.08 — far below the 0.5 unsafe threshold** — because the greedy recipes
  spread small weights over several generators and never approach the ceiling.
  The epistasis-aware selector returns identical recipes at zero cosine cost, so
  the additivity assumption is **validated for these targets**, not merely
  asserted (`results/atlas_additivity_risk.csv`,
  `norman_additivity_calibration.csv`; §8.2.2–§8.2.3).

### 5.8 External corroboration

Two independent checks that the nominations are real biology, not fitting noise
(both read off `results/modality_intervention_map.csv`):

- **The green-light nominations are already clinically pursued.** The method's
  top drugged nominations coincide with heavily-worked immune targets: JAK2 (31
  drug/candidate records; JAK inhibitors), MAPK14 (28; p38 inhibitors), CD3D
  (16; anti-CD3), RARA (11). The reachability cone is not surfacing obscure
  genes — where it says green-light, the pharmacopoeia agrees.
- **Literature positive control.** A hand-curated panel of 15 canonical T-cell
  activation / polarization regulators (ZAP70, ICOS, CD3D, CD5, IL7R, RELA,
  NFAT5, MAPK14, IRF1, IRF9, RARA, IFNGR1, TET2, JAK2, KEAP1) is recovered
  **15/15** as atlas nominations — an independent confirmation that the greedy
  knockdown recipe recovers established T-cell biology.

*(The Open Targets `knownDrugs` GraphQL field schema-mismatched this API version;
drug **counts** from the working `drugAndClinicalCandidates` field are used
instead. External Perturb-seq/CRISPRa datasets, e.g. Schmidt 2022, were left as
future corroboration — the sandbox reaches Open Targets but not those hosts.)*

---

## 6. Generalizability — does the method transfer, and what is it good for?

*Findings and application map for the convex-cone reachability oracle (absorbing the former `GENERALIZABILITY.md`). The live cross-dataset demo is reproduced end-to-end by `notebooks/03_generalizability_and_impact.ipynb` (K562 CRISPRa); the cross-cell-type transfer by `notebooks/07_cross_celltype_transfer.ipynb` (K562/RPE1). Companion data: `results/dataset_catalog.csv`, `results/tractability_grounding.csv`, `results/norman_table1..5_*.csv`, `results/generalizability_survey.md`, `analysis_cache/czi_data/per_perturbation_transfer.csv`, `analysis_cache/czi_data/CROSS_CELLTYPE_TRANSFER.md`.*


### 6.1 The scope is an input contract, not an assay

The method is **not confined to the repo's CD4⁺ T-cell CRISPRi dataset.** Its scope is fixed by
an *input contract* — a measured perturbation-effect matrix **E** (P × G) and a target direction
**d** in the same G-gene space — not by any assay. We proved transfer by running the **identical
`reachability.py`, unchanged**, on a screen that differs on three axes at once:

| | Repo dataset (notebooks 01–02) | Cross-dataset demo (notebook 03) |
|---|---|---|
| Source | `GWCD4i.DE_stats.h5ad` | **GSE133344** — Norman, Weissman et al. 2019, *Science* (PMID 31395745) |
| Cell system | Primary human CD4⁺ T cells | K562 (CML line) |
| Modality | CRISPRi (knockdown, loss-of-function) | **CRISPRa (activation, gain-of-function)** |
| Design | Single perturbations | **Combinatorial** (single + double guides) |

On the K562 CRISPRa screen the oracle returned a certified, null-calibrated, honest verdict — and,
because the data is combinatorial, it also let us **measure** the one assumption the T-cell screen
could never test.

### 6.2 The live demo (K562 CRISPRa, held-out CEBPA state)

**Reachability verdict — the CEBPA master-TF state is PARTIALLY-REACHABLE without CEBPA.**
Asking whether the CEBPA-overexpression state can be reproduced from the *other* 233 perturbations:

- reachable cosine **0.878**, residual **0.479**, support 86 perturbations
- KKT/Farkas max violation **1.7 × 10⁻¹²** → the cone fit is certified optimal
- shuffled-target null p95 = 0.411; observed at the 100th percentile → **null-z ≈ 37**
- held-out-**gene** validation cosine 0.856 (null +0.20), **z ≈ 23.5** → the fit is not overfitting

**The method discriminates — it is not a rubber stamp.** The same held-out construction across
master TFs spans **null-z ≈ 3 (IRF1, at the outside cutoff) → 37 (CEBPA) → 39 (ETS2)**, all via the
identical `reachability.py` null. A positive verdict is informative because the oracle can and does
say *outside*.

**Minimal recipe is biologically coherent.** Greedy selection picks **CEBPE** — CEBPA's own
C/EBP-family paralog — as the single best surrogate (cos 0.817); a 2-perturbation recipe reaches
**96 % of the full-cone fit** (knee at k = 2; 0.843 / 0.878).

**The certificate names what is missing.** The residual is a specific CEBPA-driven
myeloid-differentiation program no combination supplies: **MNDA, HP, VSIG4, ALOX5AP, PILRA, JAML,
NCF1, SIGLEC14**. This — *what is missing*, not just *how close* — is the method's distinctive output.

**Additivity is a bounded approximation, measured directly.** Across 126 testable doubles, the
measured double vs the sum of its two singles has **median cosine 0.71** (18 % well-predicted at
cos > 0.8; 16 % strongly non-additive at cos < 0.6). The most non-additive pairs are coherent
genetic interactions (mitotic kinesins KIF18B+KIF2C; apoptotic BAK1+BCL2L11). The single-perturbation
T-cell screen *structurally could not* run this validation; the combinatorial screen turns the
additivity assumption from an article of faith into a reported number with flagged exceptions.

---


### 6.3 Cross-cell-type transfer: is reachability a property of the biology or of one cell type's basis?

*Notebook `07_cross_celltype_transfer.ipynb`; full writeup `analysis_cache/czi_data/CROSS_CELLTYPE_TRANSFER.md`.*

The K562 demo above shows the *code* transfers to a new assay. A sharper question is whether the
**cone geometry** is a property of the target biology or of the specific cell type's measured basis.
The CZI Virtual Cell Models platform re-hosts the **Replogle 2022** genome-scale essential-gene
CRISPRi screens in two human cell types — **K562** and **RPE1** — which lets us test this directly.
We built row/column-aligned effect matrices `E_K562` and `E_RPE1` over the **843 single-gene
perturbations and 2,832 readout genes shared by both**, then ran the unchanged `reachability.py`. On-target
self-effects are negative for **100 %** of perturbations in both cell types, confirming the bases are
correctly built. The result is a three-level answer:

| Level | Question | Verdict | Evidence |
|---|---|---|---|
| **Effect direction** | Does a single knockdown do the same thing in both cell types? | **Transfers (moderately)** | Matched cross-type cosine median **+0.35** vs ≈0 shuffled-gene null; survives deflating the shared essential-stress direction (96 % stay positive, 69 % gene-specific at p<0.05) |
| **Reachability verdict** | Is a target reachable in both? | **Transfers (~100 %), but for a subtle reason** | Cross-type reachable cosines 0.50–0.73 stay above the null, so the binary verdict agrees — but within-type residuals (K562 0.61, RPE1 0.37) beat cross-type (0.87, 0.68), so the honest, discriminating signal is the **graded residual**, which does carry a cell-type penalty |
| **Minimal recipe** | Are the *same knockdowns* the answer in both? | **Does NOT transfer** | Same-target recipes overlap at median Jaccard **0.11** — ≈20× above a random-subset null of **0.006** (so not random, but far from identity); and reaching the *same* target from the *other* cell type's basis collapses overlap to **0.05** |

**What it means for the Th2→Th1 headline.** This is a robustness result with a sharp boundary. It
**defends the method** — the convex cone produces coherent, above-null, partly cell-type-invariant
structure in *three* human cell types (CD4⁺ T, K562, RPE1), so it is not an artifact of one dataset.
It **bounds the prescription** — the *direction* toward a target state transfers, but the *specific
minimal recipe* does not, so the GATA3↓/TBX21↓-style recipe is validated **for CD4⁺ T cells** and
must be re-fit on another cell type's basis before it is trusted there. This is exactly what the
project's philosophy (trust null-calibrated, held-out metrics; never the raw in-sample cosine) predicts.
*Caveat: K562/RPE1 are two non-T lines, so this is evidence of general cell-type-invariance, not a
re-measurement of the Th2→Th1 recipe in a second T-cell system (no second genome-scale CD4⁺ T screen exists).*

---


### 6.4 Three assumptions, and where each holds or breaks

1. **Non-negative weights (w ≥ 0).** "You can apply a perturbation or not, not its negative."
   Exactly right for single-modality screens (all-knockdown or all-activation); it is what makes the
   reachable set a *cone* and yields the Farkas certificate. The `signed_reachability` extension
   admits the sign-flip as a separate generator to quantify how much of a target needs the *opposite*
   modality.
2. **Additivity of co-perturbation (double ≈ e_A + e_B).** A standard linear prior, but an
   approximation — and, per above, one this dataset lets us bound at median cosine 0.71 rather than
   assume.
3. **Linearity in a shared gene space.** Effects are fixed vectors in one coordinate system; target
   and generators must share the gene axis. Cross-context transfer (resting vs stimulated, dataset to
   dataset) is folded into the residual, which the shuffled-target null and held-out-gene validation
   are designed to expose rather than hide.

---


### 6.5 Opportunity: other public datasets (13 real accessions, 5 modalities)

Every accession was retrieved live from NCBI GEO this session — none asserted from memory. Full
table with organisms, sizes, and access notes in `results/dataset_catalog.csv`.

| Modality family | Representative accessions | Reachability question it enables |
|---|---|---|
| CRISPR Perturb-seq (RNA) | GSE314342 (repo), GSE90546 (Adamson UPR), GSE90063 (Dixit), GSE124703 (iPSC-neurons) | Is a stimulus/disease state reachable by a knockdown/knockout mix; which genes must instead be activated? |
| CRISPR + protein readout | GSE153056 (ECCITE mixscape), GSE278572 (Treg/Teff Perturb-CITE) | Same, with a joint RNA+protein target direction. |
| CRISPRa / combinatorial | **GSE133344 (demo)**, GSE146194 (Replogle paired guides) | Additivity of doubles; reachability of a state on an activation basis. |
| ORF overexpression | GSE216463 (Joung TF atlas, ~1,800 TFs) | Which TF-overexpression combination reaches a target fate (directed differentiation). |
| Chemical / L1000 | GSE139944 (sci-Plex), GSE92742 / GSE70138 (LINCS L1000) | Minimal drug combination reaching a state; reverse of a disease signature. |
| Cytokine / ligand | GSE202186 (Immune Dictionary, 86 cytokines) | Which cytokine combination reaches a target immune-cell state. |

---


### 6.6 Application map — four directions, one contract

1. **Cross-cell-type screens** — input: any pooled knockdown/knockout effect matrix; target: a
   polarization/stimulus shift; decision: reachable by knockdown in *this* cell type, the minimal
   guide set, and the activation-only residual if outside. (Proves cell-type and chemistry
   independence.)
2. **Drug / small-molecule signatures** — input: drug effect vectors (sci-Plex, LINCS L1000);
   target: a desired state or **d = −(disease signature)**; decision: minimal non-negative drug
   combination, or a certificate of genes needing a non-drug modality.
3. **Cell reprogramming / fate control** — input: a gain-of-function basis (TF-overexpression atlas,
   cytokine panel); target: a developmental/therapeutic fate; decision: the directed-differentiation
   recipe, and the certificate for fates no measured activator spans.
4. **Disease reversal / target discovery** — input: any measured basis in the disease-relevant cell
   type; target: **d = −(disease signature)**; decision: the reversal recipe, plus a certificate that
   — filtered by tractability — routes each required gene to its intervention class.

**Tractability grounding (Open Targets / ChEMBL, live).** The certificate's genes sort cleanly by
how you would act on them: TF-activation targets **TBX21, STAT4, FOXP3, GATA3** have **zero
small-molecule tractability and zero approved drugs** → cell-engineering / CRISPRa problems; enzyme
targets **JAK1 (25 drugs/candidates)** and **DGKA (ligandable)** → drug-repurposing / medicinal-chemistry
problems. Detail in `results/tractability_grounding.csv`.

---


### 6.7 Why the honesty ports

The shuffled-target null (a per-target bar for *achievable-by-chance* reachability when P ≪ G) and
held-out-gene validation (the guard against a dense non-negative fit) are properties of the
**algorithm**, not of the T-cell data. So the falsifiable reachable/outside verdict — and the
constructive certificate when a state is out of reach — travel to every dataset above. The
demonstrated instances — Th2→Th1 in CD4⁺ T cells, CEBPA in K562 CRISPRa, and the K562↔RPE1
cross-cell-type transfer across three human cell types — are the existence proof; the catalog is
the opportunity.

---

## 7. The experimental-design toolkit — from a reachability verdict to a screen you can run

*Absorbs the former `DESIGN.md`.*

**From a reachability verdict to a screen you can run.**

Notebooks 01–03 built and stress-tested a convex-cone *reachability oracle*: given a dictionary
of single-perturbation effect vectors and a desired cell-state shift, it asks whether a
non-negative combination of available perturbations can reproduce the target direction, and
returns a signed decomposition of the answer. Notebook 04 turns that oracle into a **design
tool** a screen planner can actually use, and — the point that matters most — it **quantifies
how much the method's headline number can be trusted**.

`design_experiment(E, d)` takes a current→target transition and returns a **design card**:

1. a null-calibrated reachability **verdict** (four levels),
2. ranked **knockdown** and **activation** recipes from the signed decomposition,
3. a per-move **delivery-technology call** grounded in Open Targets tractability, and
4. an optimal next-screen **library** (the knee of the sparsity-vs-reach curve).

The toolkit is exercised on the genome-wide CD4⁺ T-cell CRISPRi screen across **12 transitions**
(4 target states — Th1, Th2, younger, older — × 3 culture conditions — Rest, Stim8hr, Stim48hr).

---

### 7.1 Why this is the right next step for the field

A reachability score on its own is a diagnostic. A screen planner needs three further things,
and none of them are answered by a single cosine:

- **Which perturbations, and in which direction?** A target shift generally decomposes into a
  part reachable by *knockdown* and a part that requires *induction*. These are different
  experiments (CRISPRi vs CRISPRa/ORF). Reporting them together as one "recipe" hides the
  experiment you cannot run.
- **How is each move actually delivered?** "Perturb gene X" is not a protocol. Whether X is a
  small-molecule target, an antibody-addressable surface protein, or a genetic-only handle
  decides the technology and the cost.
- **How big a screen, and how much should I believe it?** The greedy sparsity curve gives an
  optimal library size; the calibration analysis gives the honest expectation for what that
  library will achieve.

The toolkit answers all three, and refuses to let the optimistic in-sample number stand
unqualified.

---


### 7.2 The three scientific results

### 1 · Every T-cell transition needs a modality the knockdown screen cannot supply

The signed decomposition splits each target's norm into knockdown-reachable (LOF),
activation-only (GOF), and neither. Across all 12 transitions, **22–26 % of the total shift is
activation-only** — a direction a CRISPRi (knockdown) screen *structurally cannot reach*, no
matter how many guides it includes. In every transition the activation-support gene set is
larger than the knockdown-support set. This is the central design message: a knockdown-only
screen is, geometrically, the wrong instrument for a large share of these state changes.

![Modality triage across all 12 transitions]({{artifact:art_d8a5ea7a-b82e-4e78-89c5-837a6b14a937}})

### 2 · Grounding the recipe in real delivery technology

The unique recipe genes across all transitions (**270 genes**, **360 moves**) were resolved
against Open Targets. Two findings shape a screen:

- **14 knockdown targets sit in the Open Targets "clinical drug" tractability tier** (an
  approved or clinical-stage drug exists) — immediate repurposing / chemical-genetics
  candidates: ACACA, CALM1, CCR4, DPYD, ICOS, JAK2, KIF5B, MAP3K10, MAPK14, MEN1, PDE4B, PPP5C,
  RARA, VKORC1. (Two further genes, IFNGR1 and IL10RB, have `n_drugs > 0` but fall in weaker
  tiers — SM-ligandable and AB-surface respectively — so they are *not* counted here; the
  distinction is preserved in `design_modality_tractability.csv`.)
- **68 % of activation moves are not small-molecule-druggable** — the induction direction
  resolves to CRISPRa or ORF overexpression, not a compound. Open Targets' small-molecule and
  antibody handles enable *blocking/degradation*, not *induction*, so tractability annotations
  must be read direction-aware. The toolkit does this automatically.

### 3 · The in-sample reachable cosine systematically overstates real reach

This is the honesty result, and it is the reason the notebook exists in this form.

![Reliability diagram]({{artifact:art_be72897c-a66c-46aa-9185-206be61eb765}})

Across all 12 transitions the in-sample reachable cosine (mean **0.580**) overstates the
held-out-gene realized cosine (mean **0.355**) by a mean gap of **0.225** (range 0.18–0.27).
Every point sits below the identity line — the overstatement is **systematic bias, not
estimation noise**, and the gene-panel bootstrap (below) confirms it is not a resampling
artifact. A screen planned on the in-sample number would be planned on an expectation ~0.22
cosine-units too high.

Each card therefore reports the **held-out** number as the honest estimate and carries a
**confidence label**: **2 HIGH** (toward_Th1_Rest, toward_Th2_Rest), **6 MEDIUM**, **4 LOW**.

![Design-card summary]({{artifact:art_922fcad1-b4dc-43a7-8d8c-4e1706c6dfe8}})

---


### 7.3 The calibration caveat (read this before trusting a verdict)

**True leave-one-donor-out (LODO) validation is impossible on these inputs**, and it is
important to be explicit about why. The local knockdown effect vectors `E` are
**donor-collapsed**: each perturbation is a single pooled effect vector, not a set of per-donor
vectors. There is no donor axis left locally to hold out, so we cannot measure whether a recipe
fit on some donors predicts the effect in a held-out donor — the validation that would most
directly speak to reproducibility in a new experiment.

In its place the toolkit uses three **honest substitutes**, each measuring something real but
none equal to LODO:

- **Held-out-*gene* reliability** — fit the NNLS weights on half the readout genes, score the
  cosine on the other half. This is the load-bearing number; it catches a dense fit that reaches
  the target in-sample by overfitting correlated generators. It is *not* a held-out *donor* and
  *not* an experimentally realized outcome.
- **Gene-panel bootstrap CIs** — 85 %-without-replacement resamples of the signature genes
  (N = 12) on the 4 Rest transitions give tight intervals on the reach cosine and the
  LOF/GOF/neither split, showing the point estimate is stable under gene-panel perturbation.
  Note these intervals bracket the **bootstrap resample mean**, which sits slightly *above* the
  full-panel point estimate (resampling 85 % of genes fits marginally tighter), so the CI is a
  *robustness* statement about the resampled estimator, not an error bar on the full-panel
  cosine itself — the design cards label it as such and never juxtapose the two as "point ± CI".
- **Nulls** — the held-out-gene shuffled-target null (primary significance) plus a
  random-perturbation null. The latter gives cosines *above* the observed (negative z), which is
  itself a **structure diagnostic**: permuting gene labels within a generator column destroys the
  biological covariance the real dictionary carries, so a shuffled dictionary can reach a
  direction-biased target *more* easily — confirming the observed reach rests on real structure,
  not label coincidence.

A further honesty fix worth recording: **cross-donor QC is sparse.** The screen's
`donor_correlation_all_mean` statistic is populated for only ~14 % of perturbation-conditions,
and it is **neither necessary nor sufficient** for on-target significance (81 % of
on-target-significant perturbations have no donor-correlation value). It is therefore used only
as a partial, best-available reproducibility weight on recipe genes — with the covered fraction
always reported — and never as a reproducibility *flag*.

**Bottom line:** the confidence labels and CIs are the best available substitute for LODO on
donor-collapsed data. They are not a claim that a HIGH-confidence card will reproduce in a new
donor cohort — only that its reach is stable to gene-panel resampling, significant against the
target-specific null, and less overstated than a LOW card. Plan accordingly.

---


### 7.4 What's in the toolkit

#### API (in `reachability.py`)

```python
from reachability import design_experiment
card = design_experiment(
    E,                       # (P, G) single-perturbation effect dictionary
    d,                       # (G,)   desired current->target state shift
    perturbation_names=...,  # length-P gene names (for the recipe)
    readout_names=...,       # length-G gene names (for the certificate)
    hvg_mask=(d != 0),       # readout genes that define the target
    k_max=12, top=15, n_shuffles=20, seed=0,
)
card.summary()               # one-line headline
card.verdict                 # outside | weakly | partially | reachable
card.knockdown_recipe        # [{gene, weight, rank}, ...]   CRISPRi / inhibit / degrade
card.activation_recipe       # [{gene, weight, rank}, ...]   CRISPRa / ORF / agonize
card.library                 # [{gene, k, cosine_at_k, marginal_gain}, ...] up to k_max
card.optimal_k               # knee of the greedy spectrum
```

`design_experiment` composes the existing `signed_reachability`, `held_out_gene_validation`,
`reachability_spectrum`, and `activation_certificate` primitives, adds a dependency-free knee
finder (`_knee`) and a four-level verdict grader (`_grade_verdict`), and returns a single
`DesignResult` dataclass. It is covered by the module self-test (`python reachability.py`).

#### Notebook (`notebooks/04_experimental_design_toolkit.ipynb`)

22 cells; runs in ~6 s. It runs one **live** `design_experiment()` on a small synthetic in-cone
target (to show the API end-to-end and fast), then loads the 12 precomputed T-cell cards (each
takes 1–4 min live because of the held-out loop). Includes `render_design_card()`,
`design_card_markdown()`, and an **interactive ipywidgets picker** (`build_picker()`) with
one-click Markdown export. The picker degrades to a static card when `ipywidgets` is unavailable,
so the notebook always runs.

#### Result tables (`results/`)

| file | rows | what |
|---|---|---|
| `design_cards.csv` | 12 | one flat design card per transition (verdict, confidence, reach, gap, recipes, library) |
| `design_modality_tractability.csv` | 360 | per-(transition, gene, direction) move with OT buckets, drug count, delivery call |
| `design_modality_summary.csv` | 12 | per-transition modality split + druggable counts |
| `design_calibration.csv` | 12 | predicted vs realized cosine, gap, bootstrap CIs, nulls, donor-repro weight, confidence |

#### Figures (`notebooks/figures/`)

- `nb04_fig1_design_summary.png` — in-sample vs held-out reach for all 12 transitions (confidence-colored)
- `nb04_fig2_modality_triage.png` — LOF/GOF/neither split per transition
- `nb04_fig3_reliability.png` — reliability diagram (every point below identity)
- `nb04_fig4_library_curve.png` — cumulative reach vs library size with knees marked

#### Exported cards (`notebooks/cache/cards_export/`)

All 12 transitions as standalone Markdown design cards.

---


### 7.5 Reproducing

```bash
# from the repo root, in the `cellreach` environment
python reachability.py                                   # self-test (includes design_experiment)
cd notebooks && python _nb04_script.py                   # runs every notebook cell (EXIT 0)
```

The 12 cards, calibration table, and modality table are precomputed (the held-out loop and the
NNLS bootstrap are the expensive parts). Their build logic is the code in `_nb04_script.py`
together with the two analysis tracks that produced `design_calibration.csv` and
`design_modality_tractability.csv`; the notebook itself only loads and presents them.

---

## 8. Post-hackathon method advances

*This section merges the former `RESULTS.md §6` (the methodological deepenings validated against data) with the whole of `METHOD_IMPROVEMENTS.md` (the "novel algorithm vs. novel formulation" answer, the nine post-hackathon in-silico results, and the methodological + experimental agendas). Deduplicated: the positive-control, epistasis, anisotropy-null, and DEG-weighted results live once, in §8.2; the nine in-silico results in §8.3 reference them rather than restating them.*

### 8.1 Is this a novel algorithm or a novel formulation?

*Companion to `limitations_and_reinforcement_plan.tex` (the L1–L8 self-critique),
`NOVELTY.md` (positioning), and `RESULTS.md` (delivered results). That trio answers "what
is new" and "how could it be wrong." This section answers the two questions an expert asks
next: **(1) is the pipeline a novel algorithm or a novel formulation, and how do we make
the method itself stronger; (2) if we were the pharma team deciding whether to act on a
nomination, what evidence would we demand.** It ships **nine new in-silico results** that
were run to close the sharpest of those gaps, and lays out the rest as a ranked agenda.*

> **Update — full methodological + experimental agenda now executed.** The first pass shipped
> two results (§8.3.1 non-negativity ablation, §8.3.2 collateral specificity). A second pass then
> ran the remaining computationally-tractable items end to end: **generator-uncertainty
> propagation (§8.3.3), dictionary effective-rank/conditioning (§8.3.4), a group-sparse cone that
> unifies fraction and recipe (§8.3.5), a two-part activation-certificate validation (§8.3.6), a
> directional-genetics cross-check (§8.3.7), ChEMBL mechanism-of-action grounding (§8.3.8), and a
> forward-predictor head-to-head (§8.3.9).** The only items that remain open require a wet lab
> or a GPU-hosted foundation model; they are marked as such in §8.4–§8.5. Every result below is
> backed by a saved CSV and a publication-grade figure.

---

**The honest one-paragraph answer.**

The core machinery is **classical**: non-negative least squares, convex-cone membership,
Farkas/KKT duality, separating hyperplanes. None of that is new mathematics, and the paper
must not sell it as a new optimizer. What *is* new is the **reduction** — recognising that
*measured* CRISPRi loss-of-function vectors plus the biological fact that you cannot apply a
negative knockdown make cell-state reachability *exactly* convex-cone membership, which
yields a **feasibility verdict with a constructive infeasibility certificate** that none of
the 91 prior surveyed methods produce. So the correct classification is **a novel formulation
and pipeline, carrying one genuinely new algorithmic contribution** — the closed-form
anisotropy-corrected null (`analytic_anisotropy_null()`, `§8.2.4`), which replaces
~1000 shuffle refits with a validated closed form (Pearson 0.995 vs the empirical null).
Everything below either sharpens that formulation or states what it still owes a reviewer.

The nine results in §8.3 now let us say something stronger than "the formulation is sound":
**the formulation survives the checks a critic reaches for first.** The verdict is robust to
generator measurement noise (§8.3.3); the near-degenerate dictionary is a *stated* geometric
property, not an apology (§8.3.4); the two-object dense/greedy split collapses into one
group-sparse fit that keeps ≥95% of the reachable cosine (§8.3.5); the most novel output — the
activation certificate — has a reproducible ranking and is partially corroborated by
literature and human genetics, with an honest caveat about which genes it names (§8.3.6); and
the pharmacological grounding surfaces exactly the kind of *negative* result (wrong-direction
drugs, discordant genetics) that a target-selection team must see before acting (§8.3.7–§8.3.8).
Where the oracle is compared head-to-head on prediction, the only thing that beats it is a
model that violates the non-negativity physics (§8.3.9).

---


**Results index (the nine in-silico results of §8.3).**

| §    | Result                              | Headline finding                                                                 | Gap closed        |
|------|-------------------------------------|----------------------------------------------------------------------------------|-------------------|
| 1.1  | Non-negativity ablation             | Sign constraint costs ≤0.04 cosine, buys the whole certificate; OLS edge is 50% impossible activations | L4 (method)       |
| 1.2  | Collateral specificity              | Th1 recipes age the transcriptome (+0.09, z = 7.1); mean 26% off-target leak     | new (pharma)      |
| 1.3  | Generator-uncertainty propagation   | Verdict never flips under dictionary noise (flip rate 0.0 / 12 cells)            | unlisted (method) |
| 1.4  | Effective rank & conditioning       | Stable rank ~10–15, participation rank ~3,400; cond. 379–706                     | new (method)      |
| 1.5  | Group-sparse cone                   | One fit gives fraction + recipe, ≥95% cosine retained, GATA3 resolved            | method            |
| 1.6  | Certificate validation              | Ranking reproducible (z ≈ 35); 4 SUPPORTS / 7 PLAUSIBLE / 9 UNKNOWN; top genes are markers | L1 in-silico half |
| 1.7  | Directional genetics                | 2 concordant, 7 indeterminate, IRF1 wrong-direction — count ≠ direction          | new (pharma)      |
| 1.8  | ChEMBL MoA grounding                | 6/10 LOF-direction drug confirmed; IFNGR1 & RARA agonist-only (wrong direction)  | new (pharma)      |
| 1.9  | Forward-predictor head-to-head      | Cone leads realizable methods; ridge wins only via 50% negative weights          | L3 (pharma)       |

---

### 8.2 The four data-validated methodological deepenings

Four items from the "what to improve" list (§4), the novelty roadmap, and the evaluation
literature were taken from sketch to validated method after the hackathon build. Each is
calibrated or validated against data, and each is wired into `reachability.py` behind its
self-test.

#### 8.2.1 Positive control, broadened from 2 genes to a regulator panel

The original positive control rested on two genes (GATA3↓, TBX21↓). It now runs against a
curated panel of Th1/Th2 regulators. The result has structure worth stating honestly:

- **Master transcription factors behave exactly as biology predicts.** Splitting the panel into
  master-TF drivers (Th2: GATA3, STAT6, BATF, MAF; Th1: TBX21, STAT4, STAT1, RUNX3) versus
  induced cytokine receptors/markers (IL4R, CCR4, IL12RB2, IL18R1, IL18RAP), **every Th2-driver
  knockdown outranks every Th1-driver knockdown on the toward-Th1 axis (AUROC = 1.00,
  rank-sum p = 0.014; 7/8 signed-concordant with biology, binomial p = 0.035).**
- **The full 13-gene panel is only weakly enriched** (Th2-drivers-at-top AUROC = 0.69,
  p = 0.052), because the cytokine *receptors* are induced markers of a state, not drivers of
  it — knocking them down does not move the cell along the polarization axis, and several sit
  near the distribution median. This is a real, interpretable limitation of a marker-based
  panel, not a failure of the method.
- **Caveat (stated plainly):** the master-TF-vs-marker split was constructed *after* seeing the
  weak full-panel result — it is an exploratory refinement, not a pre-registered test. The
  clean AUROC = 1.00 should be confirmed on an independent axis or dataset before being treated
  as confirmatory. Both the full-panel and TF-subset statistics are reported side by side so the
  reader can judge.

*Figure `fig_positive_control_enrichment.png`; tables `results/positive_control_enrichment.csv`
(13-gene per-gene alignment/rank/percentile) and `results/positive_control_stats.csv` (5 tests).*

#### 8.2.2 An epistasis/additivity penalty, calibrated on measured double perturbations

The cone composes single-knockdown effect vectors **additively**. That assumption was calibrated
against the one public screen that measures both singles and their doubles (Norman et al. 2019,
K562 CRISPRa; 126 doubles with both singles present). Two candidate epistasis mechanisms were
tested against the *measured* non-additivity of each double:

- The intuitive prior — **effect-vector collinearity** ("hitting the same program twice should be
  sub-additive") — is **refuted**: Spearman ρ ≈ −0.16 (n.s.) with directional non-additivity; if
  anything, collinear pairs are slightly *more* additive.
- The mechanism the data support is **magnitude saturation**: combined effect magnitude falls
  below the additive sum, and the deficit grows with the combined magnitude (Spearman ρ ≈ +0.58,
  p < 0.01). Fitting `achieved = a / (1 + a/(M*·s))` gives a scale-free ceiling **M\* = 13.9**
  (in units of the dictionary's median single-effect norm), R² = 0.57, ≈ 12 % mean magnitude
  loss. Dividing by the median single-effect norm makes the law transfer from Norman
  log-fold-change to a z-score dictionary.
- Directional (angular) non-additivity, by contrast, **improves** with effect size
  (median cos(measured, additive) 0.64 → 0.81 from the weakest to strongest magnitude quartile) —
  i.e. it is largely low-SNR measurement noise, so the calibrated penalty corrects the
  **magnitude channel only** and leaves the fit direction (hence the reachable cosine and the
  verdict) intact.

This is shipped as `additivity_risk()` (a per-recipe risk score in [0,1)) and an
`epistasis_penalty` option on `reachability_spectrum` (0.0 reproduces the additive selection
bit-identically). The risk score validates out of the box against measured deficits:
**Spearman ρ = +0.46 (p = 5.6e-8)**; high-risk pairs show a median magnitude deficit of +0.18
vs −0.03 for low-risk pairs (MWU p = 4.5e-6).

#### 8.2.3 Re-annotating the atlas recipes with additivity risk — and finding they are safe

Applying the calibrated risk score to all 12 atlas recipes (k = 7) gives a reassuring **null
result**: every recipe carries **additivity-risk 0.04–0.08 — far below the 0.5 unsafe threshold**
— because the greedy recipes spread small non-negative weights over several generators and never
approach the saturation ceiling. The epistasis-aware selector therefore returns **identical
recipes at zero cosine cost** for the atlas; the penalty only changes selection under a
deliberately tightened ceiling (stress test). The scientific statement is now *earned* rather than
assumed: **for these targets, additivity is a calibrated, validated approximation, not a
liability.**

*Figure `fig_additivity_risk.png` (saturation-law fit + per-recipe risk); tables
`results/atlas_additivity_risk.csv` (12 recipes with risk + gene list) and
`norman_additivity_calibration.csv` (calibration record).*

#### 8.2.4 A closed-form, anisotropy-corrected null

The shuffled-target null was, until now, purely empirical (~1000 refits per target), and the
aging axis exposed why a naive z-against-zero is wrong (§3.4): a target whose values are mostly
one-signed keeps a large uniform ("DC") component under shuffling, which a non-negative cone
reaches for free. That baseline now has a **closed form**. Decomposing a shuffled target into a
uniform component (fraction `a² = G·mean(d)²/‖d‖²`, preserved by shuffling) plus a mean-zero
residual, reached through two orthogonal channels, gives

> **E[null cosine] ≈ √( (a·ρ₁)² + (1 − a²)·κ )**

where `ρ₁` is the dictionary's reachable cosine of the uniform direction and `κ` is the chance
reachable cosine² for a random mean-zero target (the isotropic "AC floor"). Two wrong hypotheses
were ruled out first (DC-fraction alone — Th1 has null 0.35 with a² = 0.002; and the
participation ratio — effective generator count tracks P, not PR). The law is validated against
empirical shuffled nulls across synthetic dictionaries and the real CD4 Rest dictionary
(four axes), in both the in-sample and held-out gene-split regimes: Pearson(pred, empirical) =
0.995 in-sample / 0.998 held-out, with a real-data max abs error of 0.047 in-sample / 0.070
held-out. It reproduces both the ≈ 0 null of the sign-balanced Th1/Th2 axes and
the elevated ≈ 0.26–0.34 null of the up-dominated aging axis. It is shipped as
`analytic_anisotropy_null()`, returning an **anisotropy-corrected z** from ~10–20 fits instead of
~1000 shuffles, and gated by the module self-test.

*Figure `fig_anisotropy_null.png`; table `results/anisotropy_null_validation.csv` (20
synthetic + real validation rows, both regimes: RMSE 0.029 in-sample / 0.045 held-out).*

#### 8.2.5 DEG-weighted evaluation — the verdict is not a signal-dilution artifact

The reachability verdict is scored with a cosine between the target shift and its closest
reachable point, taken **over all genes**. Mejia et al. (*Needles in the Haystack*, ICML 2026;
expanded preprint bioRxiv 10.1101/2025.10.20.683304) show that unweighted transcriptome-wide metrics are prone to
**signal dilution**: when a perturbation moves only a handful of genes, the score is dominated
by the many unchanged background genes, which can flatter a fit and let an uninformative
baseline look competitive. Their fix — DEG-aware metrics (weighted MSE / weighted ΔR²) reported
against explicit negative *and* positive controls — was adopted here as a direct robustness test
of our own metric. Three primitives were added to `reachability.py`, all non-breaking (the
default, unweighted path reproduces every number in this document bit-for-bit; §1 of notebook
08 and the module self-test both verify this on synthetic and real data):

1. **A DEG-weighted cosine** — the reachability score in the weighted inner product
   `⟨x,y⟩_w = Σ w_j x_j y_j` with `w_j = |d_j|` (the WMSE analog); optionally the fit itself is a
   weighted NNLS. This scores agreement where the perturbation actually acts and cannot be
   inflated by the quiet background — at the metric level, an identical-magnitude error scores
   ~12× higher when placed on the DEGs than on the background under weighting, whereas the
   unweighted metric cannot tell the two apart (notebook 08 §2).
2. **A positive control** — the *interpolated duplicate* (a known-reachable target built from
   the generators themselves): the metric awards it a near-maximal cosine (0.972 unweighted /
   0.987 weighted on Norman; 0.974 / 0.988 on Tier-2), confirming the metric *can* reward a
   truly reachable target, so a mid-range observed cosine is a statement about the target, not a
   ceiling of the metric.
3. **The dynamic-range fraction** `(observed − floor) / (ceiling − floor)` — a scale-free
   placement of the verdict between the shuffled-target floor and the positive-control ceiling.

**Result — the headline holds and strengthens.** For the Th2→Th1 target the DEG-weighted
reachable cosine *rises* in every condition (Rest 0.627 → 0.803, Stim8hr 0.524 → 0.736,
Stim48hr 0.533 → 0.737), and the held-out-gene z stays far above the z = 3 significance floor
throughout (28.3 → 14.1, 19.9 → 9.8, 19.2 → 8.6). The same pattern holds for six held-out
Norman double perturbations (cosine +0.07 to +0.11, every held-out z > 10). The reachability
signal therefore lives in the genes the perturbations genuinely move — it is not an artifact of
scoring the unchanged background.

**A subtlety the calibration surfaces.** For Tier-2 Rest the *raw* weighted cosine jumps
(0.627 → 0.803), but the **dynamic-range fraction barely moves** (0.447 → 0.499) — because the
shuffled-target floor *also* rises under weighting (0.346 → 0.619). In other words, weighting
lifts the reachable and the random-baseline scores together; the calibrated placement of the
verdict between its controls is stable. This is exactly the failure mode Mejia et al. warn about
— a raw metric shift that vanishes once controls are applied — and it is why we report the
calibrated fraction alongside the raw cosine rather than the cosine alone.

*Figures `notebooks/figures/fig5_deg_weighted_verdicts.png` (unweighted → DEG-weighted cosine
and held-out z per target) and `fig6_calibration.png` (floor/observed/ceiling placement +
dynamic-range fraction). Tables `results/deg_weighted_verdict_comparison.csv`,
`results/positive_control_ceiling.csv`, `results/calibration_dynamic_range.csv`. Full
reproducible analysis in `notebooks/08_deg_weighted_evaluation.ipynb`.*

---

### 8.3 Nine in-silico results closing the sharpest methodological and experimental gaps

#### 8.3.1–8.3.2 — the first two results (non-negativity ablation, collateral specificity)

The limitations plan lists a constraint ablation (L4) as priority #2 and, separately, the
appraisal that prompted this document flagged a *collateral-specificity* test that was on no
prior list. Both are computational, both are now done, and both are decision-relevant.

#### 8.3.1 The non-negativity constraint earns its place — and it is what makes a certificate possible (resolves L4)

The conceptual core is that knockdown combinations are **non-negative**, so reachability is
convex-*cone* membership rather than ordinary regression. L4 asked the fair question: does
the sign constraint change the answer, or is the cone decoration on a least-squares fit? We
ran, for all 12 atlas cells on the identical held-out-gene split the headline uses, three
fits — the non-negative cone (NNLS), unconstrained least-squares (OLS), and a
nearest-single-knockdown baseline — and added one diagnostic: the OLS solution **clipped to
the non-negative orthant**, i.e. the physically realisable part of what OLS prescribes.

![Non-negativity ablation across all 12 atlas cells. (a) Held-out-gene cosine to target for four fits. Unconstrained least-squares (purple squares) edges out the cone (green) by at most 0.04 — but that edge is bought entirely with pseudo-activations: projecting the OLS solution onto the physically realisable non-negative orthant (purple triangles) collapses it to ~0.13 and below zero in 4 of 12 cells. The nearest single knockdown (grey) is far worse, so the multi-gene combination does real work. (b) In-sample residual: the constrained fit always leaves a large residual (mean 0.75) — that residual is the infeasibility certificate — while the underdetermined OLS fits every target to machine zero and can therefore never declare a target unreachable.]({{artifact:9319d411-bde8-45ce-bd9d-7b1e02caf2b5}})

**What the ablation shows, in numbers:**

- **OLS's apparent edge is an illusion of unphysical coefficients.** Unconstrained
  least-squares beats the cone on held-out cosine in 10/12 cells, but by at most **0.044**
  (headline Th1/Rest: 0.489 vs 0.446). That margin is carried entirely by **negative
  coefficients — a mean of 50% of the OLS weight mass** — which prescribe *activating* genes,
  something CRISPRi physically cannot do. Clip OLS to the realisable non-negative orthant and
  its held-out cosine **collapses from a mean of 0.37 to 0.068**, going *negative* in 4 of 12
  cells (headline drops 0.489 → **0.131**). The cone's 0.446 is the honest ceiling of what a
  knockdown screen can actually reach; OLS's 0.489 is a number no wet-lab arm can realise.
- **The constraint is what makes infeasibility provable.** Because the dictionary is
  underdetermined (P > held-out genes), OLS fits **every** target to an in-sample residual of
  ~4×10⁻⁶ — it can never say "unreachable." The cone leaves a large residual (mean **0.755**),
  and that residual *is* the Farkas/separating-hyperplane certificate. The certificate is a
  direct consequence of the sign constraint, not an add-on.
- **The multi-gene cone does real work.** The nearest-single-knockdown baseline reaches only
  **0.177** on average vs the cone's 0.355 — the combination roughly doubles a single
  knockdown's reach, so the cone is not a dressed-up "pick the best gene."

The one-line version for the paper: *the sign constraint costs at most 0.04 of held-out
cosine and buys the entire ability to declare a target unreachable; the unconstrained fit's
apparent advantage is 50%-composed of activations that cannot be built.* Table:
[`ablation_nonnegativity.csv`]({{artifact:14743c9d-3d59-4439-a078-aeebc221027f}}).

#### 8.3.2 Recipes carry a small but systematic off-target movement — a specificity readout the verdict alone lacks (new)

A reachability verdict says a recipe moves the cell *toward* a target. It says nothing about
whether the same recipe **also drags the cell along an unwanted axis** — exactly the
specificity question a discovery team asks at nomination, not in Phase I. The atlas gives a
clean internal test: its four target axes form two near-orthogonal pairs — `toward_Th1` is
the exact sign-flip of `toward_Th2`, `toward_younger` the exact sign-flip of `toward_older`,
and polarization is orthogonal to aging (cosine 0.012). So *cross-pair* movement is a pure
off-target readout: a polarization recipe should not move the aging axis at all. We projected
every atlas recipe's achieved transcriptome shift onto all four axes.

![Collateral specificity of the atlas recipes. (a) Achieved movement (cosine of the recipe's realised shift) of each recipe (rows) along each axis (columns), averaged over the three culture conditions. The boxed diagonal is on-target movement; the off-diagonal cross-pair cells (polarization recipe → aging axis, and vice-versa) should be zero if recipes were specific. They are not: a toward-Th1 recipe carries a systematic +0.09 move along the *older* axis, and a toward-Th2 recipe a +0.07 move toward *younger*. (b) Collateral ratio (|dominant off-target| / |on-target|) per recipe. Polarization recipes leak onto the aging axis at a mean 33% of their on-target magnitude; aging recipes leak onto polarization at 19%.]({{artifact:47f725a0-0551-4ca7-9a28-042509523d12}})

**What the specificity test shows:**

- **The off-target movement is real, not a chance projection.** On the headline Th1/Rest
  recipe, movement along the orthogonal *older* axis is **+0.100**, against a
  shuffled-axis null of mean +0.037 (SD 0.009) — **z = +7.1, p < 0.0005**. The recipe
  provably ages the transcriptome while polarizing it.
- **It is directional and reproducible.** Across all three conditions, toward-Th1 recipes
  move the cell toward *older* (mean cosine **+0.094**) and toward-Th2 recipes toward
  *younger* (**+0.074**). A pro-inflammatory (Th1) knockdown program carries a small
  pro-ageing signature; a Th2 program carries a slightly youthful one. This is a falsifiable
  biological hypothesis the verdict alone would never surface.
- **The magnitude is material for triage.** The dominant off-target move averages **26% of
  the on-target move** across the atlas — **33% for polarization recipes**, 19% for aging
  recipes. A team told "this recipe is 45%-reachable toward Th1" should also be told "and it
  moves ~⅓ as far along the ageing axis," because that collateral may be the deciding safety
  fact.

This costs nothing new to compute — it reuses the existing recipes and effect matrices — and
it converts the oracle from a single-axis verdict into a **multi-axis specificity profile**.
Tables: [`collateral_movement.csv`]({{artifact:1a0eea70-4d2c-498b-a776-ef0444d921f3}}),
[`collateral_specificity_summary.csv`]({{artifact:ac417152-195e-490f-b913-e450cecbf795}}).

#### 8.3.3 The verdict is robust to generator measurement noise — the biggest *unlisted* gap, now closed (new)

The cone treats each effect vector as exact, but every generator is a noisy DESeq2 estimate
and the source `h5ad` ships its per-effect standard error (`lfcSE`) to prove it. A verdict of
"partially outside the cone" could in principle be an artifact of noisy generators sitting
near the boundary. The prior bootstrap resampled only the *target* gene panel; it never
resampled the *dictionary*. We closed that gap directly.

The noise model is not assumed — it is read off the data. On the atlas z-scale the screen
reports `zscore = log_fc / lfcSE` **exactly** (verified: max abs difference 0.0 across all
21,221 significant effects), so each dictionary entry carries **unit** measurement noise on
the z-scale. The dictionary bootstrap is therefore the parameter-free `E_boot = E + N(0,1)`
per element. We drew **B = 200** noisy dictionaries per cell (FISTA-accelerated NNLS, 12
cells) and re-ran the full reachability verdict on each.

![Generator-uncertainty propagation across all 12 atlas cells. Each cell's reachable fraction (point estimate and bootstrap 95% CI), held-out cosine, and verdict-flip rates under B=200 noisy dictionary draws with E_boot = E + N(0,1) on the z-scale. The reachable-fraction CIs are narrow and sit well below the 0.5 reachable/infeasible threshold for every cell; no draw crosses it (flip rate 0.0 everywhere).]({{artifact:0a77895a-83f6-4602-8962-3a96d3dbb9af}})

**What the propagation shows:**

- **The load-bearing verdict never flips.** The reachable/infeasible threshold (reachable
  fraction 0.5) is crossed by **zero** bootstrap draws in **all 12 cells**
  (`flip_rate_reachable_0p5 = 0.0` throughout). The reachable-fraction 95% CIs are narrow and
  far below 0.5 — headline Th1/Rest point **0.393**, CI **[0.399, 0.416]**. The claim "these
  states are only partially reachable by knockdown" is not a boundary artifact.
- **Significance survives the noise.** The held-out cosine carries a small, expected *negative*
  bias when the dictionary is degraded (Th1/Rest **0.446 → 0.349** bootstrap mean), but its
  z-score against the shuffled-target null stays strongly significant (**27.1 → 20.4**). Noise
  makes the number honest-smaller, not non-significant.
- **One honest caveat — the fine grade is not robust on the aging axis.** The *coarse* verdict
  is stable, but the *fine* grade (partially- vs weakly-reachable) flips in nearly every draw
  for four aging cells (`flip_rate_module_grade = 1.0` for younger/Stim8hr, younger/Stim48hr,
  older/Rest, older/Stim8hr; 0.04 for older/Stim48hr; 0.0 for all six Th1/Th2 cells). Those
  aging cells sit exactly on the partially/weakly grade boundary, so noise moves them across
  the sub-threshold — but never across the reachable/infeasible line. The paper should report
  the coarse verdict as robust and flag the aging-axis fine grade as noise-sensitive.

Table: [`generator_uncertainty_verdicts.csv`]({{artifact:9b37cb80-e41e-49c1-b470-f6263ca396ca}}).
This is the single most important remaining robustness check in the appraisal, and it now has
a distributional answer rather than a point estimate + null-z.

#### 8.3.4 The dictionary's near-degeneracy is a stated geometric property, not an apology (new)

With ~6,900 correlated generators acting in a ~6,200-gene signature subspace, the cone is
nearly degenerate — and that is *why* the greedy top genes (LAT2, RARA, ATF7IP2) are not the
canonical master TFs (`§3.2`). We quantified the geometry rather than footnoting it.

![Singular spectrum of the effect dictionary on the Th1/Th2 signature subspace, for all three culture conditions. Left: the spectrum decays fast — a few dominant axes carry most of the L2 energy (stable numerical rank ≈13). Right: cumulative variance shows a long low-variance tail, so the participation-ratio effective rank is ≈3,435. The dictionary is anisotropic (energy concentrated in a few directions), not globally collinear.]({{artifact:7d886da7-d22b-443c-8807-c07028aed15c}})

**What the spectrum shows:**

- **Two effective-rank measures, both reported honestly.** The **stable numerical rank**
  (σ_max² / Σσ²)⁻¹ is only **~10–15** (Rest 10.4, Stim8hr 15.5, Stim48hr 12.7): a handful of
  directions carry almost all the energy. The **participation-ratio** effective rank is
  **~3,400** (Rest 3388, Stim8hr 3430, Stim48hr 3488): a long tail of low-variance axes still
  carries independent variance. Both are true and the figure shows both.
- **Conditioning is moderate, not pathological.** The full condition number is **379–706**;
  restricted to the 99%-variance subspace it is **61–75**. The cone is workable, not rank-1.
- **This grounds the "greedy genes aren't the TFs" observation.** The best *single* generator
  toward Th1 is **LAT2** (cosine 0.226), not GATA3 or TBX21 — a direct consequence of the
  anisotropy (mean pairwise generator cosine is only 0.006, so the dictionary is anisotropic,
  not uniformly collinear). What looked like a wart is a measured property of the geometry.

Table: [`effective_rank_report.csv`]({{artifact:2fe4cd86-7525-4fa5-b39c-32a61e4ec3f2}}).

#### 8.3.5 One group-sparse cone unifies the reachable fraction and the recipe (new)

The reachable *fraction* came from a dense ~700-generator NNLS while the *recipe* came from a
separate greedy path — two objects, because the dictionary is collinear. A reviewer will ask
why. We built a **non-negative group-sparse cone**: collinear generators are grouped, and one
fit yields both the fraction and a group-sparse recipe. Because a group atom is the mean of
its members, a non-negative atom weight maps to non-negative weights on the original
generators, so the group-sparse solution is an *exactly valid* point in the original cone.

![Group-sparse cone versus the dense cone, all 12 atlas cells. Left: the one group-sparse fit reproduces the dense reachable cosine, clustering tightly along the y = x line (retention ≥95% in every cell). Right: the Th1/Rest recipe with the greedy picks and the canonical Th2 master TF GATA3 resolved into distinct explicit groups — GATA3 lands in its own group at rank 11 of 59, so the collinearity structure is made visible rather than worked around.]({{artifact:563c0968-d1f2-43d9-9c98-7799a14ce395}})

**What the unified object shows:**

- **≥95% of the reachable cosine is retained from one fit.** Across the 12 cells the
  group-sparse cosine is **95.2–98.9%** of the dense value (Th1/Rest **0.627 → 0.606**, 96.7%),
  using **46–91** active groups. The single object gives the fraction *and* the recipe; the
  two-object split is no longer necessary.
- **The collinearity becomes explicit and interpretable.** GATA3 (canonical Th2 master TF)
  falls in its own group at rank **11/59** for Th1/Rest, and **87% (13/15)** of the greedy
  recipe genes reappear among the active groups. The grouping *names* the redundancy the dense
  cone hides.

Table: [`group_sparse_atlas.csv`]({{artifact:0d0693ec-7f36-49c5-b557-ce12a813e9f6}}).

#### 8.3.6 The activation certificate is reproducible and partially corroborated — with an honest caveat about which genes it names (resolves L1's in-silico half)

The certificate — the genes that must be *activated* when a state is unreachable by knockdown
— is the project's most novel and least-tested output. The machine-precision residual proves
the convex program is optimal; it says nothing about biological correctness. The primary CD4
screen is CRISPRi (no activation arm), so a within-screen modality hold-out is impossible. We
did the two tests that *are* possible without a wet lab.

![Activation-certificate validation. Left: the certificate gene ranking is stable across a random gene-axis split — top genes hug the y = x line, and cross-half Spearman rho is 0.65 (half-vs-full 0.84), z ≈ 35 against the shuffle null of ~0. Right: literature and Open-Targets cross-reference of the top genes — the canonical Th1 master TFs (IFNG, TBX21, STAT4, STAT1, green stars) are in the certificate but rank deep, while the top-ranked genes are state-markers and negative regulators.]({{artifact:a9f7ac89-c55c-42c4-849c-8654b57994ca}})

**(A) The ranking is reproducible, not an artifact of which genes you fit on.** Rebuilding the
certificate from a random half of the signature genes and scoring on the other half, the top
genes keep their rank — LYAR (published rank 1 → mean half-fit rank **3.0**), IKZF3 (2 →
3.6), CRTAM (4 → 2.5), all with `frac_in_top30 = 1.0`. Cross-half Spearman **ρ = 0.65**,
half-vs-published **ρ = 0.84**, both **z ≈ 35** against a shuffle null of ~0.

**(B) Partial corroboration — and the important caveat.** Of the top 20 certificate genes,
the literature/genetics cross-reference classifies **4 SUPPORTS, 7 PLAUSIBLE, 9 UNKNOWN**, and
**8/20 have a strong (≥0.5) Open-Targets genetic link** to a Th1-driven autoimmune disease.
But the scientifically honest finding is *which* genes rank where: the **canonical Th1 master
regulators are present but ranked deep** (IFNG #35, TBX21 #144, STAT4 #364, STAT1 #473 — deep
precisely because they are *already* highly expressed, so their unmet upward demand is small),
while the **top-ranked certificate genes are state-markers and negative regulators** (IKZF3/
Aiolos, CBLB, LAG3 — genes whose *deletion* boosts T-cell function). The certificate names the
genes with the largest *unmet upward demand*, which is not the same as the causal Th1
activators — a caveat the paper must state plainly. The honest ceiling without a wet lab is
**"interpretable and reproducible, partially corroborated"** — not "validated."

![Certificate cross-reference detail: per-gene certificate score with literature direction-of-effect classification (left) and independent Open-Targets genetic corroboration against the top Th1-driven autoimmune disease (right).]({{artifact:a094179e-53bb-4cd7-9fa3-907bd5953788}})

Tables: [`heldout_modality_certificate.csv`]({{artifact:8baa1960-0b60-4bc1-8b76-7f9f2e2becaa}}),
[`certificate_literature_crossref.csv`]({{artifact:5faeb0a6-5e59-44ff-b3b3-f13a81faad7a}}).

#### 8.3.7 Directional genetics — association *count* hides conflicting *direction* (new)

"IRF1 has 17 immune-disease associations" is a count; a drug team needs the *direction*. The
oracle prescribes **knockdown** (lowering), so genetic support is *concordant* only if the
disease-risk allele raises expression (then lowering is protective). We cross-checked the
green-light nominations against eQTL direction (eQTL Catalogue / GTEx whole-blood / OneK1K
CD4-naive) and GWAS risk-allele direction.

**What the directional check shows** (10 genes): **2 CONCORDANT** — IFNGR1 (lead eQTL
β = −0.062) and CASP7 (β = +0.376) have risk/expression directions consistent with a
protective knockdown; **1 MIXED with a wrong-direction flag** — **IRF1**, whose 17-association
headline masks a discordant asthma variant (OneK1K CD4-naive lead β = −0.365), so a knockdown
could be *risk-increasing* for that trait even though other IRF1 variants are concordant; and
**7 INDETERMINATE** — an eGene but no disease-variant colocalisation, or not an eGene in the
queried blood/T-cell datasets. The takeaway for the paper: **association count is not
directional evidence**, and for at least one high-count nomination the direction is
partly wrong — exactly the check that stops a portfolio-enrichment signal from being
mistaken for a validated direction. (eQTL direction is necessary, not sufficient, for causal
direction — colocalisation / MR is the next step.)

#### 8.3.8 ChEMBL mechanism-of-action grounding — can a drug even reproduce the knockdown? (new)

A CRISPRi knockdown is genetic loss-of-function; a drug must *reproduce that loss*. There is
no LINCS/CMap L1000 connector available, so ChEMBL mechanism-of-action is the feasible proxy:
for each green-light nomination with a clinical-grade compound, does an approved or clinical
drug act in the **LOF direction** (inhibitor / antagonist / blocker) the knockdown implies?

![Genetics and pharmacology grounding of the green-light nominations. Left: directional-genetics concordance — the IRF1 discordant (wrong-direction) segment is flagged in red. Right: ChEMBL mechanism-of-action — LOF-consistent drug records per target, with the two wrong-direction targets (IFNGR1, RARA; only agonist drugs exist) in red.]({{artifact:d76ad46f-8ae5-4a89-b5b5-6bb19edc15f8}})

**What the MoA grounding shows** (10 genes): **4 APPROVED_LOF_MATCH** — JAK2 (ruxolitinib;
25/25 mechanisms are inhibitors), CD3D (muromonab-CD3), KEAP1 (dimethyl fumarate), VKORC1
(phenprocoumon); **2 CLINICAL_LOF_MATCH** — MAPK14 (ARRY-797, phase 3), ICOS (rozibafusp alfa,
phase 2); and — the decision-relevant negatives — **2 WRONG_DIRECTION** — IFNGR1 and RARA have
*only agonist/GOF* drugs (interferon gamma-1B, tamibarotene), so no approved compound can
reproduce the required knockdown; **CD5 AMBIGUOUS** (only a binding/cross-linking agent);
**CASP7 NO_CLINICAL_DRUG**. So the modality triage's "clinical drug" label is confirmed at the
*mechanism* level for six nominations and **contradicted for two** — the pharmacology, not the
druggability tier, is what tells you whether "knockdown works" can become "drug works."

Tables: [`directional_genetics_crosscheck.csv`]({{artifact:856c3ea9-4489-4806-905e-ed9fbc61dde7}}),
[`chembl_moa_grounding.csv`]({{artifact:3bad3dc6-9e59-421a-a07c-0b24f8da2503}}).

#### 8.3.9 Forward-predictor head-to-head — the only thing that beats the cone cheats physics (resolves L3 / pharma item 8)

The escape hatch a reviewer reaches for is "maybe the cone just can't predict." We ran the
head-to-head on a **common held-out-gene split** across all 12 cells: the cone (dense NNLS and
the k = 7 actionable recipe) against an unconstrained ridge **forward predictor** (the
GEARS/CPA-class task, done linearly — no GPU is configured, so scGPT/GEARS-class models are
remote-dependent and out of scope here) and a plain **DE-ranking** baseline.

![Forward-predictor head-to-head on the reachable task. Left: held-out-gene cosine for four methods across all 12 atlas cells. The unconstrained forward ridge (grey) wins on raw prediction, but the physically-realizable methods are the cone family (green); among realizable methods the cone leads and its k=7 recipe beats DE ranking in 11/12 cells. Right: the forward model's entire edge over the cone rests on ~50% negative weight — unrealizable "anti-knockdowns."]({{artifact:f22d4cac-7bba-412c-ac3e-4f6c060e9657}})

**What the head-to-head shows:**

- **The forward ridge wins on raw prediction — via impossible coefficients.** Its mean
  held-out cosine is **0.508** vs the cone's **0.355** (edge +0.153, winning 12/12), but that
  edge rests entirely on a **mean 50% negative weight mass** — coefficients that prescribe
  *activating* a gene, which CRISPRi cannot do. This is the §8.3.1 lesson at the whole-atlas
  scale: the unconstrained advantage is unrealizable.
- **Among physically realizable methods, the cone leads.** The cone's k = 7 recipe (mean
  **0.305**) beats DE ranking (mean **0.265**) in **11/12 cells**. The oracle's product is the
  feasibility verdict + a buildable recipe, not a prediction score — and on prediction, nothing
  that respects the non-negativity physics does better.

Table: [`forward_predictor_headtohead.csv`]({{artifact:1a3d84e4-89ca-446c-b64b-8c39d8ef800d}}).

---

### 8.4 Methodological agenda — making the method itself stronger

Ranked by value-per-effort. Items marked ✅ are delivered (above or in `RESULTS.md`); the
rest are open, with the cheapest in-silico wins first.

1. **✅ Non-negativity ablation (§8.3.1)** — L4 resolved.
2. **✅ Propagate generator uncertainty into the verdict (§8.3.3)** — the biggest *unlisted*
   gap, now closed. The dictionary is bootstrapped with the data-derived unit-noise model
   (`zscore = log_fc/lfcSE` exactly), B = 200 per cell; the reachable/infeasible verdict never
   flips (flip rate 0.0 in all 12 cells), with the aging-axis fine-grade caveat noted.
3. **✅ Held-out / stability certificate test (§8.3.6)** — L1's in-silico half resolved. A
   within-screen activation hold-out is impossible (CRISPRi has no activation arm), so the
   certificate ranking's *reproducibility* (cross-half ρ = 0.65, z ≈ 35) plus a
   literature/genetics cross-reference stand in. The remaining half — a prospective activation
   arm — is wet-lab (W1 in the limitations plan). The honest caveat (top genes are markers/
   negative regulators, not canonical drivers) is now documented.
4. **✅ Unify the dense-vs-sparse split with one group-sparse object (§8.3.5).** The non-negative
   group-sparse cone gives the fraction and the recipe from one fit, retaining ≥95% of the
   dense reachable cosine, and makes the collinearity explicit (GATA3 resolved into its own
   group).
5. **✅ Report the dictionary's effective rank / conditioning (§8.3.4).** Stable numerical rank
   ~10–15, participation-ratio rank ~3,400, condition number 379–706 — the near-degeneracy is
   now a stated geometric property that grounds why the greedy genes (LAT2, RARA) are not the
   canonical regulators.
6. **Move additivity from calibrated-borrowed to modeled (L2) — open.** The magnitude-
   saturation law (`M* = 13.9`, `§8.2.2`) is a real improvement but fit on K562
   CRISPRa. The bounded second-order interaction term (L2 option 3) converts additivity from an
   assumption into a modeled, bounded residual. This needs the in-domain double-knockdown data
   (W2 in the limitations plan), so it is gated on a wet-lab measurement.

**On the "novel algorithm" question specifically:** foreground the anisotropy-corrected null
(`§8.2.4`) as *the* methodological result. It is currently an §8.2.4 subsection; it is
the one piece of genuinely new derivation here (a closed form with a testable prediction,
validated at Pearson 0.995), and it is what lets the significance test be both cheap and
interpretable. Name it as a contribution, not a footnote.

---

### 8.5 Experimental agenda — what a pharma reviewer would require before acting

Putting on the target-selection hat: the framing (measured-not-inferred, genetic support,
modality triage) is the right register, but as evidence for *acting* on a nomination the
following are what a discovery team would demand. Ordered by how load-bearing they are; items
already tracked in the limitations plan are marked.

1. **A functional endpoint, not a transcriptional cosine (new).** The whole verdict is
   signature-matching, and the project's own guardrail concedes a signature match is not
   proof of functional rescue. A held-out cosine of 0.45 says nothing about whether the cell
   *behaves* Th1. Require the reachable recipe to move a **functional readout** — IFN-γ
   secretion, proliferation, persistence — not just the transcriptome. This transcriptome →
   phenotype gap is the first thing a team asks and is on no current list.
2. **Prospective test of *both* sides of one verdict (L1).** Reachable side: the minimal
   knockdown set actually shifts cells toward Th1. Infeasible side: the certificate's CRISPRa
   genes (LYAR, IKZF3, CRTAM, CBLB…) actually install the Th1 program. The infeasible-side
   experiment is the one that makes the paper *right* rather than merely hard to dismiss.
3. **◑ Genetic-to-pharmacological effect translation (new) — mechanism-level check delivered
   (§8.3.8), signature-level check still open.** A CRISPRi knockdown is chronic, near-complete
   genetic LOF; a drug gives partial, reversible, dose-dependent, off-target LOF. The ChEMBL
   mechanism-of-action grounding now confirms an approved/clinical **LOF-direction** compound
   exists for 6 of 10 green-light nominations and, decisively, flags **2 wrong-direction**
   targets (IFNGR1, RARA — agonist-only). What remains is the *signature-level* connectivity
   check — does a tool compound's transcriptome match the CRISPRi vector? — which needs a
   LINCS/CMap L1000 connector (not available here) or a dose-response.
4. **✅ Collateral-specificity screen (§8.3.2, new).** Delivered above for the four atlas axes;
   the natural extension is projecting each recipe onto off-target states a team actually
   fears (exhaustion, Treg→Teff, viability) once those signatures are in hand.
5. **In-domain epistasis (L2).** ~50–100 double knockdowns measured in CD4⁺ T cells, spanning
   the recipe sizes the oracle proposes. No team will trust an additivity ceiling borrowed
   from a different cell type *and* a different modality (K562 CRISPRa → CD4 CRISPRi).
6. **Robustness across donors / genotypes (L6).** The effect matrix is donor-collapsed, so
   leave-one-donor-out is impossible on the primary data. Per-donor effect vectors + LODO, or
   the external-validity claim stays "demonstrated once."
7. **✅ Directional genetics, not just association count (§8.3.7, new).** "IRF1 has 17
   immune-disease associations" is a portfolio-enrichment signal; the directional cross-check
   found 2 concordant nominations, 7 indeterminate, and — importantly — **IRF1 flagged
   MIXED/wrong-direction** (a discordant asthma eQTL), demonstrating that count is not
   direction. The remaining rigor step is formal colocalisation / Mendelian randomization to
   move "indeterminate" to a verdict.
8. **✅ Head-to-head on the overlapping task (§8.3.9, L3).** The recipe was benchmarked against a
   linear forward predictor and DE ranking on a common held-out split. The oracle does not win
   on raw prediction (the unconstrained ridge does — via 50% impossible negative weights) but
   leads among physically realizable methods, closing the "maybe it just cannot predict"
   escape. The GEARS/CPA/scGPT-class comparison is remote-GPU-dependent and remains open.

The genuinely *new* asks from the pharma lens — not in the original limitations plan — were the
**functional endpoint (1, wet-lab, still open)**, the **genetic→pharmacological translation
(3, mechanism-level now delivered)**, the **collateral-specificity profile (4, delivered for
the atlas axes)**, and **directional genetics (7, delivered)**. Of these, only the functional
endpoint now strictly requires a wet lab.

---

### 8.6 One documentation fix (done)

The positioning docs were already reconciled to the survey's true size — **91 prior methods**
(the `method_comparison_matrix.csv` has 92 rows because it includes this work's own row) — in
`NOVELTY.md` (and the former `IMPACT.md`, now merged into it), consistent with `RELATED_WORK.md` and
`limitations_and_reinforcement_plan.tex` (L3). One residual slipped through: a single prose
line in `NOVELTY.md §2e` still described the survey as "34 methods across three research
communities" while the quantitative counts two paragraphs below already read "91 prior methods
/ 14 measured." That stray "34" is now corrected to "91 prior methods," so the survey size is
stated consistently throughout. `RELATED_WORK.md` states the headline counts over the 91 prior
methods (this work's row excluded); nothing further remains to reconcile.

---

### 8.7 Bottom line — a novel formulation with one novel algorithm

It is a **novel formulation with one novel algorithm**, not a new optimizer — and it should be
sold as exactly that. The nine results shipped here strengthen it on every axis a critic
reaches for first, and — just as important — they surface the honest negatives that make the
work credible:

- **The cone earns its place.** The sign constraint is what makes an infeasibility certificate
  possible (§8.3.1); the unconstrained alternative's apparent edge is 50%-composed of impossible
  activations, at the single-cell scale (§8.3.1) and at the whole-atlas prediction scale (§8.3.9).
- **The verdict is robust where it matters.** Propagating generator measurement noise leaves
  the reachable/infeasible verdict unflipped in all 12 cells (§8.3.3); the near-degenerate
  geometry is now a *stated* property (§8.3.4); and the two-object dense/greedy split collapses
  into one group-sparse cone that keeps ≥95% of the reachable cosine (§8.3.5).
- **The novel output is reproducible but bounded.** The activation certificate's ranking is
  reproducible (z ≈ 35) and partially corroborated, with the honest caveat that its top genes
  are markers/negative-regulators rather than the canonical drivers (§8.3.6).
- **The pharmacology check does its job — by saying no.** Directional genetics flags a
  high-count nomination (IRF1) as wrong-direction (§8.3.7); ChEMBL grounding confirms a
  LOF-direction drug for six nominations and flags two (IFNGR1, RARA) as agonist-only (§8.3.8).
  A grounding step that only ever confirmed would not be worth running.

What genuinely remains needs a bench or a GPU, not more convex analysis: a **functional
endpoint** (does the recipe make the cell *behave* Th1, not just match its transcriptome —
W1), **in-domain double-knockdown epistasis** (W2), **per-donor leave-one-donor-out** (W3), a
**signature-level LINCS/CMap connectivity** check, and a **GEARS/CPA/scGPT-class** predictor
comparison. The convex core has now been pushed as far as in-silico work can take it.

---

## 9. Bottom line for the hackathon

The method does what the project claims: it converts a screen into a **falsifiable reachability
decision** with a **constructive certificate**, and on real data that decision is both nontrivial
and biologically coherent — *knockdown removes Th2, cannot install Th1, and here are the genes
that force the difference.* The honest framing is not "we polarize cells by knockdown" but
"we can tell you, with a certificate, exactly how far a knockdown screen can get and where it
provably stops." That is the decision-useful claim, and it is supported end-to-end.

### Artifacts produced

*Method module: `reachability.py`. Cached inputs: `analysis_cache/atlas_work/inputs.npz` (all
three per-condition effect matrices + four transition targets) and
`notebooks/cache/E_{Rest,Stim8hr,Stim48hr}.npz`.*

**Notebooks**
- `notebooks/01_exploratory_data_analysis.ipynb` — Tier-1 QC + data audit
- `notebooks/02_reachability_on_tier2.ipynb` — the Th2→Th1 headline pipeline (validated end-to-end)
- `notebooks/03_generalizability_and_impact.ipynb` — cross-dataset K562 CRISPRa demo
- `notebooks/04_experimental_design_toolkit.ipynb` — `design_experiment()`: verdict → calibrated recipe → modality triage (see §7)
- `notebooks/05_target_id_showcase.ipynb` — "From screen to shortlist" pharma target-ID walkthrough
- `notebooks/08_deg_weighted_evaluation.ipynb` — DEG-weighted evaluation & metric calibration (Needles-in-the-Haystack robustness test; §8.2.5)
- `notebooks/06_reinforcement_analyses.ipynb` — L1/L2/L4/L5 reinforcement battery (see `analysis_cache/nb_out/REINFORCEMENT_RESULTS.md`)
- `notebooks/07_cross_celltype_transfer.ipynb` — K562/RPE1 cross-cell-type transfer (see `analysis_cache/czi_data/CROSS_CELLTYPE_TRANSFER.md`)

**Headline figures** (`notebooks/figures/`)
- `fig1_reachability_spectrum.png` — reachable cosine vs k with null band
- `fig2_decomposition_certificate.png` — reachable/activation split + CRISPRa candidates
- `fig3_condition_comparison.png` — reachability across culture conditions
- `fig4_positive_control.png` — GATA3↓ / TBX21↓ placement in the alignment distribution
- `fig5_heldout_null.png` — held-out cosine vs shuffled-target null (z ≈ 24, 60 shuffles)

**Expansion figures** (`notebooks/figures/`)
- `fig_atlas_decomposition.png` — LOF/GOF/neither across the 12 atlas cells
- `fig_modality_triage.png` — genetics × druggability triage of the 102 knockdown nodes
- `nb03_fig1_norman_spectrum.png`, `nb03_fig2_norman_decomposition.png` — K562 CRISPRa demo

**Post-hackathon advance figures (§8.2)**
- `fig_positive_control_enrichment.png` — regulator-panel positive control (TF AUROC = 1.00)
- `fig_additivity_risk.png` — Norman saturation-law fit + per-recipe additivity risk (all < 0.1)
- `fig_anisotropy_null.png` — closed-form anisotropy-corrected null vs empirical shuffled null

**Tables** (`results/`)
- Headline (nb02): `table1_verdict.csv`, `table2_minimal_recipe.csv`,
  `table3_activation_certificate.csv`, `condition_comparison.csv`, `table5_null_summary.csv`
- Expansion: `atlas_reachability.csv` (12 cells), `modality_intervention_map.csv` (102 nodes),
  `genetics_crossverification.csv` (52 genetically-supported nominations),
  `tractability_grounding.csv`, `dataset_catalog.csv` (13 transfer datasets)
- K562 demo: `norman_table1_verdict.csv` … `norman_table5_null_summary.csv`
- Post-hackathon (§8.2): `positive_control_enrichment.csv`, `positive_control_stats.csv`,
  `atlas_additivity_risk.csv`, `norman_additivity_calibration.csv`,
  `anisotropy_null_validation.csv`
- Reinforcement (nb06, `analysis_cache/nb_out/`): `L4_constraint_ablation.csv` (NNLS-vs-unconstrained),
  `L5_reachable_cosine_ceiling.csv` (achievable ceiling), `L2_magnitude_capped_recipes.csv` +
  `L2_recipe_reliability_detail.csv` (per-recipe additivity reliability),
  `L1_certificate_test_scaffold.json` (synthetic dual-modality certificate test)
- Cross-cell-type (nb07, `analysis_cache/czi_data/`): `cross_celltype_effects.npz` (aligned K562/RPE1 basis),
  `per_perturbation_transfer.csv` (843 perturbations), `transfer_summary.json`
