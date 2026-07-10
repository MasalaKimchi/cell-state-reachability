# Response to Reviewer 2 — dataset limitations

*Companion to the adversarial dataset appraisal, now consolidated into the "Reviewer 2" section of [`CAUSAL.md`](../CAUSAL.md). Each limitation raised
there is addressed here by one of four remedies — **reframe**, **new analysis**, **state**, or
**augment** — with the supporting number and a one-line residual-risk verdict. Every number is
recomputed from the local data; result tables are named inline and shipped in `results/`.*

**Reviewer posture retained.** I have tried to argue *against* the method wherever the data allowed
it, and to report the results that cut against us (the DEG-survival gap, the concordance-filter
inflation) as plainly as the ones that vindicate us. Two findings genuinely weaken specific claims and
are now stated as limitations rather than defended.

---

## At a glance

| # | Limitation | Remedy | Result / statement | Residual risk |
|---|---|---|---|---|
| **G1** | LOF-only; GOF budget "unreachable" | reframe | 25–31 % GOF is *modality-relative* to CRISPRi | Low once reframed |
| **G2** | Partial / heterogeneous knockdown | state + augment | 73 % compliance; verdict rescale-invariant (prior IV layer) | Moderate (recipe rank) |
| **G3** | Heteroscedastic generator noise | state | median 539 cells, p5=123; verdict robust, noise model homoscedastic | Moderate |
| **G4** | Off-target / neighbor contamination | **new analysis** | drop 808 flagged atoms → Δcos **−0.006**, top-6 recipe identical | **Low — resolved** |
| **D1** | 4 young donors, narrow demographics | state | structural; scope claim narrowed | Moderate (generalization) |
| **D2** | Aging axis has no aging span | **reframe + new analysis** | aging axes **fail** strict null (kSE\*≈0.001 vs 0.033) | **Low once demoted** |
| **D3** | Donor collapsed; no true LODO | augment (future) | per-donor h5ads accessible but 118–172 GB each | Moderate |
| **D4** | Stim48hr confounded with run | **new analysis** | headline is at run-balanced Rest; Stim48hr is *weaker* | **Low — resolved** |
| **T1** | Targets cross-study / cross-assay | state | structural cross-platform assumption | Moderate |
| **T2** | Two disagreeing source studies merged | **new analysis** | verdict stable (z 29–37); headline cos inflated by concordance filter | Low–moderate |
| **T3** | Silent truncation by gene intersection | **new analysis** | only **76 %** of top-50 target DEGs are measured | **Moderate — real caveat** |
| **P1** | Preprint-stage author-derived DE | state | disclose provenance | Low |
| **P2** | No local arrayed validation | state (already disclosed) | — | Low |

Consolidated evidence figure:

![Reviewer 2 trust checks]({{artifact:d5a86da6-010f-4e89-9987-4a1dcd7b0d75}})

*(a) Th2→Th1 verdict is stable whether the target is built from both source studies or either one alone.
(b) Removing all 808 off-target-flagged generators barely moves the reachable cosine. (c) Only ~76 % of
the strongest target DEGs are measured in the screen. (d) The aging axes fail the strict anisotropy null
while Th1/Th2 clear it.*

---

## New analyses (the substantive additions)

### G4 — Off-target specificity: **resolved**
`results/reviewer2_recipe_specificity.csv`. 808 of 6,871 Rest generators (11.8 %) carry an
`offtarget_flag`; three headline recipe genes (MUTYH, PIGS, VTI1A) are among them. Re-fitting the dense
NNLS and greedy recipe with **all flagged generators removed**:
- reachable cosine **0.6266 → 0.6203** (−0.006, ~1 % relative);
- the **top-6 recipe is identical** (LAT2, APPBP2, RARA, ICOS, SNAP23, SNX4); 11/15 overlap;
- the flagged genes drop out and are replaced by unflagged genes deeper in the spectrum.

**Verdict:** the headline verdict and leading nominations are *not* off-target artifacts. This concern is
retired. *(Remedy: new analysis.)*

### D4 — Stim48hr run confound: **resolved for the headline**
`results/reviewer2_stim48_confound.csv`. All four Stim48hr libraries are on sequencing run `CD4i_R2`
(Rest and Stim8hr are balanced 2/2). But the headline Th1/Th2 verdict is reported at **Rest**
(reachable 0.63/0.64, held-out z 45/21), which is run-balanced, and the confounded **Stim48hr condition
gives *weaker* cosines (0.53)** — so it cannot be inflating the headline.

**Verdict:** the headline is safe; any *Stim48hr-specific* claim is flagged as run-confounded. *(Remedy:
new analysis + state.)*

### T2 — Ota vs Höllbacher: **verdict stable, magnitude partly inflated by the concordance filter**
`results/reviewer2_ota_hollbacher_split.csv`. The two source studies share 11,616 genes at a
between-study cosine of only **0.689** — they genuinely disagree on ~31 % of genes. The pipeline already
mitigates this by keeping only the **68.5 %** of genes whose sign agrees across both studies. Recomputing
the verdict three ways on E_Rest:

| target | reachable cos | held-out cos | held-out z |
|---|---|---|---|
| merged (headline, sign-concordant) | 0.627 | 0.446 | 28.7 |
| Ota 2021 only | 0.560 | 0.438 | 37.2 |
| Höllbacher 2021 only | 0.511 | 0.342 | 30.4 |

**Verdict:** the *reachability verdict does not depend on source choice* — all three are strongly
reachable and far above null. But the headline reachable cosine (0.63) sits **above either study alone**,
because the sign-concordance filter keeps only mutually-agreed genes, which raises the cosine. This is
worth one honest sentence: the exact headline number benefits from a filtering step; the qualitative
verdict does not. *(Remedy: new analysis; disclose the filter.)*

### T3 — Silent truncation: **a real caveat**
`results/reviewer2_deg_survival.csv`. Of the **top-50 strongest DEGs** in the Th2/Th1 target signature,
only **38 (76 %) are measured** in the screen's 10,282-gene axis; the rest are silently dropped before
the cone sees them. The fraction is stable at 74–80 % from top-50 through top-1000 (aging signature is
identical at 76 %).

**Verdict:** the reachable cosine is computed on a target **missing roughly a quarter of its most-defining
genes** — a data-*coverage* limitation, not a method flaw, but one that must be stated. It caps how
completely any verdict on this screen can speak to the full biological target. *(Remedy: state, with the
number.)*

---

## Reframes

### G1 — the LOF/GOF split is modality-relative
The signed decomposition gives LOF/GOF/neither = 39/25/35 (signed) for Th2→Th1: loss-of-function alone
covers well under half the shift. This is currently read as a biological fact about the target. It is
more precisely a **fact about the assay** — CRISPRi can only remove function, so the size of the
"gain-of-function-locked" budget is set by the modality of the one screen used. A CRISPRa arm would
convert an unknown part of that budget into reachable territory.

**Reframe:** present 25–31 % GOF not as "unreachable" but as "unreachable *by knockdown alone*, and here
is the falsifiable CRISPRa prediction that follows." This turns a limitation into the headline next
experiment. *(The certificate already emits the specific must-activate genes to test.)*

### D2 — the aging axis as a negative demonstration
`a1_sensitivity_radius.csv`. At Stim48hr the aging axes **fail the strict anisotropy null**: younger
baseline cosine 0.597 vs null-p99 0.616 (z ≈ 1.90, `boot_clears_null = False`), older 0.566 (z ≈ 1.20),
both with a coordinated-bias radius **kSE\* ≈ 0.001** — essentially zero margin — versus kSE\* ≈ 0.033 for
Th1. This is exactly what should happen: the generators come from donors aged 22–34, so the dictionary
*cannot* legitimately reach an aging signature, and the trust layer says so.

**Reframe:** demote the aging axis from a first-class result to an explicit **negative demonstration** —
"here is the oracle correctly refusing a target the data cannot support." As-is it invites a reader to
distrust the Th1/Th2 axes that *are* solid; reframed, it becomes evidence the calibration works. *(Remedy:
reframe + the null-failure number.)*

---

## Statements (structural limitations to disclose)

- **D1 — donor panel.** Four donors, ages 22–34, 3 F / 1 M, all blood type O+. Every generator and verdict
  is estimated from this panel. Claim scope: *young, healthy, predominantly female primary human CD4⁺ T
  cells* — no claim beyond it.
- **T1 — cross-assay targets.** Generators are single-cell CRISPRi pseudobulk z-scores; targets are
  bulk/sorted-population DE from external studies (Ota, Höllbacher, Yaza). They are reconciled by
  gene-name intersection + cosine. This cross-platform commensurability is an assumption, now stated.
- **G2 / G3 — imperfect, heteroscedastic instruments.** 73 % of guides reach significant knockdown;
  among those, median residual expression 0.115 (p90 = 0.49). Per-perturbation support is
  right-skewed (median 539 cells, p5 = 123, min 17). The existing IV/compliance layer shows the verdict
  is invariant to per-generator positive rescaling and the EIV bootstrap shows it is robust to unit noise;
  the honest residual is that the noise model is homoscedastic and the recipe *ranking* is not
  rescale-invariant.
- **P1 — provenance.** The dictionary is the authors' preprint-stage DE output (Zhu 2025, bioRxiv), not
  recomputed from raw counts; its pseudobulk/shrinkage choices propagate into every generator.

---

## Augmentation — what CZI-accessible data can and cannot fix
`results/czi_augmentation_map.csv`. I assessed every accessible CZI/GEO dataset against the limitations.
The honest bottom line: **no accessible dataset fixes the core donor-diversity limitation (D1/D2).**

- **Cell-line screens** (Replogle K562/RPE1, cz-benchmarks S3) — already used for cross-system transfer
  (notebook 07). They test generalization but are cell lines, LOF-only, and add no donor diversity.
- **Per-donor T-cell h5ads** (marson2025_data S3) — would enable a true leave-one-donor-out (D3) and
  per-atom noise from raw cells, but each file is **118–172 GB** (infeasible on an 18 GB-RAM laptop this
  session) and still covers only the same 4 young donors — so even these do not add an aging span or new
  donors.
- **Joung 2023 TF Atlas** (GSE216481, GEO) — the intended CRISPRa/GOF test for G1, but in hESC, not a
  matched activation arm of this screen.
- **CELLxGENE (connector, used this session)** — provides no perturbation vectors, so it cannot fix any
  dictionary limitation, but it does give an independent **biological-coherence** check
  (`results/reviewer2_marker_coverage.csv`): **93–95 % of canonical Th1/Th2 marker genes are measured in
  the screen** and 71–84 % appear in the target signature, including GATA3 (the Th2 master TF) in the Th2
  markers. This confirms the target axes are biologically real, not artifacts of the two stitched studies.

---

## Ready-to-paste Limitations paragraph

> **Dataset scope and its consequences.** Our dictionary is a single genome-scale CRISPRi screen in
> primary human CD4⁺ T cells from four young donors (ages 22–34), and this bounds every claim in three
> ways. First, because CRISPRi only removes function, the loss-of-function/gain-of-function split we
> report (≈39 %/25 % for Th2→Th1) is relative to the assay: the gain-of-function budget we label
> unreachable is unreachable *by knockdown*, and its size would change under a CRISPRa arm — which the
> infeasibility certificate turns into a specific, falsifiable prediction. Second, the donor panel is
> demographically narrow, so we make no claim beyond young, healthy, predominantly female CD4⁺ T cells;
> in particular the aging-axis targets, whose signature comes from a separate older cohort, correctly
> *fail* our strict anisotropy null (coordinated-bias radius ≈0.001 vs ≈0.03 for the polarization axes),
> and we present them as a negative demonstration that the calibration refuses targets the data cannot
> support. Third, the target vectors are external, cross-assay differential-expression signatures
> reconciled to the screen by gene intersection: only ~76 % of a target's strongest DE genes are measured
> in the screen, so verdicts are computed on partially-observed targets. We verified that the polarization
> verdict is stable to the choice of source study (reachable cosine 0.51–0.63 across Ota-only,
> Höllbacher-only, and their sign-concordant merge; held-out z = 29–37) and to the removal of all
> off-target-flagged generators (Δ reachable cosine −0.006, identical leading recipe), and that the
> headline is reported at the run-balanced Rest condition rather than the run-confounded Stim48hr
> condition. Extending the dictionary to additional donors, an activation arm, and in-domain
> combinatorial perturbations is the primary direction for future work; we note that while per-donor
> data for this screen is publicly accessible, no currently accessible perturbation dataset adds the donor
> diversity or activation modality these limitations call for.

---

## Result files
- `results/reviewer2_ota_hollbacher_split.csv` — T2 verdict on each source study
- `results/reviewer2_recipe_specificity.csv` — G4 off-target re-fit
- `results/reviewer2_deg_survival.csv` — T3 top-DEG survival fractions
- `results/reviewer2_stim48_confound.csv` — D4 run-confound cross-check
- `results/reviewer2_marker_coverage.csv` — CELLxGENE marker coherence check
- `results/czi_augmentation_map.csv` — dataset-by-limitation augmentation map
- `reviewer2_trust_checks.png` — consolidated evidence figure
