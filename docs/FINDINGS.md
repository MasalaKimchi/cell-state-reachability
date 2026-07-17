# Updated Findings

This is the human-readable companion to [`results/findings.json`](../results/findings.json).
The JSON ledger controls displayed values; [`results/manifest.json`](../results/manifest.json)
controls artifact identity.

## Claim boundary

The current evidence is a retrospective, screen-relative test of transcriptional
direction under a non-negative linear-combination model. It does not demonstrate cell-state
conversion, functional rescue, intervention efficacy, biological necessity, or target
validation.

The dictionary contains donor-collapsed primary-CD4 CRISPRi differential-expression
z-score profiles measured in resting cells. The target is constructed from two external
Th1-vs-Th2 source contrasts. The perturbations were not applied to established polarized
Th2 cells.

## 1. Split-stable directional alignment

Across 12 fixed random-gene splits, the held-out cosine is **0.446 ± 0.010** (mean ± SD;
range 0.433–0.459). The older fixed split gives **0.448**. All 60 diagnostic target
shuffles fall below that fixed-split value, yielding a conservative plus-one empirical
**p = 1/61**.

These tests show stability to these particular random gene partitions and shuffles. Genes
are correlated, so neither test substitutes for module holdout, donor replication, or a
structured biological null.

Source: [`headline_heldout_split_stability.csv`](../results/evidence/headline_heldout_split_stability.csv).

## 2. The cone is interpretable, not the best predictor

On the registered split, held-out cosines are:

| Method | Held-out cosine |
|---|---:|
| Mean direction | -0.022 |
| PCA-1 | 0.006 |
| PCA-5 | 0.225 |
| Matched random-cone mean | 0.400 |
| Non-negative cone | **0.446** |
| Unconstrained least squares | 0.489 |
| PCA-20 | 0.504 |
| Ridge (λ=100) | 0.507 |

The cone is useful when non-negativity and measured-effect orientation are part of the
question. It is not supported as a generally superior predictor.

Source: [`baseline_comparison.csv`](../results/evidence/baseline_comparison.csv).

## 3. A circular calibration claim was retired

The previous “71% of modality ceiling” statement divided held-out cosine by a quantity
that was numerically the in-sample cosine. It was not an independent biological ceiling.
The retained descriptive replacement is a floor-to-synthetic-reachable-reference dynamic-range
fraction: **0.463** under uniform weighting and **0.615** under transition-specificity
weighting.

The reference comprises 40 reachable-by-construction noisy synthetic targets; it is not a
biological positive control. This is a sensitivity description, not a probability that the
target is attainable.

Sources: [`metric_calibration.csv`](../results/evidence/metric_calibration.csv) and
[`metric_calibration_provenance.json`](../results/evidence/metric_calibration_provenance.json).

## 4. Rest is the strongest observed context

Held-out cosine is **0.446** in Rest, **0.293** at Stim8hr, and **0.299** at Stim48hr.
Stim48hr is run-confounded because all four relevant libraries are on the same sequencing
run. Rest remains the case-study context; the result is not a polarization time-course
claim.

Sources: [`context_condition_comparison.csv`](../results/evidence/context_condition_comparison.csv)
and [`context_runbalance_caveat.json`](../results/evidence/context_runbalance_caveat.json).

## 5. Registered confounder challenges

- Removing 46 cell-cycle genes changes held-out cosine by about **-0.002**.
- The registered cone exceeds a norm-matched null (descriptive z ≈ **2.85**).
- Positive per-generator rescaling changes fitted cone geometry by at most numerical
  roundoff in the tested cases.

These checks address specific mechanisms. They do not establish absence of batch,
cell-state-composition, donor, or other biological confounding.

Sources: [`confounder_robustness.csv`](../results/evidence/confounder_robustness.csv) and
[`confounder_robustness_summary.json`](../results/evidence/confounder_robustness_summary.json).

## 6. The generator filter was already active

All **6,871** Rest generators in the frozen dictionary were source-flagged significant.
A stricter top-half effect-size filter improves held-out cosine from 0.4464 to about
**0.4580**. This addresses generator admission under the source test; it does not turn
the dictionary into noise-free causal effects.

Sources: [`generator_significance.csv`](../results/evidence/generator_significance.csv) and
[`generator_significance_summary.json`](../results/evidence/generator_significance_summary.json).

## 7. Directional ranking has a narrow curated benchmark

On a curated 13-gene polarization panel, the signed cone ranking has directional AUROC
**1.0**, versus **0.5** for an effect-magnitude baseline. The panel is small and curated;
this does not replace comparison with a frozen GRN-derived ranking or prospective assay.

Sources: [`ranking_validation.csv`](../results/evidence/ranking_validation.csv) and
[`ranking_validation_summary.json`](../results/evidence/ranking_validation_summary.json).

## 8. Combination semantics are unstable

Across 126 Norman K562 CRISPRa doubles, one retired legacy-v0 threshold label flips when
measured doubles replace additive approximations, and **100/126 staged proxy labels flip**.
Those labels used 200 target shuffles and a sequential non-negative `E` then `-E`
dominant-component heuristic; they are not mechanistic modality assignments. Median
cosine between measured and additive double effects is 0.712.

This is a useful failure demonstration: a coarse threshold can look stable while the
retired staged proxy changes. It does not validate primary-CD4 CRISPRi combinations.

Sources: [`additivity_verdict_flip.csv`](../results/evidence/additivity_verdict_flip.csv),
[`additivity_verdict_flip_summary.json`](../results/evidence/additivity_verdict_flip_summary.json),
and [`norman_table4_additivity.csv`](../results/evidence/norman_table4_additivity.csv).

## 9. Target construction matters

The two external source contrasts have cosine **0.689** across 11,616 shared genes, with
68.5% sign concordance. This supports a shared direction while also showing that exact
magnitude depends on source and filtering choices.

Sources: [`reviewer2_ota_hollbacher_meta.json`](../results/evidence/reviewer2_ota_hollbacher_meta.json)
and [`reviewer2_ota_hollbacher_split.csv`](../results/evidence/reviewer2_ota_hollbacher_split.csv).

Only **9,831/25,672** target-table genes intersect the screen, and **38/50** strongest
target DE genes survive that intersection. Source-only held-out cosines are **0.438** for
Ota and **0.342** for Höllbacher, versus 0.446 for the sign-concordant merge. The headline
therefore applies only to the measured, jointly signed subset; it is not evidence about the
unobserved target coordinates.

Source: [`reviewer2_deg_survival.csv`](../results/evidence/reviewer2_deg_survival.csv).

## Input provenance

The target contrasts trace to Ota et al., *Cell* 2021,
doi:10.1016/j.cell.2021.03.056 (NBDC E-GEAD-397), and Höllbacher et al.,
*ImmunoHorizons* 2020, doi:10.4049/immunohorizons.2000037 (GEO GSE149090). The
dictionary is Zhu et al. 2025, doi:10.64898/2025.12.23.696273 (SRP643211 / GSE314342).
The cross-system combination diagnostic is Norman et al., *Science* 2019,
doi:10.1126/science.aax4438 (GSE133344). Exact transformations are in
[`METHODS.md`](METHODS.md).

## What remains unknown

Claim-bearing biological evidence still requires:

- grouped/module holdouts and structured null calibration;
- true donor-held-out effects;
- measured in-domain combinations;
- paired measured CRISPRi and CRISPRa dictionaries;
- chromatin, protein/cytokine, viability, durability, and phenotype readouts;
- prospective validation from an established polarized starting state.

Until then, the appropriate output is a model-relative direction test and a reasoned next
measurement—not a reachability verdict or intervention recommendation.
