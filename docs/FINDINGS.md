# Updated Findings

This narrative is controlled by [`results/findings.json`](../results/findings.json); the
artifact identities are controlled by [`results/manifest.json`](../results/manifest.json).

## Claim boundary

The project tests directional support under a non-negative linear-combination model. The
primary dictionary contains donor-collapsed primary-CD4 CRISPRi differential-expression
profiles from the source study's post-expansion `Rest` condition. The target is a
source-study-reused, selectively constructed cross-sectional population contrast oriented
toward the reported Th1 centroid. It is not independent, an observed trajectory, or a
measurement from established polarized Th2 cells.

## 1. Source-bound directional alignment is variable across gene splits

With the gene universe ordered lexicographically, split IDs hashed, and NumPy
`default_rng(seed)` frozen for seeds 0–11, held-out cosine is **0.444 ± 0.018** (mean ±
SD; range **0.417–0.473**) over 12 half-gene splits. This is split variability, not donor
uncertainty. Correlated genes make these splits neither independent replicates nor a
substitute for module holdout.

The earlier **0.446 ± 0.010** table is retired: its script consumed an unhashed, deleted
`analysis_cache/atlas_work/inputs.npz`, and that object's gene ordering was not preserved.
The separately archived fixed split (**0.448154**) agrees with the new source-bound seed-0
split within `3e-10`, confirming source alignment but not reconstructing the old multisplit
protocol. Its 60 target shuffles were not regenerated and **p = 1/61 is no longer a current
finding**.

Source: [`source_reconstruction.json`](../results/source_reconstruction.json).

## 2. Target construction materially narrows the estimand

The target table has 25,672 genes in the source union, 11,616 shared by Ota and
Höllbacher, 7,960 with concordant effect signs, and 6,188 remaining after screen
intersection. Across the 8,950 shared, screen-measured genes used for leakage-safe source
transfer, between-source cosine is **0.791** for log fold change and **0.698** for z-score.
These sources differ in cohort, assay, and preprocessing; agreement is not independence
from the source study.

## 3. Cross-source transfer is directional, not magnitude accurate

The cone is fit to one target source and evaluated against the other without selecting
genes by held-out-source sign. On log fold change, cone cosine exceeds the better of the
mean-response ray and best single atom in every one of six correlated random-gene splits:

| Direction | Mean cone cosine | Mean cosine gain | Mean cone nRMSE | Best-single nRMSE |
|---|---:|---:|---:|---:|
| Ota → Höllbacher | 0.255 | +0.087 | 1.185 | 1.003 |
| Höllbacher → Ota | 0.290 | +0.090 | 1.020 | 0.981 |

The directional gain is consistent across these splits, but magnitude error is worse than
the best single and sometimes worse than zero prediction. Baselines remain incomplete
(no nested PCA, ridge, unconstrained least squares, or matched random cone), and selecting
the better baseline on the test coordinates is disclosed. This is a construction-sensitivity
diagnostic, not predictive superiority, donor generalization, or biological validation.
Z-score results are sensitivity-only.

## 4. Independent Arce transfer is modest and context dependent

The compact Arce benchmark compares the frozen negative Zhu IL2RA log-fold-change score
with an independent CRISPR-KO CD25/IL2RA screen. Selection uses only four-guide coverage
in every Arce context and source-side Rest admission, leaving 480 targets. No Arce outcome
is used for selection or choosing the sign; orientation is fixed from assay semantics.

| Context | Signed Spearman | Signed Kendall | Direction agreement | Magnitude top-25 overlap |
|---|---:|---:|---:|---:|
| Resting Teff | 0.148 | 0.102 | 0.538 | 8 |
| Stimulated Teff | 0.084 | 0.057 | 0.521 | 9 |
| Resting Treg | 0.088 | 0.061 | 0.535 | 5 |

Resting-Teff rank association is strongest; the other contexts are weaker. Top-k overlap
ranks absolute magnitudes, while global Spearman/Kendall retain the prespecified sign.
Permutations exchange target labels, not donors. All 18 context-by-metric permutation
p-values are unadjusted exploratory diagnostics across correlated tests and support no
FWER- or FDR-controlled inference. The result supports only modest cross-study,
cross-modality ranking alignment, not state reachability, donor generalization, or causal
treatment validation.

Sources: [`arce_external_validation_meta.json`](../results/evidence/arce_external_validation_meta.json)
and [`arce_il2ra_context_predictions.csv`](../results/evidence/arce_il2ra_context_predictions.csv).

## 5. The systemic harness exposes nonspecific success modes

The data-free harness independently checks the NNLS solver and fails closed on axis and
provenance corruption. Its structured-specificity scenario produces high raw cosine from
a nuisance common response, shows random-gene optimism relative to module holdout, and
shows large sign-selection inflation while retaining power for a sparse true alternative.
The fast maxT check has 24/500 false family rejections; its exact one-sided 95% upper bound
is 0.067 under the declared 0.075 gate. This certifies software and synthetic statistical
contracts only.

## What remains unknown

Claim-bearing biological evidence still requires donor- and guide-held-out effects,
module/pathway holdouts, stronger nested baselines, calibrated structured nulls, measured
in-domain combinations, paired CRISPRi/CRISPRa, chromatin and protein/function endpoints,
durability and fitness, and prospective testing from an established polarized starting
state. Until then the appropriate output is a model-relative direction test and the next
measurement—not a reachability verdict or intervention recommendation.
