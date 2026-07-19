# Findings

This narrative is controlled by [`results/findings.json`](../results/findings.json); the
artifact identities are controlled by [`results/manifest.json`](../results/manifest.json).

**How to read this document.** The project's contribution is fail-closed,
model-relative point-estimate cone geometry plus an adversarial harness that exposes known
failure modes. Sections are tagged by tier: **[MAIN]** is the numerical method and harness;
**[SUPPORTING]** is the Zhu Th2→Th1 proving ground; **[STRESS]** maps empirical limits; and
**[EXTENSION]** applies the same geometry to retrospective effect-dictionary coverage.
None of these tiers turns numerical representability into biological reachability.

## Claim boundary

The project tests directional support under a non-negative linear-combination model. The
primary dictionary contains donor-collapsed primary-CD4 CRISPRi differential-expression
profiles from the source study's post-expansion `Rest` condition. The target is a
source-study-reused, selectively constructed cross-sectional population contrast oriented
toward the reported Th1 centroid. It is not independent, an observed trajectory, or a
measurement from established polarized Th2 cells.

## 1. [SUPPORTING] Source-bound directional alignment is variable across gene splits

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

## 2. [SUPPORTING] Target construction materially narrows the estimand

The target table has 25,672 genes in the source union, 11,616 shared by Ota and
Höllbacher, 7,960 with concordant effect signs, and 6,188 remaining after screen
intersection. Across the 8,950 shared, screen-measured genes used for leakage-safe source
transfer, between-source cosine is **0.791** for log fold change and **0.698** for z-score.
These sources differ in cohort, assay, and preprocessing; agreement is not independence
from the source study.

## 3. [SUPPORTING] Cross-source transfer is directional, not magnitude accurate

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

## 4. [STRESS] Donor-pair transfer retains weak direction but fails magnitude

The released donor H5MU contains six two-donor modalities. All six share 1,584 complete
`Rest` perturbation atoms, and 8,949 genes are common to both target sources and the
donor object. For each of three complementary 2-vs-2 donor partitions, the runner fits
NNLS weights on one donor-pair dictionary, one target source, and one random half of the
genes. It applies those weights unchanged to the complementary donor pair, opposite
target source, and held-out genes. The common-ray and best-single identity/scalar are
also selected only from the training side and frozen before testing.

The two mixed-run partitions contribute 24 run-balanced correlated challenges across
both donor directions, both target-source directions, and three fixed gene splits. Median
cone cosine is **0.031**, a **+0.032** gain over the training-selected best single; 18/24
challenges have positive cosine gain. Magnitude fails: median normalized RMSE is **1.153**
for the cone versus **1.018** for the best single, with improvement in only 1/24
challenges and none over the zero/common-response baselines. The 12 fully run-confounded
donors-1–2 versus donors-3–4 challenges are reported separately and not headlined.

This is a published-pipeline sensitivity, not leakage-free donor generalization. The
modalities are two-donor summaries, released presence is conditioned on DE eligibility,
gene splits are correlated, and four fixed donors cannot support donor-population
inference at 0.05. The result therefore narrows the utility claim: some directional
structure transfers, but predictive magnitude does not.

Source: [`donor_pair_transfer.json`](../results/donor_pair_transfer.json).

## 5. [STRESS] Released guide-rank transfer does not demonstrate reciprocal robustness

The official VCP v1.0 dataset card (accessed 2026-07-19) defines `guide_1` and `guide_2`
as the first and second sgRNA ID after alphanumeric sorting within target-condition, and
public pseudobulk/guide-library artifacts carry IDs. The exact 29.4 GB ranked H5MU does not
itself embed sgRNA ID or sequence, and this repository has not reconstructed and
hash-cross-verified the exact rank-to-sgRNA crosswalk. `guide_2` is an exact subset of
`guide_1` in this object.

The analysis retains all 8,323 common category-labeled `Rest` atoms, reports 2,752
`guide_1`-only category-labeled `Rest` atoms as excluded, and withholds another 35 nominal
`guide_1` `Rest` keys whose categorical metadata are missing. It uses 8,950 target genes
selected by identifier presence only, fits 12 models (two training ranks × two target
sources × three gene splits), freezes each coefficient vector and training-selected
baselines, and scores both ranks against the same and opposite target sources, yielding 24
correlated rows.

The same-source result fails reciprocal-position robustness:

| Frozen evaluation | Median cone cosine | Median paired cone gain over training best single | Positive cosine gain | Median cone / best-single nRMSE |
|---|---:|---:|---:|---:|
| Within training rank, held genes | 0.251227 | +0.132000 | 12/12 | 1.035181 / 0.996362 |
| Reciprocal held rank, same target source | −0.019197 | −0.020345 | 3/12 | 1.386696 / 1.017276 |
| Reciprocal held rank + opposite target source | −0.034383 | −0.041518 | 2/12 | 1.347783 / 1.015038 |

The same-source median paired change from within-rank to reciprocal-rank cone
cosine is **−0.291089**. Normalized-RMSE gain over the frozen best single is positive in
**0/12** reciprocal rows. The within-rank-versus-reciprocal-rank cone predictions have median cosine
**0.510527**, so partial prediction alignment does not preserve target alignment. One
additional reverse-atom-order fit changes the canonical solution by at most **1.76×10⁻¹⁵**
relative L2; this is a narrow representation check, not proof of coefficient
identifiability.

The source object is full-file-hash-bound to
`964eeafb3356a7322a1d5b1121802c6a1433456f3591e2d5797817df3bf9c2f6`. Its 114
`guide_1` categorical-missing rows have suffixes `Rest` 35 / `Stim8hr` 18 / `Stim48hr`
61; all 38 `guide_2` missing rows have suffix `Stim48hr`. They are recorded without
imputation, and the selected category-labeled common-`Rest` target mapping is complete.
More importantly, the pinned author pipeline first selected
`keep_min_cells & keep_effective_guides & keep_total_counts`, required at least three
passing replicates per guide-condition and at least five cells per guide in each
condition/sample, constructed `cond_targets` from the unreleased `for_DE_by_guide.csv`,
kept exactly two testable guides per target-condition, applied final `keep_test_genes`, and
fit `~ log10_n_cells + target` without a donor term. This fixed descriptive reciprocal-rank
benchmark fails to demonstrate robustness across the released rank summaries. It is not
named-sgRNA replication or leakage-safe physical-guide generalization.

Source: [`guide_pair_transfer.json`](../results/guide_pair_transfer.json).

## 6. [STRESS] Independent Arce transfer is modest and context dependent

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

## 7. [STRESS] Arce donor/guide strata expose both robustness and heterogeneity

The same hash-bound archive contains 100,087 S14 singlet cells across two donors, two
guides for each member of the authors' preselected 28-regulator panel, nine Non-Targeting guides, and four Teff/Treg ×
resting/stimulated contexts. All 520 donor×guide×context strata are present. The runner
first takes each guide's median supplied `activation.score`, then subtracts the median of
the nine Non-Targeting guide medians within donor and context; cells are never counted as
biological replicates.

All 116 published S8 target×context means and medians reproduce from S14 with maximum
absolute error **5.7×10⁻¹³**. Across 28 targets, donor A-versus-B rank concordance is
**0.73–0.93** across contexts, but only **50–64%** of targets have the same nonzero effect
sign across both guides and both donors. One retained stratum has 8 cells; every other
stratum has at least 20. These are descriptive robustness diagnostics. S8 is derived from
S14, its pooled-cell p-values are unused, and two donors cannot estimate population
uncertainty. The score's gene set, formula, normalization, and independence are not frozen,
so it is not interpreted as CD25 protein, functional activation, Th1/Th2 identity, or
state conversion. Panel selection used prior screen/state-specific evidence, so these
metrics are conditional within-panel concordance, not genome-wide generality or independent
validation. The endpoint-mismatched Zhu IL2RA-to-`activation.score` correlation is
intentionally not reported.

Source: [`arce_activation_guide_effects.csv`](../results/evidence/arce_activation_guide_effects.csv).

## 8. [SUPPORTING] Source-selected arrayed follow-up is transcriptomically specific and cytokine-consistent

The source authors arrayed nine screen-selected perturbations and released bulk RNA-seq
plus IL-10/IL-21 flow percentages for six additional donor labels (Donor5–Donor10). The
runner binds both tables to author commit `848d62f`, verifies every input byte, uses the
`Stim8hr` screen profiles, intersects 8,976 Ensembl genes, and masks all nine panel target
genes from every profile before transcriptome scoring (8,967 coordinates per profile).

Every arrayed bulk profile retrieves its matching screen perturbation first among all nine
(top-1 **9/9**, MRR **1.0**). Median raw cosine is **0.549**, versus a median **0.283**
gain over the common source-panel response. After subtracting each assay's across-panel
gene mean—a Systema-style specificity sensitivity—median cosine is **0.580** and retrieval
remains 9/9. Median normalized RMSE is **1.012** raw and **1.052** centered, so magnitude
calibration is not supported despite strong identity and directional structure.

Flow effects are `log2(percent positive / donor mean NTC percent positive)` within donor;
targets then give every observed donor equal weight. Donor coverage is uneven (three to
six). Across all nine targets, Spearman correlation between donor-median flow effect and
screen/bulk RNA cytokine log fold change is **0.817/0.717** for IL-10 and **0.833/0.850**
for IL-21. Synchronized exhaustive target-label diagnostics over these four correlations
are conditional on the selected panel; maximum-statistic tail fractions are 0.020, 0.062,
0.015, and 0.011, respectively. Unequal donor count/composition means target
exchangeability is not established. These fractions are not inferential p-values or
multiplicity-adjusted inference, do not undo upstream selection, and are not
donor-population inference.

This supports same-study cross-platform replication and direct cytokine consistency for
the measured panel. It does not establish held-out target discovery, guide robustness,
population generalization, durable state conversion, function, fitness, chromatin
remodeling, or intervention efficacy.

Sources: [`zhu_arrayed_validation_meta.json`](../results/evidence/zhu_arrayed_validation_meta.json),
[`zhu_arrayed_profile_metrics.csv`](../results/evidence/zhu_arrayed_profile_metrics.csv), and
[`zhu_arrayed_flow_effects.csv`](../results/evidence/zhu_arrayed_flow_effects.csv).

## 9. [MAIN] The systemic harness exposes nonspecific success modes

The data-free harness independently checks the NNLS solver and fails closed on axis and
provenance corruption. Its structured-specificity scenario produces high raw cosine from
a nuisance common response, shows random-gene optimism relative to module holdout, and
shows large sign-selection inflation while retaining power for a sparse true alternative.
The fast maxT check has 24/500 false family rejections; its exact one-sided 95% upper bound
is 0.067 under the declared 0.075 gate. This certifies software and synthetic statistical
contracts only.

## 10. [EXTENSION] Retrospective effect-dictionary coverage is distinct from prospective design

`library_coverage.py` lifts the projection to a catalog while keeping strict cone membership
separate from an application-chosen cosine bar. It also measures leave-one-out atom
redundancy and scores supplied candidate effect vectors against certified dual-separator
normals. The optional realized comparator appends each already-measured candidate and
recomputes mean catalog cosine; it is not an experimental ground truth.

The v3 audit assigns library/catalog/candidate roles before feature selection. Zhu admits
all 6,871 source-eligible effects rather than selecting by benchmark magnitude. Norman
groups canonical genes before role assignment, freezes one measured representative per
gene, and proves zero library/candidate gene overlap. Twelve deterministic partitions add
best-single, common-response, and signed-span comparators. They quantify algorithmic
sensitivity over correlated measured rows—not biological uncertainty or generalization.

The three portable caches now have hash-gated deterministic builders. This closes the
computational reconstruction gap, but not durable availability: the upstream
Norman/Replogle URLs are mutable and historical object versions cannot be retrieved
anonymously. A release archive remains necessary.

<!-- BEGIN VALIDATED LIBRARY COVERAGE SUMMARY -->
| Audit | Cone mean cosine | Strict membership | Cosine ≥0.5 | Simple comparator | Certificate vs realized gain |
|---|---:|---:|---:|---|---|
| Zhu CD4 CRISPRi | 0.436 | 0/150 | 32/150 | Best atom: 0.241; 0/150 | ρ=0.861; top-1 differs |
| Norman K562 CRISPRa | 0.958 | 0/131 | 131/131 | Constituent sum: 0.914; two-atom cone: 0.924 | ρ=0.921; top-1 differs |
| Replogle K562 CRISPRi | 0.776 | 0/150 | 143/150 | Best atom: 0.674; 116/150 | ρ=0.873; top-1 differs |

| 12-split sensitivity | Cone cosine median [range] | Cosine ≥0.5 median [range] | Strict-positive splits | Cone − best atom | Signed span − cone |
|---|---:|---:|---:|---:|---:|
| Zhu CD4 CRISPRi | 0.427 [0.409, 0.443] | 0.200 [0.150, 0.283] | 0/12 | +0.188 | +0.084 |
| Norman K562 CRISPRa | 0.865 [0.807, 0.878] | 0.983 [0.900, 1.000] | 0/12 | +0.106 | +0.033 |
| Replogle K562 CRISPRi | 0.774 [0.745, 0.807] | 0.950 [0.900, 0.967] | 0/12 | +0.095 | +0.063 |

Strict membership is absent from every reference catalog and every sensitivity split, so a
soft directional bar cannot be renamed reachability. The cone beats best-single and
common-response rays in all 36 splits, while the signed span is uniformly better: useful
alignment exists, but non-negativity is a real capacity constraint. Rank correlations are
strong, yet certificate and realized top candidates disagree in 3/3 audits. Norman's
0.958 reference uses all 152 single rows; its 0.865 sensitivity median uses 40 measured
canonical-gene representatives and is the appropriate partition-robustness result. Even in
the reference design, constituent-only baselines already clear the soft bar for 131/131
doubles, so the full cone adds alignment rather than establishing a biological manifold.
<!-- END VALIDATED LIBRARY COVERAGE SUMMARY -->

Source: [`library_coverage_crossdataset.json`](../results/library_coverage_crossdataset.json).

## 11. [STRESS] The Goudy triple is a negative, confounded cross-experiment result

The registered GSE306915 run has execution status **`PASS`**, geometric status
**`FAILS_DECLARED_GEOMETRIC_MODEL`**, and interpretation status
**`INCONCLUSIVE_CROSS_EXPERIMENT_CONFOUNDING_LOW_RELIABILITY`**. At the primary
control-only CPM ≥1 filter, every donor-specific component cone remains strictly outside
the declared cone:

| Model | Selection/scoring contract | Median cosine | Median nRMSE | Strict inside |
|---|---|---:|---:|---:|
| Component cone | Same-donor in-sample oracle | 0.0949 | 0.9952 | 0/4 |
| Equal sum | Three unit coefficients | 0.0646 | 2.3084 | — |
| LODO component cone | Three training donors → held donor | 0.0881 | 0.9982 | — |
| LODO best single | Identity/scale selected on training donors | 0.0957 | 0.9968 | — |
| On-target-frozen score | Three target coordinates → disjoint transcriptome | 0.0690 | 1.8478 | — |

Reliability is weak at the effect level. Median pairwise-donor cosine is **0.0480** for the
triple and **0.0037 / 0.0809 / 0.0419** for the FAS / RC3H1 / SUV39H1 singles. Although
the two AAVS1 expression profiles have median cosine 0.9942, their difference norm is
**0.890–1.634×** the measured triple-effect norm. Across the prespecified control-only CPM
thresholds, component-cone median cosine ranges from **0.0696 to 0.1999**; no threshold
rescues the model.

The near-exact three-coordinate on-target fit is a saturated 3×3 calibration diagnostic;
only its frozen, disjoint-transcriptome score appears above. The additive residual is
defined, but the experiment/control/guide-burden confounding prevents it from identifying
a statistical interaction. A module-level error is unavailable because no module set was
preregistered. Together with unresolved triple guide identity and incomplete author-key
coverage for donor-3/4 multiplex roles, these results support only a descriptive negative
cross-experiment stress conclusion.

Source: [`goudy_combination_validation.json`](../results/goudy_combination_validation.json).

## 12. [STRESS] Schmidt top-effect transfer is conditional; genome-wide agreement is limited

The independent Schmidt genome-wide functional screens add two fixed donors, CRISPRa and
CRISPRi, and IL2/CD4 and IFNG/CD8 contexts. The runner verifies the exact Zenodo archive,
all used members and schemas, and the deposited author transformation before analysis.
Eligibility uses only common gene identity and guide count. At the primary threshold of at
least three paired guides per gene/donor, 18,568 genes are complete in all four screens.

Whole-universe same-reagent donor concordance is limited:

| Screen | Signed Spearman | Absolute-effect Spearman | Sign agreement |
|---|---:|---:|---:|
| CRISPRa IFNG | 0.312 | 0.121 | 0.599 |
| CRISPRa IL2 | 0.332 | 0.129 | 0.603 |
| CRISPRi IFNG | 0.137 | 0.079 | 0.538 |
| CRISPRi IL2 | 0.135 | 0.088 | 0.533 |

After orienting CRISPRi as in the author script, whole-universe CRISPRa-versus-CRISPRi
signed Spearman is only **0.036** for IFNG and **0.020** for IL2. Those comparisons also
change Calabrese versus Dolcetto guide libraries. IL2-versus-IFNG signed Spearman is 0.337
for CRISPRa and 0.145 for CRISPRi, but cytokine changes jointly with CD4/CD8 cell type and
screen context.

Conditioning on the training donor's 200 largest absolute effects produces stronger—but
selection-conditional—held-donor results across eight correlated directions per class:

| Frozen source top-200 evaluation | Median signed Spearman | Median absolute-effect Spearman | Sign agreement | Global held-target top-200 overlap |
|---|---:|---:|---:|---:|
| Same screen, held donor | 0.887 | 0.608 | 0.945 | 0.458 |
| Other donor + modality + guide library, same context | 0.300 | 0.183 | 0.597 | 0.098 |
| Other donor + cytokine + cell type, same modality | 0.749 | 0.395 | 0.863 | 0.340 |

The last two rows jointly change donor, so neither difference is attributable to the
named screen axis alone.

The full ≥1/3/6-guide × top-50/100/200/500 grid is retained and explicitly labeled
exploratory/post-hoc. The same guides are reused within each modality, so the positive
control is same-reagent reproducibility rather than guide-held-out replication. Together,
the result shows that conditional extreme-effect concordance can coexist with weak
whole-genome agreement. It does not supply paired transcriptomic effect dictionaries,
modality or library equivalence, cytokine/cell-type generality, donor-population inference,
functional efficacy, state reachability, or a target recommendation.

Source: [`schmidt_external_validation.json`](../results/schmidt_external_validation.json).

## What remains unknown

All seven release requirements remain open:

1. leakage-safe pseudobulk donor- and guide-held-out evaluation;
2. module/pathway holdouts and calibrated structured nulls;
3. nested PCA, ridge, unconstrained and capacity-matched baselines;
4. same-experiment, guide-burden- and control-matched measured in-domain perturbation
   combinations;
5. paired CRISPRi and CRISPRa dictionaries;
6. independent whole-state protein, chromatin, function, fitness and durability; and
7. prospective established-state validation.

Here guide-held-out means an exact, hash-cross-verified rank-to-sgRNA ID mapping—not rank
labels alone. That crosswalk plus structural-QC-only pseudobulk physical-guide evaluation
is the next gate. Until all applicable gates pass, the appropriate output is model-relative
geometry—not a reachability verdict or intervention recommendation.
