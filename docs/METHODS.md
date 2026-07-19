# Method

*This project's contribution is the method described here: fail-closed, model-relative
point-estimate cone geometry with a numerical infeasibility certificate, backed by the
adversarial harness in `scripts/run_validation_harness.py`. The case-study sections below
exercise the method without turning numerical representability into a biological verdict.*

## Question and geometry

Given measured perturbation profiles and a target transcriptional direction in a shared
gene space, how closely can a non-negative linear combination of those profiles align with
the target — and, when it cannot, can we certify that the target lies outside the reachable
cone? The method tests one dictionary, target, and metric; it is narrower than
perturbation-response prediction. It reports both the best fit and whether the declared
numerical cone-membership contract was established, and it fails closed when its numerical
diagnostics do not certify that contract.

Let `E` have shape perturbations × genes, `A = E.T`, and `d` be the target. Coefficients
solve

```text
minimize_w>=0  1/2 ||sqrt(W) (d - A w)||²
```

for diagonal gene metric `W`. Zero-weight coordinates are excluded and the largest
retained weight is normalized to one. The fitted direction is `A w`; the residual is
`rho = d - A w`; and a metric-relative separator is proportional to `W rho`. The separator
describes this measured cone only. It does not prove biological impossibility or identify
a gene that must be activated.

## Numerical contract

The core reports coefficients, fitted direction, weighted cosine, residual fraction,
relative objective, atom-scale-invariant KKT diagnostics, numerical geometry status, and
a separator only when residual energy clears a declared tolerance. KKT, polarity, and
orthogonality certify at `1e-8` or the projection fails closed. Exact-zero, non-finite,
misaligned, all-zero-weight, and unrepresentable-scale inputs fail before a result is
returned. No biological verdict threshold exists. In noisy high-dimensional data, exact
point-estimate cone membership is generically stringent; the separator certifies the
declared numerical problem, not measurement-stable biological infeasibility.

## Retrospective catalog coverage and candidate scoring

`library_coverage.py` applies the same projection to rows of a target catalog. It keeps two
different summaries explicit:

- **strict membership** is the core's `inside_tolerance` status;
- **thresholded coverage** is the fraction whose best cosine clears a declared,
  application-chosen bar.

Those quantities are not interchangeable. `atom_redundancy` recomputes catalog summaries
after removing each library atom. `rank_acquisitions` scores supplied candidate effect
vectors against the unit dual-separator normals of strict-outside targets, weighted by
residual fraction. When gene weights are present, the certified separator—not the raw
residual—is the relevant normal. The optional realized comparator appends each candidate
and recomputes mean catalog cosine. It is a deterministic retrospective comparator, not an
experimental ground truth; the candidate effect vectors are already measured inputs.

The v3 cross-dataset runner applies this layer to Zhu, Norman, and Replogle effect
dictionaries. Zhu and Replogle split measured rows before selecting the 400 highest-variance
features from the current 50-atom library; catalog and candidate outcomes cannot select
features. Norman's descriptive all-single comparison uses 152 position-specific singles and
131 doubles. Its acquisition audit first groups 105 canonical genes, freezes one measured
representative per gene (`GENE+ctrl`, with one declared fallback), and then assigns 40 genes
to the library and 30 to the supplied-candidate pool. Canonical genes cannot cross roles.

Twelve deterministic partitions report the cone, best-single ray, common-response ray, and
unconstrained signed-span capacity ceiling at cosine thresholds 0.3, 0.5, 0.7, and 0.9.
Each split must satisfy signed span ≥ cone ≥ either ray within numerical tolerance. Norman
uses the same measured representatives and 40-gene library size as its acquisition audit,
against one fixed 60-double catalog. These partitions measure algorithmic sensitivity over
correlated observed rows, not sampling uncertainty. Candidate labels are unique and one-to-one
with canonical labels, so exact-row and canonical top-1 agreement are identical and reported
once.

Zhu and Replogle remain same-screen compression audits; Norman is a measured
single-to-double additivity/alignment diagnostic. None is a prospective library-design
benchmark. The three pickle-free NPZ inputs are rebuilt deterministically from registered,
hash-matching raw bytes by `scripts/build_library_coverage_caches.py`. Public redistribution
URLs remain mutable, so reconstructability does not replace a durable release archive.

## Source-bound held-out challenge

`held_out_alignment` fits on one gene subset and scores the frozen coefficients on a
disjoint subset. Source-bound case-study genes are ordered lexicographically by symbol;
NumPy `default_rng(seed)` generates half-gene splits for seeds 0–11; and each split's
fit/score gene identifiers are SHA-256-bound in the report. Random-gene splits treat
correlated genes as exchangeable, so module/pathway holdout remains required.

The former 60-shuffle target diagnostic was retired because its deleted intermediate did
not preserve gene order. No claim-bearing p-value is attached to the current case-study
cosine. The systemic harness separately tests structured common-response, module,
sign-selection, and multiplicity failure modes.

## Case-study inputs and target construction

- Dictionary: donor-collapsed primary-human-CD4 CRISPRi differential-expression profiles
  in the source study's post-engineering/post-expansion `Rest` condition.
- Target: source-study-reused, selectively constructed cross-sectional population
  contrast oriented toward the reported Th1 centroid; not independent or a trajectory.
- Unit: perturbation-condition profile, not donor effect.
- Primary source-transfer scale: log fold change; z-score is sensitivity-only.

The registered merged target retains genes present in both Ota and Höllbacher target
contrasts, requires concordant effect signs, averages their z-scores, reverses the reported
Th2-versus-Th1 orientation, and intersects with screen genes. Counts are 25,672 source
union, 11,616 shared, 7,960 sign-concordant, and 6,188 registered coordinates. Every
availability, sign, and screen filter is a selection step.

The target sources are Ota et al., *Cell* 2021, DOI
`10.1016/j.cell.2021.03.056` (E-GEAD-397), and Höllbacher et al.,
*ImmunoHorizons* 2020, DOI `10.4049/immunohorizons.2000037` (GSE149090). The
dictionary is Zhu et al. 2025, DOI `10.64898/2025.12.23.696273` (SRP643211 / GSE314342).
The exact target CSV and 16.8 GB H5AD bytes are verified before source-bound analysis.

## Cross-source transfer

Ota and Höllbacher are evaluated on the 8,950 genes present in both sources and the
screen. This universe uses identifiers/presence only; held-out-source signs and magnitudes
do not select coordinates. For each direction and split, coefficients are fit to one
source on fit genes and scored against the other on score genes. Mean-response and
best-single baselines are trained on the fit source. Improvement is reported against the
larger test cosine of those comparators, a disclosed test-selected descriptive baseline.
Cosine, normalized RMSE, norm ratio, sign agreement, support size, and split hashes are
retained.

## Independent Arce benchmark

Arce S1 is read directly from its hash-verified Zenodo archive. Genes require exactly
four guides in all three contexts. Intersecting with Zhu `Rest` availability and
source-side `ontarget_significant` admission—without reading Arce outcomes—leaves 480
targets. S1 positive log fold change denotes enrichment of knockouts lowering CD25/IL2RA,
so the regulator score is frozen as negative Zhu IL2RA transcript log fold change.

Per context, the runner reports signed Spearman/Kendall after the prespecified orientation,
direction agreement, absolute-magnitude top-25/50/100 overlap, and target-label permutation
diagnostics. Directional agreement uses a null center fixed analytically from the two sign
margins. Targets, not cells or donors, are the exchangeable units. All 18 permutation
p-values are unadjusted exploratory diagnostics across correlated contexts/metrics and
provide no multiplicity-controlled inference or donor uncertainty. The claim ceiling is
cross-study/cross-modality ranking alignment, not whole-state reachability or causal
validation.

S14 is streamed from the same verified archive with its complete 32-column schema checked,
while retaining only registered QC and score fields. Every cell must have one sgRNA,
`has_sgrna=true`, a unique cell ID, finite supplied `activation.score`, and global HTO
Singlet status. For target *t*, guide *g*, donor *d*, and context *c*, the guide contrast is
the guide-cell median minus the median of the nine Non-Targeting guide-cell medians within
the same donor/context. The two guide contrasts are equally weighted for a target/donor
summary. All strata are retained and their cell counts disclosed.

The runner exactly reproduces S8 pooled target/context means and medians as an archival
provenance gate. It does not use S8 pooled-cell p-values, bootstrap cells as donors, or
emit donor-population inference from two donors. Guide-pair and donor rank/sign concordance
are descriptive. Because the local tables do not freeze the supplied score's genes,
formula, normalization, or independence, it is not equated with functional activation or
with the one-gene Zhu IL2RA predictor.
The 28 regulators were preselected by the source authors using prior screen/state-specific
evidence; target-rank and sign concordance are therefore conditional panel diagnostics,
not genome-wide estimates or independent validation.

## Schmidt two-donor functional-screen transfer stress

The Schmidt archive is opened only after its 26,152,593 bytes, SHA-256, upstream MD5, and
exact ZIP member allow-list match the frozen config. The runner then verifies the author R
script and each of four complete MAGeCK sgRNA tables: CRISPRa/CRISPRi × IL2/IFNG. It rejects
schema drift, unexpected members, unsafe paths, duplicate sgRNA keys, nonfinite numeric
fields, malformed donor suffixes, unpaired guide identities, and row-count drift. The
paper states that the two modalities use the same two donors; the deposited script maps
`r0`/`r1` to Donor1/Donor2 and defines each donor gene effect as the median guide LFC.

Orientation is prespecified from the author script: CRISPRa LFC is retained, while CRISPRi
LFC is multiplied by −1 so that positive values denote positive cytokine regulators. This
supports signed rank/sign comparisons only; it does not make the modalities scale matched.
The outcome-independent universe excludes `NO-TARGET`, requires gene identity in all four
screens and both donors, and applies only a unique-guide-count threshold. LFC, sign, FDR,
rank, and hit calls never enter eligibility. The primary ≥3-guide universe has 18,568 genes;
thresholds 1/3/6 are retained as an explicitly exploratory sensitivity.

Whole-universe metrics report signed Spearman/Kendall, sign agreement, absolute-effect
Spearman, and cosine for donor reproducibility within screen, CRISPRa-versus-CRISPRi within
context, and IL2-versus-IFNG within modality. The second comparison is also a
Calabrese-versus-Dolcetto guide-library comparison; the third is inseparable from CD4/CD8
cell type and screen context.

For conditional top-effect transfer, each source screen and donor direction selects the
top *K* genes by absolute oriented effect using the training donor only. That exact gene
set is frozen before scoring the other donor in the same screen, the other donor plus
modality/library in the same context, or the other donor plus cytokine/cell-type context
in the same modality. The cross-screen rows therefore estimate joint shifts, not isolated
modality/library or context effects. Global
held-target top-*K* overlap is computed against the complete frozen universe, not within
the source-selected subset. The primary *K*=200 result is shown with 50/100/200/500
sensitivity. The eight directions in each class are correlated descriptive challenges;
no p-values, confidence intervals, or donor-population tests are emitted. Guides are reused
within a modality, so the ceiling is same-reagent reproducibility in two fixed donors—not
guide-held-out replication, transcriptomic reachability, modality equivalence, context
generality, efficacy, or target recommendation.

## Source-selected arrayed RNA/protein follow-up

The Zhu author repository contributes arrayed bulk-RNA differential-expression summaries
for nine source-selected perturbations and IL-10/IL-21 flow percentages in follow-up donor
labels Donor5–Donor10. The tables and aggregate H5AD are full-file SHA-256-bound. The
screen side uses `Stim8hr` log fold change without consulting aggregate admission,
guide-correlation, donor-correlation, significance, or follow-up outcomes.

Bulk and screen profiles intersect on 8,976 Ensembl IDs. All nine panel target genes are
masked from every profile before computing cosine, Pearson, Spearman, normalized RMSE, and
retrieval rank among all nine source profiles. Raw profiles are compared with the common
source-panel response. A specificity sensitivity subtracts the across-nine-perturbation
mean separately for every gene and assay before rescoring; this panel centering is not a
replacement biological estimand.

Within each donor and cytokine, the flow baseline is the arithmetic mean of every
available NTC-labeled measurement. A perturbation effect is
`log2(percent_positive / donor_NTC_mean)`. Donor measurements are never pooled as cells;
target summaries are the equal-donor median, with mean/range and donor count retained.
Coverage is unbalanced and no missing target/donor is imputed.

All `9! = 362,880` target-label assignments are enumerated. Retrieval top-1/MRR and four
primary Spearman associations (screen/bulk RNA × IL-10/IL-21 donor-median flow) are
recomputed synchronously; the correlation family records a one-sided maximum-statistic
tail fraction. Because donor count, composition, and NTC replication differ by target,
target exchangeability is not established. These exhaustive fractions are conditional
diagnostics, not inferential p-values or multiplicity-adjusted inference. They do not
replay upstream panel selection and provide neither held-out-discovery nor donor-population
inference. The claim ceiling is measured same-study cross-platform replication and
cytokine consistency.

## Staged donor/guide validation contract

The released donor-pair H5MU stage is implemented on `Rest` log fold change. Its six
modalities represent every two-donor combination among four donors. The analysis requires
the 1,584 perturbation atoms present in every modality and constructs an 8,949-gene target
universe from both target sources plus donor-H5MU gene symbols without using held-out
source signs. For each of three complementary 2-vs-2 partitions, both donor directions,
both target-source directions, and seeds 0–2, the runner fits NNLS on one donor pair,
source, and half-gene split. It applies identical coefficients to the complementary donor
pair, other source, and held-out genes.

The common-response ray and best-single atom/scalar are selected solely on the training
dictionary, target source, and fit genes; their identity and scale remain frozen on test.
The two mixed-run donor partitions yield 24 primary descriptive challenges. The donors
1–2 versus 3–4 partition yields 12 fully run-confounded sensitivities reported separately.
All challenges are correlated and no donor-population p-value or confidence interval is
emitted. The released modalities omit DE-ineligible targets, so complete-case presence is
potentially effectiveness-selected.

The released guide-rank H5MU stage is also implemented on `Rest` log fold change. The
29,424,424,894-byte input is bound to SHA-256
`964eeafb3356a7322a1d5b1121802c6a1433456f3591e2d5797817df3bf9c2f6`; root, modality,
observation, variable, layer, and categorical encodings are checked before analysis.

The official VCP v1.0 dataset card (accessed 2026-07-19) defines `guide_1` and `guide_2`
as the first and second sgRNA ID after alphanumeric sorting within target-condition, and
public pseudobulk/guide-library artifacts carry sgRNA IDs. The ranked H5MU modalities do
not themselves embed sgRNA ID or sequence. This repository has not reconstructed and
hash-cross-verified their exact rank-to-sgRNA crosswalk, so report rows retain opaque H5MU
keys and rank labels rather than naming sgRNAs. `guide_2` is an exact subset of `guide_1`
in the H5MU.

The reciprocal universe retains all 8,323 common category-labeled `Rest` opaque keys
without outcome-ranked top-*K* selection and reports the 2,752 `guide_1`-only
category-labeled `Rest` keys as excluded. Another 35 nominal `guide_1` `Rest` keys are
withheld because their categorical metadata are missing. Across the full object, the 114
`guide_1` categorical-missing rows have suffixes `Rest` 35 / `Stim8hr` 18 / `Stim48hr`
61; the 38 `guide_2` rows have `Rest` 0 / `Stim8hr` 0 / `Stim48hr` 38. Missing categories
are recorded, not imputed, and every selected category-labeled common-`Rest` row has an
exact target mapping. The 8,950 target coordinates are selected only by shared
Ota/Höllbacher and H5MU identifiers; no held-source sign or magnitude enters the universe.

For each of two training ranks, two target sources, and seeds 0–2, NNLS is fit once on
half the genes. Its coefficients, common-response scale, and best-single opaque key/scale
are then frozen for both the within-rank and reciprocal-rank dictionaries on the disjoint
genes. Each of the 12 fits is scored against both the training target source and the
opposite source, yielding 24 correlated rows: 12 isolate released guide-rank transfer and
12 jointly change guide rank plus target source. A thirteenth,
non-claim-bearing fit reverses the atom order; its maximum relative L2 difference from the
canonical representation is `1.76e-15`. This checks one solver representation, not
coefficient identifiability.

At author commit `848d62f`, the guide-DE pipeline first applied
`keep_min_cells & keep_effective_guides & keep_total_counts`, requiring at least three
passing replicates per guide-condition and at least five cells per guide in each
condition/sample. It constructed `cond_targets` from the unreleased
`for_DE_by_guide.csv`, retained exactly two testable guides per target-condition, applied
the final `keep_test_genes` filter, and fit `~ log10_n_cells + target` without a donor
term. Released presence can therefore select on effectiveness.

A leakage-safe donor/guide claim therefore requires structural-QC re-estimation from the
pseudobulk count object after an exact rank-to-sgRNA ID crosswalk. Universe construction
must exclude `keep_effective_guides`, aggregate `ontarget_significant`, held-out
DE/sign/correlation fields, and any filter whose outcome independence is not demonstrated.
The fixed reciprocal-rank result does not establish named-sgRNA replication or failure.
With four donors, all donor partitions remain fixed-cohort descriptive challenges; exact
donor-population significance at 0.05 is mathematically unavailable.

## Goudy triple cross-experiment stress contract

`scripts/run_goudy_combination_validation.py` is a negative **STRESS** benchmark on one
measured FAS + RC3H1 + SUV39H1 triple across four fixed donor labels. Parsing begins only
after exact SHA-256 and MD5 verification of three registered inputs: the deposited
GSE306915 normalized-count archive, the GEO family SOFT file, and the author sample key
pinned to commit `53155c9`. The SOFT and author key are crosswalked to the count columns,
GSMs, titles, experimental day, guide identity, and experiment identifier. The key confirms
36/40 declared analysis roles; donor-3/4 multiplex controls and triples are absent from it,
and the triple guide identity is unresolved. Constituent singles and AAVS1 controls are
author-key-confirmed in Co065, while donor-1/2 multiplex controls and triples are confirmed
in Co066; the corresponding donor-3/4 experiment identifiers remain unresolved.

Deposited non-negative values are used as author-normalized linear CPM without a second
library-size normalization, then transformed sample-wise as `log2(CPM + 1)`. Constituent
and unrelated-single effects subtract the within-donor mean of two AAVS1-guide controls;
the triple effect subtracts that donor's multiplex NTC. The primary universe requires mean
linear CPM ≥1 across the 12 declared controls only. Sensitivity repeats the complete score
at thresholds 0.5, 1, 2, 5, and 10; perturbation outcomes never select the mask. FAS,
RC3H1, and SUV39H1 rows are excluded from transcriptome scoring. At thresholds 5 and 10,
SUV39H1 has already failed the control-only filter, so only the two retained target rows
are removed; the scoring coordinates remain disjoint from every retained target coordinate.

For each donor, the component cone is an in-sample NNLS projection of the measured triple
onto the three target-matched constituent singles. The equal-sum comparator fixes all three
coefficients to one, and the in-sample best-single oracle selects and scales one constituent
on that donor's triple. Leave-one-donor-out (LODO) scoring instead builds the expression
mask from training-donor controls, fits one shared coefficient vector by stacking only the
three training donors, and freezes it on the held donor. Its best-single baseline likewise
selects identity and non-negative scale from training donors only. A separate on-target
diagnostic fits three coefficients on the three target coordinates, freezes them, and scores
the disjoint target-masked transcriptome. The 3×3 fit is saturated and is retained only as a
calibration diagnostic.

The additive residual is the measured triple effect minus the unit-coefficient sum of the
constituent effects. It cannot separate a statistical interaction from experiment,
control-type, or guide-burden differences. Reliability diagnostics report pairwise-donor
effect cosines, AAVS1 guide-to-guide expression agreement and difference norms relative to
the measured effects, and control-filter sensitivity. Module error is unavailable because
no module set was preregistered; choosing one after seeing these outcomes would not be a
valid holdout. The claim ceiling is therefore a descriptive negative test of whether this
single measured triple is represented by target-matched singles under inseparable
cross-experiment and control designs—not a general additivity, population, prospective,
state-conversion, or functional claim.

## Removed legacy behavior

The maintained surface intentionally excludes categorical reachability labels, staged
`E`/`-E` modality proxies, activation certificates, greedy recipes, uncalibrated analytic
nulls, automatic experiment recommendations, and evidence tables derived from deleted,
unhashed intermediates. Git history preserves provenance without presenting those outputs
as current evidence.

## Reproducibility

`./reproduce.sh` runs exact numerical tests and the demo, checks the frozen systemic
harness, validates every canonical scalar against the source-bound and Arce reports,
rejects any failed report, and verifies artifact hashes. External runners regenerate their
reports from gitignored registered inputs; acquisition is documented in
[`data/README.md`](../data/README.md). The full cross-dataset acquisition audit performs
thousands of cone projections and typically takes about 5–10 minutes on a laptop.
