# Scientific Validation Plan

**Review status (2026-07-19): strengthened measured follow-up and repaired the first
catalog-extension audit; still not ready for inferential or biological reachability
claims.** The repository supports model-relative point-estimate geometry, retrospective
donor-collapsed diagnostics, source-selected arrayed RNA/cytokine replication, and
same-screen effect-dictionary compression, plus bounded negative released-guide-position
and cross-experiment CRISPRoff stress reports. It does not yet support a calibrated
reachability rule, prospective library design, state conversion, intervention utility,
donor generalization, target-specific significance, or coefficient-based prioritization.

This document is the execution contract for advancing the project. It separates work
that can be completed from public data from claims that require a prospective experiment.

## Corrections established during adversarial review

1. The registered merged target contains **6,188 genes**, not 9,831. The complete lineage
   is 25,672 genes in the union of target-source tables, 11,616 shared by both sources,
   7,960 with concordant signs, and 6,188 after screen intersection. The separate 9,831
   value is only the raw target-table/screen overlap.
2. The Höllbacher and Ota targets are cross-sectional population contrasts, not observed
   Th2-to-Th1 trajectories. A line between population centroids is not a demonstrated
   dynamical path.
3. Zhu `Rest` is a post-engineering/post-expansion assay condition from initially naïve
   primary CD4 cells, not an established memory-Th2 starting state.
4. Zhu already used these perturbation signatures and its
   `4_polarization_signatures` workflow to nominate Th1/Th2 regulators. This project is a
   methods reanalysis and audit, not an independent biological discovery.
5. Z-score atoms combine effect and precision. Their non-negative coefficients are
   geometric weights, not perturbation doses or additive intervention amounts.

Primary anchors: [Zhu et al.](https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1),
[source analysis](https://github.com/emdann/GWT_perturbseq_analysis_2025),
[Höllbacher et al.](https://doi.org/10.4049/immunohorizons.2000037), and
[Ota et al.](https://doi.org/10.1016/j.cell.2021.03.056).

## Current data-quality assessment

The intended analytical grain is one perturbation-condition differential-expression
profile by gene. It is not a donor replicate or a cell.

| Dimension | Verified state | Consequence |
|---|---|---|
| Main matrix | 33,983 perturbation-condition rows × 10,282 unique Ensembl genes | Dense summary is suitable for geometry, not donor inference |
| Conditions | Rest 11,287; Stim8hr 11,415; Stim48hr 11,281 profiles | Context comparison is available |
| Admitted generators | Rest 6,871; Stim8hr 7,155; Stim48hr 7,195 source-significant targets | Admission is upstream-selected, not independently tested here |
| Biological samples | 4 donors × 3 conditions = 12 libraries | Donor-level files exist but the maintained result collapses them |
| Run structure | Rest/Stim8hr donors 1–2 on R1 and 3–4 on R2; all Stim48hr on R2 | The 1–2 vs 3–4 Rest/Stim8 split is run-confounded; mixed-run donor partitions are run-balanced; cross-condition Stim48 support is limited |
| Released guide positions | 8,323 common `Rest` atoms; `guide_2` is a subset of `guide_1`; 2,752 `guide_1`-only rows | Positional labels only; physical guide IDs are absent and released presence is effectiveness-selected |
| Effect layers | log FC, z-score, p, adjusted p, base mean, and LFC SE | Log-FC primary analysis with uncertainty is feasible |
| Target input | 6,188 final coordinates after source-sharing, sign, and screen filters | Selection leakage and excluded-coordinate scope must be explicit |
| Missing inference unit | Current result has no donor-level sampling distribution | Random-gene split SD must never be reported as biological uncertainty |

The donor-pair and reciprocal guide-position published-pipeline sensitivities and compact
arrayed bulk-RNA/IL-10/IL-21 follow-up are now executable. The next primary-model gate is
leakage-safe pseudobulk re-estimation with physical guide identities.

## Execution program

### Workstream 0 — software and lineage certification

**Status:** implemented in this revision.

Actions:

1. Normalize target scale before solving and avoid overflow-prone squared diagnostics.
2. Refuse uncertified KKT or separator outputs.
3. Add an exhaustive active-set oracle for small NNLS problems.
4. Add labelled axes, source hashes, namespace, units, modality, context, and timepoint.
5. Add deterministic grouped splits, maxT p-values, fault injection, and a frozen JSON
   harness report.
6. Cross-assert the frozen 25,672 → 11,616 → 7,960 → 6,188 lineage counts. Source-level
   recomputation remains part of Workstream 1.

Exit gate: all input corruptions fail closed; oracle fitted-point error ≤1e-8; degenerate
fit error ≤1e-8; no group leakage; maxT simulation familywise error ≤0.10 in the fast CI
check; deterministic report matches its frozen artifact.

### Workstream 1 — reconstruct real-data lineage

**Lead experts:** computational biologist/data engineer, with statistical-genetics review.

Inputs from the public Zhu release:

- `GWCD4i.pseudobulk_merged.h5ad`;
- `GWCD4i.DE_stats.by_donors.h5mu`;
- `GWCD4i.DE_stats.by_guide.h5mu`;
- `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv`;
- the existing `GWCD4i.DE_stats.h5ad` and polarization table.

Actions:

1. Create a dataset card recording accession, exact object key, bytes, SHA-256, license,
   gene namespace, orientation, units, donors, guides, conditions, runs, and missingness.
2. Build one versioned `source summary → target construction → alignment → split → fit →
   null → report` command. Never silently substitute a collapsed average for a missing
   donor or guide estimate.
3. Make log fold change the primary effect scale; use LFC SE for declared precision
   weighting/uncertainty. Keep z-score geometry as sensitivity-only.
4. Reproduce both source-specific targets without sign-based cross-source selection.
5. Reproduce Zhu's `pert2state_model` ranking as the required same-question baseline.

Current status: the aggregate H5AD/target phase and released donor-pair and guide-position
H5MU sensitivities are implemented with full input hashes; the unreconstructable archived
0.446 multisplit result was withdrawn and replaced by the source-bound 0.444 result. The
donor test is not leakage-safe because released presence is eligibility-selected. The guide
test is likewise effectiveness-selected and exposes only positional modality labels, not
physical guide identities. Pseudobulk re-estimation and the same-question baseline remain
exit blockers.

### Workstream 2 — independent transfer and structured inference

**Lead experts:** statistical geneticist, perturbation-ML scientist, and immunologist.

Actions:

1. **Source transfer:** fit/tune with Ota, then score the frozen direction against
   Höllbacher; reverse the roles. Define the shared screen-measurable universe without
   using held-out source signs or magnitudes.
2. **Donor-pair transfer:** fit on one donor pair and score on the disjoint pair across
   every complementary 2-vs-2 partition. Report donors 1–2 vs 3–4 as fully
   run-confounded, and donors 1–3 vs 2–4 and 1–4 vs 2–3 as run-balanced partitions. Treat
   the released H5MU as published-eligibility sensitivity because DE-ineligible targets
   are absent.
3. **Released guide-position transfer:** the implemented stress fits each positional
   modality against both target sources and evaluates the frozen coefficients within and
   across positions. It retains all 8,323 common `Rest` atoms, reports 2,752
   `guide_1`-only rows, and emits no guide-generalization inference. The labels do not map
   to released physical guide IDs, and modality presence reflects the authors'
   `keep_effective_guides` filter.
4. **Leakage-safe re-estimation — next gate:** rebuild donor and physical-guide effects
   from pseudobulk counts using only declared coverage/depth QC. Exclude
   `keep_effective_guides`, aggregate admission, held-out DE/sign/correlation fields, and
   any filter not proven outcome-independent. Fit on one identified guide, evaluate its
   identified mate, reverse the roles, and report excluded single-guide targets.
5. **Grouped holdouts:** use pathway/module, source, donor, guide, context, and run blocks.
6. **Structured nulls:** preserve gene covariance using rotations or module-level
   randomization. Use ≥9,999 resamples for a primary family and maxT/Westfall–Young control
   across contexts, target directions, metrics, and preprocessing variants.
7. **Nested baselines:** tune only inside training folds. Include best single atom,
   common-effect mean, PCA/ridge, unconstrained least squares, sparse/non-negative ridge,
   and capacity-matched random dictionaries.
8. **Metrics:** cosine, normalized RMSE, norm ratio, DE-sign recovery, top-DE precision and
   recall, pathway-direction agreement, target retrieval rank, and improvement over the
   strongest frozen baseline.
9. **Identifiability:** report active-matrix rank, singular values, strict-complementarity
   margin, fitted-point stability, and coefficient equivalence sets. Never interpret one
   coefficient vector when alternatives fit equivalently.

Exit gates:

- the multiplicity-adjusted lower confidence bound for paired improvement over the frozen
  best nested baseline is above zero, plus a preregistered practical margin for any utility claim;
- source-, donor-, and guide-held-out effects agree in direction;
- structured-null familywise error is calibrated and primary adjusted p ≤0.05;
- 95% interval coverage is 0.925–0.975 in simulation;
- calibration cannot pass by universal abstention;
- no headline depends on sign-concordance selection using its evaluation source.

Current donor-pair result: all six modalities share 1,584 complete `Rest` atoms and an
8,949-gene source-defined universe. Across 24 run-balanced correlated challenges, median
cosine is 0.031 with +0.032 gain over the training-selected frozen best single, but median
normalized RMSE is worse (1.153 versus 1.018) and improves in only 1/24 challenges. This
supports neither predictive utility nor donor-population generalization. It establishes
the execution contract and makes the guide/pseudobulk stages more urgent, not less.

Current released guide-position result: 12 fits yield 24 correlated rows. Same-source
held-gene median cone cosine falls from 0.251227 within position to −0.019197 reciprocally
(paired change −0.291089); gain over the training-selected best single is positive in only
3/12 rows, and normalized-RMSE gain is positive in 0/12. Joint guide-position plus
target-source transfer is also negative (median cosine −0.034383; gain −0.041518). The result
falsifies robustness across released positional summaries, not replication between known
physical guides.

### Workstream 3 — source-selected arrayed RNA/cytokine follow-up

**Status:** implemented for the author-released IL-10/IL-21 follow-up.

The pinned author repository supplies arrayed bulk-RNA profiles for nine source-selected
perturbations and IL-10/IL-21 flow percentages in Donor5–Donor10. The maintained runner
masks all nine panel target genes, evaluates 8,967 transcript coordinates, adds
Systema-style panel centering, performs nine-way retrieval, normalizes flow within donor
against all available NTC measurements, and exhaustively enumerates all 9! synchronized
target-label assignments as conditional diagnostics rather than inferential p-values.

Current result: matching perturbation retrieval is 9/9 raw and panel-centered; median
centered cosine is 0.580; RNA-to-flow cytokine rank associations are 0.717–0.850. Median
normalized RMSE remains at or above one, donor coverage is uneven, and all targets were
selected upstream. The exit claim is therefore same-study cross-platform replication and
cytokine consistency only—not held-out discovery, population generalization, or state
conversion. The separately documented `Th1Th2_validation_summary.suppl_table.csv` remains
optional until a working immutable source is verified; it is no longer an execution blocker.

### Workstream 4 — combination and modality stress tests

**Lead experts:** perturbation biologist, causal/ML evaluator, and data curator.

1. **Implemented negative cross-experiment stress:** [Goudy et al. GSE306915](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306915)
   compares one measured `FAS + SUV39H1 + RC3H1` triple with target-matched constituent
   singles across four fixed donor labels under different experiment, control, and
   guide-burden designs.
   Execution is `PASS`; geometry is `FAILS_DECLARED_GEOMETRIC_MODEL`; biological status is
   `INCONCLUSIVE_CROSS_EXPERIMENT_CONFOUNDING_LOW_RELIABILITY`. Component-cone median
   cosine/nRMSE is 0.0949/0.9952 with strict membership 0/4; equal-sum cosine/nRMSE is
   0.0646/2.3084, LODO cosine is 0.0881, and triple pairwise-label cosine is 0.0480. No
   declared control-only expression threshold rescues the model. Triple guide identity and
   the D3/D4 triple/NTC experiment IDs remain unresolved. Use
   [GSE306917](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306917) only as a
   two-donor CD55/CD81 CRISPRoff mechanism check, not mechanism evidence for the triple.
2. **Implemented compact Schmidt stress:** [Schmidt et al.](https://doi.org/10.1126/science.abj4008)
   GSE174292/GSE174255 contributes two-donor genome-wide CRISPRa/CRISPRi IL2/IFNG
   functional screens. The hash-gated runner freezes author orientation, uses guide
   coverage only for eligibility, co-reports weak whole-genome agreement with conditional
   top-effect transfer, and exposes all 1/3/6-guide × top-50/100/200/500 sensitivities.
   CRISPRa/i also changes Calabrese/Dolcetto libraries; IL2/IFNG also changes CD4/CD8;
   every cross-screen top-effect row also changes donor, and guides are reused within
   modality. The result is therefore not an estimate of any isolated axis, scale equivalence,
   guide-held-out replication, context generality, donor-population inference, or the open
   paired-transcriptomic-dictionary requirement. GSE174284 and related same-paper
   GSE190604/GSE190846 remain potential orthogonal follow-up rather than inputs here.
3. Use [Arce et al.](https://doi.org/10.1038/s41586-024-08314-y) GSE271090/GSE278572 for
   context sign reversals and RNA/protein agreement. The first S1 IL2RA/CD25 ranking
   benchmark and S14/S8 donor/guide supplied-score robustness are implemented; orthogonal
   S20/S4 RNA/protein endpoints and independent donors remain.
4. Retain Norman/Replogle only as out-of-domain numerical and combination benchmarks.
   Norman is K562 in the original study; any redistributed object with conflicting cell-type
   metadata must be identity-checked before promotion.
5. Treat the catalog layer as a retrospective effect-dictionary compression/oracle audit:
   candidate vectors are already measured, and Zhu/Replogle library, catalog, and candidates
   come from the same screen. Never call this prospective acquisition.
6. Freeze strict `inside_tolerance` separately from cosine-threshold coverage. For Norman,
   compare the all-single cone with constituent sum, constituent-only cone, best component,
   unrelated singles, and common response. For Zhu/Replogle, report gains over best-single
   and common-response rays rather than raw threshold fractions alone.
7. Split atoms before any benchmark-only feature or magnitude selection; select features
   from the declared training/library side. Sweep frozen seeds, library sizes, candidate
   sizes, and thresholds, and report uncertainty over configurations.
8. Store input bytes/hashes, source versions, licenses, schema, preprocessing differences,
   split IDs, both candidate orders, and per-candidate scores in a frozen config/report.

Exit gate: reported conclusions name the assay, cell system, modality, donor unit, and
claim ceiling. No cancer-cell-line result is used to validate primary-CD4 biology.
The Goudy stress report meets only this descriptive gate; the open primary-T-cell
requirement is a same-experiment measured design with verified guide identities, every
constituent single, and control- and guide-burden matching.

### Workstream 5 — prospective state-transition test

**Lead experts:** human T-cell experimentalist, cytometry/multiome specialist,
biostatistician, and independent preregistration reviewer.

#### Stage A: feasibility

Use three donors and analyze two starting systems separately:

- FACS-purified established memory Th2 cells: CD3+ CD4+ CD45RA− CD25low CD127high
  CXCR3− CCR6− CCR4+, with CRTH2+ and CRTH2− strata separated or analyzed explicitly;
- naïve-derived in-vitro polarized Th2 cells as a distinct mechanistic model.

Controls: non-targeting, mock, vector-burden, Th2 maintenance, Th1-polarizing cytokines,
GATA3/STAT6 CRISPRi, TBX21 or IL12RB2/STAT4 CRISPRa, and direction-opposite controls.

Feasibility gates: ≥95% starting purity; verified editing with two concordant guides;
≥80% control viability/recovery; positive and negative controls separate on the frozen
state axis. Verify baseline GATA3 and IL-4/IL-5/IL-13 in a matched aliquot and show that
engineering/activation preserves the starting program. Stop if the positive control fails.

#### Stage B: single perturbations and mechanism

Pilot six donors for variance estimation, then power an independent single-perturbation
confirmation cohort from the donor-level effect. Test only perturbations admitted by a
preregistered outcome-independent rule that is stable across epsilon-optimal coefficient
equivalence sets and held-out effect estimates; a single NNLS coefficient vector is never
a selection score.
The current dictionary cannot nominate CRISPRa interventions; any CRISPRa candidate
requires a separately measured and frozen CRISPRa dictionary/model. Otherwise CRISPRa is
an assay-control arm only.
Measure 24 h, 72 h, day 7, and after seven-day washout:

- perturbation-linked scRNA-seq;
- CITE-seq/flow for CXCR3, CCR4, CRTH2, CD25, CD127, activation, and exhaustion;
- intracellular/secreted IFN-γ, IL-2, TNF, IL-4, IL-5, IL-13, and IL-10;
- editing, proliferation, apoptosis, recovery, and cell cycle;
- targeted ATAC or multiome for the TBX21/IFNG/IL12RB2 and GATA3/Th2-cytokine loci.

Pilot go gate: point estimate ≥25% of a matched positive-control state distance, lower 95%
donor-level confidence bound above zero, two-guide concordance, coordinated RNA and
protein/cytokine movement, and no comparable off-axis movement. Confirmation requires the
lower 95% confidence bound to exceed the frozen 25% practical margin. The margin is an
operational value to freeze, not a biological constant.

#### Stage C: combinations, durability, and function

Only after independent confirmation of frozen singles, optionally freeze at most three
pairs and one triple. A strong single may be the final result without a combination stage.
Include all constituent singles and matched guide load. Determine sample size by simulation
from Stage B; roughly 12–15 donors is only an orientation for large paired effects and is
not the final power calculation.

Require improvement over the best single, directly measured interaction versus additive
prediction, multiplicity-controlled donor-level uncertainty, ≥80% viability/recovery,
and sign agreement in held-out donors. At days 7 and 14 after washout, require persistence,
within-lineage evidence using heritable lineage barcodes or clonally split/index-sorted
cells; scTCR tracking is supportive only. Add a mechanism-blocked functional assay. This
stage is necessary but not sufficient for cautious "partial reprogramming" language:
same-cell identity change, alternative-state exclusion, durable orthogonal function, and
held-out-donor replication must all succeed.

## Systemic harness

The fast harness now covers independent NNLS oracle comparison, atom scaling and
degeneracy, labelled-axis corruption, group leakage, exact uncertainty around maxT error,
and structured common-response/module/sign-selection specificity. It is intentionally
data-free and runs in pull-request CI.

The next statistical/nightly expansion must sweep block-correlation levels and add full
adaptive-pipeline replay, partial alternatives, heteroskedasticity, donor heterogeneity,
donor-specific sign reversal, guide efficiency/missingness/off-target mixtures, run
confounding, common-response strength, combination synergy/antagonism/saturation,
interval coverage, power, and abstention risk-coverage. Use at least 2,000 simulations per
scenario, while primary real-data inference uses ≥9,999 structured resamples. Require a
one-sided 95% familywise-error upper bound ≤0.075 at nominal 0.05, interval coverage
0.925–0.975, ≥0.80 power at the frozen practical margin, and ≥0.80 retained coverage so
universal abstention cannot pass.

Recent evaluation frameworks motivate perturbation-specific references, multiple metrics,
retrieval/rank tests, and strong simple baselines: [Systema](https://doi.org/10.1038/s41587-025-02777-8)
and [PerturBench](https://papers.nips.cc/paper_files/paper/2025/hash/8aee537279a66ced96319dfca3c00002-Abstract-Datasets_and_Benchmarks_Track.html).

## Expert review and handoff design

| Expert | Owns | Required handoff |
|---|---|---|
| Data curator/engineer | Dataset cards, hashes, schemas, missingness, ingestion | Immutable manifest and source-to-table run |
| Statistical geneticist | Donor/source units, uncertainty, structured nulls, multiplicity | Frozen estimand and inferential analysis plan |
| Perturbation ML scientist | Nested splits, baselines, capacity controls, metrics | Benchmark config and calibration report |
| Numerical optimization reviewer | Oracle, KKT, conditioning, equivalence sets | Certificate/failure audit |
| T-cell biologist | Construct validity, starting state, controls, off-axis states | Biological assay specification |
| Cytometry/multiome specialist | Protein/chromatin panels and QC | Blinded endpoint and batch plan |
| Reproducibility auditor | Clean-room rerun and claim-to-artifact trace | Signed release checklist |
| Independent manuscript reviewer | Adversarial review after results are frozen | Blocker/caveat/readiness rating |

No workstream approves its own claim upgrade. Data and code review precede statistical
review; statistical freezing precedes outcome access; biological and reproducibility
reviews precede manuscript wording.

## Execution order and release gates

1. **Source-bound reconstruction.** The target lineage, aggregate-screen diagnostics, and
   bidirectional source transfer, donor-pair and released guide-position sensitivities, and
   source-selected arrayed RNA/cytokine follow-up now regenerate from full-hash-bound
   inputs. The next gate is leakage-safe pseudobulk re-estimation with physical guide
   identities, keeping log fold change primary and z-score sensitivity-only. Exit only
   when a clean external-data run reproduces every claim-bearing table in its declared
   scope.
2. **Transfer and calibration.** Run leakage-safe Ota-to-Höllbacher and reverse transfer,
   followed by donor-, guide-, module-, pathway-, context-, and run-held-out evaluation.
   Add covariance-preserving nulls, nested baselines, multiple metrics, multiplicity,
   interval coverage, power, abstention, and coefficient-equivalence analysis. Exit only
   if target-specific performance clears the best frozen baseline with calibrated
   uncertainty.
3. **Independent compact benchmarks.** The aggregate Arce S1 IL2RA/CD25 transfer,
   S14/S8 donor/guide supplied-score robustness, Schmidt two-donor functional-screen
   transfer, and same-study Zhu arrayed bulk-RNA/IL-10/IL-21 follow-up are implemented
   with separate claim ceilings. Next add orthogonal independent RNA/protein strata and
   matched transcriptomic CRISPRa/i dictionaries. The Goudy negative cross-experiment
   stress is implemented and remains biologically inconclusive; next require same-experiment,
   control- and guide-burden-matched measured combinations. Exit only with
   dataset-specific grain, confounding, and claim ceilings attached to every result.
4. **Retrospective catalog calibration.** Complete training-side feature selection,
   multi-seed/threshold/library-size sensitivity, simple and regularized comparators, and
   uncertainty for the corrected catalog audit. A prospective design claim additionally
   requires an externally frozen target catalog and candidate-effect information available
   without measuring the candidate in the evaluation screen.
5. **Prospective biology.** Execute the preregistered established-memory-Th2 and
   naïve-derived-Th2 stages above with independent donors, matched modalities, direct
   combinations, fitness, durability, lineage, and functional endpoints. Only coordinated,
   durable, within-lineage molecular and functional movement can motivate cautious
   partial-reprogramming language.

## Manuscript deliverables

Required tables:

1. dataset cards and missingness;
2. estimands, split units, nulls, metrics, and multiplicity family;
3. numerical oracle and fault injection;
4. synthetic FPR, power, coverage, and abstention;
5. source/donor/guide/context transfer with intervals;
6. coefficient equivalence and combination regret;
7. source-selected arrayed profile retrieval and equal-donor cytokine effects.

Required figures:

1. leakage-safe evaluation design;
2. structured-null calibration and power;
3. cross-source/donor/context forest plot;
4. multi-metric baseline comparison;
5. fitted-point and equivalence-set stability;
6. additive prediction versus measured combinations;
7. arrayed RNA retrieval and RNA-to-cytokine concordance with donor coverage;
8. prospective RNA/protein/chromatin/function concordance if generated.

## Claim ceiling

Until prospective stages succeed, do not claim biological attainability, conversion,
reprogramming, rescue, efficacy, a minimal intervention, dose, synergy, causal regulator,
donor/population generality, predictive superiority, or independent discovery. Do not
interpret 0.444 as a converted fraction, treat the retired `1/61` as current inference,
use a coefficient as candidate importance, or read a separator as biological impossibility.

The strongest current statement remains: **a retrospective, donor-collapsed CRISPRi
effect dictionary contains model-relative directional alignment with a source-study-reused,
selectively constructed cross-sectional population contrast oriented toward the reported
Th1 centroid; the target is neither independent nor a trajectory.**
