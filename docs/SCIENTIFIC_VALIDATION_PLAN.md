# Scientific Validation Plan

**Review status (2026-07-17): Needs revision before inferential or biological
reachability claims.** The current repository supports a retrospective, donor-collapsed
directional-geometry result. It does not yet support state conversion, intervention
utility, donor generalization, target-specific significance, or coefficient-based
prioritization.

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
| Effect layers | log FC, z-score, p, adjusted p, base mean, and LFC SE | Log-FC primary analysis with uncertainty is feasible |
| Target input | 6,188 final coordinates after source-sharing, sign, and screen filters | Selection leakage and excluded-coordinate scope must be explicit |
| Missing inference unit | Current result has no donor-level sampling distribution | Random-gene split SD must never be reported as biological uncertainty |

The smallest remediation is to ingest the accessible donor-, guide-, pseudobulk-, and
target-construction Zhu files, while separately resolving the documented validation table.

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
- `Th1Th2_validation_summary.suppl_table.csv` only after a working source route is verified;
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

Current status: the aggregate H5AD/target phase is implemented with full input hashes;
the unreconstructable archived 0.446 multisplit result was withdrawn and replaced by the
source-bound 0.444 result. Donor-, guide-, and pseudobulk-level reconstruction and the
same-question baseline remain exit blockers.

### Workstream 2 — independent transfer and structured inference

**Lead experts:** statistical geneticist, perturbation-ML scientist, and immunologist.

Actions:

1. **Source transfer:** fit/tune with Ota, then score the frozen direction against
   Höllbacher; reverse the roles. Define the shared screen-measurable universe without
   using held-out source signs or magnitudes.
2. **Donor-pair transfer:** fit on one donor pair and score on the disjoint pair across
   every complementary 2-vs-2 partition. Report donors 1–2 vs 3–4 as fully
   run-confounded, and donors 1–3 vs 2–4 and 1–4 vs 2–3 as run-balanced partitions.
3. **Guide transfer:** define the target universe without aggregate/held-out-guide
   effectiveness, restrict reciprocal testing to targets with two usable guides, select/fit
   using one guide and evaluate another, then reverse, and report excluded single-guide
   targets.
4. **Grouped holdouts:** use pathway/module, source, donor, guide, context, and run blocks.
5. **Structured nulls:** preserve gene covariance using rotations or module-level
   randomization. Use ≥9,999 resamples for a primary family and maxT/Westfall–Young control
   across contexts, target directions, metrics, and preprocessing variants.
6. **Nested baselines:** tune only inside training folds. Include best single atom,
   common-effect mean, PCA/ridge, unconstrained least squares, sparse/non-negative ridge,
   and capacity-matched random dictionaries.
7. **Metrics:** cosine, normalized RMSE, norm ratio, DE-sign recovery, top-DE precision and
   recall, pathway-direction agreement, target retrieval rank, and improvement over the
   strongest frozen baseline.
8. **Identifiability:** report active-matrix rank, singular values, strict-complementarity
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

### Workstream 3 — documented measured follow-up (retrieval blocked)

**Lead experts:** T-cell biologist and assay statistician; analysis must be frozen before
opening outcome columns.

First resolve retrieval of Zhu's documented but currently unavailable
`Th1Th2_validation_summary.suppl_table.csv`. The VCP schema describes arrayed CRISPRi
under non-polarized, Th1-polarized, and Th2-polarized culture with bulk-RNA state scores
and flow-protein outcomes, but no working public object key has been verified.

Actions:

1. Freeze the candidate universe, score, sign, and exclusions before reading validation
   outcomes.
2. Evaluate continuous Th1 and Th2 score direction, rank, and error without refitting.
3. Evaluate IFN-γ protein direction/rank and culture-context reversals.
4. Report that the arrayed targets were selected after the source discovery screen; this
   is source-selected measured follow-up, not held-out or prospective validation.

Exit gate: source-selected candidates show RNA-score and IFN-γ flow concordance in the
preregistered direction across donors. RNA can challenge stress, apoptosis, and cell-cycle
signatures if gene-level data are accessible; viability selection remains untested unless
an additional direct readout is located.

### Workstream 4 — combination and modality stress tests

**Lead experts:** perturbation biologist, causal/ML evaluator, and data curator.

1. Use [Goudy et al. GSE306915](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306915)
   to compare the measured four-donor `FAS + SUV39H1 + RC3H1` triple against its matched
   singles. Report additive prediction error, interaction residual, norm ratio, and module
   error. One triple cannot establish general selection regret. Use
   [GSE306917](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306917) only as a
   two-donor CD55/CD81 CRISPRoff mechanism check, not mechanism evidence for the triple.
2. Use [Schmidt et al.](https://doi.org/10.1126/science.abj4008) primary-T-cell CRISPRa
   resources (GSE190604/GSE174255/GSE190846) to test gain-of-function orientation and
   cytokine agreement. Never treat CRISPRa as scale-matched negative CRISPRi.
3. Use [Arce et al.](https://doi.org/10.1038/s41586-024-08314-y) GSE271090/GSE278572 for
   context sign reversals and RNA/protein agreement. The first S1 IL2RA/CD25 ranking
   benchmark is implemented; donor/guide-grain S14 and orthogonal S20/S4 endpoints remain.
4. Retain Norman/Replogle only as out-of-domain numerical and combination benchmarks.

Exit gate: reported conclusions name the assay, cell system, modality, donor unit, and
claim ceiling. No cancer-cell-line result is used to validate primary-CD4 biology.

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
run confounding, combination synergy/antagonism/saturation, interval coverage, power, and
abstention. Use at least
2,000 simulations per scenario, while primary real-data inference uses ≥9,999 structured
resamples.

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
   bidirectional source transfer now regenerate from full-hash-bound H5AD/CSV inputs.
   Extend the same contract to donor-, guide-, and pseudobulk-level objects, keeping log
   fold change primary and z-score sensitivity-only. Exit only when a clean external-data
   run reproduces every claim-bearing table in its declared scope.
2. **Transfer and calibration.** Run leakage-safe Ota-to-Höllbacher and reverse transfer,
   followed by donor-, guide-, module-, pathway-, context-, and run-held-out evaluation.
   Add covariance-preserving nulls, nested baselines, multiple metrics, multiplicity,
   interval coverage, power, abstention, and coefficient-equivalence analysis. Exit only
   if target-specific performance clears the best frozen baseline with calibrated
   uncertainty.
3. **Independent compact benchmarks.** The aggregate Arce S1 IL2RA/CD25 transfer is now
   implemented. Next add its donor/guide and RNA/protein strata, then Schmidt CRISPRa/i
   orientation and cytokine transfer. Treat the Goudy triple as descriptive until sample
   mapping and reuse terms are resolved. Exit only with dataset-specific grain,
   confounding, and claim ceilings attached to every result.
4. **Prospective biology.** Execute the preregistered established-memory-Th2 and
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
6. coefficient equivalence and combination regret.

Required figures:

1. leakage-safe evaluation design;
2. structured-null calibration and power;
3. cross-source/donor/context forest plot;
4. multi-metric baseline comparison;
5. fitted-point and equivalence-set stability;
6. additive prediction versus measured combinations;
7. prospective RNA/protein/chromatin/function concordance if generated.

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
