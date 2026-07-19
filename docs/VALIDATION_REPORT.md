# Validation Report

**Scope:** consolidated maintained surface, 2026-07-19
**Canonical facts:** [`results/findings.json`](../results/findings.json)  
**Artifact identity:** [`results/manifest.json`](../results/manifest.json)

The numerical contract certifies the declared point-estimate optimization problem; the
adversarial harness demonstrates known ways evaluation can fail. Neither establishes that a
particular real-data score is artifact-free or that exact cone membership is a calibrated
biological decision rule. The source-reconstructed and independently evaluated sections
state their narrower empirical ceilings below.

## What is technically certified

The numerical tests cover exact inside/boundary/outside projections, weighted separator
polarity and orthogonality, scaling extremes, zero/duplicate/near-duplicate atoms, frozen
held-out coefficients, conservative empirical p-values, and fail-closed malformed inputs.
KKT and separator diagnostics certify at `1e-8`; no biological verdict threshold exists.

The six-scenario systemic harness adds an independent active-set oracle, eight injected
axis/provenance faults, scale/degeneracy challenges, grouped splits, exact uncertainty
around maxT familywise error, and structured specificity. It demonstrates that common
response can yield high raw cosine, random-gene splits can be optimistic relative to
module holdout, and sign selection can inflate scores. The maxT fast check has 24/500
false families and an exact one-sided 95% upper bound of 0.067 under the 0.075 gate. This
is synthetic contract certification, not biological validation.

## What is source reconstructed

`scripts/run_source_reconstruction.py` verifies the full bytes of the 16.8 GB Zhu H5AD
and target CSV, validates grain/schema, rebuilds target lineage, hashes every split's
fit/score genes, and performs bidirectional Ota/Höllbacher transfer on log fold change and
z-score. Any hash or frozen-split mismatch makes the report `FAIL`; the canonical validator
rejects failed reports.

The audit retired the previous 0.446 ± 0.010 table and p=1/61 because their deleted
intermediate did not preserve gene order. The source-bound result is 0.444 ± 0.018. A
separate fixed-split value reproduces within `3e-10`, supporting source alignment without
reconstructing the retired pipeline.

## What is independently evaluated

`scripts/run_donor_pair_transfer.py` verifies the 16.9 GB donor H5MU and target-table
bytes, validates all six two-donor modalities, freezes gene splits, and applies NNLS
weights plus training-selected baselines unchanged to the complementary donor pair and
opposite target source. Tests cover H5 schema corruption, axis alignment, training-only
baseline selection, run-balance classification, and deterministic output. The 24
run-balanced challenges show weak directional gain but worse magnitude error. Released
presence is DE-eligibility-selected, so this certifies a published-pipeline fixed-cohort
sensitivity—not leakage-free donor generality.

`scripts/run_guide_pair_transfer.py` verifies all 29,424,424,894 bytes of the guide H5MU
(SHA-256 `964eeafb3356a7322a1d5b1121802c6a1433456f3591e2d5797817df3bf9c2f6`) and the
target table before checking the complete MuData/AnnData schema. It requires `guide_2` to
be an exact subset of `guide_1`, records rather than imputes the 114/38 rows with missing
categorical metadata, and verifies that every selected common-`Rest` row has an exact
target mapping. Focused fixtures cover hash/schema drift, target-mapping integrity,
training-only baseline selection, frozen coefficient reuse, summary reconstruction, and
byte-exact check mode.

The runner retains all 8,323 common `Rest` atoms, reports 2,752 `guide_1`-only exclusions,
and performs 12 fits scored in 24 rows. One additional reverse-atom-order fit has maximum
relative L2 difference `1.76e-15`. The frozen report SHA-256 is
`0fd291002f55c7d6a505d3d19c15f53747bc4cf5b6473540a06978de53b4e721`. This certifies
the released positional-summary computation only. Physical guide IDs are absent,
`guide_2` is a selected subset, the author rule includes `keep_effective_guides`, the
intermediate `for_DE_by_guide.csv` is unavailable, and the upstream DE formula omits a
donor term. The negative transfer therefore does not certify failure between physical
guides or leakage-safe guide generalization.

`scripts/run_arce_external_validation.py` verifies Arce archive/member hashes and schema,
enforces outcome-independent eligibility, extracts the Zhu `Rest` IL2RA predictor, and
evaluates 480 perturbations in three CRISPR-KO CD25 screen contexts. Tiny workbook/H5
fixtures test failure modes, guide eligibility, selection isolation, orientation,
determinism, and check-mode drift. The same runner streams S14, verifies all 520
donor/guide/context strata, and reproduces 116 S8 aggregates to <5.7×10⁻¹³. Unequal-size
fixtures prove that controls are weighted by guide rather than pooled cell count. This is
modest transfer plus descriptive supplied-score robustness, not donor-population or
whole-state validation.

`scripts/run_schmidt_external_validation.py` verifies the complete 26.15 MB Zenodo archive
before ZIP parsing, then enforces an exact member allow-list, author-script identity, four
sgRNA-member hashes, schemas, row counts, finite numeric values, paired donor suffixes, and
unique guide keys. Nine focused tests cover deterministic fixture execution and check mode,
the pre-parse hash gate, schema/nonfinite/duplicate/suffix corruption, outcome-independent
eligibility, and source-only top-*K* selection. The canonical ≥3-guide universe contains
18,568 complete genes; all 288 rows in the 3-threshold × 4-top-*K* × 3-transfer-class ×
8-direction grid are retained and cross-bound to the ledger. The frozen report SHA-256 is
`87745d31bcd08bf70d2a6c16287db9abe5fa836d04ed41fda632d7cc09da8810`.

This certifies the computation, not a broad biological claim. Whole-universe same-reagent
donor Spearman is 0.135–0.332 and modality-plus-library Spearman is 0.020–0.036. The
top-200 signed-rank medians (0.887 same-screen donor, 0.300 donor-plus-modality/library,
0.749 donor-plus-context) are conditional on training-donor extreme-effect selection and are shown
with absolute-rank, sign, global-overlap, and full sensitivity diagnostics. Because guides
are reused, libraries change with modality, cytokine changes with CD4/CD8 context, and only
two fixed donors exist, no guide-held-out, context-general, population, efficacy, state, or
target-recommendation claim is certified.

`scripts/run_zhu_arrayed_validation.py` verifies the aggregate screen and two compact
author-table hashes, masks all nine panel target genes from every profile, and evaluates source-to-arrayed
bulk-RNA identity both raw and after panel centering. It normalizes IL-10/IL-21 flow within
donor against all available NTC measurements and enumerates all 362,880 target-label
permutations with synchronized maxT for four RNA-to-flow rank associations. Tiny fixtures
cover hash drift, outcome-independent screen selection, duplicate bulk keys, donor-control
normalization, all-panel-target masking, retrieval, and exhaustive enumeration determinism. The panel
was source-selected and donor coverage is unbalanced, so this certifies cross-platform
follow-up only.

## What is certified in the Goudy cross-experiment stress test

`scripts/run_goudy_combination_validation.py` verifies the normalized-count, family-SOFT,
and pinned author-key bytes before parsing. It crosswalks the declared roles across all
three sources, enforces donor- and role-specific controls, selects genes from controls only,
masks retained target coordinates, and keeps the on-target fit disjoint from transcriptome
scoring. The report exposes the four author-key-unresolved donor-3/4 multiplex roles rather
than filling them by assumption, and it records the single/triple experiment, control-type,
and guide-burden confounding.

Fifteen focused tests cover pre-parse hash gates, count/SOFT/author-key schema failures,
the exact crosswalk omissions, role-specific controls, target masking, training-only LODO
filtering and fitting, saturated-calibration labeling, metadata-conflict fail-closure, and
byte-exact check mode. The frozen report reproduces with SHA-256
`f73651b7fd0eee376eda2238466cbde527b6dd28dbb27e919615291752f33933`.
An independent provenance and mathematical re-audit found no blocker- or high-severity
finding in this frozen artifact.

Execution is technically reproducible, but that does not rescue the empirical model:
component-cone cosine is 0.0949 in sample and 0.0881 LODO, donor-effect reproducibility is
low, and the result changes materially with the control-only expression threshold. The
certified conclusion is a negative, confounded cross-experiment stress result—not evidence
for a general additive or prospective model.

## What is certified in the retrospective catalog extension

`library_coverage.py` adds catalog bookkeeping over the same projection core. Its data-free
contract tests distinguish strict membership from a soft cosine bar, detect a planted
duplicate, use the certified weighted dual separator, keep certificate order separate from
the optional realized comparator, and enforce the nested signed-span ≥ cone ≥ ray capacity
ordering. Realized gains must be non-negative, and candidate labels must be one-to-one with
canonical perturbation labels.

When the exact registered local caches are present,
`scripts/run_library_coverage_crossdataset.py --check` verifies their identities before
reproducing the Zhu, Norman, and Replogle retrospective audits. A separate builder now
reconstructs all three pickle-free cache artifacts byte-for-byte from hash-matching raw
sources; synthetic fixtures exercise the same contracts in CI. The v3 report proves
split-first feature selection, Norman canonical-gene role disjointness, a fixed Norman
double catalog, and 12 deterministic sensitivity partitions. It certifies deterministic
software behavior on already-observed effect dictionaries—not prospective selection of an
unmeasured perturbation, experimental outcome, sampling uncertainty, or biological
reachability. Durable public availability of the mutable upstream bytes remains open.

The corrected frozen audit retains its negative top-pick result: certificate and realized
top candidates differ in all three audits. Exact-row and canonical-label agreement are
mathematically identical because each candidate pool is required to contain one unique row
per canonical label. The validated numerical summary is maintained once in
[Findings](FINDINGS.md#10-extension-retrospective-effect-dictionary-coverage-is-distinct-from-prospective-design).

## Canonical evidence status

| Question | Evidence | Status |
|---|---|---|
| Source-bound random-gene alignment? | 0.444 ± 0.018 across 12 hash-frozen splits | Descriptive; correlated splits |
| Target-source directional transfer? | Positive cosine gain over mean/best-single in 6/6 splits both ways | Directional only; nRMSE not improved |
| Published-eligibility donor-pair transfer? | Run-balanced median cosine gain +0.032; nRMSE 1.153 vs 1.018 | Weak direction; magnitude fails; four fixed donors |
| Released guide-position transfer? | Same-source held-gene cosine 0.251227 within position and −0.019197 reciprocal; paired gain over best single positive in 3/12, nRMSE gain in 0/12 | Negative positional-summary stress; physical guide IDs absent; effectiveness-selected universe |
| Independent functional-screen transfer? | Arce Spearman 0.148 / 0.084 / 0.088 | Modest, context dependent |
| Arce selected-panel concordance? | Within-dataset A-vs-B target-rank concordance 0.73–0.93; four-stratum sign agreement 50–64% | Descriptive; preselected 28 regulators, two donors, supplied score |
| Schmidt functional-screen concordance? | Whole-genome donor 0.135–0.332; top-200 signed-rank medians 0.887 donor / 0.300 donor+modality/library / 0.749 donor+context | Conditional two-fixed-donor stress; joint shifts, same reagents, no population or guide generality |
| Zhu arrayed transcriptome replication? | Matching target retrieves 9/9; median panel-centered cosine 0.580 | Source-selected same-study follow-up; nRMSE 1.052 |
| Zhu cytokine consistency? | Screen/bulk RNA versus donor-median flow Spearman 0.717–0.850 | Six follow-up donor labels with unequal target coverage; conditional panel diagnostic |
| Software/statistical contracts? | All six systemic scenarios pass | Synthetic certification only |
| Goudy cross-experiment triple stress? | Component cone 0.0949 (strict 0/4); LODO 0.0881 vs training-selected best single 0.0957; triple donor cosine 0.0480 | Declared geometric model fails; confounded and low reliability |
| Retrospective dictionary coverage? | 12 split-first partitions; strict-positive splits 0/12 for each dataset; cone beats simple rays in all 36 | Same-screen compression/additivity audit; not prospective design or sampling uncertainty |
| Norman single-to-double alignment? | Full-single reference 0.958; 40-gene representative sensitivity median 0.865 [0.807, 0.878] | Retrospective additivity diagnostic; constituent baselines already strong |
| Leakage-safe primary-model donor/physical-guide generality? | Donor-pair and positional-guide released-object sensitivities only | Not tested; pseudobulk re-estimation is next |
| Functional state conversion? | No established-state prospective assay | Not tested |

## Artifact consistency

`scripts/validate_findings.py` cross-checks the ledger against all external reports and
the systemic harness, requires report `PASS`, requires full source hashes, checks split-ID
hashes, and validates SHA-256/bytes/executable bits for every canonical file. The central
figure reads only the ledger and is rendered twice to confirm deterministic PNG/PDF output.

## Reproduce

```bash
python -m pip install -r requirements.txt
./reproduce.sh
```

External regeneration additionally requires the gitignored registered inputs and
`requirements-external.txt`. Scientific readiness remains **Needs revision** until the
seven open gates in the [scientific validation plan](SCIENTIFIC_VALIDATION_PLAN.md) pass:
leakage-safe pseudobulk donor/physical-guide evaluation; module/pathway holdouts and
structured nulls; nested baselines; matched in-domain combinations; paired CRISPRi/a
dictionaries; independent whole-state endpoints and durability; and prospective
established-state validation.
