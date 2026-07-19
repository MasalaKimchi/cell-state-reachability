# External data

Scientific source data stay local and are ignored by Git. The repository ships a compact
frozen evidence bundle and a data-free validation harness, not the large source matrices.

## Primary source: Zhu et al.

Genome-scale CRISPRi Perturb-seq in primary human CD4 cells from four donors and three
assay conditions:

- preprint: https://doi.org/10.64898/2025.12.23.696273
- GEO/SRA: GSE314342 / SRP643211
- source code: https://github.com/emdann/GWT_perturbseq_analysis_2025
- public data prefix: `s3://genome-scale-tcell-perturb-seq/marson2025_data/`
- VCP release: [v1.0 dataset page](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq?access_dataset=true)
- licensing: VCP labels the dataset MIT and applies its Acceptable Use Policy; the
  analysis code is separately MIT; mirror terms must be recorded independently

List the public objects without credentials:

```bash
aws s3 ls --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/
```

Highest-priority objects:

| Object | Bytes | Use |
|---|---:|---|
| `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` | 6,155,771 | Exact target-construction input; fetched from the author repository |
| `IL10IL21bulkRNAseq_DESeq2_results.csv` | 13,952,871 | Arrayed bulk-RNA follow-up for nine source-selected perturbations |
| `IL10_IL21_arrayed_validation.csv` | 2,200 | IL-10/IL-21 flow percentages for Donor5–Donor10 and matched NTC controls |
| `GWCD4i.DE_stats.by_donors.h5mu` | 16,866,278,447 | Disjoint donor-pair transfer |
| `GWCD4i.DE_stats.by_guide.h5mu` | 29,424,424,894 | Negative reciprocal transfer across released guide-rank summaries |
| `GWCD4i.pseudobulk_merged.h5ad` | 44,566,657,140 | Controlled re-aggregation and alternate DE models |
| `GWCD4i.DE_stats.h5ad` | 16,786,240,107 | 33,983 perturbation-condition profiles × 10,282 genes; log FC, z-score, p, adjusted p, base mean, and LFC SE |

Frozen Zhu source identities used by the maintained configs and runners:

| Object | SHA-256 |
|---|---|
| `GWCD4i.DE_stats.h5ad` | `c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62` |
| `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` | `c47d2df21414ca85e7aa255f4148904eec700fbcd9debc2f734ec97049698444` |
| `IL10IL21bulkRNAseq_DESeq2_results.csv` | `c20418a9285b10104dbae362b825971f86f97425800a92269e4433ce780e666d` |
| `IL10_IL21_arrayed_validation.csv` | `f60cdda392d6f29d10a539727ff7324b04d17e35c0512c889b733e00380b83dc` |
| `GWCD4i.DE_stats.by_donors.h5mu` | `2ee3cf90925600eb044619021da2bdd47d661f306a204586652256facf17af64` |
| `GWCD4i.DE_stats.by_guide.h5mu` | `964eeafb3356a7322a1d5b1121802c6a1433456f3591e2d5797817df3bf9c2f6` |

The target CSV is pinned to author commit
[`848d62f`](https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/848d62fc2b7027f7218d6fc5f5b0c37255dc94af),
not the mutable default branch.

The VCP schema also documents `Th1Th2_validation_summary.suppl_table.csv`, but no working
public S3 or author-repository object was verified on 2026-07-17. Treat retrieval as
unresolved; do not present this table as an executable benchmark until a source route is
confirmed.

Download one allow-listed object with:

```bash
./data/fetch_de_stats.sh Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv
./data/fetch_de_stats.sh IL10IL21bulkRNAseq_DESeq2_results.csv
./data/fetch_de_stats.sh IL10_IL21_arrayed_validation.csv
./data/fetch_de_stats.sh GWCD4i.DE_stats.by_donors.h5mu
./data/fetch_de_stats.sh GWCD4i.DE_stats.by_guide.h5mu
```

The default remains `GWCD4i.DE_stats.h5ad`. HTTPS downloads resume from a `.part` path,
verify registered byte length, and then rename atomically. Record SHA-256 and retrieval
date in the dataset card before analysis.

With both frozen files present, run the source-bound reconstruction separately from the
small CI suite:

```bash
python -m pip install -r requirements-external.txt
python scripts/run_source_reconstruction.py --profile
python scripts/run_source_reconstruction.py --check results/source_reconstruction.json
```

This reconstructs target lineage and aggregate-screen source-transfer diagnostics. It
does not substitute for donor- or guide-level objects.

With the two compact follow-up tables present, the source-selected arrayed benchmark adds
8,967 transcript coordinates per perturbation after masking all nine panel target genes,
plus donor-normalized IL-10 and
IL-21 flow readouts:

```bash
python scripts/run_zhu_arrayed_validation.py --check
```

The nine targets were chosen upstream using the same source study. Donor coverage ranges
from three to six, and target-label permutations are conditional on that panel. This is
same-study cross-platform follow-up in six additional donor labels—not held-out discovery,
donor-population inference, or state conversion.

The donor-pair H5MU stage is implemented:

```bash
python scripts/run_donor_pair_transfer.py --check
```

It fits on one two-donor modality and target source, then freezes coefficients and
training-selected baselines for the complementary donor pair, opposite target source,
and held-out genes. Across 24 run-balanced correlated challenges, median cosine gain over
the frozen best single is +0.032, while normalized RMSE is worse (1.153 versus 1.018).
This is a published-pipeline sensitivity, not predictive utility. The next stage is an
exact guide-ID crosswalk plus structural-QC-only re-estimation from pseudobulk counts.
Released H5MU modalities omit
DE-ineligible targets, so presence itself can select on effectiveness. Only pseudobulk can
rebuild a leakage-safe universe without `keep_effective_guides`, held-out DE, sign, or
correlation fields. Four donors support fixed-cohort robustness, not donor-population
significance.

The released guide-position H5MU stage is also implemented:

```bash
python scripts/run_guide_pair_transfer.py --check
```

The registered object is
`s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.by_guide.h5mu`
(29,424,424,894 bytes; SHA-256
`964eeafb3356a7322a1d5b1121802c6a1433456f3591e2d5797817df3bf9c2f6`). The public
object metadata record S3 version ID `SQHf_ZhmdbCteM9f4HG1zs.2k25CL3gb`, last modified
`2026-05-28T23:20:10Z`, and multipart ETag
`"2e6705636ebaa276c7bc7c5a148ad096-3508"`. The fetcher uses that ETag as an `If-Match`
version guard; it is not treated as a content hash. Anonymous version-addressed retrieval
was unavailable, so the full SHA-256 is the canonical content identity. The downloaded
H5MU has no embedded license field, so the VCP MIT metadata and Acceptable Use Policy
remain the recorded data terms.

The [official VCP v1.0 dataset card](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq?access_dataset=true)
(accessed 2026-07-19) defines `guide_1` and `guide_2` as the first and second sgRNA IDs
after alphanumeric sorting within target-condition; public pseudobulk and guide-library
artifacts carry sgRNA IDs. The ranked H5MU modalities themselves do not embed sgRNA ID or
sequence, and this repository has not reconstructed and hash-cross-verified their exact
rank-to-sgRNA crosswalk. `guide_2` is an exact subset of `guide_1` in the H5MU.

Reciprocal testing keeps all 8,323 common category-labeled `Rest` atoms, reports the 2,752
`guide_1`-only category-labeled `Rest` atoms as excluded, and uses 8,950
identifier-selected target genes. Another 35 nominal `guide_1` `Rest` keys are withheld
because their required categorical metadata are missing. Across the full object, the 114
`guide_1` categorical-missing rows have key suffixes `Rest` 35 / `Stim8hr` 18 /
`Stim48hr` 61; the 38 `guide_2` rows are `Rest` 0 / `Stim8hr` 0 / `Stim48hr` 38. These
rows are recorded rather than imputed, and every selected category-labeled common-`Rest`
row has an intact target mapping.

The pinned author pipeline first selected
`keep_min_cells & keep_effective_guides & keep_total_counts`, requiring at least three
passing replicates per guide-condition and at least five cells per guide in each
condition/sample. It built `cond_targets` from the unreleased `for_DE_by_guide.csv`, kept
exactly two testable guides per target-condition, applied the final `keep_test_genes`
filter, and fit `~ log10_n_cells + target` without a donor term. Accordingly, the observed
drop from median same-rank held-gene cosine 0.251227 to reciprocal-rank cosine −0.019197
does not demonstrate robustness in this fixed descriptive reciprocal-rank benchmark. It
is neither named-sgRNA replication nor leakage-safe physical-guide generalization. Exact
rank-to-ID crosswalk plus pseudobulk re-estimation using only structural coverage/depth QC
is the next gate.

## Target sources

- Ota et al., *Cell* 2021, DOI 10.1016/j.cell.2021.03.056, NBDC E-GEAD-397.
- Höllbacher et al., *ImmunoHorizons* 2020, DOI
  10.4049/immunohorizons.2000037, GEO GSE149090.

The current target uses the Zhu supplementary polarization table. Its registered merged
analysis has 6,188 genes after requiring both sources, concordant signs, and screen
measurement. Source-transfer validation must return to the independent Ota and
Höllbacher inputs and must not select coordinates using the held-out source.

## Compact external benchmarks

These archives are small enough to exercise the external-data harness before downloading
the 17–45 GB donor/guide/pseudobulk objects:

| Dataset | Exact archive | Bytes | Terms and allowed use |
|---|---|---:|---|
| Arce Perturb-CITE, Zenodo 13924126 / GSE278572 | [`data_tables.zip`](https://zenodo.org/api/records/13924126/files/data_tables.zip/content), MD5 `886ed0fea0b9dc0625355c2e4928077c`, SHA-256 `dc9e2efb04d24f1a6d4b8db6a8b1d5cd01c935777c3740088be339de5b5062b4` | 57,967,623 | **Implemented**: S1 CRISPRi-transcript → CRISPR-KO CD25 transfer and S14/S8 donor/guide supplied-score robustness; not state or donor-population validation |
| Schmidt screens, Zenodo 5784651 / GSE174292 | [`Genome-wide-screens.zip`](https://zenodo.org/api/records/5784651/files/Genome-wide-screens.zip/content), MD5 `e0392eb7b2512720bb8cbf705ce9854f`, SHA-256 `15571c41d76b2462d15f167f8920b0ec335f685b1582d18b0264f65f21b2fefd` | 26,152,593 | **Implemented**: two-fixed-donor CRISPRa/CRISPRi functional-screen concordance; not paired transcriptomic dictionaries or donor/guide generality |

The Arce runner verifies the archive plus S1/S8/S14 workbook bytes, requires four guides in every
screen context, freezes Zhu `Rest` admission before reading Arce outcomes, and checks the
current deterministic evidence:

```bash
python scripts/run_arce_external_validation.py --check
```

The retained S1 benchmark has 1,347 four-guide genes, 1,259 present in the Zhu `Rest`
dictionary, and 480 source-admitted analysis targets. S1 is aggregate and cannot supply
donor uncertainty. Arce S9 is significant-only (absence is censored, never zero), and S4
contains unlabelled technical replicates that must be aggregated within biological key.
Schmidt's cell-level marker tables are likewise not donor-level inference.

S14 contributes 100,087 singlet cells, 520 complete guide×donor×context strata, and four
contexts. S8's 116 pooled summaries are exactly reproduced but are not independent data;
its pooled-cell tests are excluded. The supplied `activation.score` lacks a frozen local
formula/gene set and is used only for descriptive within-object robustness.

The Schmidt runner reads the four MAGeCK sgRNA summaries directly from the exact Zenodo
archive after verifying its full bytes, SHA-256, upstream MD5, complete member allow-list,
used-member hashes, 15-column schema, row counts, paired donor suffixes, and finite values.
The deposited author script fixes `r0`/`r1` as Donor1/Donor2, aggregates each gene/donor by
median guide LFC, and reverses CRISPRi orientation. Paper Methods state that CRISPRa and
CRISPRi used the same donors. Eligibility uses only common gene identity and guide coverage;
FDR, hit calls, signs, and magnitudes do not enter the universe.

```bash
./data/fetch_de_stats.sh Genome-wide-screens.zip
python scripts/run_schmidt_external_validation.py --check
```

At the primary ≥3-guide contract, 18,568 genes are complete across all four screens and
both donors. Whole-genome same-reagent donor Spearman is 0.135–0.332, and CRISPRa-versus-
oriented-CRISPRi agreement is 0.020–0.036. Conditional on the training donor's 200 largest
absolute effects, median signed Spearman is 0.887 in the same screen/held donor, 0.300 for
the joint donor-plus-modality/library transfer, and 0.749 for the joint donor-plus-
cytokine/cell-type context transfer. The latter two do not isolate a single axis. The
1/3/6-guide × top-50/100/200/500 grid is explicitly exploratory and post-hoc. Because the
same guides are reused within modality, Calabrese/Dolcetto changes with modality, IL2/CD4
changes jointly with IFNG/CD8, and there are only two donors, this is not guide-held-out,
modality-equivalence, cell-type/cytokine-generalization, or donor-population evidence.

## Goudy GSE306915 negative cross-experiment stress

The implemented Goudy benchmark is a **STRESS**-tier bulk-RNA CRISPRoff diagnostic in
primary human bulk T cells, sourced from [Goudy et al.](https://doi.org/10.1038/s41587-025-02856-w)
and GEO GSE306915. It asks only whether one measured `FAS + RC3H1 + SUV39H1` triple is
represented by target-matched constituent singles under the deposited design. All three
registered inputs remain local and gitignored:

| Registered input | Bytes | SHA-256 | Source |
|---|---:|---|---|
| `GSE306915_normalized_counts_CO065.csv.gz` | 9,358,021 | `02307f1019429530fae91d8da3d808a1c8e04241fe4657832205c94d01f43d42` | [GEO supplementary normalized counts](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306915/suppl/GSE306915_normalized_counts_CO065.csv.gz) |
| `GSE306915_family.soft.gz` | 10,412 | `9059377ff91eee08ba71b52c787d4166baa0b2e29a9a3b02ba29566c63bbe5c4` | [GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306915/soft/GSE306915_family.soft.gz) |
| `rna_seq_meta_key.csv` | 3,671 | `ba27f8502a517dab3c25c7c8001e85e303659e42d84e96c0e48b20d51fbe3e2f` | [Author key at commit `53155c9`](https://raw.githubusercontent.com/GilbertLabUCSF/T_Cell_CRISPRoff/53155c9207d8b4f70f0ae5d60e1f4c0513d41bd7/Cas9_KO_vs_CRISPRoff_KD/rna_seq_meta_key.csv) |

The author key is pinned to commit
`53155c9207d8b4f70f0ae5d60e1f4c0513d41bd7`. The paper is CC-BY-4.0; GEO does not
specify a data license for the deposited files. Raw GEO bytes and the pinned author key
are not redistributed by this repository. Fetch and verify each registered object, then
check the frozen report:

```bash
./data/fetch_de_stats.sh GSE306915_normalized_counts_CO065.csv.gz
./data/fetch_de_stats.sh GSE306915_family.soft.gz
./data/fetch_de_stats.sh rna_seq_meta_key.csv
python scripts/run_goudy_combination_validation.py --check
```

The pinned key has 71 rows, 67 of which overlap the count matrix, and confirms 36/40
declared analysis roles. Constituent singles and AAVS1 controls are author-key-confirmed
in Co065; D1/D2 multiplex NTC and triple samples are confirmed in Co066. The D3/D4
multiplex NTC and triple rows are absent from that key, so their experiment IDs remain
unresolved, and the triple guide identity is unresolved for all donors. Same-guide
matching is therefore not claimed. GEO labels eight D1 single/control characteristics as
Day 18 while titles, the paper, and pinned key indicate Day 7; the runner uses Day 7 and
makes no duration claim. Single-versus-triple status is inseparable from experiment,
control type, and guide burden.

Execution is `PASS`, but geometry is `FAILS_DECLARED_GEOMETRIC_MODEL` and biological
interpretation is `INCONCLUSIVE_CROSS_EXPERIMENT_CONFOUNDING_LOW_RELIABILITY`. At the
canonical control-only filter, the component cone has median cosine 0.0949, median nRMSE
0.9952, and strict membership 0/4; equal-sum cosine/nRMSE is 0.0646/2.3084, LODO cosine is
0.0881, and triple pairwise-donor cosine is 0.0480. The declared filter sweep changes the
point estimate but does not rescue the model. This is a bounded negative stress result,
not evidence for a general model or a causal interpretation of the residual.

## Retrospective effect-dictionary caches

The cross-dataset coverage audit uses three local, gitignored NPZ caches. Their exact byte
lengths, SHA-256 identities, array schemas, source links, split designs, and claim ceiling
are registered in
[`configs/library_coverage_crossdataset.json`](../configs/library_coverage_crossdataset.json).
The frozen report refuses a cache whose registered identity or schema differs.

| Audit | Cell system / modality | Registered cache | Scope |
|---|---|---|---|
| Zhu | Primary human CD4+ T cells / CRISPRi | `slice/zhu_rest_effects_portable_v1.npz` | Same-screen dictionary compression |
| Norman | K562 / CRISPRa | `slice/norman_effects_portable_v1.npz` | Measured single-to-double additivity/alignment |
| Replogle | K562 / CRISPRi | `slice/replogle_effects_portable_v1.npz` | Same-screen dictionary compression |

All three caches are pickle-free and byte-stable. Rebuild them from hash-matching source
bytes with:

```bash
./data/fetch_de_stats.sh norman_perturbation.h5ad
./data/fetch_de_stats.sh replogle_k562_essential_perturbation.h5ad
python scripts/build_library_coverage_caches.py --dataset all
python scripts/build_library_coverage_caches.py --dataset all --check
```

Zhu additionally uses the H5AD and target table registered in
`configs/source_reconstruction.json`. The builder verifies every raw byte identity before
loading, uses upstream processed expression as-is for Norman and Replogle, and writes only
`E`, `perts`, unique `genes`, and cell counts where available. Replogle's display symbol
`TBCE` occurs twice, so its unique upstream `gene_id` is the transport axis. Norman's paper
identifies K562 whereas the redistribution's `obs.cell_type` says A549; the benchmark flags
that discrepancy and never uses the conflicting field computationally.

This makes the computation reconstructable from matching raw bytes, but not durably
obtainable: the Norman/Replogle public URLs are mutable and historical S3 version IDs are
not anonymously retrievable. A release archive or DOI for the approximately 275 MB
canonical cache bundle remains the manuscript-grade availability step, subject to source
redistribution terms. Every candidate effect is already measured; the audit does not rank
unmeasured genes or validate prospective experiments.

## Next dataset priorities

- Paired primary-T-cell CRISPRa/CRISPRi transcriptomic effect dictionaries; the implemented
  Schmidt scalar functional screens do not provide matched whole-state dictionaries.
- Same-experiment primary-T-cell measured combinations with verified guide identities,
  all constituent singles, and control- and guide-burden-matched designs.
- Zhu donor/guide objects above, after compact benchmark calibration is stable.

See [`docs/SCIENTIFIC_VALIDATION_PLAN.md`](../docs/SCIENTIFIC_VALIDATION_PLAN.md) for
the exact evaluation and claim ceiling for each resource.

## Data policy

- Never force-add H5AD, H5MU, CSV, NPZ, or raw-count payloads under `data/`.
- Every analysis input needs an accession, exact object key, version/retrieval date,
  byte length, SHA-256, license/terms status, gene namespace, orientation, units,
  donor/guide/context/batch fields, and missingness profile.
- A missing donor or guide estimate remains missing; it is never silently replaced by a
  collapsed average.
- Claim-bearing tables must be regenerated by the external-data runner, not copied from
  an interactive notebook.
