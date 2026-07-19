# Results

This directory contains only the current canonical findings and directly supporting
artifacts.

| Path | Role |
|---|---|
| `findings.json` | Machine-readable values, interpretations, and open requirements |
| `manifest.json` | SHA-256, byte length, and executable bit for every maintained artifact |
| `validation_harness.json` | Frozen data-free systemic stress report |
| `source_reconstruction.json` | Full-file-hash-bound target lineage, frozen splits, and source-transfer report |
| `donor_pair_transfer.json` | Frozen-weight complementary donor-pair and target-source transfer sensitivity |
| `guide_pair_transfer.json` | Negative reciprocal transfer across released positional guide summaries; physical guide IDs absent |
| `goudy_combination_validation.json` | Negative GSE306915 cross-experiment CRISPRoff stress report with source crosswalk, reliability, filter sensitivity, and fail-closed statuses |
| `schmidt_external_validation.json` | Hash-gated two-fixed-donor CRISPRa/CRISPRi functional-screen concordance and complete guide/top-*K* sensitivity grid |
| `library_coverage_crossdataset.json` | Split-first Zhu/Norman/Replogle audit with strict/soft metrics, 12-partition sensitivity, nested comparators, and hash-bound portable-cache provenance |
| `evidence/arce_il2ra_context_predictions.csv` | Per-target independent cross-modality benchmark |
| `evidence/arce_activation_guide_effects.csv` | Donor×guide contrasts for the supplied Arce score |
| `evidence/arce_external_validation_meta.json` | Arce provenance, transfer, robustness, and claim ceilings |
| `evidence/zhu_arrayed_profile_metrics.csv` | Raw and panel-centered screen-to-arrayed bulk-RNA profile replication |
| `evidence/zhu_arrayed_flow_effects.csv` | Within-donor NTC-normalized IL-10/IL-21 follow-up effects |
| `evidence/zhu_arrayed_validation_meta.json` | Pinned provenance, exact panel-label diagnostics, summaries, and claim ceiling |

Legacy tables derived from an unhashed, deleted `inputs.npz` were removed rather than
presented as reproducible evidence. Their Git history remains available for provenance;
the retired summary is recorded in `source_reconstruction.json`.

`python scripts/validate_findings.py` fails closed on report status, canonical values,
paths, hashes, and executable bits.
