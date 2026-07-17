# Results

This directory contains only the current canonical findings and directly supporting
artifacts.

| Path | Role |
|---|---|
| `findings.json` | Machine-readable values, interpretations, and open requirements |
| `manifest.json` | SHA-256, byte length, and executable bit for every maintained artifact |
| `validation_harness.json` | Frozen data-free systemic stress report |
| `source_reconstruction.json` | Full-file-hash-bound target lineage, frozen splits, and source-transfer report |
| `evidence/arce_il2ra_context_predictions.csv` | Per-target independent cross-modality benchmark |
| `evidence/arce_external_validation_meta.json` | Arce provenance, attrition, metrics, nulls, and claim ceiling |

Legacy tables derived from an unhashed, deleted `inputs.npz` were removed rather than
presented as reproducible evidence. Their Git history remains available for provenance;
the retired summary is recorded in `source_reconstruction.json`.

`python scripts/validate_findings.py` fails closed on report status, canonical values,
paths, hashes, and executable bits.
