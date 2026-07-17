# Results

This directory is the complete maintained result bundle.

| Path | Role |
|---|---|
| `findings.json` | Canonical machine-readable values, interpretations, and open requirements |
| `manifest.json` | SHA-256 and byte length for every maintained artifact |
| `evidence/` | Selected tables that directly support the updated findings |

The public narrative is [`docs/FINDINGS.md`](../docs/FINDINGS.md). Values displayed in
the README and central figure must first appear in `findings.json`.

## Evidence groups

### Headline and baseline

- `headline_heldout_split_stability.csv`
- `historical_fixed_split_null.csv`
- `baseline_comparison.csv`
- `metric_calibration.csv`
- `metric_calibration_provenance.json`

### Context and source construction

- `context_condition_comparison.csv`
- `context_regulator_strength.csv`
- `context_runbalance_caveat.json`
- `reviewer2_ota_hollbacher_meta.json`
- `reviewer2_ota_hollbacher_split.csv`
- `reviewer2_marker_coverage.csv`
- `reviewer2_stim48_confound.csv`

### Specific robustness challenges

- `confounder_robustness.csv`
- `confounder_robustness_summary.json`
- `generator_significance.csv`
- `generator_significance_summary.json`
- `ranking_validation.csv`
- `ranking_validation_summary.json`
- `reviewer2_deg_survival.csv`
- `reviewer2_recipe_specificity.csv`
- `reviewer2_recipe_specificity_meta.json`
- `positive_control_enrichment.csv`
- `positive_control_stats.csv`

### Combination and transport diagnostics

- `additivity_verdict_flip.csv`
- `additivity_verdict_flip_summary.json`
- `norman_table4_additivity.csv`
- `generality_second_dataset_summary.csv`

## Policy

- Evidence tables may preserve historical column names such as `verdict` for provenance;
  those labels are not public API or current claims.
- “Positive control” in historical evidence names a synthetic reachable-by-construction
  calibration reference, not a biological control.
- No recipe, druggability, causal-oracle, clinical, or activation-certificate artifacts are
  maintained.
- Large data and intermediate caches are external to the repository.
- `python scripts/validate_findings.py` fails if values, paths, or hashes drift.
