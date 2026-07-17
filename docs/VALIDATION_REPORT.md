# Validation Report

**Scope:** consolidated maintained surface, 2026-07-17  
**Canonical facts:** [`results/findings.json`](../results/findings.json)  
**Artifact identity:** [`results/manifest.json`](../results/manifest.json)

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

`scripts/run_arce_external_validation.py` verifies Arce archive/member hashes and schema,
enforces outcome-independent eligibility, extracts the Zhu `Rest` IL2RA predictor, and
evaluates 480 perturbations in three CRISPR-KO CD25 screen contexts. Tiny workbook/H5
fixtures test failure modes, guide eligibility, selection isolation, orientation,
determinism, and check-mode drift. This is modest cross-study/cross-modality ranking
alignment, not donor-level or whole-state validation.

## Canonical evidence status

| Question | Evidence | Status |
|---|---|---|
| Source-bound random-gene alignment? | 0.444 ± 0.018 across 12 hash-frozen splits | Descriptive; correlated splits |
| Target-source directional transfer? | Positive cosine gain over mean/best-single in 6/6 splits both ways | Directional only; nRMSE not improved |
| Independent functional-screen transfer? | Arce Spearman 0.148 / 0.084 / 0.088 | Modest, context dependent |
| Software/statistical contracts? | All six systemic scenarios pass | Synthetic certification only |
| Donor/guide generality? | Donor-collapsed primary case study; aggregate Arce S1 | Not tested |
| Functional state conversion? | No established-state prospective assay | Not tested |

## Artifact consistency

`scripts/validate_findings.py` cross-checks the ledger against both external reports and
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
donor/source/guide, stronger-baseline, structured-null, measured-outcome, and prospective
gates in the [scientific validation plan](SCIENTIFIC_VALIDATION_PLAN.md) pass.
