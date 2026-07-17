# Validation Report

**Scope:** consolidated maintained surface, 2026-07-17  
**Canonical facts:** [`results/findings.json`](../results/findings.json)  
**Artifact identity:** [`results/manifest.json`](../results/manifest.json)

## What is validated

### Numerical core

The test suite covers:

- exact inside, boundary, and outside projections in an orthant;
- the metric-correct weighted separator (`W rho`), polarity, and orthogonality;
- zero-weight coordinate exclusion;
- target scaling from `1e-6` to `1e6`;
- positive atom rescaling through `1e16`, duplicate atoms, and zero atoms;
- frozen held-out coefficients;
- conservative plus-one empirical p-values with ties;
- fail-closed handling of zero, non-finite, misaligned, fractional/boolean-index,
  non-boolean-mask, and invalid-weight inputs;
- absence of removed legacy decision APIs.

The maintained implementation emits geometry and diagnostics only. It does not emit a
biological reachability verdict.

### Frozen findings

`scripts/validate_findings.py` recomputes the 12-split mean, SD, range, and values; checks
the fixed-split null count, exceedances, and plus-one p-value; cross-checks every scalar in
the findings ledger against named evidence; verifies every declared source path; and
validates SHA-256, byte length, and executable status for every canonical artifact.

### Figure

`docs/figures/make_at_a_glance.py` reads only `results/findings.json`. It uses a headless
backend and fixed PDF metadata epoch. Consecutive renders must produce identical PNG and
PDF hashes before the manifest is updated.

## Updated evidence status

| Question | Current evidence | Status |
|---|---|---|
| Stable to registered random-gene splits? | 0.446 ± 0.010 across 12 fixed splits | Supported for these splits |
| Above the 60 diagnostic target shuffles? | 60/60 below historical fixed split; p=1/61 | Descriptive only |
| Better than simple baselines? | Better than mean/PCA-1/PCA-5/random-cone mean; worse than PCA-20/LS/ridge | Mixed; no superiority claim |
| Driven by cell-cycle genes? | Ablation changes held-out cosine by about -0.002 | Not load-bearing in this check |
| Driven only by generator magnitude? | Exceeds registered norm-matched null | Not explained by this check |
| Generator significance filter absent? | All 6,871 Rest generators source-flagged significant | Concern addressed for admission |
| Robust combination semantics? | Retired-v0 diagnostics: 1/126 threshold flips, 100/126 staged proxy flips | Removed heuristics are unstable |
| Target fully observed? | 9,831/25,672 genes; 38/50 strongest DE genes | No; measured subset only |
| Donor-general? | Donor-collapsed inputs | Not tested |
| Functional conversion/rescue? | No prospective functional assay | Not tested |

## Known limitations

1. Random-gene partitions can leak correlated modules.
2. Sixty shuffles are insufficient for a strong inferential tail claim.
3. Z-score profiles are not calibrated physical doses or additive intervention units.
4. Donor-collapsed effects do not permit donor-held-out evaluation.
5. The target and perturbation dictionary come from different studies and contexts; the
   sign-concordance and screen-intersection filters select a restricted subset.
6. Direct combinations, matched activation effects, chromatin, and functional endpoints
   are absent.
7. Selected supporting tables are frozen evidence; most require large external inputs to
   regenerate.

## Reproduce

```bash
python -m pip install -r requirements.txt
./reproduce.sh
```

Exit code zero requires tests, demo, findings validation, and artifact hashes to pass.
