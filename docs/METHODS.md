# Method

## Question and geometry

Given measured perturbation profiles and a target transcriptional direction in a shared
gene space, how closely can a non-negative linear combination of those profiles align with
the target? The method tests one dictionary, target, and metric; it is narrower than
perturbation-response prediction.

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
returned. No biological verdict threshold exists.

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
[`data/README.md`](../data/README.md).
