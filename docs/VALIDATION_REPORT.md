# Validation Report

**Audit date:** 2026-07-12  
**Scope:** core method, tests, real-data artifacts, manuscript source/PDF, public pitch,
interactive explorers, fact ledgers, and reproduction instructions.

## Overall assessment

**Ready to share as a carefully bounded method demonstration. Not ready to present as
validated cell-state control or target validation.**

The core convex optimization is sound for the question it actually answers: directional
membership/projection relative to a measured perturbation dictionary under non-negative,
additive composition. The primary numerical outputs reconcile with committed tables, the
KKT conditions are satisfied, the flagship held-out result is stable across gene splits,
and the compiled manuscript is visually clean.

The original public narrative overstated several conclusions. It conflated
perturbation–condition profiles with unique knockdowns, described a greedy panel as globally
minimal, treated a full-vector dual certificate as per-gene activation proof, described the
GOF proxy incorrectly, and interpreted small permutation samples too strongly. Those issues
have been corrected in the manuscript, README, submission pitch, demo script, explorer copy,
fact ledgers, and method documentation.

## What was independently verified

| Check | Evidence | Result |
|---|---|---|
| Tier-2 grain | Direct HDF5 inspection of `GWCD4i.DE_stats.h5ad` | 33,983 rows × 10,282 genes; 11,526 distinct targeted genes; Rest/8h/48h row counts 11,287/11,415/11,281 |
| Working dictionaries | `analysis_cache/atlas_work/inputs.npz` | 6,871/7,155/7,195 generators by condition |
| Target direction | Notebook construction and marker signs | `toward_Th1 = Th1 − Th2`; GATA3 down, TBX21/IFNG/STAT4/STAT1 up |
| Headline geometry | `table1_verdict.csv`, `atlas_reachability.csv` | in-sample 0.627; held-out 0.448; LOF 0.393; residual 0.779; KKT ≈ 1.1e−11 |
| Dense null | `table5_null_summary.csv` | observed 0.448 exceeds all 60 shuffles; plus-one empirical p = 1/61; z = 24.44 is descriptive |
| Split stability | New 12-seed rerun | mean held-out 0.4464, SD 0.0101, range 0.4326–0.4587 |
| Signed proxy | `atlas_reachability.csv`, code inspection | staged 39% LOF / 25% sign-flipped GOF proxy / 35% neither; shares sum to one within rounding |
| Atlas aggregate | Recomputed from CSV | 12 cases; mean LOF 0.3377; maximum LOF 0.4132; Rest has highest held-out cosine on all four axes |
| Norman additivity | Recomputed from `norman_table4_additivity.csv` | 126 doubles; median measured-vs-additive cosine 0.712 |
| Candidate triage | `modality_intervention_map.csv` | 102 unique greedy LOF candidates; 45 hard-to-drug; 10 clinical-drug tier; 9 genetics+tractability highlights |
| Literature matrix | `method_comparison_matrix.csv` | 92 rows including this work; 91 prior; 14 prior rows rated fully measured; no prior row rated fully measured + full feasibility verdict |
| Core software | `pytest tests/test_reachability.py -vv -s` | 11/11 passed in 154.89 s on the final audit run |
| Manuscript render | all 25 pages rendered and inspected | no clipping, overlap, broken figures, or unreadable glyphs |
| Technical dossier render | all 101 pages rendered to a contact sheet; first pages inspected at full resolution | no clipping, overlap, missing figures, or broken glyphs observed |

## Material issues corrected

### 1. Target sign was reversed in the Methods prose — high severity

The code correctly constructs `toward_Th1` by sign-flipping a Th2-vs-Th1 contrast, but the
Methods described the target as source minus destination. The manuscript now states
`destination − source`, matching the code and verified marker directions.

### 2. The matrix grain was mislabeled — high severity

“33,983 knockdowns” implied 33,983 unique gene interventions. The H5AD rows are
perturbation–condition profiles across 11,526 distinct targeted genes. The manuscript now
reports both the source matrix and the quality-filtered per-condition generator counts.

### 3. Global minimality was not established — high severity

NNLS solves the best cone projection; the spectrum uses greedy forward selection. Neither
solves the combinatorial minimum-cardinality problem. “Minimal recipe” is retained only in
legacy filenames; prose now says “ranked mixture” and “greedy sparse candidate panel.”

### 4. The certificate was biologically overinterpreted — high severity

The complete residual is a valid separating direction. A positive coordinate means the
closest cone projection under-delivers that readout; it does not prove that no other cone
point can raise the individual gene, or that activating the gene causes the target state.
The biological output is now framed as a CRISPRa/de-repression hypothesis ranking.

### 5. The GOF decomposition was described incorrectly — high severity

The Methods said the residual was projected onto the positive-expression orthant and invoked
mutually polar cones. The implementation instead performs a staged projection onto a cone
generated by `-E`. The manuscript now describes that sign-symmetry proxy and its limitation:
real CRISPRa effects need not mirror CRISPRi.

### 6. Permutation z-scores were overclaimed — high severity

A z-score estimated from 60 shuffles is not a literal 24-sigma tail probability, and an
eight-shuffle atlas cannot provide a conventional empirical tail p-value. The headline is now
reported as “above all 60 shuffles, plus-one empirical p = 1/61,” with z retained only as a
descriptive standardized separation. Atlas z values are labelled screening estimates.

### 7. State-level decisions and candidate-level annotations were conflated — medium severity

GO/REDIRECT/STOP applies to a target direction and dictionary. The 102-item table is instead
the union of top-10 greedy LOF candidate panels. The revised narrative keeps those levels
separate and no longer assigns each candidate an individual reachability verdict.

### 8. The reproduction command was misdescribed — medium severity

`reproduce.sh` verifies software on synthetic fixtures; it does not regenerate the real-data
headline from a clean clone because the 16.8 GB matrix is untracked. README and script output
now state this explicitly and provide the separate real-analysis route.

### 9. The public hook implied more clinical leverage than the evidence supports — medium severity

Clinical attrition remains useful context, but this transcriptomic test does not explain or
predict clinical failure. The new pitch leads with the route-map metaphor and promises only
an early, screen-relative modality check.

## Calculation and methodology notes

- For Euclidean projection onto a closed convex cone, fitted and residual components are
  orthogonal, so the KKT/separating result is mathematically appropriate.
- The KKT residual certifies optimizer convergence and the model-relative projection. It does
  not validate additivity, causal transport, intervention magnitude, or phenotype.
- The cone is unbounded. Because the principal metric is cosine, the claim is about direction,
  not an achievable biological dose or complete state transition.
- The staged LOF/GOF/neither shares form a nested Pythagorean partition. LOF and GOF fitted
  vectors need not be pairwise orthogonal.
- The one-sided 39/31/30 and staged 39/25/35 decompositions are distinct and remain correctly
  separated.

## External factual checks

- The CZI dataset card reports ~22 million cells, four donors, and three stimulation
  conditions: [Primary Human CD4+ T Cell Perturb-seq](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq?access_dataset=true).
- The simple-baseline motivation matches the published paper:
  [Ahlmann-Eltze, Huber & Anders, Nature Methods (2025)](https://www.nature.com/articles/s41592-025-02772-6).
- The 2.6× genetics figure matches
  [Minikel et al., Nature (2024)](https://www.nature.com/articles/s41586-024-07316-0).
- The historical 30.7% Phase-II transition rate matches the
  [BIO/Biomedtracker/Amplion 2006–2015 report](https://go.bio.org/rs/490-EHZ-999/images/Clinical%20Development%20Success%20Rates%202006-2015%20-%20BIO%2C%20Biomedtracker%2C%20Amplion%202016.pdf).

## Remaining caveats and required next experiments

1. **No direct certificate validation.** The strongest distinctive biological claim still
   needs a CRISPRa/de-repression experiment guided by the unmet-readout ranking.
2. **No in-domain combination calibration.** Additivity is calibrated on K562 CRISPRa
   doubles, not primary CD4⁺ CRISPRi combinations.
3. **No leave-one-donor-out effect matrix.** Published effects are donor-collapsed.
4. **Transcriptome is a proxy endpoint.** Cytokines, chromatin, viability, and in-vivo
   function remain untested.
5. **Atlas null depth is shallow.** Eight shuffles per case are adequate for screening but
   not formal multiple-testing control.
6. **The novelty survey is structured, not exhaustive proof.** “No prior method” must remain
   qualified as “none found in the 91-method survey.”
7. **Open Targets is a snapshot.** Tractability and drug counts can change and require an
   explicit refresh date before external publication.
8. **A2 conditional analysis is scaffolded, not run.** It requires raw single-cell counts.

## Sharing recommendation

Lead with this sentence:

> Given a measured perturbation dictionary, we test whether non-negative combinations can
> point toward a target transcriptomic direction—and return the full separating residual
> when they cannot.

Avoid “impossible,” “smallest recipe,” “must activate,” “clinically validated,” and bare
atlas z-scores. Use “outside this measured cone under the additive model,” “greedy sparse
panel,” “positive unmet readout,” “wet-lab hypothesis,” and “eight-shuffle screening
estimate.”
