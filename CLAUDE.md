# CLAUDE.md — operating manual for this repo

*Orientation for an autonomous coding agent (or a new human contributor) working in
`cell-state-reachability`. Read this before touching code. It records the enduring facts,
guardrails, and literature anchors.*

*This file deliberately does **not** repeat the repo map, the reproduce commands, or the doc
index — those live once, in [`README.md`](./README.md). For **what the method found**, read
[`docs/RESULTS.md`](./docs/RESULTS.md) first.*

---

## 0. What this project is (one paragraph)

Genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells gives a *measured* causal
effect vector for each gene knockdown. This project treats those vectors as a dictionary
and asks a geometric question: **is a target transcriptomic state (e.g. Th2→Th1
polarization) reachable by some non-negative combination of knockdowns, and if so, what
is the smallest ranked set that reaches it?** Because knockdowns are loss-of-function and
weights are non-negative, the reachable shifts form a **convex cone**; the headline
output is a falsifiable **reachable / provably-outside-the-cone** verdict with a residual
and — when outside — a constructive activation certificate, each nomination carrying a
decomposable confidence score. This is *not* differential expression and *not* a
similarity ranking — it is a reachability decision.

---

## 1. Current state

The method is implemented in **`reachability.py`** (~1,200 LOC: cone fit, signed
LOF/GOF decomposition, Farkas/KKT certificate, honest nulls, reachability spectrum,
additivity risk). The 16.8 GB Tier-2 effect matrix is local, and the full analysis has run.

Everything runs CPU-only on an 18 GB laptop; compute is not a constraint (full solver +
1,000-iteration null at 34k×2k is ~0.57 s, 1.4 GB RAM).

Commands and the notebook order are in [`README.md`](./README.md#reproduce).

**Public API of `reachability.py`** (verdict layer): `reachability(E, d)` →
reachable_cosine / residual_norm / verdict / non-negative weights;
`signed_reachability(...)` → LOF/GOF/neither split; `reachability_spectrum(...)` (with
`epistasis_penalty`); `additivity_risk(...)`; the shuffled-target and
random-perturbation nulls; `analytic_anisotropy_null(...)` (closed-form). Reconstruction
helpers (`greedy_minimal_set`, `omp_minimal_set`), target-state construction, confidence
decomposition, and the PubMed/Open Targets/Enrichr evidence layer round out the module.

---

## 2. Verified quantitative facts (trust these; reproduced from the data)

| fact | value | note |
|---|---|---|
| Local CSVs | 7 files, 37.7 MB | 6/7 row-exact; sgRNA table is **26,504** rows (not 31,110) |
| Sign convention | `toward_Th1` → GATA3 ↓, TBX21/IFNG/STAT4/STAT1 ↑, IL4/IL5/IL13 ↓ | both source contrasts |
| Cross-source concordance | 11,616 shared genes · 68.5% sign-concordant · ρ=0.562 (z-score) | ρ=0.533 on log-fc |
| GATA3 rank (`toward_Th1`) | #1 Rest, #2 Stim8hr, #2 Stim48hr | positive control, NOT a discovery |
| Headline Th2→Th1 verdict | held-out cosine 0.448 (in-sample 0.627), null z ≈ 24, KKT residual 1.1×10⁻¹¹ | partially reachable |
| Signed split (Th2→Th1) | LOF 39% / GOF 25% / neither 35% | knockdown never the majority modality (atlas mean LOF 0.34) |
| Solver speed @ 34k×2k | 0.57 s greedy + 1,000-null; 1.4 GB RAM | compute is not a constraint |
| Autoimmune enrichment | 17 diseases, 185 hits FDR<0.05, top OR 58 (Crohn's) | disease-impact layer |
| Conditions / donors | Rest, Stim8hr, Stim48hr / 4 donors (age 22–34) | in DE_stats `culture_condition` |

Data grain: `DE_stats` is per **perturbation × condition** (33,983 rows × 10,282 genes),
key columns `target_contrast_gene_name`, `culture_condition`, `ontarget_effect_size`,
`ontarget_significant`, `offtarget_flag`, `crossdonor_correlation_mean/min`,
`crossguide_correlation`. The aging join uses `gene_name`; polarization uses `variable`.

> **Known inconsistency to resolve.** The signed split above (39/25/35) matches
> `docs/figures/fig_central_illustration.png`, but `docs/figures/fig2_decomposition_certificate.png`
> shows 39/31/30 for what reads as the same quantity. Establish which is right before
> either number reaches the manuscript.

---

## 3. GUARDRAILS (non-negotiable — this project's credibility depends on them)

- **Knockdown-only.** CRISPRi is loss-of-function. Never invent activation effects; keep
  weights non-negative (`knockdown_only=True`). If a target needs a gene *up*, the honest
  answer is "not reachable by knockdown; here is what a CRISPRa arm would test" — the
  certificate, not a fabricated recipe.
- **Additivity is a calibrated assumption, not a free one.** Multi-gene sets assume no
  epistasis; the `additivity_risk()` score (calibrated on Norman doubles — a magnitude-
  saturation law, not collinearity) bounds the extrapolation. Flag every multi-gene
  nomination as an extrapolation to be tested.
- **Nominations are hypotheses, not validated targets.** Always framed for wet-lab test.
- **Claim novelty for the METHOD, not the regulator lists.** The source paper (Zhu et al.)
  already reports Th1/Th2 + aging regulators; recovering them = validation, not discovery.
  Novelty = the measured-effect convex-cone reachability verdict + certificate +
  held-out-gene validation.
- **Nulls before claims.** No reachability/alignment number is reported without its null
  and CI. Report the held-out number as the headline, never the in-sample cosine.
  *This applies to documentation too: a result with no committed figure or table is not a
  reported result. See the §8.3 status table in `docs/RESULTS.md`.*
- **No invented judging rubric.** The event page publishes no scored criteria; do not
  optimize for or cite an imagined rubric.
- **Reproducibility hygiene.** Fixed seed everywhere; `fig, ax = plt.subplots()` then
  `fig.savefig()` (never bare `plt.savefig()`); fetches in their own step; pin the env.
- **Do not commit the raw data** (`data/*` gitignored; never commit the Tier-2 h5ad).

---

## 4. Literature anchors (verified — cite these exactly)

- **Baseline-first justification:** Ahlmann-Eltze, Huber & Anders, *Nat Methods*
  22:1657–1661 (2025), doi:10.1038/s41592-025-02772-6. Five foundation + two other DL
  models vs deliberately simple baselines; none beat them. For double perturbations, all
  models had prediction error *substantially higher* than the additive baseline (Fig. 1a,b).
  (Verified against full text.) This is the field-level result `NOVELTY.md` uses to motivate
  the inverse-feasibility reframe, not merely to license the linear model.
- **Prior art — Mogrify:** Rackham et al., *Nat Genet* 48:331–335 (2016). Minimal TF set →
  cell conversion (gain-of-function + network). *Delta:* measured LOF vectors.
- **Prior art — CellOracle:** Kamimoto et al., *Nature* 614:742–751 (2023). Inferred-GRN
  perturbation → cell-identity transition vector + randomized null. *Delta:* measured
  effect matrix + convex-cone reachability verdict.
- **Source dataset:** Zhu et al., bioRxiv, doi:10.64898/2025.12.23.696273.
- **Cross-dataset generalization:** Norman et al. 2019 (K562 CRISPRa doubles, additivity
  calibration); Replogle et al. 2022, *Cell*, doi:10.1016/j.cell.2022.05.013, PMID 35688146
  (K562/RPE1 essential-gene screens, cross-cell-type transfer).
- **External validation arms:** Shifrut et al., *Cell* (2018, PMID 30449619); Schmidt et al.,
  *Science* (2022, PMID 35113687) — CRISPRa arm tests up-regulation hypotheses CRISPRi can't.

*The full 91-method survey with every DOI/PMID is in `docs/RELATED_WORK.md` and `results/references.csv`.*
