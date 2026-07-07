# Roadmap v2 — Cell-State Reachability (7 days, solo, CPU-only)

This is the working roadmap. It supersedes v1: same dataset and honest-scope stance,
but the method is reframed from "sparse reconstruction scored by cosine" to a
**reachability-of-cell-states** treatment, a **cross-dataset external-validation**
layer was added, and — per the no-auth mandate — the **graded core now runs entirely
on the open, no-login Tier-1 CSVs** (no CZI/Synapse/Wiley account required anywhere on
the critical path).

> **Repo state: planning / brainstorming only.** All runnable code (`src/`, `app/`,
> `tests/`, build config) has been removed; this repo now holds the design docs plus
> the local Tier-1 CSVs. The "Done" items below are **validated findings from the
> earlier prototype** — preserved as ground truth to rebuild from, not live code.

## Status checklist

**Done (validated in the prototype before code was stripped)**
- [x] Schema gate (`data_loader.validate_local_data`) — all 7 CSVs pass
- [x] Fixed target sign (`toward_Th1` = −Th2_vs_Th1); validated: GATA3 down, TBX21/IFNG up
- [x] Cross-source concordance computed (Ota vs Höllbacher: 68.5% sign-concordant, ρ=0.56)
- [x] Fixed confidence QC columns (real DE_stats names + off-target penalty fires)
- [x] Confidence `coverage` field (flags scores built from thin evidence)
- [x] **No-auth Tier-1 directional nominations** (`tier1_directional_nominations`) — GATA3 ranks #3/2135
- [x] Consolidated `pathways.py` → `evidence.py`; Open Targets now returns real disease associations
- [x] External-evidence layer confirmed 100% no-auth (PubMed / Open Targets / Enrichr)
- [x] Folder renamed → `cell-state-reachability`; core smoke tests pass on CPU

**Pending (needs the user or is optional)**
- [ ] Commit the code-removal + doc updates and push (run locally; sandbox can't push)
- [ ] Rebuild the pipeline from this plan when implementation resumes
- [ ] Decide whether to track the 7 Tier-1 CSVs (currently gitignored despite README)
- [ ] *(optional)* Tier-2 gene-level matrix for full reachability — CZI login; not on critical path
- [ ] Reachability cone + spectrum, held-out-gene eval (Tier-2, once/if matrix is fetched)
- [ ] Streamlit cards wired to live no-auth evidence

## 0. Feasibility verdict — is this possible on a CPU-only MacBook?

**Yes, comfortably. Compute was never the binding constraint here.**

| Workload | Size | Laptop reality |
|---|---|---|
| Tier-1 CSVs (`data/`) | ~36 MB | Seconds. |
| Tier-2 `GWCD4i.DE_stats.h5ad`, one layer float32, subset to sig-on-target perts + ~2k HVGs | ~34k×2k ≈ 0.25 GB (full layer ~1.4 GB) | Fits in RAM on a 16 GB Mac; solvers are seconds, bootstraps are minutes. |
| Raw 22M-cell matrix | 100s of GB | Out of scope — never touched. |

The real constraints are **statistical and logistical**, not FLOPs: only 4 donors,
single-gene CRISPRi (no combinatorial ground truth), and one external download as a
single point of failure. The way to push the frontier on this hardware is therefore a
**better-posed problem, not a bigger model** — a stance the 2025 *Nature Methods*
linear-baseline result explicitly licenses.

## 0a. Framing decision (why we pivoted)

The featured Marson/Pritchard dataset is a **CD4+ T-cell** genome-scale CRISPRi
Perturb-seq screen, not AML/HSPC data. We keep the *counterfactual method* and apply
it to a transition the dataset can support **and independently sanity-check**:

- **Primary axis — Th1 ↔ Th2 polarization.** The dataset ships a Th2-vs-Th1 target
  signature from **two independent source contrasts** (Ota 2021 and Höllbacher 2021),
  letting us cross-check regulator nominations for consistency across signature sources.
- **Secondary axis — CD4+ T-cell aging (demoted to appendix).** The local donors are
  all young (ages 22–34), so there is no in-sample aged reference; this axis is
  exploratory only and no longer part of the graded core.

## 0b. The reframe: reachability, not just reconstruction

CRISPRi yields **loss-of-function** effect vectors combined with **non-negative**
weights, so the set of transcriptome shifts achievable by knockdowns is a **convex
cone**. The scientific question sharpens to:

> *Is the target state transition inside — or how far outside — the cone spanned by the
> available knockdown effects?*

This is a small non-negative least-squares / LP, milliseconds on CPU, and it yields
what a bare cosine cannot: a falsifiable verdict — **reachable** (here is the closest
reachable point and its residual) vs **provably unreachable by any knockdown
combination** (because the move requires gene activation). *Unreachable is a real
result, not a failure.*

**Headline deliverable — the reachability spectrum.** Best achievable alignment vs
sparsity *k*, plotted for the true target against **shuffled-target** and
**random-perturbation** nulls. One figure that (a) reframes the output from fragile
per-gene hits to the geometry of what knockdowns can and cannot do to a CD4+ T cell,
and (b) neutralizes the overfitting critique because the null band shows the
achievable-by-chance floor at every *k*.

## 0c. What data we actually have locally (ground truth for this plan)

Seven derived supplementary CSV tables (~36 MB) drive everything below. Verified
shapes/columns (these are now enforced by a **schema gate**, see §1c):

| File | Rows | Grain | Key fields we use |
|---|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | perturbation × condition | `ontarget_effect_size/significant`, `offtarget_flag`, `crossdonor_correlation_mean/min`, `crossguide_correlation`, DE-gene counts |
| `Th2_Th1_polarization_signature…csv` | 37,288 | gene × source-contrast | `log_fc`, `zscore`; gene col is **`variable`**; two contrasts stacked (`contrast`) |
| `CD4T_aging_signature…csv` | 10,000 | gene | `log_fc`, `zscore`; gene col is **`gene_name`** |
| `guide_kd_efficiency.suppl_table.csv` | 73,765 | guide × condition | `signif_knockdown`, `rank`, t-stat |
| `sgrna_library_metadata.suppl_table.csv` | 31,110 | guide | TSS distance, off-target / nearby non-target genes |
| `cluster_autoimmune_enrichment…csv` | 5,236 | cluster × disease × gene-set | `odds_ratio`, `p_adj_fdr`, `intersecting_genes`, 17 autoimmune diseases |
| `sample_metadata.suppl_table.csv` | 12 | sample | 4 donors × 3 conditions; age 22–34, sex, ethnicity |

Three culture conditions everywhere: **Rest, Stim8hr, Stim48hr** — every effect and
target vector is condition-specific.

**Honesty correction (unchanged from v1):** `DE_stats` is a per-perturbation *summary*
table, **not** the gene-level logFC matrix. The reachability/reconstruction solver
needs the gene-level matrix in `GWCD4i.DE_stats.h5ad` (fetched via `vcp-cli`).

## 0d. Two tiers of analysis (no-auth first)

**Tier 1 — CSV-only, NO AUTH, THE GRADED CORE.** Signatures, per-perturbation
QC/effect summaries, guide QC, off-target design, disease enrichment, donor metadata —
**plus directional knockdown nominations** via `tier1_directional_nominations`. This is
an honest 1-D proxy (not full gene-space reachability): it scores each knockdown by
whether the target wants that gene's transcript lower, weighted by on-target effect and
cross-donor reproducibility. It demonstrates the counterfactual *idea* end-to-end with
zero login and zero download. Sanity-checked: for `toward_Th1`, the Th2 master
regulator **GATA3 ranks #3 of 2,135**.

**Tier 2 — adds `GWCD4i.DE_stats.h5ad` (OPTIONAL).** Upgrades the 1-D proxy to the full
gene-space **reachability cone + spectrum** and held-out-gene evaluation. Requires a
(free) CZI login and a non-trivial download, so it is explicitly **off the critical
path**: a graded, defensible submission ships on Tier 1 alone. Note from a real run —
`vcp data search --exact` returns only the 12 raw per-donor datasets; the derived
matrix is not a `--file`-selectable hit, so Tier 2 is "nice to have," not required.

## 1. Success criteria (falsifiable, not vibes)

A submission is "done" when:

1. **Cross-source robustness (Tier 1).** The polarization target is built from the
   **sign-concordant core** of Ota vs Höllbacher, and the concordance is *reported*
   (currently: 11,616 shared genes, **68.5% sign-concordant**, Spearman 0.56). Nominated
   regulators rank consistently across the two sources.
2. **Reachability with honest nulls (Tier 2).** The reachability spectrum sits above
   BOTH a shuffled-target null and a random-perturbation null, with effect size +
   bootstrap CI, and holds under leave-one-donor-out (n=4, reported as robustness).
3. **Held-out-gene generalization (Tier 2).** Weights fit on a random half of the target
   signature genes; alignment is scored on the held-out half. This is the primary
   validity claim — it needs no extra donors and is not inflatable by overfitting.
4. **External corroboration (Tier 1/2).** Top nominations are cross-checked against
   ≥1 independent public dataset (see §5b) and a known-biology sanity gate passes
   (e.g. `toward_Th1` moves GATA3 down and TBX21/IFNG up — verified).
5. Every hypothesis card renders a decomposable confidence score, an
   off-target/knockdown-quality flag, disease linkage, and citations.
6. The whole pipeline reproduces from a fixed seed + pinned environment on CPU in
   under ~15 min (Tier 1 in seconds).

## 1b. Analysis catalog (grouped by tier; ★ = high value-for-effort)

### A. Target-state & signature analysis — *Tier 1*
1. ★ Build both target vectors in z-score space, with **correct sign conventions**
   (`toward_Th1` = negative of the Th2_vs_Th1 table; `toward_young` = negative of
   aged-vs-young). *[fixed — was inverted in v1]*
2. ★ **Cross-source concordance** — pivot polarization by `contrast`, keep the
   sign-concordant core, report agreement. *[fixed — v1 silently kept one contrast]*
3. Signature anatomy; overlap between aging and polarization signatures.
4. Condition-resolved targets across Rest / Stim8hr / Stim48hr.

### B. Reachability & directionality — *Tier 2 (gene-level)*
5. ★ **Reachability cone** — NNLS/greedy fit of the target inside the non-negative
   knockdown cone; report reachable-vs-unreachable + residual.
6. ★ **Reachability spectrum** — alignment vs k, with shuffled-target + random nulls.
7. **Additivity-aware selection** — penalize co-selecting genes with collinear effect
   vectors (shared downstream program → predicted sub-additivity). Turns the epistasis
   caveat into an explicit model feature.
8. Effect-magnitude landscape / hub vs narrow perturbations (Tier 1 summary).
9. Condition-dependent regulators (sharp Rest↔Stim change).

### C. Confidence, reproducibility & QC — *Tier 1*
10. ★ Per-nomination confidence from the **correct** DE_stats columns
    (`crossdonor_correlation_mean`, `crossguide_correlation`, `ontarget_significant`,
    `offtarget_flag`) + per-guide `signif_knockdown`. *[fixed — v1 used wrong names →
    silent 0.5 defaults + off-target penalty never fired]*
11. Knockdown-quality gate; off-target audit from `sgrna_library_metadata`.
12. Bootstrap stability selection; leave-one-donor-out (n=4).

### D. Disease-relevance layer — *Tier 1*
13. ★ Autoimmune enrichment linkage (17 diseases) → interpretable "why this matters".
14. Disease-specific target shortlists; negative-control check.

### E. External validation — *Tier 1/2, NEW (see §5b)*
15. ★ **Independent-screen overlap** — do our nominations recur as T-cell state
    regulators in Shifrut 2018 / Schmidt 2022?
16. ★ **Independent-expression check** — do nominated regulators show the expected
    Th1-vs-Th2 direction in DICE sorted-cell RNA-seq?
17. **Cross-modal coherence** — does a `toward_Th1` set de-enrich Th2-driven disease
    (e.g. asthma) genes in the autoimmune enrichment table?

### F. Evidence, interpretation, packaging — *Tier 1*
18. ★ Literature/Open Targets/Consensus evidence per gene (support, not proof).
19. Pathway / gene-set enrichment.
20. ★ Hypothesis cards + Streamlit explorer; hard-coded limitations banner.
21. **Experiment-design queue** — rank the single most informative *next* wet-lab
    experiment (pair whose additive-vs-epistatic prediction is most decision-relevant),
    reframing "combinations untestable" as the product itself.

## 1c. Engineering hardening applied (bug fixes shipped)

- **Schema gate** (`data_loader.validate_local_data`): asserts every CSV's required
  columns + row counts; `python -m src.data_loader --check` fails loudly on drift.
  Highest-ROI change — it catches all bugs below.
- **`target_states.py`**: reads the correct gene column (`variable` vs `gene_name`);
  builds the target from the concordant core instead of silently collapsing the two
  contrasts; corrects the `toward_Th1` sign (validated by canonical biology).
- **`confidence.py`**: reads the real DE_stats reproducibility columns (with alias
  fallbacks), so the score is no longer stuck at neutral defaults and the off-target
  penalty actually applies.
- Concordance stats are numpy-only (no scipy dependency in the hot path).

## 2. Day-by-day plan (revised)

**Day 0 — De-risk the download + lock the data contract.** Register on CZI VCP,
install `vcp-cli`, pull `GWCD4i.DE_stats.h5ad`. Run the schema gate. (Front-loaded
because Tier 1 alone can't demonstrate the method.)

**Day 1 — Data load + known-biology sanity.** Load 7 CSVs; assert shapes. Build target
vectors; confirm the sanity gate (`toward_Th1` → GATA3 down, TBX21/IFNG up). Compute
and record Ota-vs-Höllbacher concordance.

**Day 2 — Confidence + QC wiring.** Correct-column reproducibility, off-target audit,
per-guide knockdown gate.

**Day 3 — Reachability engine.** NNLS cone fit + reachability spectrum with
shuffled-target and random-perturbation nulls; additivity-aware redundancy penalty.
Enforce CRISPRi = knockdown-only (non-negative weights); surface up-regulation needs
as clearly-labeled non-testable hypotheses.

**Day 4 — Honest benchmarking.** Held-out-gene generalization (primary validity claim);
bootstrap stability selection; LODO (n=4). Additive linear baseline is the primary
model; optional scGen-style comparison only.

**Day 5 — Disease linkage + external validation + evidence.** Autoimmune enrichment
join; **independent-screen and DICE checks (§5b)**; cross-modal coherence; PubMed/Open
Targets/Consensus citations; pathway enrichment.

**Day 6 — Explorer UI + hypothesis cards.** Streamlit: pick axis + condition → ranked
sets → cards (confidence breakdown + reachability + external-validation badges +
pathways + citations + limitations banner). Plus the experiment-design queue.

**Day 7 — Reproducibility, write-up, demo.** Pin exact versions, fix seeds, expand
`tests/`, record demo, finalize README/slides/one-page limitations.

## 3. Risks & mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| h5ad download/registration friction | Medium | **Front-loaded to Day 0**; schema gate catches partial pulls |
| Cosine/alignment inflated by overfitting in G≫constraints | High impact | Reachability spectrum + shuffled-target null + **held-out-gene** eval |
| Confidence score silently degraded (v1 bug) | Was High | **Fixed** + schema gate prevents recurrence |
| Sign error in target vector (v1 bug) | Was High | **Fixed**; known-biology sanity gate is a permanent test |
| Additivity/epistasis oversold | Medium | Additivity-aware penalty; flag multi-gene sets; prioritize k=1–2; experiment-design queue |
| Small n (4 donors) oversold | Medium | Validity rests on held-out-**genes**, not held-out-donors; LODO framed as robustness |
| No local arrayed wet-lab validation | Medium | Replaced by cross-source concordance + **independent public datasets (§5b)** + literature |
| Overclaiming causality/rescue | High impact | Hard-coded limitations banner; hypotheses framed as falsifiable |

## 4. Explicitly out of scope

- Wet-lab validation or any clinical claim.
- Retraining a foundation model or GPU-scale DL.
- Reprocessing the 22M-cell raw matrix from FASTQ/counts.

## 5. Key references

- Zhu et al. *Genome-scale perturb-seq in primary human CD4+ T cells.* bioRxiv 2025.
- Ahlmann-Eltze, Huber & Anders. *Deep-learning-based gene perturbation effect
  prediction does not yet outperform simple linear baselines.* Nat Methods 22,
  1657–1661 (2025).
- Lotfollahi et al. scGen (Nat Methods 2019); CPA (Mol Syst Biol 2023);
  Roohani et al. GEARS (Nat Biotechnol 2023).

## 5b. External validation datasets (independent, public, open) — NEW

The local screen has no arrayed wet-lab validation table, so we triangulate against
independent public data. Each row states the concrete check.

| Resource | Type | Independent check it enables | Access |
|---|---|---|---|
| **Shifrut et al. 2018, *Cell*** — genome-wide CRISPR screens in primary human T cells (PMID 30449619) | Loss-of-function screen | Do our `toward_Th1`/state nominations recur as T-cell function/proliferation regulators in an independent knockout screen? | Open; supplementary hit tables |
| **Schmidt et al. 2022, *Science*** — CRISPRa/CRISPRi screens decoding cytokine regulation (PMID 35113687; PMC9307090) | Gain- **and** loss-of-function screen | Direct sign check: e.g. FOXQ1 selectively dampens Th2 cytokines — do our Th1-pushing nominations agree in direction? A CRISPRa arm also tests the up-regulation hypotheses our CRISPRi assay cannot. | Open; Addgene / BioStudies S-EPMC9307090 |
| **DICE** — Database of Immune Cell Expression/eQTLs (91 donors, 13 sorted immune subsets incl. Th1 and Th2 CD4) | Sorted-cell RNA-seq + eQTL | Independent expression direction: do nominated regulators show the expected Th1-vs-Th2 TPM difference? eQTL layer links them to genetic variation. | Free — dice-database.org/downloads |
| **Open Targets** (already wired in `evidence.py`) | Target–disease genetics | Do `toward_Th1` nominations carry autoimmune/allergic genetic associations consistent with the disease-linkage layer? | Open GraphQL API |
| **Canonical master-regulator biology** (GATA3 = Th2 master TF; T-bet/TBX21 = Th1) | Prior knowledge | Permanent sanity gate: `toward_Th1` must move GATA3 down, TBX21/IFNG/STAT4 up. **Verified.** | Literature (e.g. Nawijn 2001 J Immunol) |

**Why this matters scientifically.** Two orthogonal validations — an independent
perturbation screen (Shifrut/Schmidt) agreeing on *which* genes are regulators, and an
independent expression atlas (DICE) agreeing on *the direction* — would convert the
nominations from "internally reproducible" to "reproducible across labs, assays, and
data modalities." The Schmidt CRISPRa arm is especially valuable: it is the only way,
short of new wet-lab work, to put evidence behind the gain-of-function hypotheses that
CRISPRi structurally cannot test.
