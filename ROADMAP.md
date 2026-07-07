# Roadmap — Counterfactual Biology Explorer (7 days, solo, CPU-only)

This is the working roadmap. The polished versions live in the accompanying Word
document and slide deck. This copy is the source of truth for the repo.

## 0. Framing decision (why we pivoted)

The featured Marson/Pritchard dataset is a **CD4+ T-cell** genome-scale CRISPRi
Perturb-seq screen, not AML/HSPC data. Forcing an AML narrative onto it would
require sourcing separate AML scRNA-seq and treating the Perturb-seq only as a weak
prior — more moving parts, weaker validation, higher risk in a one-week solo build.

Instead we keep the *counterfactual method* exactly as proposed and apply it to a
transition the dataset can support **and independently sanity-check**:

- **Primary axis — Th1 ↔ Th2 polarization.** Chosen because the dataset ships a
  Th2-vs-Th1 target signature derived from **two independent source contrasts**
  (Ota 2021 and Höllbacher 2021), letting us cross-check regulator nominations for
  consistency across signature sources.
- **Secondary axis — CD4+ T-cell aging (aged → young-like).** The dataset ships an
  aging signature (Yaza 2022) we can build the reverse-aging target vector from.

The counterfactual question becomes: *what is the minimal set of gene knockdowns
whose measured effects best move the CD4+ T-cell transcriptome from state A toward
state B?*

## 0b. What data we actually have locally (ground truth for this plan)

Seven derived supplementary CSV tables (~36 MB) are checked into `data/` and drive
everything below. Verified shapes/contents:

| File | Rows | Grain | Key fields we use |
|---|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | perturbation × condition | on/off-target flags, DE-gene counts, `crossdonor_correlation_mean/min`, `crossguide_correlation`, `ontarget_effect_size/significant` |
| `Th2_Th1_polarization_signature...csv` | 37,288 | gene × source-contrast | `log_fc`, `zscore` for Th2-vs-Th1 (Ota 2021 **and** Höllbacher 2021) |
| `CD4T_aging_signature...csv` | 10,000 | gene | `log_fc`, `zscore` for aged-vs-young (Yaza 2022) |
| `guide_kd_efficiency.suppl_table.csv` | 73,765 | guide × condition | `signif_knockdown`, guide vs NTC expression, `rank`, t-stat |
| `sgrna_library_metadata.suppl_table.csv` | 31,110 | guide | TSS distance, off-target / nearby non-target genes, alt alignments |
| `cluster_autoimmune_enrichment...csv` | 5,236 | cluster × disease × gene-set | `odds_ratio`, `p_adj_fdr`, `intersecting_genes` across 17 autoimmune diseases |
| `sample_metadata.suppl_table.csv` | 12 | sample | 4 donors × 3 conditions; age (22–34), sex, ethnicity |

Descriptive totals worth knowing: 11,526 unique target genes in `DE_stats`;
21,216 perturbation×condition rows with a *significant on-target knockdown*; 2,837
flagged for off-target activity; 54,094 of 73,765 guides pass `signif_knockdown`.

**Two conditions matter for scope:**
- Three culture conditions everywhere: **Rest, Stim8hr, Stim48hr** (~11.3k
  perturbations each). Every effect and target vector is condition-specific, so
  state transitions are defined *within* a condition (or explicitly compared across).
- Only **4 donors** — "held-out donor" means leave-one-donor-out (LODO) with n=4,
  so report it as a small-n robustness check, not a formal generalization test.

**Honesty correction from the original plan:** `DE_stats.suppl_table.csv` is a
per-perturbation **summary** table (counts, reproducibility, flags), **not** the
full gene×perturbation logFC matrix. The sparse reconstruction solver (`E ∈
R^{P×G}`) needs the gene-level matrix in `GWCD4i.DE_stats.h5ad`, which is **not
local yet**. Likewise, no arrayed wet-lab validation table is present locally, so
our polarization "ground truth" is internal reproducibility + literature, not an
arrayed-hit overlap. This splits the work into two tiers.

## 0c. Two tiers of analysis (runnable *now* vs. after one download)

**Tier 1 — CSV-only, runnable today** (no auth, no h5ad, seconds on CPU). Uses
target signatures, per-perturbation summaries, guide/off-target QC, disease
enrichment, and donor metadata. This alone is a complete, defensible submission.

**Tier 2 — adds `GWCD4i.DE_stats.h5ad`** (one `vcp-cli` download; see
`data/README.md`). Unlocks the full gene-level counterfactual reconstruction solver
and quantitative reconstruction metrics. Strictly additive — Tier 1 stays the
graded core so a download failure never sinks the submission.

## 1. Success criteria (falsifiable, not vibes)

A submission is "done" when:

1. **(Tier 1)** Nominated regulators for an axis are reproducible: they rank
   consistently across the two polarization source contrasts (Ota vs Höllbacher),
   and their perturbations carry high cross-donor/cross-guide correlation and
   significant on-target knockdown. Report rank-overlap + a null.
2. **(Tier 2)** The sparse solver reconstructs the target direction with cosine
   similarity meaningfully above a random-perturbation null (report effect size +
   bootstrap CI), and holds under leave-one-donor-out (n=4).
3. Every hypothesis card renders a confidence score with its component breakdown,
   an off-target / knockdown-quality flag, and ≥1 literature/Open Targets citation.
4. The whole pipeline reproduces from a fixed seed + pinned environment on CPU in
   under ~15 min (Tier 1 in seconds).

## 1b. Extensive analysis catalog (what we can actually do with this data)

Grouped by tier. Each item names the file(s) it needs so nothing is aspirational.
"★" marks high-value-for-effort items worth demoing.

### A. Target-state & signature analysis — *Tier 1*
1. ★ **Build both target vectors** — `d_polarization` (Th2→Th1) and `d_aging`
   (aged→young-like, sign-flipped) from the signature CSVs, in z-score space.
2. ★ **Cross-source signature agreement** — correlate Ota 2021 vs Höllbacher 2021
   Th2-vs-Th1 logFCs; keep the concordant core as a robust target and report the
   discordant tail as lower-confidence. (Unique to polarization; free robustness.)
3. **Signature anatomy** — top up/down genes, volcano plots, overlap between the
   aging and polarization signatures (shared vs axis-specific biology).
4. **Condition-resolved targets** — quantify how the achievable transition differs
   across Rest / Stim8hr / Stim48hr.

### B. Perturbation-effect & directionality analysis — *Tier 1 (summary) / Tier 2 (gene-level)*
5. ★ **Directional-alignment ranking (Tier 2)** — score every perturbation by
   cosine/dot alignment of its effect vector with each target direction; rank the
   knockdowns that individually push toward Th1 (or young-like).
6. **Effect-magnitude landscape (Tier 1)** — using `n_up_genes`/`n_down_genes`/
   `n_downstream` from `DE_stats`, map which perturbations are "hubs" (broad
   transcriptional impact) vs. narrow, per condition.
7. **Condition-dependent perturbations (Tier 1)** — flag genes whose DE-gene count
   or on-target effect changes sharply between Rest and Stim (context-specific
   regulators are the interesting ones).
8. **Sparse minimal-set solver (Tier 2)** — LASSO / OMP / greedy forward selection
   over the gene-level dictionary; ranked minimal sets of size k=1…N with
   reconstruction quality and a random-set null.

### C. Confidence, reproducibility & QC — *Tier 1*
9. ★ **Per-nomination confidence score** — combine `crossdonor_correlation_mean`,
   `crossguide_correlation`, `ontarget_significant`, `offtarget_flag`
   (`DE_stats`) with per-guide `signif_knockdown` and `rank` (`guide_kd_efficiency`).
10. **Knockdown-quality gate** — drop or down-weight nominations whose guides don't
    achieve significant knockdown; report the achievable-vs-nominal gap.
11. **Off-target auditing** — use `sgrna_library_metadata` (TSS distance, nearby
    non-target genes, alternate alignments) to flag nominations whose phenotype may
    be confounded by off-target cutting.
12. **Reproducibility atlas** — distribution of cross-donor vs cross-guide
    correlation; show that high-confidence perturbations cluster where both agree.

### D. Disease-relevance layer — *Tier 1* (this is the "healthier state" bridge)
13. ★ **Autoimmune enrichment linkage** — `cluster_autoimmune_enrichment` connects
    perturbation clusters to 17 autoimmune diseases (asthma, Crohn's, RA, SLE, MS,
    T1D, psoriasis, …). Prioritize nominations whose downstream programs are
    enriched for disease GWAS genes → an interpretable "why this matters" per card.
14. **Disease-specific target shortlists** — invert the enrichment: for a chosen
    disease, list the enriched clusters/`intersecting_genes` as candidate axes.
15. **Negative-control check** — the table's `negative_control_disease` column lets
    us confirm enrichment is disease-specific, not generic.

### E. Donor & covariate analysis — *Tier 1*
16. **Leave-one-donor-out robustness** — with 4 donors, re-derive rankings dropping
    each donor; report rank stability (honest small-n caveat).
17. **Covariate sanity** — check age/sex balance in `sample_metadata` before making
    any aging-axis claim (n=4, ages 22–34 → aging axis is exploratory only).

### F. Evidence & interpretation — *Tier 1, uses connected MCP tools*
18. ★ **Literature/target evidence** — for each nominated gene, query **PubMed**,
    **Open Targets**, and **Consensus** for known roles in Th polarization / T-cell
    aging / the linked autoimmune disease; attach citations (support, not proof).
19. **Pathway / gene-set enrichment** on nominated sets for interpretability.
20. **Druggability annotation** — optional **ChEMBL** / **Open Targets** lookup for
    whether a nominated regulator is tractable (framed as future-direction only).

### G. Packaging — *Tier 1*
21. ★ **Hypothesis cards** — one card per nomination: confidence breakdown +
    knockdown quality + off-target flag + disease linkage + citations + a hard-coded
    limitations banner.
22. **Streamlit explorer** — pick axis + condition → ranked nominations → cards.

## 2. Day-by-day plan

### Day 1 — Data load + sanity (Tier 1 data is already local)
- Load the 7 CSVs via `data_loader.py`; assert shapes against §0b as a correctness
  gate. Subset `DE_stats` by condition; index guide QC by `perturbed_gene_id`.
- **Reproduce a known result** (e.g. a canonical Th regulator appearing with the
  expected sign in the polarization signature) before building anything new.
- *In parallel, kick off the Tier-2 download:* register on CZI Virtual Cells
  Platform, install `vcp-cli`, pull `GWCD4i.DE_stats.h5ad`. If it stalls, Tier 1
  proceeds unaffected.

### Day 2 — Target-state vectors + QC/confidence wiring (Tier 1)
- Build `d_polarization` (Th2→Th1) and `d_aging` (aged→young-like) in z-score space
  (target_states.py); compute the Ota-vs-Höllbacher concordance core (catalog #2).
- Wire per-perturbation reproducibility columns and per-guide `signif_knockdown`
  into `confidence.py` (catalog #9–12).

### Day 3 — Counterfactual engine
- **Tier 1 first:** directional alignment ranking + effect-magnitude landscape
  (catalog #5–7) so there is a working ranked nomination list with *only* CSVs.
- **Tier 2 (if h5ad landed):** sparse minimal-set solver — LASSO, OMP, greedy
  forward selection over the gene-level dictionary (catalog #8). Enforce **CRISPRi =
  knockdown-only** as a hard constraint; surface upregulation needs as clearly
  labeled non-testable hypotheses. Output ranked sets k=1…N with reconstruction
  quality + random-set null.

### Day 4 — Confidence + honest benchmarking
- Confidence module: dataset reproducibility + off-target audit
  (`sgrna_library_metadata`, catalog #11) + **bootstrap stability selection** +
  **leave-one-donor-out** (n=4) robustness.
- Benchmarks: additive linear baseline (primary), random-perturbation null, and
  (optional/stretch) a small scGen-style latent-arithmetic model — reported as a
  comparison only, per the Nature Methods 2025 finding.

### Day 5 — Disease linkage + evidence + interpretation
- **Disease-relevance layer** (catalog #13–15): join nominations to
  `cluster_autoimmune_enrichment`; attach the enriched autoimmune disease(s) and
  `intersecting_genes` to each card, with the negative-control check.
- `evidence.py`: for each nominated gene, query **PubMed**, **Open Targets**, and
  **Consensus** for known roles; attach citations. Evidence is *support*, never proof.
- `pathways.py`: gene-set / GO enrichment on the nominated set for interpretability.

### Day 6 — Explorer UI + hypothesis cards
- Lightweight Streamlit app (`app/explorer.py`): pick an axis and condition → view
  ranked minimal sets → expandable hypothesis cards (confidence breakdown + pathways
  + citations + limitations banner).
- CPU-only friendly; no GPU, no heavy model serving.

### Day 7 — Reproducibility, write-up, demo
- Pin environment, fix seeds, add `tests/`, record a short demo.
- Finalize README, slides, and a one-page limitations statement.

## 3. Risks & mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| `GWCD4i.DE_stats.h5ad` download/registration friction | Medium | **Tier 1 needs no download** — the graded core runs on the local CSVs; h5ad is additive |
| RAM blow-up loading the h5ad (Tier 2) | Medium | Load one layer (z-score) as float32; subset to significant on-target perturbations + HVGs; cache to `.npz` |
| No arrayed wet-lab validation table locally | Medium | Replace "arrayed-hit overlap" with cross-source-signature concordance + internal reproducibility + literature; state the swap openly |
| Overclaiming causality/rescue | High impact | Hard-coded limitations banner on every output; hypotheses framed as falsifiable |
| Combinatorial predictions untestable | Medium | Flag multi-gene sets as extrapolations; prioritize k=1–2 nominations |
| Small n (4 donors) oversold | Medium | LODO framed as robustness, not generalization; aging axis flagged exploratory |
| Scope creep across both axes | Medium | Polarization is the graded deliverable; aging is a stretch |
| MCP auth for some literature tools | Low | PubMed, Open Targets, bioRxiv, ChEMBL, ClinicalTrials, Consensus are connected no-auth; Wiley/Synapse/BioRender need OAuth (optional) |

## 4. Explicitly out of scope

- Wet-lab validation or any clinical claim.
- Retraining a foundation model or GPU-scale DL.
- Reprocessing the 22M-cell raw matrix from FASTQ/counts.

## 5. Key references

- Zhu et al. *Genome-scale perturb-seq in primary human CD4+ T cells* bioRxiv 2025.
- Ahlmann-Eltze, Huber & Anders. *Deep-learning-based gene perturbation effect
  prediction does not yet outperform simple linear baselines.* Nat Methods 22,
  1657–1661 (2025).
- Lotfollahi et al. scGen (Nat Methods 2019); CPA (Mol Syst Biol 2023);
  Roohani et al. GEARS (Nat Biotechnol 2023).
