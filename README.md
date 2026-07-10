# Counterfactual Biology Explorer

**Built with Claude — Life Sciences Hackathon (Research Track)**

*What is the minimal set of gene perturbations that shifts a primary human CD4+ T cell from one transcriptional state toward a target "healthier" state — and how confident should we be?*

---

## TL;DR

Rather than reporting *what is different* between two cell states (differential
expression), this project generates *interpretable, falsifiable counterfactual
hypotheses* about **what changes might move a cell toward a target state**. It does
this by treating each genome-scale CRISPRi perturbation as a measured causal
"effect vector" and asking a geometric question: **is a target state reachable by some
combination of knockdowns**, and if so, what is the **smallest set** that gets closest?
The headline output is a falsifiable **reachable / provably-outside-the-cone** verdict,
not just a similarity score. Every hypothesis ships with an explicit **confidence
score** built from the dataset's own reproducibility metrics, held-out validation, and
orthogonal literature/database evidence.

## Why this framing (and an honest scope note)

The hackathon abstract originally targeted Acute Myeloid Leukemia (AML) using the
Marson/Pritchard Perturb-seq dataset. **That dataset is a genome-scale CRISPRi
Perturb-seq screen in primary human CD4+ T cells — not AML/HSPC data.** We
therefore reframe the identical *method* onto a state transition that this dataset
can actually support:

- **Primary axis:** Th1 ↔ Th2 helper-T polarization balance
- **Secondary axis:** CD4+ T-cell aging signature (aged → young-like)

Both target-state signatures are provided directly by the dataset authors. The
polarization signature ships from **two independent source contrasts** (Ota 2021,
Höllbacher 2021), which we use for a cross-source robustness check in lieu of an
arrayed wet-lab validation table (not present in our local data). See
[`ROADMAP.md`](./docs/ROADMAP.md) for the full rationale and the analysis catalog.

## How this maps to the hackathon (from the event page)

*Verified against the [event details page](https://cerebralvalley.ai/e/built-with-claude-life-sciences)
on 2025-07-07. The page publishes **no scored judging rubric** — only that finalists
are chosen by a panel (Anthropic + our partner, Gladstone Institutes) and that those
selections are final. The priorities below are read off the **stated deliverable and
prize structure**, not an official criteria list.*

- **Track — Research ("Build From the Bench").** The mandate is to *start from a
  biological question, use Claude Science to answer it, and submit something discrete —
  "a finding, a trained model, an analysis others can reproduce."* A reproducible
  analysis with figures is exactly the expected shape; this repo is built to that.
- **Featured dataset, verbatim prompt.** This CD4+ T-cell Perturb-seq screen (Marson +
  Pritchard/Stanford) is a highlighted Research-track dataset, framed as **"Find new
  drug targets in this T cell Perturb-seq dataset."** Our knockdown nominations *are*
  target hypotheses — align the framing accordingly.
- **Gladstone Special Prize.** A dedicated award for the entry with **"most potential
  to advance science that can overcome disease"** ($10k). Our autoimmune-enrichment +
  Open Targets disease layer is the direct hook for it — make disease relevance
  prominent, not an appendix.
- **Logistics.** Fully virtual, 7-day build, team size ≤ 2; each participant gets one
  month of Max 20x + $200 API credits. Prize pool $100k in credits across both tracks.
- **What this implies for us.** With no published rubric, we optimize for what the page
  *does* reward: a **discrete, reproducible deliverable** (fixed seed + pinned env), a
  **new-drug-target** reading of the output, and a **disease-impact** narrative — all
  wrapped in the honest-limitations stance below so the claims survive expert scrutiny.

## Data

Marson & Pritchard labs, *Genome-scale perturb-seq in primary human CD4+ T cells*
(Zhu et al., bioRxiv 2025). Hosted on the CZI Virtual Cells Platform under an MIT
license.

- Dataset card: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- Analysis repo & supplementary tables: https://github.com/emdann/GWT_perturbseq_analysis_2025

**Feasibility on a laptop (CPU-only):** the raw dataset is ~22M cells and is *not*
tractable on a laptop. We build on the authors' **precomputed derived artifacts**,
and the work is split into two tiers by what's needed.

**Tier 1 — local (checked into `data/`, ~36 MB of CSVs):** the **warm-up layer** —
target signatures, per-perturbation QC/effect summaries, guide QC, off-target design,
autoimmune-disease enrichment, and donor metadata. Tier 1 supports a 1-D directional
ranking (used as a sanity check in `notebooks/01`); the headline convex-cone method runs
on the Tier-2 matrix below.

| Local file | Rows | What it is |
|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | Per **perturbation × condition** summary (DE-gene counts, on/off-target flags, cross-donor & cross-guide reproducibility). *Summary, not the gene-level matrix.* |
| `Th2_Th1_polarization_signature…csv` | 37,288 | Th2→Th1 target signature, two source contrasts |
| `CD4T_aging_signature…csv` | 10,000 | Aged→young target signature |
| `guide_kd_efficiency.suppl_table.csv` | 73,765 | Per-guide knockdown QC |
| `sgrna_library_metadata.suppl_table.csv` | 26,504 | Guide design / off-target annotation |
| `cluster_autoimmune_enrichment…csv` | 5,236 | Perturbation-cluster × 17 autoimmune diseases |
| `sample_metadata.suppl_table.csv` | 12 | 4 donors × 3 conditions, demographics |

**Tier 2 — local (`data/GWCD4i.DE_stats.h5ad`, 16.8 GB):** the full gene-level
effect matrix that the **reachability cone runs on** — i.e. the headline method.
This is where every result in [`RESULTS.md`](./docs/RESULTS.md) — the headline verdict, the
modality triage, and the generalizability transfers — is computed.

| Artifact | Shape | What it is |
|---|---|---|
| `GWCD4i.DE_stats.h5ad` | 33,983 pert×cond × 10,282 genes | Per-perturbation effect matrix (logFC, z-score, p) — **the reachability dictionary `E`** |

The matrix is 16.8 GB on disk (not the ~1.4 GB the early data card estimated); on an
18 GB-RAM laptop it is read selectively — subset to significant on-target perturbations
+ HVGs per condition and cached to `analysis_cache/atlas_work/inputs.npz` and `notebooks/cache/`, never
loaded whole. See [`data/README.md`](./data/README.md).

## Method (baseline-first, honestly benchmarked)

1. **Perturbation dictionary** `E ∈ R^{P×G}`: each row is one perturbation's measured
   causal effect on the transcriptome (z-scored logFC), per stimulation condition.
2. **Target direction** `d ∈ R^G`: the desired transcriptomic shift (e.g. the Th2→Th1
   signature, or the reverse-aging vector).
3. **Reachability first, then minimal set.** Because CRISPRi is loss-of-function and
   weights are **non-negative**, the reachable transcriptome shifts form a **convex
   cone**. We first solve a non-negative least-squares fit of `d` inside that cone to get
   the *reachable-vs-outside* verdict + closest point + residual; then apply sparse
   selection (LASSO / OMP / greedy forward) *within* the cone to find the **smallest**
   perturbation set that reconstructs the reachable component. Sparsity operationalizes
   "minimal"; non-negativity operationalizes "achievable by knockdown."
4. **Confidence** = combination of (a) dataset-native reproducibility (cross-guide and
   cross-donor correlation, on-target knockdown significance, off-target flags),
   (b) bootstrap **stability selection** frequency, (c) held-out-donor generalization,
   and (d) orthogonal **literature/Open Targets** evidence for the nominated gene.

A deliberately simple **linear/additive baseline is the primary model, not a
fallback.** A 2025 *Nature Methods* benchmark (Ahlmann-Eltze, Huber & Anders) shows
current deep-learning perturbation predictors do **not** yet beat simple linear
baselines — so any DL component here is an explicitly-optional comparison, never an
unbenchmarked claim.

## Related work — and our specific novelty delta

"Find a minimal set of perturbations that moves a cell toward a target state" is a
**named, established problem**, so our contribution has to be stated as a precise
*delta*, not a category. The two closest prior methods:

- **Mogrify** (Rackham et al., *Nat Genet* 48:331–335, 2016) "combines gene expression
  data with regulatory network information to predict the reprogramming factors
  necessary to induce cell conversion," applied across 173 human cell types, and
  validated two new transdifferentiations. That is the minimal-set-to-target idea — but
  with **transcription factors / gain-of-function** and a **regulatory-network** basis.
- **CellOracle** (Kamimoto et al., *Nature* 614:742–751, 2023) uses **inferred GRNs** to
  simulate KO/overexpression effects, "converted into a vector map of transitions in
  cell identity," benchmarked against a **randomized-model null**. That is the
  perturbation-to-state-transition idea, with null-model hygiene close to ours.

**Our delta (what is genuinely new here):**

1. **Measured, not inferred, effect vectors.** Our perturbation dictionary is the
   *empirically measured* genome-scale CRISPRi effect matrix — not a GRN inferred from
   wild-type data (CellOracle) nor a network-influence heuristic (Mogrify).
2. **A convex-cone reachability verdict.** Loss-of-function effects + non-negative
   weights span a convex cone; we return a formal **reachable vs. (provably) outside the
   knockdown cone** answer with the closest reachable point and a residual — not just a
   similarity score. *Precise statement: "unreachable" means the target has a component
   outside the non-negative span of the measured knockdown effect vectors (which are
   themselves mixed-sign); this often — not always — corresponds to needing activation.*
3. **Held-out-gene validity + a screen-native confidence decomposition** (see Method §4).

Note also that the biology our pipeline surfaces (Th1/Th2 and aging regulators) is what
the **source paper itself reports** — the Zhu et al. abstract states perturbation
signatures let them nominate "regulators of Th1 and Th2 polarization and of age-related
T cell phenotypes." We therefore treat regulator recovery as **validation that the
method works**, and claim novelty only for the method + decision layer above — never for
the regulator lists themselves.

## Honest limitations (read before trusting any output)

- **CRISPRi is loss-of-function only.** We can directly nominate *knockdowns*.
  Gain-of-function hypotheses are extrapolations that this assay cannot test.
- **Additivity is an assumption.** The screen perturbs single genes; any multi-gene
  "minimal set" is an *untested combinatorial extrapolation* that ignores epistasis.
- **Transcriptome ≠ phenotype.** Matching a transcriptional signature does not prove
  functional rescue.
- **The target state is a proxy.** "Healthier" is operationalized as a transcriptomic
  signature, not a clinical outcome.
- Outputs are **ranked, falsifiable hypotheses for future experimental validation** —
  not conclusions.

## Status: built and validated

The method is implemented, the Tier-2 matrix is local, and the full analysis has run.
The headline Th2→Th1 verdict, the 12-cell reachability atlas, the modality triage of 102
knockdown nodes, and a cross-dataset K562 CRISPRa transfer demo are all reproduced in the
notebooks and written up in the results docs. Two later reinforcement passes harden the
claims against their own stated limitations:

- **Reinforcement battery** (`06_reinforcement_analyses.ipynb`, `docs/REINFORCEMENT_RESULTS.md`) —
  the non-negativity constraint costs only +0.018 held-out cosine but is the sole source of the
  certificate (**L4**); the modest 0.448 cosine is **71% of the achievable knockdown ceiling** (**L5**);
  every recommended recipe is additive-safe 12/12 with a 5.6× margin to the cap (**L2**); and a runnable
  dual-modality certificate test verifies at AUROC 0.999 on synthetic ground truth (**L1**).
- **Cross-cell-type transfer** (`07_cross_celltype_transfer.ipynb`, `docs/CROSS_CELLTYPE_TRANSFER.md`) —
  running the unchanged `reachability.py` on the Replogle 2022 K562 and RPE1 essential-gene screens
  (via CZI) shows effect *direction* transfers across cell types while the specific minimal *recipe*
  does not — a robustness result with a sharp, honest boundary.

**Where to start reading:**

| Doc | What it covers |
|---|---|
| [`RESULTS.md`](./docs/RESULTS.md) | The headline Th2→Th1 verdict + the full expansion: signed reachability, the 12-cell atlas, the genetics × druggability **modality triage**, **cross-dataset & cross-cell-type generalizability**, the **experimental-design toolkit** (`design_experiment()`), and the nine method-improvement in-silico results. **Start here.** |
| [`NOVELTY.md`](./docs/NOVELTY.md) | What is scientifically new (the convex-cone-reachability delta), the **real-world / disease-impact** case, and the **field positioning** (why the inverse-feasibility question matters when forward-prediction models lose to linear baselines). |
| [`RELATED_WORK.md`](./docs/RELATED_WORK.md) | The citation-grounded survey of 91 prior methods across three research communities, with the capability matrix and landscape map. |
| [`CAUSAL.md`](./docs/CAUSAL.md) | The design-based causal-inference reframe, the IV/compliance trust layer, the A1–A6/B1–B4 research agenda, the assumption-by-assumption **validation ledger**, and an adversarial dataset critique. |
| [`ROADMAP.md`](./docs/ROADMAP.md) | The 3-day hackathon build plan (packaging + cross-dataset replication of existing results). |
| [`CLAUDE.md`](./CLAUDE.md) | Operating manual: verified facts, guardrails, literature anchors. |

Supporting writeups, all under [`docs/`](./docs):

| Doc | What it covers |
|---|---|
| [`REINFORCEMENT_RESULTS.md`](./docs/REINFORCEMENT_RESULTS.md) | The L1/L2/L4/L5 reinforcement battery (notebook 06). |
| [`CROSS_CELLTYPE_TRANSFER.md`](./docs/CROSS_CELLTYPE_TRANSFER.md) | Replogle 2022 K562/RPE1 transfer, and its independent reproduction (notebook 07). |
| [`GENERALIZABILITY_SURVEY.md`](./docs/GENERALIZABILITY_SURVEY.md) | The 13-candidate-dataset survey and application map. |
| [`REVIEWER_RESPONSE.md`](./docs/REVIEWER_RESPONSE.md) | Dataset-limitations response: the G4/D4/T2/T3 analyses and a ready-to-paste limitations paragraph. |

> **Caveat on `RESULTS.md` §8.3.** Seven of the nine "in-silico results" in that section have no
> committed figure or table — they were rendered in a working session and never written to the repo.
> The section opens with a status table saying exactly which. Regenerate before publication.

**Reproduce:**

```bash
bash reproduce.sh                     # pytest (11 tests) + reachability._selftest()

# the method module + the batch drivers produce every headline output:
python scripts/run_atlas.py           # 12-cell atlas → analysis_cache/atlas_work/*.json, results/atlas_reachability.csv
python scripts/run_nulls.py           # held-out-gene significance per cell (lean pass)
python scripts/run_bootstrap.py       # gene-panel subsampling CI on the headline verdict
python scripts/run_a1_sensitivity.py  # A1 verdict sensitivity radius  (feeds notebook 09)
python scripts/run_iv_compliance.py   # IV / compliance layer          (feeds notebook 09)
```

Or step through the notebooks — [`notebooks/README.md`](./notebooks/README.md) has the reading
route. Build order: `01` EDA → `02` headline → `03` generalizability (K562 CRISPRa) → `04` design
toolkit → `05` target-ID showcase → `06` reinforcement battery → `07` cross-cell-type transfer →
`08` DEG-weighted evaluation → `09` causal-validation dossier, plus `bring_your_own_target` for an
arbitrary target signature.

**Repository layout:**

```
cell-state-reachability/
├── README.md                # this file — framing, method, related work, hackathon fit
├── CLAUDE.md                # operating manual (verified facts, guardrails, literature anchors)
├── reachability.py          # the method: cone fit, signed decomposition, certificate, nulls, spectrum
├── reproduce.sh             # one-command reproduce (pytest + in-module self-test)
├── environment.yml          # conda environment
├── requirements.txt         # pip pins (numpy/scipy/pandas/…)
├── tests/                   # test_reachability.py — 11 tests, run by reproduce.sh
├── docs/                    # every narrative writeup (start with RESULTS.md)
│   ├── RESULTS.md                    #   headline verdict + atlas + modality triage ← primary
│   ├── NOVELTY.md                    #   method delta (vs prior art) + disease impact + positioning
│   ├── RELATED_WORK.md               #   the 91-prior-method survey (citation-grounded)
│   ├── CAUSAL.md                     #   causal reframe + IV/compliance + validation ledger + critique
│   ├── ROADMAP.md                    #   the 3-day hackathon build plan
│   ├── REINFORCEMENT_RESULTS.md      #   L1/L2/L4/L5 battery (nb06)
│   ├── CROSS_CELLTYPE_TRANSFER.md    #   K562/RPE1 transfer + reproduction (nb07)
│   ├── GENERALIZABILITY_SURVEY.md    #   13-dataset survey + application map
│   ├── REVIEWER_RESPONSE.md          #   dataset-limitations response + limitations paragraph
│   └── figures/                      #   figures embedded by the docs above (tracked)
├── scripts/                 # analysis drivers: run_atlas/run_nulls/run_bootstrap/run_iv_compliance,
│                            #   run_a1_sensitivity, run_deg_weighted_eval, build_effect_matrices, build_nbB
├── notebooks/               # 01–09 + bring_your_own_target + README.md (reading route)
│   ├── figures/             #   notebook figures (gitignored; doc-embedded ones copied to docs/figures/)
│   └── cache/              #   small cached bundles + exported design cards
├── app/                     # interactive explorer — 8 self-contained HTML views + data + DEPLOY.md
│   └── previews/            #   PNG previews of each view
├── results/                 # atlas + modality + K562 tables, a-series outputs, references.csv
├── manuscript/              # LaTeX manuscript (sections/ + figures/)
├── data/
│   ├── README.md            # data provenance + tiers
│   ├── GWCD4i.DE_stats.h5ad # Tier-2 effect matrix (16.8 GB, local, gitignored)
│   └── *.suppl_table.csv    # Tier-1 supplementary tables (local, gitignored)
└── analysis_cache/          # cached intermediates — heavy .npz gitignored, small tables tracked
    ├── atlas_work/          #   cached inputs.npz + per-cell atlas JSONs + bootstrap_ci.json
    ├── nb_out/              #   reinforcement outputs (L1–L5) + figR1
    ├── czi_data/            #   K562/RPE1 aligned effects + per-perturbation transfer table
    └── czi_fig/             #   cross-cell-type transfer figures (nb07)
```

The method delta, disease impact, and field positioning live in [`NOVELTY.md`](./docs/NOVELTY.md);
the data provenance lives in [`data/README.md`](./data/README.md); the build plan and its
risk log are in [`ROADMAP.md`](./docs/ROADMAP.md).

## License

MIT. Data: MIT (CZI Virtual Cells Platform). Please cite Zhu et al. 2025.
