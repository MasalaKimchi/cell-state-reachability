# Cell-State Reachability

**Built with Claude — Life Sciences Hackathon · Research / Lab Track**

*Given a genome-scale screen of gene perturbations, can we reach a desired cell state —
and if not, what is provably stopping us?* This project answers that with a **convex-cone
reachability oracle**: it returns a falsifiable **reachable / provably-outside-the-cone**
verdict for a real CRISPRi Perturb-seq screen in primary human CD4⁺ T cells, plus the
minimal knockdown recipe, a numerical activation certificate, and a confidence score.

> **New here?** [`SUBMISSION.md`](./SUBMISSION.md) is the 2-minute elevator view.
> This README is the **full repository map** — what every file is and where to start.

---

## Start here — the five things worth opening

| If you want… | Open |
|---|---|
| **The 2-minute pitch** | [`SUBMISSION.md`](./SUBMISSION.md) |
| **The paper** (preprint PDF) | [`manuscript/main.pdf`](./manuscript/main.pdf) |
| **All technical write-ups, merged** — method, novelty, related work, causal foundations + 5 appendices | [`docs/Technical_Dossier.pdf`](./docs/Technical_Dossier.pdf) (162 pp) · editable source [`docs/Technical_Dossier.md`](./docs/Technical_Dossier.md) |
| **Interactive explorers** (7 single-file, no server) | [`app/index.html`](./app/index.html) |
| **Reproduce it** | [`reproduce.sh`](./reproduce.sh) · [`notebooks/`](./notebooks/) |

---

## What this is

Each measured CRISPRi effect vector is a direction in expression space. Because knockdown
only ever *subtracts*, the states reachable by any non-negative mix of knockdowns form a
**convex cone**, and reachability becomes one geometric question — *does the target vector
lie inside that cone?* — solved by non-negative least squares on the real effect matrix of
**33,983 knockdowns × 10,282 genes**. The output is not a similarity score but a
**falsifiable verdict** carrying (1) a minimal ranked knockdown recipe, (2) a Farkas/KKT
activation certificate naming the genes no knockdown mix can deliver, and (3) a confidence
score from the dataset's own reproducibility, a permutation null, and orthogonal
literature/genetics evidence.

For the flagship **Th2 → Th1** switch the verdict is *partially reachable*: **39% of the
shift is reachable by knockdown (LOF), 25% provably requires activation (GOF), 35% is
neither** — held-out cosine 0.448, KKT/Farkas residual 1.1 × 10⁻¹¹. Master-regulator
controls land correctly (GATA3 recovered at rank 155/6871; TBX21 correctly anti-aligned
under knockdown), the operator transfers unchanged to a K562 CRISPRa screen (Norman 2019;
held-out CEBPA at cosine 0.878, z = 36.97), and across a 12-cell atlas knockdown is never
the majority modality. Full numbers, tables, and figures are in the
[Technical Dossier](./docs/Technical_Dossier.pdf) and the [manuscript](./manuscript/main.pdf).

## Scope, stated up front

Nominations are **hypotheses for wet-lab testing**, not validated targets. CRISPRi is
loss-of-function only (activation hypotheses need a separate CRISPRa arm); multi-gene
recipes assume bounded additivity; the system is one primary-cell screen across four
donors; matching a transcriptional signature is not functional rescue. This is an
experiment-triage instrument, not a target-validation engine. The manuscript's
[limitations PDF](./manuscript/limitations_and_reinforcement_plan.pdf) and the dossier's
Appendix E (Response to Reviewer 2) treat these limits head-on.

## Data

Marson & Pritchard labs, *Genome-scale Perturb-seq in primary human CD4⁺ T cells*
(Zhu et al., bioRxiv 2025), on the CZI Virtual Cells Platform (MIT license).

- Dataset card: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1

The raw dataset (~22M cells) is not laptop-tractable; the work builds on the authors'
precomputed derived artifacts, in two tiers (both local, both gitignored — see
[`data/README.md`](./data/README.md)):

- **Tier 1** (~36 MB of CSVs): target signatures, per-perturbation QC/effect summaries,
  guide QC, off-target design, autoimmune-disease enrichment, donor metadata. Supports the
  1-D directional sanity check in `notebooks/01`.
- **Tier 2** (`data/GWCD4i.DE_stats.h5ad`, 16.8 GB): the full gene-level effect matrix —
  33,983 pert×cond × 10,282 genes — that the reachability cone runs on. Read selectively
  (subset to significant on-target perturbations + HVGs per condition, cached), never
  loaded whole, on an 18 GB-RAM laptop.

## Method, in four steps

1. **Perturbation dictionary** `E ∈ ℝ^{P×G}`: each row is one perturbation's measured
   effect on the transcriptome (z-scored logFC), per stimulation condition.
2. **Target direction** `d ∈ ℝ^G`: the desired shift (Th2→Th1 signature, or reverse-aging).
3. **Reachability, then minimal set.** Non-negative weights → the reachable shifts form a
   convex cone; a non-negative least-squares fit gives the reachable-vs-outside verdict,
   closest point, and residual, after which sparse selection *within* the cone finds the
   smallest perturbation set. Sparsity = "minimal"; non-negativity = "achievable by knockdown."
4. **Confidence** = dataset-native reproducibility + bootstrap stability + held-out-donor
   generalization + orthogonal literature/Open Targets evidence.

A linear/additive baseline is the **primary** model, not a fallback: a 2025 *Nature Methods*
benchmark shows current deep-learning perturbation predictors do not yet beat simple linear
baselines, so any DL comparison here is explicitly optional and benchmarked, never assumed.

## Novelty delta

"Minimal set of perturbations that moves a cell toward a target state" is a named problem
(Mogrify 2016; CellOracle 2023). The delta here: **measured, not inferred** effect vectors
(the empirical CRISPRi matrix, not a GRN inferred from wild-type data); a **convex-cone
reachability verdict** with a provable outside-the-cone certificate rather than a similarity
score; and a held-out-gene validity test plus a screen-native confidence decomposition. A
survey of 91 prior methods (2011–2026) finds none returning a feasibility verdict on
measured effects. The regulator biology the pipeline recovers (Th1/Th2, aging) is what the
source paper itself reports — so regulator recovery is treated as *validation the method
works*, and novelty is claimed only for the method + decision layer. Full argument and the
91-method survey: [Technical Dossier](./docs/Technical_Dossier.pdf), Parts 2–3.

## Repository map

```
cell-state-reachability/
├── SUBMISSION.md            # ← 2-minute pitch (problem → method → result → impact)
├── README.md                # ← this file: full repo map + navigation
├── CLAUDE.md                # operating manual: verified facts, guardrails, literature anchors
├── DEMO_VIDEO_SCRIPT.md     # walkthrough script for the demo video
├── reachability.py          # THE METHOD: cone fit, signed decomposition, certificate, nulls, spectrum
├── reproduce.sh             # one-command reproduce (pytest + in-module self-test)
├── environment.yml          # conda environment
├── requirements.txt         # pip pins (numpy / scipy / pandas / …)
├── LICENSE                  # MIT
│
├── docs/
│   ├── Technical_Dossier.pdf    # 162-pp merged dossier — all technical write-ups (READ THIS)
│   ├── Technical_Dossier.md     #   editable source of the dossier
│   └── figures/                 #   13 figures embedded by the dossier
│
├── manuscript/              # LaTeX preprint (ICML style)
│   ├── main.pdf                 #   the paper
│   ├── main.tex + sections/     #   editable LaTeX source (00_abstract … 90_supplement)
│   ├── figures/                 #   fig1–5, figS1–9 (png + pdf)
│   ├── references.bib, main.bbl #   bibliography
│   ├── limitations_and_reinforcement_plan.pdf   # honest-limitations companion
│   └── manuscript_source.tar.gz #   self-contained source bundle
│
├── notebooks/               # 01–09 + bring_your_own_target (+ README.md reading route)
│   ├── figstyle.py, make_deg_figures.py   #   shared plotting helpers
│   └── cache/                   #   small cached bundles + exported design cards
│
├── app/                     # interactive explorer — 8 self-contained HTML views + DEPLOY.md
│   ├── index.html               #   landing page linking all 7 explorers
│   └── previews/                #   PNG preview of each view
│
├── scripts/                 # analysis drivers (run_atlas / run_nulls / run_bootstrap /
│                            #   run_a1_sensitivity / run_iv_compliance / run_deg_weighted_eval /
│                            #   build_effect_matrices / a2_conditional_reachability_scaffold)
├── results/                 # atlas + modality + K562 tables, a-series & reviewer-2 outputs,
│                            #   references.csv (+ README.md cataloguing every file)
├── tests/                   # test_reachability.py — 11 tests, run by reproduce.sh
├── data/                    # README.md + Tier-1 CSVs + Tier-2 h5ad (data local & gitignored)
└── analysis_cache/          # cached intermediates (heavy .npz gitignored, small tables tracked)
```

## Reproduce

```bash
bash reproduce.sh                     # pytest (11 tests) + reachability._selftest()

# the method module + batch drivers produce every headline output:
python scripts/run_atlas.py           # 12-cell atlas → results/atlas_reachability.csv
python scripts/run_nulls.py           # held-out-gene significance per cell
python scripts/run_bootstrap.py       # gene-panel subsampling CI on the headline verdict
python scripts/run_a1_sensitivity.py  # A1 verdict sensitivity radius (feeds notebook 09)
python scripts/run_iv_compliance.py   # IV / compliance layer         (feeds notebook 09)
```

Or step through the notebooks — [`notebooks/README.md`](./notebooks/README.md) has the reading
route. Build order: `01` EDA → `02` headline → `03` generalizability (K562 CRISPRa) →
`04` design toolkit → `05` target-ID showcase → `06` reinforcement battery →
`07` cross-cell-type transfer → `08` DEG-weighted evaluation → `09` causal-validation dossier,
plus `bring_your_own_target` for an arbitrary target signature.

## License

MIT. Data: MIT (CZI Virtual Cells Platform). Please cite Zhu et al. 2025.
