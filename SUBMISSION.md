# Cell-State Reachability

### An oracle that tells you what's *impossible* — before you run the experiment.

**Built with Claude — Life Sciences Hackathon · Research / Lab Track**

> Every AI model in biology predicts what a perturbation *will* do. We built the first
> one that returns a falsifiable **reachable / provably-outside-the-cone** verdict — plus
> the minimal knockdown recipe and a numerical certificate — for a real genome-scale
> CRISPRi screen in primary human CD4⁺ T cells.

---

## Start here

| If you want… | Open |
|---|---|
| **The paper** (preprint PDF) | [`manuscript/main.pdf`](manuscript/main.pdf) |
| **Honest limitations + reinforcement plan** | [`manuscript/limitations_and_reinforcement_plan.pdf`](manuscript/limitations_and_reinforcement_plan.pdf) |
| **Technical Dossier** — all technical write-ups merged (method, novelty, related work, causal foundations, appendices) | [`docs/Technical_Dossier.pdf`](docs/Technical_Dossier.pdf) · editable source [`docs/Technical_Dossier.md`](docs/Technical_Dossier.md) |
| **Interactive explorers** (7 single-file, no server) | [`app/index.html`](app/index.html) |
| **Reproducible notebooks** (`01`–`09` + `bring_your_own_target`) | [`notebooks/`](notebooks/) |
| **Core library** | [`reachability.py`](reachability.py) |
| **One-command reproduction** | [`reproduce.sh`](reproduce.sh) |
| **Headline results & tables** | [`docs/Technical_Dossier.pdf`](docs/Technical_Dossier.pdf) · [`results/`](results/) |

*New here? This page is the elevator view. [`README.md`](README.md) is the full repo map.*

---

## Problem

Nine of ten clinical drug programs never reach approval, and **Phase II is the lowest
transition of the pipeline (30.7%)** — where roughly **half of late-stage failures are
lack of efficacy**, i.e. the target was never a way to move the biology. Every AI tool
in this space is a *forward* predictor: it scores what a perturbation would do. **None
can tell a team, from measured data, that a desired cell-state change is simply out of
reach for the chosen modality** — so programs discover unreachability only after the
screen, or after the trial.

## Method

We treat each measured CRISPRi effect vector as a direction in expression space. Because
knockdown only ever *subtracts*, the states reachable by any non-negative mix of
knockdowns form a **convex cone**. Reachability becomes one geometric question — *does
the target vector lie inside that cone?* — solved by **non-negative least squares
(Lawson–Hanson)** on the real effect matrix of **33,983 knockdowns × 10,282 genes**. The
output is not a similarity score but a **falsifiable verdict** carrying (1) a minimal
ranked knockdown recipe, (2) a **Farkas'/KKT activation certificate** naming the genes no
knockdown mix can deliver, and (3) a confidence score built from the dataset's own
reproducibility, a permutation null, and orthogonal literature/genetics evidence.

## Key result

For the flagship **Th2 → Th1** polarization switch (resting CD4⁺ T cells), the verdict is
**partially reachable**: **39% of the target shift is reachable by knockdown (loss-of-function), 25%
provably requires gene *activation* (gain-of-function), 35% is neither** — with a held-out cosine of
**0.448** and a KKT/Farkas optimality residual of **1.1 × 10⁻¹¹** as numerical proof.
Master-regulator positive controls land correctly — the Th2 driver **GATA3** (must go
down) is recovered at **rank 155/6871 (97.7th percentile)**, and the Th1 driver **TBX21**
(must go up) is correctly anti-aligned under knockdown (**rank 6775/6871**). The operator
generalizes unchanged: transferred to a different assay, cell type, and direction (Norman
et al. 2019 K562 CRISPRa), it recovers a held-out **CEBPA** master-TF state at cosine
**0.878** (permutation *z* = 36.97). Across a **12-cell atlas** (4 transitions × 3
activation conditions) every cell clears its null (held-out *z* = 4.7–45.0), and knockdown
is *never* the majority modality — activation and irreducible components always dominate.

## Impact

The verdict converts a gene list into a **GO / STOP / REDIRECT** decision *before* the
expensive step. On the flagship axis we take the **102 knockdown nominations** the oracle
says are actually needed and cross them against Open Targets druggability × immune-disease
genetics: **45/102 (44%) are hard-to-drug** and only **10 of 102 are clinical-grade
today**. **IRF1** is the headline collision — a top-genetics node with no conventional
drug handle; **JAK2** and **ICOS** are green-lit. A survey of **91 prior methods
(2011–2026)** finds **zero** that return a feasibility verdict on measured effects — this
work is the sole occupant of the *measured × achievability* quadrant. Since **human
genetic support carries 2.6× approval odds** (Minikel 2024), a feasibility-plus-genetics
triage attacks failure exactly where the funnel leaks most.

*Scope, stated up front: nominations are **hypotheses for wet-lab testing**, not validated
targets; CRISPRi is loss-of-function only (activation hypotheses need a separate CRISPRa
arm); multi-gene recipes assume bounded additivity; one primary-cell system, four donors.
This is an experiment-triage instrument, not a target-validation engine.*

---

**Dataset** · Genome-scale CRISPRi Perturb-seq in primary human CD4⁺ T cells — Zhu et al.
2025 (Marson & Pritchard labs, CZI Virtual Cells Platform).
**Built with Claude Science** · method, 91-method survey, manuscript, and 7 explorers in a
one-week build.
