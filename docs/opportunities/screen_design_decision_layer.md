# A3 — The Screen-Design Decision Layer

**What CombiCone does that a forward predictor cannot: it tells you *what to run
next* and *which axis your library can never reach*, and it certifies both.**

A forward predictor answers "what will combination X do?" The screen-design
questions that actually gate a campaign are different:

1. **Triage** — of the thousands of combinations I *could* run, which are worth
   the well? (rank by predicted emergence)
2. **Certify** — of the ones I *did* run, which are genuinely emergent above
   measurement noise, not just large? (two-bar infeasibility certificate)
3. **Recommend** — what is the next diverse batch to measure? (certificate-guided
   acquisition)
4. **Name the missing axis** — which perturbation is my library structurally
   unable to represent, so that no amount of *combining* what I have will reach
   it? (aggregate the infeasibility certificates into an "unmet-demand" direction)

Question 4 is the one no forward predictor provides. A predictor reasons about
inputs it is given; the certificate is a statement about the *library itself* —
derived from the geometry of what the library **cannot** represent. This
document runs the whole flow end-to-end on the real Norman K562 CRISPRa
substrate and reproduces the falsifiable test behind question 4.

Everything below is **COMPUTED from real repo data** (`combicone_substrate.npz`,
105 single-gene effect atoms × 5045 genes, 131 measured double-perturbation
effects; Norman et al. 2019, GEARS-processed). No values are simulated. The one
re-derivation (the three statistical controls) is checked against the repo's own
published `results/screenloop_campaign.json` and matches to 6 decimals for every
deterministic quantity.

**Scope (model-relative, verbatim from the code).** *A separator certifies a
direction outside the non-negative cone of the supplied effect atoms under the
chosen metric; a nominated perturbation is a ranked hypothesis about an unmet
library axis, not a validated intervention.* Unreachable means "outside the cone
of the measured singles," not "biologically impossible."

---

## Part 1 — Reproducing the held-out-gene recovery (question 4, falsified)

### The test

`screenloop.held_out_single_recovery` is the falsifiable version of "name the
missing axis." For each single-gene atom that participates in ≥ 2 measured
combinations (53 of the 105 atoms qualify):

1. **Remove** that gene's atom from the library (the cone of the other 104).
2. Take **only the combinations that involved it** — now each falls outside the
   reduced cone and the core emits a **dual separator** (a Farkas/KKT certificate
   direction the library cannot reach).
3. **Aggregate** those separators, residual-weighted, into one "unmet-demand"
   direction (`nominate_atoms`).
4. **Rank all 105 candidate atoms** by alignment to that direction. If the
   certificate captured the hidden axis, the true held-out gene should return to
   the **top**.

Two honest controls are reported alongside:

- **Naive "average-the-combos" baseline** — rank candidates by cosine to the
  *mean* of the combos that involved the held-out gene (no cone, no separator).
  This is the obvious shortcut ("the missing gene dominates its own
  combinations"); the separator has to *beat* it to be doing real work.
- **Magnitude-only ranker** — rank candidates by effect L2 norm. If recovery
  were a size artifact, this would win.

### Result (COMPUTED from `combicone_substrate.npz`)

| Ranker | median rank | mean rank | top-1 | top-5 |
|---|---|---|---|---|
| **Aggregated separator (CombiCone)** | **1.0** | **1.62** | **0.981** (52/53) | **0.981** |
| Naive "average-the-combos" | 1.0 | 6.66 | 0.547 (29/53) | 0.830 |
| Magnitude-only (‖effect‖) | — | 49.5 | 0.019 (1/53) | — |

The aggregated separator recovers the true held-out gene at **median rank 1**,
**top-1 0.981** — 52 of 53 held-out genes return to rank 1 out of 105. The one
miss (UBASH3A) lands at rank 34. The naive baseline gets only 55% and has a long
tail (worst rank 67); the magnitude-only control is at chance.

**These numbers reproduce the repo's published
`results/screenloop_campaign.json → screens.norman.recovery` exactly** (median
1.0, mean 1.6226415…, top-1 0.9811320…), using the bundled `screenloop`
functions unchanged.

### Three controls (re-derived, checked against published JSON)

The campaign driver that produced the controls is not in the bundle, so
`scripts/reproduce_screen_design_recovery.py` re-derives them from their documented
definitions and verifies against the published values:

| Control | computed | published | verdict |
|---|---|---|---|
| magnitude-only top-1 | 0.018868 | 0.018868 | exact |
| Spearman(dominance, separator advantage) | −0.626584 | −0.626584 | exact |
| mean atom dominance | 0.695588 | 0.695588 | exact |
| separator top-1 where baseline fails | 0.958 (23/24) | 0.958 | exact |
| permutation-null z (top-1) | 71.9 | 72.5 | within 1 MC SE* |

*\*The permutation null is Monte-Carlo (2000 draws) and the driver's exact seed
is not bundled; both sit within one MC standard error of the analytic
expectation 1/105 = 0.0095, and the floored p-value (1/2001 ≈ 5×10⁻⁴) is
identical. Every **deterministic** control matches to 6 decimals — including the
decisive dominance–advantage Spearman, which confirms the re-derived definitions
are the exact ones used.*

### The decisive control: geometry, not magnitude

The **dominance-vs-advantage** control (right panel of the figure) is what rules
out the "the held-out gene just dominates its own combos" explanation. Define
*dominance* = mean cosine of the held-out atom to the combinations it appears in,
and *advantage* = (naive rank − separator rank). The Spearman correlation is
**ρ = −0.63**: the separator's edge over the naive baseline **grows** precisely
when the held-out gene is **least** dominant in its own combinations. The cone
geometry, not the effect magnitude, does the recovery.

The exemplar is **CBFA2T3** (dominance 0.34 — well below average): the naive
baseline buries it at **rank 67**, while the aggregated separator returns it to
**rank 1**.

![Held-out-gene recovery on the Norman substrate. (a) Rank of the true held-out gene under three rankers; only the aggregated cone separator recovers it at rank 1 (top-1 98%). (b) The separator's advantage over the naive baseline grows as the held-out gene becomes less dominant in its own combinations (Spearman ρ = −0.63) — the geometry does the work, not the magnitude.](figures/screen_design_recovery.png)

### Name-the-missing-axis, concretely

Running `nominate_atoms` directly on the CBFA2T3 hold-out
(`results/a3_walkthrough/name_missing_axis_demo.json`):

- Library = 104 atoms with **CBFA2T3 removed**.
- Combos that needed it: `FEV+CBFA2T3`, `POU3F2+CBFA2T3`, `PRDM1+CBFA2T3` — all
  3 fall outside the reduced cone (3/3 separating).
- Aggregated-separator nomination, **top-1 = CBFA2T3** (alignment +0.168); every
  other candidate projects ≈ 0 onto the separator by construction.
- The naive baseline instead nominates high-magnitude HOX genes
  (`HOXA13`, `HOXC13`, `HOXB9`, …) and puts the true axis at rank 67.

This is the decision-support payload: *"your library cannot reach the axis these
three combinations demand; the single perturbation that supplies it is
CBFA2T3."*

---

## Part 2 — End-to-end CLI walkthrough (real substrate)

All commands run against the bundled `combicone_substrate.npz`. Outputs are
captured in `results/a3_walkthrough/`.

### Step 1 — Ingest: parse the screen and report its structure

```bash
python combicone_cli.py ingest combicone_substrate.npz
```
```json
{
  "n_atoms": 105,
  "n_combos": 131,
  "n_genes": 5045,
  "singles_coverage": 1.0,
  "control_label": "ctrl",
  "separator": "+",
  "source": "combicone_substrate.npz",
  "from_npz": true
}
```
105 single-gene atoms, 131 measured doubles, 5045 genes, 100% single-coverage
(every double's constituents are measured as singles).

### Step 2 — Triage: rank unmeasured combinations by predicted emergence

```bash
python combicone_cli.py triage combicone_substrate.npz --measured-only --top 10 \
    -o results/a3_walkthrough/triage.csv
```
Top of the training-free `−agg_cos` ranking (what to run first):

| rank | combination | score |
|---|---|---|
| 1 | ETS2+TGFBR2 | +0.4731 |
| 2 | MAPK1+TGFBR2 | +0.4717 |
| 3 | CNN1+ETS2 | +0.4400 |
| 4 | OSR2+PTPN12 | +0.4275 |
| 5 | CNN1+MAPK1 | +0.4220 |

*Note the honest boundary reproduced elsewhere in the repo: in a full campaign
replay, this training-free triage does **not** beat the cheap `−cos` baseline at
wells-to-90%-discovery (Norman: triage 96 vs magnitude 104 vs cone-adaptive
120). Triage ranks; it does not claim to be the most sample-efficient policy.*

### Step 3 — Certify: which measured combinations are emergent above noise

```bash
python combicone_cli.py certify combicone_substrate.npz --method analytic \
    -o results/a3_walkthrough/certify.csv
```
`34/131 certified emergent (both bars)` — significance **and** a noise-floor
ratio ≥ 1.9×. The leaders are mechanism-coherent chromatin/MAPK synergies:

| combination | verdict | floor ratio | z |
|---|---|---|---|
| SET+CEBPE | certified emergent | 3.5× | 61.5 |
| IRF1+SET | certified emergent | 3.1× | 50.1 |
| MAPK1+PRTG | certified emergent | 2.6× | 36.8 |
| SET+KLF1 | certified emergent | 2.8× | 34.6 |
| TBX3+TBX2 | certified emergent | 2.6× | 34.5 |

The textbook DUSP9+MAPK1 phosphatase-on-substrate pair is flagged *emergent above
noise but modest* (1.7× floor, p = 2×10⁻⁴⁸) — real but below the 1.9× bar, which
is exactly the two-bar verdict working as designed.

*Boundary note: this two-bar count (34/131) is the CLI's default-threshold
verdict; the repo's `certificate_dossier.json → sensitivity` reports 35 certified
under a separate gamma-sweep definition. They are different computations and are
not equated here.*

### Step 4 — Recommend: the next diverse batch to run

```bash
python combicone_cli.py recommend combicone_substrate.npz --batch-size 10 \
    --strategy diversified
```
Training-free relevance (`−agg_cos`) × greedy max-marginal-relevance diversity
(`diversity_weight=0.5`), over 5329 candidate pairs:

| run order | combination | relevance | novelty |
|---|---|---|---|
| 1 | ARRDC3+SLC38A2 | 0.586 | 0.682 |
| 2 | MAP2K3+MEIS1 | 0.443 | 0.663 |
| 3 | ELMSAN1+ETS2 | 0.410 | 0.708 |
| … | … | … | … |
| 10 | MAML2+RUNX1T1 | 0.345 | 0.784 |

The recommender chooses what to **measure**; whether a discovery is genuinely
emergent is decided by `certify`, not by the acquisition score.

---

## Reproduce it yourself

```bash
# Part 1: held-out-gene recovery + controls, checked against published JSON
python scripts/reproduce_screen_design_recovery.py \
    --substrate combicone_substrate.npz \
    --published results/screenloop_campaign.json \
    --out results/a3_walkthrough/held_out_recovery.json

# Part 2: the four CLI stages
python combicone_cli.py ingest    combicone_substrate.npz
python combicone_cli.py triage    combicone_substrate.npz --measured-only --top 10 -o results/a3_walkthrough/triage.csv
python combicone_cli.py certify   combicone_substrate.npz --method analytic          -o results/a3_walkthrough/certify.csv
python combicone_cli.py recommend combicone_substrate.npz --batch-size 10 --strategy diversified
```

**Artifacts produced**
- `docs/figures/screen_design_recovery.png` — recovery rank distribution + dominance control
- `results/a3_walkthrough/held_out_recovery.json` — recovery metrics + controls (COMPUTED)
- `results/a3_walkthrough/name_missing_axis_demo.json` — the CBFA2T3 name-the-axis demo (COMPUTED)
- `results/a3_walkthrough/triage.csv`, `certify.csv` — CLI stage outputs
- `scripts/reproduce_screen_design_recovery.py` — one-command reproduction + verification harness
