# Cell-State Reachability

### Can knockdown point a cell where you want it to go?

**Built with Claude — Life Sciences Hackathon · Research / Lab Track**

![Cell-State Reachability at a glance — the reframe; the Th2→Th1 flagship verdict (held-out cosine 0.448, above all 60 shuffled targets, KKT optimality violation 1.1×10⁻¹¹); the signed 39/25/35 LOF/GOF-proxy/neither decomposition; and the four outputs.](docs/figures/fig_at_a_glance.png)

> A GPS is useful for two reasons: it finds a route, and it tells you when the road you
> need is not on the map. This project brings that second answer to cell engineering.

Most computational biology tools predict what a perturbation might do or rank genes that
look promising. This project asks an earlier question: **given the effects a screen has
actually measured, can non-negative combinations of this intervention class point the
transcriptome toward the desired state?**

The answer is deliberately narrower than “will this cure disease?” It is a falsifiable,
screen-relative test of directional feasibility—and a way to avoid asking a knockdown
screen to do a job that appears to require another modality.

---

## The idea, without the math

Each CRISPRi perturbation leaves a measured transcriptional fingerprint: some genes rise,
others fall. You can apply that knockdown effect, or combine it with others, but you cannot
reverse the entire fingerprint and call it a “negative knockdown.” Under an explicit
first-order additivity approximation, all non-negative combinations form a geometric cone.

The target state is an arrow. The method asks how closely the cone can point along it.

It returns four things:

1. **A directional verdict** calibrated on held-out genes and shuffled targets.
2. **A ranked mixture and greedy sparse panel** for follow-up—not a proof of the globally
   smallest gene set.
3. **A dual certificate** proving when the full target direction lies outside the measured
   cone.
4. **A ranked list of unmet readouts** that can motivate CRISPRa or de-repression tests.
   Those genes are hypotheses, not validated activation targets.

The source H5AD contains **33,983 perturbation–condition profiles across 10,282 readout
genes** and 11,526 distinct targeted genes. After the source screen's quality filters, the
resting-cell analysis uses **6,871 perturbation generators** and a 6,188-gene target
signature.

## What it found

For the flagship **Th2 → Th1** direction in resting primary human CD4⁺ T cells:

- The fixed-split held-out cosine is **0.448**, versus 0.627 in-sample.
- It exceeds all **60** shuffled targets (plus-one empirical **p = 1/61**; descriptive
  z ≈ 24).
- The staged modality proxy assigns **39%** of target energy to measured knockdown
  directions, **25%** to sign-flipped gain-of-function proxy directions, and **35%** to
  neither.
- The KKT optimality violation is **1.1 × 10⁻¹¹**. That certifies the numerical cone
  projection; it does not certify biological efficacy.

Known regulators behave in the expected direction: GATA3, which should be reduced when
moving away from Th2, ranks near the knockdown-aligned end; TBX21, which should rise for
Th1, is anti-aligned when knocked down. These are positive controls, not discoveries.

Across 12 transition–condition cases, knockdown is never the majority component (mean LOF
fraction **0.34**). Each case exceeds its eight shuffled screening controls, but those
eight-shuffle z values are prioritization evidence—not precise permutation p-values.

The same optimization code also runs without retuning on a held-out CEBPA state in a K562
CRISPRa screen (cosine **0.878**). More importantly, that dataset contains 126 measured
double perturbations: the median cosine between a measured double and the sum of its
singles is **0.71**. Additivity is useful, bounded, and absolutely not guaranteed.

## Why this is interesting

The novelty is not another regulator list. It is a different kind of output:

> **“Not with this measured intervention dictionary, in this direction, under these
> assumptions—and here is the separating direction that proves it.”**

A structured survey of **91 prior methods** found no prior entry combining a
target-specific feasibility certificate with directly measured perturbation effects. That
is a survey-bounded novelty claim, not a claim that convex cones, network control, or
cell-state engineering are new.

The practical workflow has two levels:

- At the **state level**, choose whether to test a focused CRISPRi panel, add a CRISPRa or
  de-repression arm, or seek a better dictionary because the result does not clear its null.
- At the **candidate level**, annotate the top greedy LOF panels for tractability and human
  genetics. Their union contains **102 unique candidate knockdown nodes**: 45 lack a
  conventional drug handle and 10 have a clinical-grade drug in the saved Open Targets
  snapshot.

This does not explain clinical attrition or validate a drug target. It can expose a
modality mismatch before a costly combination experiment.

## Scope that stays attached to every claim

- The verdict is relative to one measured dictionary and a transcriptomic target direction.
- Matching expression is not functional rescue.
- Multi-gene recipes extrapolate from single perturbations and require wet-lab validation.
- The GOF fraction uses `-E` as an activation proxy; CRISPRa is not generally the exact
  mirror of CRISPRi.
- The primary screen is one CD4⁺ T-cell system across four donors; donor-collapsed effects
  prevent true leave-one-donor-out validation.
- Every nomination is a hypothesis for testing, not a validated target.

---

## Start here

| If you want… | Open |
|---|---|
| **The 3-minute demo** | [Watch on YouTube](https://www.youtube.com/watch?v=GJbxLxYUBMo) |
| **The guided interactive story** | [Open live](https://raw.githack.com/MasalaKimchi/cell-state-reachability/5b96d535ba8e33e63c43f69b4a89b744ce7699e7/app/index.html) · [`app/index.html`](app/index.html) |
| **The paper** | [`manuscript/main.pdf`](manuscript/main.pdf) |
| **The validation report** | [`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md) |
| **The full technical record** | [`docs/Technical_Dossier.pdf`](docs/Technical_Dossier.pdf) |
| **Core method** | [`reachability.py`](reachability.py) |
| **Verify the software** | [`reproduce.sh`](reproduce.sh) |
| **Re-run the data analyses** | [`README.md`](README.md#reproduce-the-analysis) |

**Dataset:** Zhu et al. 2025, genome-scale CRISPRi Perturb-seq in primary human CD4⁺
T cells (Marson and Pritchard labs; CZI Virtual Cells Platform, MIT license).
