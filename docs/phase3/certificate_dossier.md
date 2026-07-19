# Emergence-certificate trust dossier

*A causal-inference-style trust argument for the CombiCone emergence certificate,
computed on the real Norman substrate (105 single-gene atoms, 131 measured
doubles) and its shipped certificate. Reproduce with*
`python scripts/certificate_dossier.py --certificate <emergence_certificate.csv> --ledger docs/phase2/manuscript_numbers.json`.

The certificate makes one causal-flavored claim: *a measured combination lands
outside the non-negative cone of its constituent single-gene effects by more than
measurement noise can explain*. This document stress-tests that claim along the
three axes a skeptical reviewer would attack: **specificity** (does it fire when
it shouldn't?), **sensitivity** (how fragile is the one key assumption?), and
**construct validity** (does it measure genuine non-additivity?).

---

## 1. Negative controls — specificity

**NC-additive (the headline specificity test).** We built 150 synthetic
combinations as a literal sum of two randomly chosen measured single effects,
`D = a + b`, plus measurement-scale Gaussian noise. Such a combination is inside
the cone *by construction* (coefficients `[1, 1]`), so a correct certificate must
never call it emergent.

> **Result: 0 / 150 false-certified.** 141 were called "within measurement
> noise", 9 "above noise but modest", **0** "certified emergent". Mean floor
> ratio 1.01, 95th percentile 1.16 — the noise-floor bar (≥ 1.9×) is nowhere
> near cleared. False-positive rate **0.000**.

**NC-sham.** A pure-noise "combination" (control-scale Gaussian around zero)
scrapes bar (a) — it is statistically distinguishable from a bootstrap null — but
its floor ratio is ≈ 1.0, so bar (b) rejects it. This is the concrete
demonstration of *why both bars exist*: significance alone (bar a) can be fooled
by a vector that is directionally far from the cone but tiny; the noise-floor
ratio (bar b) cannot.

**NC-reachable (a deliberately imperfect control).** We took the 13
least-unreachable *real* doubles (bottom decile by cone residual) and re-certified
them. Verdicts: **3 certified, 10 "above noise but modest", 0 within-noise.** This
control does *not* pass cleanly — and that is informative, not a failure. A small
cone residual means the pair is *nearly* additive in effect *size*, but three of
these pairs still clear both bars because their measurement noise is tiny, so even
a small residual is many noise units away from the cone. This is the size-vs-
confidence distinction (§3) made concrete: low residual (small effect) does not
imply low z (low confidence). It is a caution for users — "certified emergent" is
a noise-relative confidence statement, and a certified pair with a small residual
is a small-but-reliably-nonzero departure, not necessarily a large one. Read the
residual alongside the verdict.

---

## 2. Sensitivity radius — how fragile is the noise assumption?

The certificate's one load-bearing assumption is that the split-half estimate
`SE = |t₁ − t₂| / 2` faithfully captures measurement noise. If we *underestimated*
that noise, a pair could look more emergent than it is. The **sensitivity radius
Γ\*** is the smallest inflation factor such that treating the true noise as
`Γ\* × SE` drops a pair below the "certified" verdict — the certificate's analogue
of a Rosenbaum sensitivity bound / E-value.

> **Result (35 certified pairs): median Γ\* = 1.25×, range 1.25×–2.0×.** The
> flagship pair SET+CEBPE tolerates ~1.9× noise inflation (floor 3.6× → 1.9× at
> Γ ≈ 2.0). Only 3% of certified pairs survive a full 2× noise misestimation.

**Honest reading.** Γ\* is a *monotone re-expression of the bar-b margin*
(empirically Γ\* ≈ floor_ratio / 1.9; certified floor/1.9 spans 1.00–1.83, median
1.13), not an independent axis. It does not manufacture new robustness — it
translates "how far above the 1.9× floor is this pair" into interpretable
"tolerable noise-misestimation" units. The takeaway is calibration, not comfort:
most certified pairs sit close to the floor and would not survive a doubling of
the noise estimate, so the two-bar verdict is doing real work and should be read
as "emergent under the stated noise model", not "emergent, full stop". The
strongest pairs (SET+CEBPE and the other top-decile) carry the widest margins.

---

## 3. Construct validity — does it measure non-additivity?

Three independent lines of evidence that the geometry captures genuine synergy:

1. **Recovers classical synergy.** The cone residual (unreachable fraction)
   correlates with the raw non-additivity `‖D − (a+b)‖ / ‖D‖` at Spearman **0.74**,
   and at partial Spearman **0.62** after controlling for effect magnitude. The
   synergy signal is not a magnitude artefact.

2. **The certificate removes the magnitude confound.** The *raw* residual is
   strongly, spuriously anti-correlated with effect magnitude (Spearman **−0.56**):
   big-effect pairs look artificially unreachable. The noise-robust z is
   decorrelated from magnitude (**+0.14**) — the confound is removed by design.

3. **Top certified pairs are known interactions.** The strongest certified pairs —
   SET+CEBPE, IRF1+SET, MAPK1+PRTG, CEBPE+RUNX1T1 — are dominated by the
   myeloid/erythroid master-regulator SET and the CEBPE differentiation axis,
   consistent with Norman et al.'s reported strong genetic interactions.

### The key subtlety: size vs. confidence

The dossier surfaces a distinction the certificate's users must respect:

| quantity | what it measures | correlates with |
|---|---|---|
| **cone residual** | *effect size* of emergence | synergy (partial ρ = 0.62 \| magnitude); z only ρ = 0.06 |
| **z / floor ratio** | *noise-relative confidence* | noise floor (ρ = −0.59); magnitude ρ = 0.14 |

The z-score is **almost uncorrelated with the cone residual (ρ = 0.06)** and is
instead dominated by the per-pair noise floor (ρ = −0.59). This is exactly what
makes it decorrelate from magnitude — but it means **z must not be read as "how
large the synergy is."** Effect size lives in the residual; z answers the
separate question "how confident are we this is above noise." Reporting both is
the honest practice, and CombiCone's two-bar verdict does precisely that: bar (a)
is the confidence test, bar (b) is an effect-size floor.

---

## What this does and does not license

- **Does:** support the claim that a certified pair is unreachable from its
  single-gene cone beyond the stated measurement-noise model, with a quantified
  fragility (Γ\*) and a zero observed false-positive rate on additive controls.
- **Does not:** assert biological impossibility of additivity (scope is
  model-relative), attach a calibrated probability to any single pair, or claim
  robustness to a 2× noise misestimation for the median pair. "Certified emergent"
  means *emergent under this noise model, at this margin* — the margin is in the
  certificate, and this dossier makes it legible.

See `fig_dossier.png` for the three-panel summary and `certificate_dossier.json`
for every number.
