# I2 — A design-based causal reading of the emergence certificate

**Methods note. CombiCone: certified triage for combinatorial perturbation screens.**

*Scope of this note.* This is a **descriptive / design-only** formalization with **one small
numeric illustration computed on real substrate data** (`combicone_substrate.npz`, the Norman
combinatorial CRISPRa screen; file label A549, canonical line K562). It re-reads the existing
convex-cone emergence certificate as a design-based counterfactual query, states the identifying
assumptions the reading requires, contrasts them with the implicit assumptions of the 2025–2026
"causal" perturbation models, and bridges the certificate's residual to the uplift /
heterogeneous-treatment-effect (HTE) literature. It introduces **no new method and no new
performance claim.** Every quantity is tagged **[COMPUTED]** (run here on real Norman effect
vectors), **[DOSSIER]** (read from `results/certificate_dossier.json`), **[MANUSCRIPT]** (from
`manuscript/sections/50_discussion.tex`), or **[DESIGN]** (a definition or assumption, not a
measurement). The note stays inside the repo's stated `claim_boundary`: emergence is **relative to
the measured dictionary under an additive model**, never biological synergy in an absolute sense.

---

## 1. Setup: single-gene effects as average treatment effects

### 1.1 The estimand

A pooled perturbation screen (Perturb-seq / CRISPRa / CRISPRi) assigns each cell a genetic
perturbation and reads out a high-dimensional cell state — here a length-`G` transcriptome
(`G = 5045` highly-variable genes in the substrate **[COMPUTED]**). Write `Y ∈ ℝ^G` for the
measured cell-state vector and, for each single-gene perturbation `g ∈ {1,…,N}`, let `do(g)`
denote the intervention "deliver the guide(s) for `g`". Using potential outcomes, the **single-gene
average treatment effect (ATE)** relative to the non-targeting control arm is

```
a_g  :=  E[ Y | do(g) ]  −  E[ Y | do(∅) ] .
```

`E[Y | do(∅)]` is the control mean (`ctrl` in the substrate); `E[Y | do(g)]` is the perturbation
pseudobulk mean (`means[g]`). In the substrate this is exactly the effect atom: for the first
single, `atoms[0]` equals `means["AHR+ctrl"] − ctrl` to machine precision (max abs difference
`0.0`, correlation `1.0`) **[COMPUTED]**. So **the CombiCone "effect atoms" are, by construction,
the estimated single-gene ATEs of the screen** — `atoms` is a stack of `N = 105` ATE vectors
`{a_g}` **[COMPUTED]**.

For a measured **combination** `S = {g₁,…,g_k}` (e.g. a Norman double), define the combination ATE

```
τ_S  :=  E[ Y | do(S) ]  −  E[ Y | do(∅) ] ,
```

again the pseudobulk contrast `means[S] − ctrl`.

### 1.2 What makes these ATEs (identification)

The contrast `a_g` estimates a causal effect — rather than a mere association — only under the
screen's **design**. Pooled screens are the closest thing functional genomics has to a randomized
experiment: guide delivery is by random viral integration at low multiplicity, so which cell
receives which perturbation is (to a good approximation) independent of the cell's baseline state.
Formally we need:

- **(R) Randomization / ignorability.** Guide assignment `⊥ {Y(do(g))}_g` given the modeled
  covariates (library batch, cell-cycle, depth). Under (R), `E[Y|do(g)] − E[Y|do(∅)]` is the ATE,
  not a confounded difference of populations. This is the assumption a pooled screen is *designed*
  to satisfy; it is the reason a perturbation screen is a legitimate causal experiment and an
  observational expression atlas is not.
- **(C) Consistency / well-defined intervention.** The guide for `g` realizes a single, stable
  version of `do(g)` (one knockdown/activation level), so `Y = Y(do(g))` for treated cells.

Under (R)+(C) each `a_g` is an unbiased design-based ATE. These two assumptions are shared by any
honest analysis of the screen; §3 lists the *additional* assumptions the emergence query needs.

---

## 2. Emergence as a design-based counterfactual

### 2.1 The additive counterfactual set (the cone)

Fix the dictionary of measured single-gene ATEs `A = {a_1,…,a_N}`. The **additive
counterfactual** for a hypothetical dosing of the singles at non-negative intensities
`c = (c_1,…,c_N) ≥ 0` is the mixture

```
μ(c)  =  Σ_g c_g · a_g ,      c_g ≥ 0 .
```

Non-negativity is not a modeling convenience — it is the design constraint that **a perturbation
cannot be delivered in negative amount** (you cannot "un-knock-down" a gene) **[MANUSCRIPT, L4]**.
The set of all such additive counterfactuals is exactly the **non-negative conic hull** of the
ATEs:

```
Cone(A)  =  { Σ_g c_g a_g : c_g ≥ 0 }  ⊆  ℝ^G .
```

`Cone(A)` is precisely "every cell state additively reachable, at any non-negative dose, by the
single-gene interventions the screen measured". It is the counterfactual comparison class for the
emergence query.

### 2.2 The counterfactual query

> **Definition (design-based emergence).** A measured combination `S` is **emergent relative to
> `A`** iff its combination ATE `τ_S` lies **outside** `Cone(A)` — i.e. no non-negative mixture of
> the single-gene ATEs reproduces the combination's measured causal effect:
> ```
> τ_S  ∉  Cone(A)      ⟺      min_{c ≥ 0}  ‖ τ_S − Σ_g c_g a_g ‖_Q  >  0 .
> ```

This is a **counterfactual contrast**: it compares the *factual* combination effect `τ_S` against
the *best additive counterfactual* — the closest state the single-agent interventions could have
produced acting together additively. Emergence is the statement "the joint intervention did
something the sum of its parts, at any non-negative doses, provably could not." **[DESIGN]**

The comparison is operationalized by `reachability.project_cone`, which solves the non-negative
least-squares (NNLS) program `min_{c ≥ 0} ‖√Q (τ_S − Σ c_g a_g)‖²` and returns:

- `coefficients` `ĉ` — the optimal non-negative dosing (the **realized additive counterfactual**);
- `fitted` `μ(ĉ)` — the best additive-reachable state;
- `residual` `τ_S − μ(ĉ)` — **the part of the combination effect no additive counterfactual
  explains**;
- `residual_fraction` `= ‖residual‖ / max(‖τ_S‖, ‖fitted‖)` — the normalized uplift (§5);
- `dual_separator` — when `τ_S ∉ Cone(A)`, a **Farkas certificate**: a hyperplane `h` with
  `⟨h, a_g⟩ ≤ 0 ∀g` and `⟨h, τ_S⟩ > 0`, proving separation. This is the "certificate" — a
  checkable *witness* of the counterfactual claim, not a prediction.

The KKT conditions of this convex program are verified to `≤ 1e-8` (`kkt_violation`), so the
counterfactual is solved to optimality, not approximated **[COMPUTED; the worked example clears it
at 3.14e-15]**.

### 2.3 Why this is the right comparison class

Two properties make `Cone(A)` the honest null for "additively reachable", and both are
**[COMPUTED]** on the substrate in §5–§6:

1. **Monotone in the reference set.** Adding atoms to `A` can only enlarge the cone, so the
   residual can only shrink. Emergence is therefore a property *of a stated reference*, and the
   method forces that reference to be named — enlarging singles → singles+doubles shrinks the
   certified set from 40 to 16 **[MANUSCRIPT]**. In §6 we show the residual against the pair's own
   two ATEs is an upper bound on the residual against all 105 (true for 131/131 doubles).
2. **Noise-gated.** Geometry alone always finds *some* residual; the certificate is trusted only
   after a split-half noise-injection null shows the residual exceeds what measurement noise
   produces (§4). This is what separates "outside the cone" from "outside the cone by more than
   noise".

---

## 3. Identifying assumptions, stated explicitly

The value of reading the certificate causally is that each identifying assumption becomes a
**named, checkable** condition rather than an implicit one. Below, each assumption is paired with
*what it buys*, *how it can fail*, and *the concrete robustness check* the repo already implements
or specifies. This is the disciplined version of the "causal-inference agenda" paragraph in
`50_discussion.tex` **[MANUSCRIPT]**.

**A1 — Randomization / ignorability of guide assignment (R).** *Buys:* the atoms `a_g` are ATEs,
so the cone is a cone of causal effects, not of correlational signatures. *Fails if:* guide
assignment correlates with baseline state (selection during the screen, growth bias of certain
perturbations). *Check:* pooled low-MOI design is the mitigation by construction; residual
selection is bounded by comparing perturbation abundance pre/post and modeling it as a covariate.

**A2 — SUTVA / no interference between cells.** *Buys:* one cell's ATE does not depend on which
perturbations *other* cells received, so pseudobulk pooling estimates a well-defined `a_g`.
*Fails if:* a perturbed cell's secreted signals alter a neighbour's state in pooled culture —
paracrine spillover — which would bias the very atoms the cone is built from. *Check:* an
arrayed-vs-pooled or MOI-titration comparison bounds the spillover **[MANUSCRIPT]**. This is the
single most consequential assumption for a *secreted-signaling* screen and the least testable
without new data; we flag it as **[DESIGN, wet-lab]** and do not claim it is satisfied.

**A3 — Exclusion / no direct guide effect.** *Buys:* the guide changes `Y` only through the
intended gene, so `a_g` is attributable to `do(g)`. *Fails if:* off-target editing or a generic
dCas9/CRISPRa burden shifts `Y` independent of the target. *Check:* multi-guide concordance
(≥2 guides per gene agreeing) is the standard mitigation; the substrate's split-half stability
(A5) is a partial internal proxy.

**A4 — Compliance / errors-in-variables (guide non-compliance).** *Buys:* the *nominal*
intervention `do(g)` equals the *realized* one, so `a_g` is the full-compliance effect. *Fails if:*
incomplete knockdown/activation attenuates the measured effect — a classic errors-in-variables
problem: the observed atom is `a_g · (realized efficiency)`, a shrunk version of the true ATE. A
cone built from *attenuated* atoms is **smaller** than the true additive-reachable set, which makes
emergence calls **anti-conservative** (a combination can appear outside a shrunken cone merely
because its parts were under-delivered). *Check:* an instrumental-variables treatment using
measured guide efficiency (e.g. residual target expression) as the first stage would rescale atoms
to full-compliance ATEs **[MANUSCRIPT]**. This is stated, not performed here — **[DESIGN]**.

**A5 — Measurement model: undirected noise vs coordinated bias.** *Buys:* the split-half
noise-injection null (§4) correctly calibrates "how much residual is just noise". *Fails if:* the
noise is not per-gene independent but **coordinated** — a systematic attenuation of the
combination-specific signal — which the isotropic-per-gene null does not model. *Check:* the
dossier reports a per-pair **sensitivity radius `Γ*`** (the certificate's analogue of a Rosenbaum
bound): the multiple of the noise floor at which the verdict flips. Median `Γ* = 1.25×`, max
`2.0×`, and only 2.9% of certified pairs survive a `2×` coordinated inflation **[DOSSIER]** — an
honest statement that most verdicts are robust to modest, but not large, coordinated bias.

**A6 — Additivity of the reference model.** *Buys:* `Cone(A)` (non-negative *linear* mixtures) is
the correct "no-interaction" counterfactual set, so a residual means genuine super-additive
structure. *Fails if:* the biologically correct null is non-linear (e.g. Bliss/HSA multiplicative
independence), in which case some cone residual reflects the wrong null rather than epistasis. This
is the repo's central declared boundary: emergence is defined **under an additive model relative to
the measured dictionary** **[MANUSCRIPT, L2]**. Additivity is calibrated only out-of-domain
(median single-vs-double cosine ≈ 0.71 on K562), and an in-domain double-knockdown panel is the
specified wet-lab test.

**Construct validity (the negative-control side of identification).** A causal test is only
credible if it does *not* fire on effects that are additive by construction. The dossier's
negative-control check synthesizes additive combinations from measured single ATEs and certifies
**0** of them: false-positive rate `0.0` over `100` additive negative-control pairs (verdict
counts: 0 certified / 6 modest / 94 within-noise) **[DOSSIER]**. *(Note: `50_discussion.tex` cites
"150 additive combinations"; the dossier JSON released with the code records `n_pairs = 100`. The
computed FPR of 0.0 is identical either way; we quote the JSON as the authoritative artifact and
flag the text discrepancy for correction.)* Positive construct validity holds too: the noise-robust
residual recovers classical non-additivity at partial Spearman **0.62** after controlling for
magnitude **[DOSSIER]**, so the pairs the test certifies are the biologically synergistic ones, not
low-magnitude noise.

### 3.1 Contrast with the 2025–2026 "causal" perturbation models

The field has taken a visible "causal turn," but in most recent models the identifying assumptions
above are **implicit and unstated**, and "causal" describes the *architecture* rather than a
*design-based estimand*:

- **PDGrapher** (Gonzalez et al., *Nature Biomedical Engineering*, vol. 10, 2026; online 2025,
  DOI 10.1038/s41551-025-01481-x) predicts combinatorial therapeutic perturbations with a "causally
  inspired" neural network. It is genuinely causal in flavour — verbatim, it "uses protein–protein
  interaction (PPI) networks or gene regulatory networks (GRNs) as approximations of the causal
  graph, operating under the assumption of no unobserved confounders," and optimizes for
  interventions that shift a diseased cell toward a treated state. But the causal graph is a
  *borrowed prior* (an approximation, in the authors' own words), and the no-unobserved-confounder
  assumption is imposed on that network rather than earned from a randomized design; the output is a
  *predicted* target set, not a falsifiable verdict with a witness.
- **X-Cell** (bioRxiv 2026) is titled "Scaling **Causal** Perturbation Prediction across Diverse
  Cellular Contexts via Diffusion Language Models" — a large diffusion model. Here
  "causal" is essentially a branding of scaled *prediction*; the model emits a predicted expression
  response, and its identifying assumptions (what makes the predicted response a causal effect
  rather than a fit) are not separated from its training objective.
- **scGPT / foundation-model causal-validation** work (the 2025–2026 benchmark wave) evaluates
  whether large pretrained single-cell models predict held-out perturbation responses. The recurring
  finding is a *predictive* one — foundation models often match simple additive or mean baselines on
  held-out double effects — which is exactly the axis CombiCone argues is the wrong one for a
  *decision* about emergence **[MANUSCRIPT]**.

The distinction is not that these models are wrong; it is **where the assumptions live**. CombiCone
reads the *screen's own randomized design* as the identification source, makes A1–A6 explicit, and
returns a *certificate* (a Farkas witness plus a noise-gated verdict) that is falsifiable on both
sides. The neural "causal" models locate causality in a **learned or borrowed graph / architecture**
and return a **point prediction**, leaving the estimand and its assumptions implicit. A screening
team using CombiCone knows *exactly* what counterfactual it is testing and under what conditions the
answer holds; the same team using a forward predictor gets a number whose causal warrant is not
separable from the model fit.

---

## 4. The noise gate as a sharp null on the counterfactual

Design-based emergence (§2.2) is an idealized, noiseless statement. In practice `τ_S` and every
`a_g` are estimated with finite cells, so a strictly-positive residual is guaranteed even for a
truly additive combination. The certificate therefore tests the *sharp null*

> **H₀:** `τ_S ∈ Cone(A)` (the combination is additively reachable), and the observed residual is
> entirely measurement noise.

The null distribution is built **from the data's own noise**: `certify_emergence` projects `τ_S`
onto the cone to get the best reachable point `f₀ = μ(ĉ)`, then repeatedly adds fresh per-gene
Gaussian noise of scale `σ` (the split-half standard error `σ = |t₁ − t₂| / 2`, from `means1`,
`means2`) to `f₀` and re-projects, forming the residual distribution a *genuinely reachable* target
of the same magnitude would show. The reported statistics are:

- `p_value` — `P(noise-only residual ≥ observed)`, a conservative plus-one empirical p (or the
  closed-form generalized-χ² analytic null, which is conservative by construction);
- `z` — `(residual − null_mean) / null_sd`, a magnitude-robust confidence;
- `floor_ratio` — `residual / null_mean`, an effect-size axis; the two-bar verdict certifies only
  when `p < α` **and** `floor_ratio ≥ 1.9`.

Reading this causally: the noise null is the counterfactual query's **estimator**, and the two bars
separate *statistical* significance of the counterfactual departure (`p`) from its *effect size*
(`floor_ratio`). A point forward-predictor emits neither — it cannot supply a null, so it cannot
say whether a predicted departure is distinguishable from noise **[MANUSCRIPT]**.

---

## 5. Bridge to uplift / heterogeneous treatment effects (HTE)

### 5.1 Emergence *is* uplift beyond additive single-agent uplift

The uplift-modeling / HTE literature (Athey–Imbens causal forests; Künzel et al. meta-learners;
the marketing "uplift" tradition) studies the **incremental** effect of a treatment beyond a
baseline — the part of an outcome attributable to the treatment over and above what a reference
policy would have produced. The emergence certificate is exactly an uplift statement, lifted from a
scalar outcome to a cell-state vector:

| HTE / uplift concept | CombiCone counterpart |
|---|---|
| Treatment effect (ATE) | Single-gene effect `a_g = E[Y\|do(g)] − E[Y\|do(∅)]` |
| Baseline / reference policy | Best non-negative additive mixture `μ(ĉ) = Σ ĉ_g a_g` |
| **Uplift** (incremental effect over baseline) | **Residual** `τ_S − μ(ĉ)` — the combination's effect the additive policy cannot reproduce |
| Uplift magnitude | `residual_fraction` (normalized ‖uplift‖) |
| Significant, non-trivial uplift | Two-bar verdict: `p < α` **and** `floor_ratio ≥ 1.9` |
| Effect heterogeneity / interaction | Super-additivity: `τ_S ≠ Σ a_g` — the reason uplift is non-zero |

Concretely: additive reachability is the "no-interaction" counterfactual, so the residual is the
**combination-specific uplift** — the interaction / heterogeneous component of the joint
intervention's effect. "Emergent" = "the combination has uplift beyond the best additive
single-agent policy, by more than measurement noise." The Farkas separator is then the *direction*
of that uplift in gene space: the set of readouts the additive policy under-delivers (`⟨h, τ_S⟩ > 0`
while `⟨h, a_g⟩ ≤ 0` for every single). This is the vector-valued analogue of an uplift score's
sign, and it is what the manuscript calls the "unmet-readout" certificate.

### 5.2 What the bridge does and does not license

- **Does:** it places emergence in a well-understood inferential frame — a *contrast against a
  policy baseline*, with a null and an effect size — and clarifies that the certificate is an
  interaction/HTE object, not a prediction. It also imports the right caution: uplift estimates are
  notoriously sensitive to baseline mis-specification, which here is the additivity assumption (A6)
  and the completeness of the dictionary `A`.
- **Does not:** it does not turn the gene-level uplift *direction* into a validated activation list.
  Ranking the coordinates of the separator is a *biological hypothesis layer* the manuscript
  explicitly flags as **never bench-tested** (limitation L1); §6 reports the uplift direction's
  concentration but makes **no claim** that raising those genes produces the state **[MANUSCRIPT,
  L1]**.

---

## 6. Worked numeric example on the Norman substrate (SET + CEBPE)

All numbers in this section are **[COMPUTED]** by running the repo's own
`reachability.project_cone` and `combicone.certify_emergence` on real Norman effect vectors from
`combicone_substrate.npz`. Reproduce with `scripts/causal_formalization_example.py`. We use the dossier's
flagship certified pair, **SET + CEBPE** (`Γ* = 2.0×`, the most robust certified pair) **[DOSSIER]**.

### 6.1 The three vectors

- `a_SET   = means["SET+ctrl"]   − ctrl`, the SET ATE, `‖a_SET‖   = 5.109`
- `a_CEBPE = means["CEBPE+ctrl"] − ctrl`, the CEBPE ATE, `‖a_CEBPE‖ = 7.315`
- `τ       = means["SET+CEBPE"]  − ctrl`, the double ATE, `‖τ‖      = 4.595`

The two single ATEs are nearly orthogonal: `cos(a_SET, a_CEBPE) = +0.170` — dissimilar singles,
which is exactly the training-free triage flag for likely emergence **[COMPUTED]**.

### 6.2 The additive counterfactual (best non-negative 2-agent mixture)

Projecting `τ` onto `Cone({a_SET, a_CEBPE})` — the two-agent additive-reachable set — gives the
optimal non-negative dosing and the uplift:

```
ĉ            = (0.584, 0.288)          # non-negative doses on (SET, CEBPE)
μ(ĉ)         = 0.584·a_SET + 0.288·a_CEBPE   # best additive counterfactual
cos(τ, μ(ĉ)) = 0.856                    # the additive policy captures most, not all, of the double
uplift       = τ − μ(ĉ),   ‖uplift‖ = 2.373
uplift_fraction (residual_fraction) = 0.516
```

So **51.6% of the SET+CEBPE double effect (by normalized residual) is uplift the best additive
mixture of its own two single ATEs cannot reproduce** **[COMPUTED]**. For contrast, the naive
equal-dose Bliss-style sum `a_SET + a_CEBPE` (fixing `c = 1, 1`) aligns with the double at
`cos = 0.820` — worse than the NNLS-optimal mixture (0.856), confirming the NNLS baseline is the
*strongest* additive counterfactual, so the uplift is not an artifact of a lazy baseline.

### 6.3 The full CombiCone certificate (cone = all 105 single ATEs)

The deployed certificate compares `τ` against the whole dictionary, not just the pair's two atoms:

```
unreachable_fraction (vs all 105 ATEs) = 0.498     # ≤ the 2-agent 0.516, as it must be
n atoms with nonzero weight            = 12        # 12 of 105 singles enter the optimal mixture
geometry_status                        = outside_model_cone
KKT violation                          = 3.14e-15  # counterfactual solved to optimality
noise-aware z                          = 63.8
floor_ratio                            = 3.59      # clears the 1.9× effect-size bar
p_value (MC, split-half noise)         = 0.00498
verdict                                = certified emergent (p=0.00498, 3.6× noise floor)
```

**The decisive honesty point:** enlarging the reference set from the pair's own **2** ATEs to all
**105** barely moves the residual (0.516 → 0.498). Even given the entire measured dictionary and
free to use 12 different single-gene interventions at any non-negative doses, no additive policy
reaches within 50% of what SET+CEBPE jointly produced. Emergence here is **not** an artifact of an
impoverished two-atom basis — it survives the richest additive counterfactual the screen can build
**[COMPUTED]**.

### 6.4 The uplift direction (reported, not validated)

The uplift vector `τ − μ(ĉ)` is concentrated: **23 genes carry 50%** and **136 carry 80%** of its
squared norm (of 5045), with the single largest coordinate holding 11.6% of the residual energy
**[COMPUTED]**. This is the "unmet-readout" direction the Farkas separator certifies. Per limitation
L1, we deliberately **do not** name these genes as an activation list — the coordinate ranking is a
biological hypothesis with no bench test in this work **[MANUSCRIPT, L1]**.

### 6.5 The pattern holds across the dossier's top certified pairs

Running the same two-agent uplift decomposition on the other top-certified pairs (real Norman
doubles) shows the uplift is a general, sizeable fraction of the double effect — and tracks the
sensitivity radius `Γ*` **[COMPUTED / DOSSIER]**:

| pair | cos(A,B) | NNLS dose (c_A, c_B) | cos(double, mix) | 2-agent uplift frac | Γ* [DOSSIER] |
|---|---|---|---|---|---|
| SET+CEBPE | +0.170 | (0.584, 0.288) | 0.856 | 0.516 | 2.00× |
| IRF1+SET | +0.195 | (0.851, 0.817) | 0.929 | 0.369 | 1.75× |
| MAPK1+PRTG | -0.284 | (0.491, 0.629) | 0.819 | 0.574 | 1.50× |
| CEBPE+RUNX1T1 | +0.415 | (0.507, 0.466) | 0.954 | 0.301 | 1.50× |
| TBX3+TBX2 | +0.261 | (0.730, 0.768) | 0.971 | 0.239 | 1.50× |
| ETS2+IKZF3 | -0.164 | (0.504, 0.924) | 0.968 | 0.251 | 1.25× |
| CEBPE+KLF1 | -0.227 | (0.261, 0.771) | 0.890 | 0.457 | 1.25× |

### 6.6 The uplift interpretation is faithful to the deployed certificate

Across **all 131 measured Norman doubles**, the simple two-agent uplift fraction and the full-105-ATE
certificate residual agree tightly: **Spearman ρ = 0.880** (`p = 1.3e-43`), the full-cone residual
is `≤` the two-agent value for **131/131** pairs (the monotone-in-reference property, verified), and
enlarging the cone from 2 → 105 atoms reduces the residual by a mean of only **0.070 per pair** (the
median falls from **0.388** to **0.313**, a drop of **0.075**) **[COMPUTED]**. The uplift reading in §5 is therefore
not a simplification that loses the deployed method — it *is* the deployed method, read as an
interaction effect against an additive policy baseline.

---

## 7. Summary and boundary

- **Estimand.** Single-gene effect atoms are design-based ATEs from a randomized pooled screen;
  their non-negative conic hull `Cone(A)` is the set of additively-reachable counterfactual
  cell states.
- **Counterfactual query.** A combination is *emergent* iff its measured combination ATE lies
  outside `Cone(A)` — no non-negative mixture of the single-gene ATEs reproduces it — witnessed by
  a Farkas separator and gated by a split-half noise null (two-bar verdict).
- **Assumptions, explicit.** Randomization/ignorability (A1), SUTVA/no-interference (A2), exclusion
  (A3), compliance/errors-in-variables (A4), an undirected-noise measurement model quantified by the
  sensitivity radius `Γ*` (A5), and additivity of the reference (A6). Each is a named, checkable
  condition; the neural "causal" models (PDGrapher, X-Cell, scGPT-family validation) leave the
  analogous assumptions implicit and locate causality in a learned/borrowed graph rather than the
  design.
- **HTE bridge.** Emergence = *uplift beyond the best additive single-agent policy*: the cone
  residual is the combination-specific interaction effect, the two-bar verdict is significant
  non-trivial uplift, and the separator is the uplift's direction in gene space.
- **Worked example.** SET+CEBPE: 51.6% of the double effect is uplift the best two-agent additive
  mixture cannot reach; 49.8% against the full 105-ATE dictionary; certified emergent
  (z = 63.8, floor 3.59×, p = 0.005). Across 131 doubles the uplift reading and the deployed
  certificate agree at Spearman 0.880.

**Boundary (binding).** Everything here is *model-relative* and *additive-model-relative*:
"emergent" means "outside the non-negative cone of the **measured** single-gene ATEs under an
**additive** reference", never biological synergy in an absolute sense. The design-based reading
does **not** add evidence for the gene-level uplift *direction* (limitation L1, never bench-tested),
does **not** verify SUTVA (A2) or full compliance (A4) on this data, and inherits the additivity
calibration caveat (A6, calibrated only out-of-domain on K562). This note is a **formalization and
one illustration**, not a new result or a strengthened claim.

---

### Provenance

| Tag | Source |
|---|---|
| [COMPUTED] | `reachability.project_cone`, `combicone.certify_emergence` run on `combicone_substrate.npz` (Norman doubles; 105 single-gene ATEs, 131 doubles, 5045 genes). Reproduce: `scripts/causal_formalization_example.py`. |
| [DOSSIER] | `results/certificate_dossier.json` — `construct_validity`, `sensitivity` (`Γ*`), `nc_additive`. |
| [MANUSCRIPT] | `manuscript/sections/50_discussion.tex` (causal-inference agenda paragraph) and `manuscript/limitations_and_reinforcement_plan.tex` (L1, L2, L4). |
| [DESIGN] | Definitions and identifying assumptions (A1–A6); no measurement claimed. |

**External references (2025–2026 "causal turn," for the §3.1 contrast).** These are used only for
the qualitative contrast in §3.1, not for any computed result. PDGrapher — Gonzalez, Lin, Herath,
Veselkov, Bronstein & Zitnik, "Combinatorial prediction of therapeutic perturbations using causally
inspired neural networks," *Nature Biomedical Engineering* vol. 10 (2026; online 2025), DOI
10.1038/s41551-025-01481-x; the quoted phrases (PPI/GRN "approximations of the causal graph,"
"assumption of no unobserved confounders") are verified against the fetched full text. X-Cell —
"X-Cell: Scaling Causal Perturbation Prediction across Diverse Cellular Contexts via Diffusion
Language Models," bioRxiv 2026 (title verified via search; no internal metrics quoted here). scGPT
and the 2025–2026 single-cell foundation-model perturbation-benchmark wave — the additive-baseline-tie
finding as cited in `50_discussion.tex`.

*One flagged discrepancy for the maintainers:* `50_discussion.tex` states "150 additive
combinations" for the negative-control check, while the released `certificate_dossier.json` records
`nc_additive.n_pairs = 100` (false-positive rate `0.0` in both). Recommend reconciling the
manuscript text to the JSON artifact.
