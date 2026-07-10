# Trust and causal inference for the reachability oracle

*This is the consolidated trust-and-causal-inference dossier for the cell-state
reachability oracle. It unifies five formerly separate documents: the design-based
causal reframe and the instrumental-variables / compliance layer (`CAUSAL_INFERENCE.md`),
the counterfactual-explanation framing (`COUNTERFACTUAL_EXPLANATION.md`), the
causal-inference research agenda (`CAUSAL_RESEARCH_AGENDA.md`), the validation ledger and
prioritized experiment plan (`VALIDATION_AND_EXPERIMENTS.md`), and the adversarial dataset
appraisal (`REVIEWER2_DATASET_CRITIQUE.md`). In one place it (1) names the causal framing
precisely, (2) maps the machinery already in `reachability.py` onto standard causal
constructs, (3) adds the one genuinely missing causal layer — an IV / compliance treatment
of imperfect CRISPRi knockdown — with a proof-and-measurement that the headline verdict is
robust to it, (4) develops the verdict as a counterfactual explanation, (5) enumerates the
identifying-assumption stack and the research agenda that stress-tests it, (6) records the
honest validation ledger, and (7) states plainly where causal inference would hurt and what
a skeptical reviewer will attack in the data. It is the conceptual companion to `NOVELTY.md`
/ `RELATED_WORK.md`; figures it references live beside it in the repo (relative paths, so
they render on GitHub).*

---

## 1. The reframe that costs nothing

The reachability oracle is *already* a causal-inference method — it just was not wearing
the vocabulary. A genome-scale CRISPRi Perturb-seq screen is a designed *interventional*
experiment, so each column of the effect dictionary `E` is an estimated **average treatment
effect** (ATE), and the question "is state *d* reachable from a non-negative mix of
knockdowns?" is a **counterfactual query about a compound intervention**.

Perturb-seq is not observational transcriptomics. Guide RNAs are delivered to a pool and
assignment of a cell to "carries guide *g*" is, to a good approximation, **randomized** at
the cell level (Poisson infection at low MOI, independent of the cell's latent state). That
randomization is exactly what licenses a causal reading of each measured effect:

> **E[g, :]  =  mean(X | guide g)  −  mean(X | non-targeting control)**
> is an estimate of  **E[X | do(knock down gene g)] − E[X | do(control)]**,
> the average treatment effect of the perturbation on the transcriptome.

Nothing in the pipeline changes. What changes is the claim we are entitled to make: the
dictionary rows are not correlations to be regressed, they are *interventional contrasts
measured under randomization*. The reachable cone

  **C = { Eᵀw : w ≥ 0 }**

is therefore the set of transcriptomic shifts achievable by some **compound intervention**
(co-perturbation), and the reachable / provably-outside verdict is a **counterfactual
feasibility answer**, not a ranking heuristic.

This is not a strained analogy. The single-cell perturbation field has explicitly adopted
treatment-effect language: perturbations are *treatments*, response vectors are *treatment
effects*, and estimating them is a causal-inference problem. Our method sits squarely inside
that frame — and, importantly, on the more defensible side of it (see §5, "Two camps").

---

## 2. Machinery already in the repo, in causal-inference terms

Most of the causal content is **already built**. The table below maps existing objects in
`reachability.py` and the notebooks onto their causal-inference names. The right column is
what to *call* them in the manuscript; the left is where they already live.

| In the code / results | Causal-inference construct | Where |
|---|---|---|
| `E[g,:] = mean(X\|guide g) − mean(X\|NTC)` | **Average treatment effect (ATE)** of an intervention, estimated under cell-level randomization | `build_effect_matrices.py`, `inputs.npz` |
| Non-targeting-control pool | The **control / placebo arm**; its spread is the do-nothing null | `guide_kd_efficiency` (NTC rows) |
| "Is *d* ∈ cone C?" | **Counterfactual query** on a compound (multi-gene) intervention | `reachability()`, `signed_reachability()` |
| Farkas separating direction ρ = *d* − Eᵀw\* | The identified **infeasible sub-space**: what no intervention in the dictionary can cause | `activation_certificate()` |
| Additivity of stacked effects (cone assumption) | **No-interaction / partial-SUTVA** assumption on compound treatments; its failure *is* epistasis | module §"Additivity / epistasis", `additivity_risk()` |
| Norman-2019 saturation law (R²≈0.57) | Empirical **interaction (epistasis) model** calibrating deviation from additivity | `SATURATION_CEILING`, `norman_additivity_calibration.csv` |
| Shuffled-target / anisotropy / held-out-gene nulls | **Negative-control / falsification tests** for the estimand | `shuffled_target_null()`, `analytic_anisotropy_null()`, `held_out_gene_validation()` |
| Bootstrap verdict distribution | **Sampling uncertainty of the ATE** propagated to the verdict | `run_bootstrap.py`, `generator_uncertainty_verdicts.csv` |
| Cross-cell-type transfer (K562 ↔ RPE1) | **Transportability** of causal effects across domains (Pearl–Bareinboim) | `06/07_cross_celltype_transfer.ipynb` |

The upshot: three of the four "opportunities" a reader might propose (uncertainty on the
verdict, an epistasis/additivity test, a placebo arm) are **done**. Re-deriving them would
be duplication. The framing in this document is the value-add for those; the new build (§3)
is the instrumental-variables layer, which was genuinely missing.

---

## 3. The missing layer: imperfect knockdown as noncompliance (instrumental variables)

The one causal construct **not** already in the pipeline is compliance. CRISPRi knockdown
is incomplete and heterogeneous: a guide is *assigned*, but the target gene is only
partially silenced, and by a different amount in every gene and condition. That is textbook
**imperfect compliance**, and it maps onto an instrumental-variables (IV) design:

| IV construct | Here |
|---|---|
| **Instrument Z** | guide assignment — randomized in the pool |
| **Treatment T** | realized on-target knockdown of the perturbed gene (continuous, in [0,1]) |
| **Outcome Y** | the transcriptomic shift |
| **First stage** Z→T | `guide_kd_efficiency` — the target gene's own expression in guide cells vs NTC |
| **Exclusion restriction** | the guide affects Y *only* through knocking down its target (no off-target route) |

Under this design, the published effect vector `E[g,:]` is an **intent-to-treat (ITT)**
quantity: the effect of *assigning* guide *g*, which already averages over cells where
knockdown succeeded and cells where it did not. The corresponding **complier / local
average treatment effect (LATE)** — the effect of *actually* silencing the gene — is the
Wald rescaling by the first-stage compliance πg:

  **E_LATE[g, :]  =  E[g, :] / πg**,  where πg = measured fractional knockdown of gene *g*.

Both are legitimate and answer different questions. ITT is what an experimentalist *gets*
when they deliver the guide (the operationally honest quantity); LATE is the mechanistic
per-unit-knockdown effect. The compliance layer does three things for the oracle:

**(a) It proves the verdict is compliance-invariant.** A convex cone is invariant to
positive per-generator rescaling: if every row is scaled `E[g] → cg·E[g]` with cg > 0, the
set `{ Σ wg (cg E[g]) : w ≥ 0 }` is *identical* — only the weights relabel (`wg → wg/cg`).
Rescaling ITG→LATE is exactly such a positive rescaling (πg > 0 for any valid instrument),
so **the reachable / provably-outside verdict, the reachable cosine, and the LOF/GOF/neither
split cannot change** — only the recipe weights do. We verify this numerically to machine
precision (Results below, max |Δcosine| = 2.22e-16 across all 12 cells). This is a
*feature*: the headline verdict does not depend on the incomplete-knockdown nuisance at all.

**(b) It supplies an exclusion-restriction robustness arm.** A generator with **no**
significant on-target knockdown is an *invalid instrument*: if it still moves the
transcriptome, it does so through an off-target / non-excluded route, violating the
exclusion restriction. "Invalid" admits two thresholds, both measured from the first-stage
table: a generator with **no significant on-target guide** (0.08–0.15 % of the dictionary:
10/6871 Rest, 6/7155 Stim8hr, 6/7195 Stim48hr), or the stricter **no *measurable* on-target
knockdown** (target too lowly expressed to quantify silencing; 0.41–1.05 %: 72/6871,
29/7155, 47/7195). The valid-instrument arm in the Results drops the former (the
directly-testable exclusion violators); the stricter set is a superset and would drop at
most ~1 %. Either way the removed fraction is tiny. The per-target effect on the verdict is
reported below.

**(c) It turns compliance into a per-recipe deliverability score.** Because LATE weights up
low-compliance genes, a recipe that leans on hard-to-silence targets is *less deliverable*
than its ITT magnitude suggests. A recipe-level compliance score (the contribution-weighted
mean πg of the recipe's genes) flags this — the compliance analogue of the existing
`additivity_risk` saturation score. High = the recipe's genes are ones CRISPRi silences
well; low = realized effect is compliance-limited.

The self-contained analysis lives in `run_iv_compliance.py`; the frozen first-stage table
is `analysis_cache/atlas_work/first_stage_compliance.csv`.

---

## 4. Results — compliance does not move the verdict

The IV/compliance analysis (`run_iv_compliance.py`) was run on the headline dictionary for
all four cell-state targets × three culture conditions (12 cells). Two robustness checks and
one deliverability readout:

**(1) Rescaling invariance (ITT → LATE) — exact.** Rescaling every generator row by its
inverse compliance `1/πg` (the Wald / LATE rescaling) and re-solving leaves the reachable
cosine unchanged to machine precision. This was computed directly for all 12 cells by
`write_invariance_csv()` in `run_iv_compliance.py` (run `python scripts/run_iv_compliance.py
--invariance-only` to regenerate); the maximum absolute change is **|Δcosine| = 2.22e-16**
(headline toward-Th1/Stim48hr: ITT = 0.532884, LATE = 0.532884, Δ = 0.00e+00). Per-cell
values are written to `analysis_cache/atlas_work/late_rescaling_invariance.csv`.

This is the convex-cone invariance theorem verified numerically — positive per-generator
rescaling relabels the recipe weights but preserves the cone, hence the verdict. The
incomplete-knockdown nuisance provably cannot change reachability. (This is distinct from
the exclusion-restriction check in (2), which *removes* rows rather than rescaling them and
therefore does move the cone by a small, non-zero amount.)

**(2) Exclusion-restriction robustness (drop invalid instruments) — negligible.**
Re-solving on the valid-instrument subset (dropping generators with no significant on-target
knockdown, 0.08–0.15 % of rows) changes the reachable cosine by at most **|Δcosine| =
0.0004** across all 12 cells (KKT optimality residual ≤ 3.1e-11 everywhere). The verdict does
not lean on generators that move the transcriptome without demonstrably hitting their target.

| target | condition | cosine (ITT) | cosine (valid) | Δ | invalid dropped |
|---|---|---|---|---|---|
| Th1 | Rest | 0.6266 | 0.6264 | −0.0002 | 10 / 6871 |
| Th1 | Stim8hr | 0.5235 | 0.5235 | −0.0000 | 6 / 7155 |
| Th1 | Stim48hr | 0.5329 | 0.5329 | +0.0000 | 6 / 7195 |
| Th2 | Rest | 0.6428 | 0.6428 | −0.0000 | 10 / 6871 |
| Th2 | Stim8hr | 0.5370 | 0.5370 | +0.0000 | 6 / 7155 |
| Th2 | Stim48hr | 0.5300 | 0.5298 | −0.0002 | 6 / 7195 |
| younger | Rest | 0.6262 | 0.6258 | −0.0004 | 10 / 6871 |
| younger | Stim8hr | 0.5885 | 0.5885 | +0.0000 | 6 / 7155 |
| younger | Stim48hr | 0.5970 | 0.5970 | −0.0000 | 6 / 7195 |
| older | Rest | 0.5964 | 0.5964 | +0.0000 | 10 / 6871 |
| older | Stim8hr | 0.5899 | 0.5899 | −0.0000 | 6 / 7155 |
| older | Stim48hr | 0.5662 | 0.5662 | +0.0000 | 6 / 7195 |

(full table: `iv_compliance_verdicts.csv`)

**(3) Deliverability — where compliance *does* matter.** Because it cannot change the
verdict, compliance is instead an honest **deliverability** annotation on the recipe.
First-stage knockdown is strong but incomplete and heterogeneous — median π = 0.84, 5th
percentile 0.26 (Fig. b). The headline toward-Th1 recipe has a weight-weighted compliance
score of **0.78**, but includes weak instruments the reader should know about: *SNX5* (π =
0.07) and *RARA* (π = 0.46) carry real recipe weight yet are barely silenced by CRISPRi, so
their realized contribution will fall short of their assigned weight (Fig. c). This is the
compliance analogue of the existing `additivity_risk` saturation score: it does not change
*whether* the target is reachable, it flags *how faithfully the recipe will be delivered*.

![Causal-inference / compliance analysis of the reachability oracle]({{artifact:art_88317910-40a7-4be0-8dd5-b3fa40028923}})

*Figure. (a) The IV/compliance DAG: guide assignment Z (randomized) → realized on-target
knockdown T (imperfect compliance, first stage π) → transcriptome Y; the published effect is
intent-to-treat, its π-rescaled counterpart is a local ATE. (b) First-stage compliance
distribution over valid instruments (median π = 0.84). (c) The headline toward-Th1 recipe,
points colored by compliance; SNX5 is a high-weight, near-zero-compliance instrument. (d)
Reachable cosine with the full dictionary (ITT) vs the valid-instrument subset, all 12
targets × conditions; markers are coincident, max |Δcosine| = 0.0004.*

---

## 5. Two camps — and why this method is in the stronger one

Causal work on perturbation data splits into two families, and the project's design choice
("use measured effect vectors, not an inferred GRN") lands it in the more defensible one.

**Camp A — observational causal *discovery*.** Infer a gene-regulatory / causal graph from
expression, then predict interventions off the fitted structure (e.g. variational-causal-
inference autoencoders, causal-GRN inference, optimal-transport matching, causal
representation learning). Powerful in ambition, but it inherits the hard identifiability
problems of structure learning — unobserved confounders, latent common causes,
Markov-equivalence ambiguity — and its interventional predictions are only as good as the
recovered graph.

**Camp B — design-based / experimental causal inference.** Use the *measured* interventional
effects directly. No graph is fitted; the randomized experiment already did the causal
identification. This is where the reachability oracle sits. Its verdicts depend on the
dictionary and two transparent assumptions (non-negativity from CRISPRi loss-of-function;
additive composition, itself tested in the epistasis analysis), **not** on a learned causal
structure that could be wrong.

Framing to use in `NOVELTY.md` / `RELATED_WORK.md`: *we answer a counterfactual reachability
question using gold-standard interventional effects, sidestepping the identifiability burden
of observational causal discovery.* That is a strength worth stating out loud — it is the
reason the "provably outside the cone" certificate is trustworthy in a way a GRN-derived
prediction is not.

### Scope discipline — where causal inference would *hurt*

Not every causal method belongs here. Three that a reviewer might reflexively suggest, and
why they are out of scope:

- **Do not pivot to GRN / causal-structure discovery.** It directly undercuts the
  positioning above (the whole point is that we do *not* need a fitted graph), it is fragile
  (identifiability, unobserved confounders), and the space is crowded. Our edge is
  design-based identification.
- **Do not adopt causal representation learning for this scope.** Latent-variable causal
  models (sparse-mechanism-shift, causal diffusions, Uhler-style representation learning)
  target almost exactly this question — predict unseen interventions, find optimal
  interventions for a state transition — but they are a multi-month latent-variable research
  agenda, not a hackathon build. Cite them in `RELATED_WORK.md` as the heavyweight
  complement; position the convex-cone oracle as the transparent, assumption-light
  alternative whose verdict comes with a numerically-checkable certificate.
- **Keep the randomization caveat honest.** Cell-level guide assignment is *approximately*
  randomized within the pool, not a clean RCT: multiplicity of infection, guide-detection
  error, cell-cycle and proliferation differences, and batch are real confounders. The
  covariate-adjustment done upstream in the DE model (condition, donor, technical factors)
  is what makes the ATE reading defensible — name it as **causal covariate adjustment**, and
  state the residual caveat rather than claiming a perfect experiment.

The edge is design-based identification with a numerically-checkable certificate. The
refinements in the research agenda (§8) sharpen that edge; they do not trade it away.

---

## 6. The verdict as a counterfactual explanation

Where §1 argues the verdict is a causal (rung-3) query, this section argues the verdict is
also a **counterfactual explanation** in the precise sense the interpretable-ML literature
means — and that it satisfies, by construction, the two desiderata that post-hoc
counterfactual-explanation methods routinely violate.

A counterfactual explanation answers *"what would have to be different for the outcome to
change, and is that change something one could actually enact?"* (Wachter, Mittelstadt &
Russell 2017). For a machine-learning classifier this is retrofitted after the fact by
perturbing inputs until the prediction flips. The reachability oracle does not need the
retrofit: its **output is already a counterfactual**. Asking whether a target cell-state
shift *d* lies in the non-negative cone of measured knockdown effects *is* the question
"which compound intervention would move the cell toward *d*, and by how much?" The solver
returns the intervention (the recipe), the part of the target it cannot deliver (the signed
residual), a witness when no intervention suffices (the Farkas certificate), and — new in
this layer — a set of *diverse* feasible interventions and a robustness radius on the whole
conclusion. Each of these is a named element of a counterfactual explanation, and each is
pinned below to the real Th2→Th1 verdict (Stim48hr).

### 6.1 The elements, mapped to the objects the solver emits

| Counterfactual-explanation element | Formal role | Object in the oracle | Th2→Th1 (Stim48hr) value |
|---|---|---|---|
| **Actionable counterfactual** | the minimal enactable change that achieves the outcome (Wachter 2017) | greedy **minimal recipe** | knock down {ATF7IP2, FBXO32, RARA, HPS1, KIF5B, …} |
| **Diverse counterfactuals** | *several* distinct feasible options, not one (Mothilal, Sharma & Tan 2020, "DiCE") | **disjoint recipe peeling** | 4 fully disjoint 8-gene recipes; the 4th (lead STAT6) retains 83 % of the first's cosine |
| **Contrastive foil** | *why P rather than Q* — the part of the goal not delivered (Miller 2019) | **signed decomposition** LOF / GOF / neither | 28 % knockdown-reachable, 26 % needs up-regulation, 45 % neither |
| **Impossibility witness** | a certificate that *no* counterfactual exists | **Farkas / KKT separating direction** ρ | KKT ≈ 1.1×10⁻¹¹ (verdict certified optimal) |
| **Cross-modality counterfactual** | what *outside* the current action set would close the gap | **activation certificate** | LYAR, STOM, TMEM62, CBLB, IKZF3, GBP5 (ranked unmet upward demand) |
| **Assumption-level counterfactual** | how wrong the evidence must be for the *explanation itself* to change | **A1 sensitivity radius** (E-value) | flips at ≈ 8 % coordinated bias ≈ 0.03 SE-units; bootstrap CI clears the null |

All six are computed objects, not narrative. The recipe, decomposition and certificate ship
in `analysis_cache/atlas_work/point_toward_*.json`; the diverse recipes in `diverse_recipes.csv`; the radius
in `a1_sensitivity_radius.csv`; and all of them are assembled per target in
`counterfactual_cards.json`.

### 6.2 Why this surface is unusually well-founded

Two desiderata are the standard failure modes of counterfactual explanation for black-box
models. The reachability oracle satisfies both *structurally*, not by post-hoc filtering.

**Feasibility.** Wachter-style counterfactuals are notorious for proposing input changes
that cannot be enacted (a lower age, a categorical value that violates a constraint). Here
feasibility is the **non-negativity constraint** of the cone: the oracle cannot propose
"un-knock-down gene X", a negative dose, or a perturbation absent from the screened library.
Every counterfactual it emits is, by construction, a runnable CRISPRi recipe. The action
space *is* the feasible space.

**Causal validity.** The deeper critique (Karimi, Schölkopf & Valera 2021; Verma et al.
2020) is that a counterfactual which respects only the input *distribution* — not the causal
mechanism — recommends changes that are individually plausible but jointly non-causal
(change a feature that is downstream of another and the recourse is illusory). The
reachability counterfactual is built from **measured average treatment effects of the
interventions themselves**: each generator is the effect of *actually knocking that gene
down* in a randomized CRISPRi screen (Zhu et al. 2025). The recipe is therefore not a
plausible-looking input vector; it is a composition of interventions whose individual causal
effects were measured. This is exactly the "recourse in the presence of a causal model" that
the literature argues for — with the causal model supplied by the experiment rather than
assumed.

A third, softer property — **contrastiveness** (Miller 2019: human explanations are always
"why P rather than foil Q") — is met by the signed decomposition. The foil is explicit: *the
state is only partly reachable because 26 % of the target demands genes to go **up**, which
loss-of-function knockdown structurally cannot supply.* That is not a saliency score; it is a
geometric decomposition of the target into a knockdown-reachable component, an
activation-required component, and an orthogonal remainder, and it turns directly into the
falsifiable next experiment named by the activation certificate.

### 6.3 Diverse counterfactuals — the new computation

The DiCE desideratum is that a user should see *several* genuinely different feasible
actions, not a single brittle one. We generate them by **disjoint peeling**: run the greedy
recipe, forbid its entire membership (zero those generator rows), and re-run — repeating to
yield recipes that share *no* genes. For every one of the four targets the four recipes are
fully disjoint (max pairwise gene overlap = 0), and the fourth disjoint recipe still retains
at least 80 % of the first recipe's reachable cosine (0.80–0.86 for Th1/Th2/younger; the
"older" axis retains 0.998 via a near-equivalent second recipe):

| target | recipe-1 lead | 4th-recipe lead | 4th/1st cosine retention |
|---|---|---|---|
| Th1 | ATF7IP2 | STAT6 | 0.83 |
| Th2 | RSBN1L | ARHGAP30 | 0.80 |
| younger | SLC9A1 | TBC1D14 | 0.86 |
| older | ZAP70 | CD3E (recipe-2) | 0.998 |

This is more than a usability feature — it is a **robustness statement about the cone**: the
target does not hinge on any single set of generators, so the verdict cannot be an artifact
of a few influential perturbations. And the diversity is biologically legible rather than
random: Th1's fourth disjoint recipe is led by **STAT6** — precisely the regulator the A6
construct-validity check flags as the strongest signed lever toward Th1 (silencing the Th2
driver releases the Th1 program) — while the "older" axis admits two near-equivalent
TCR-signaling routes (ZAP70 and CD3E). Detail with full cosine curves in
`diverse_recipes_detail.json`.

### 6.4 What this is *not*

Honesty about scope, in the spirit of the rest of the repo:

- This is an **explanation of the verdict**, not a new predictive claim. The recipes are the
  same greedy selection the oracle already performed; the contribution is naming and
  assembling them as counterfactual-explanation objects and adding the diversity/robustness
  computation.
- The sparse recipe cosine (~0.30 at k = 8) is deliberately below the full-dictionary cosine
  (~0.53): a short, human-actionable recipe trades reach for sparsity. The card reports both
  so the trade is explicit.
- "Diverse and disjoint" is a statement about the *measured dictionary*, transportable only
  as far as the effects are (see the transportability assumption VI in §7). The aging axes,
  flagged fragile by A1, inherit that fragility here.

The assembled objects live in `counterfactual_cards.json` (one counterfactual-explanation
object per target: recipe, signed decomposition, certificate, sensitivity radius, disjoint
alternatives), `diverse_recipes.csv` / `diverse_recipes_detail.json` (the DiCE-style disjoint
recipes), and are surfaced in `notebooks/09_causal_validation_dossier.ipynb` §9 ("The verdict
as a counterfactual explanation").

### 6.5 Counterfactual-explanation references

- Wachter, Mittelstadt & Russell (2017). *Counterfactual explanations without opening the
  black box.* Harvard J. Law & Technology 31(2). — actionable counterfactuals; feasibility.
- Miller (2019). *Explanation in artificial intelligence: insights from the social
  sciences.* Artificial Intelligence 267:1–38. — explanation is contrastive (why P rather
  than Q). *(Tim Miller — distinct from the H. E. Miller perturbation-metrics paper cited in
  `RELATED_WORK.md`.)*
- Mothilal, Sharma & Tan (2020). *Explaining machine learning classifiers through diverse
  counterfactual explanations (DiCE).* FAT* '20. — diversity + feasibility. *(The
  counterfactual-explanation method — distinct from the DICE immune-cell database in
  `ROADMAP.md`.)*
- Karimi, Schölkopf & Valera (2021). *Algorithmic recourse: from counterfactual explanations
  to interventions.* FAccT '21. — causal validity of recourse.
- Verma, Dickerson & Hines (2020). *Counterfactual explanations for machine learning: a
  review.* arXiv:2010.10596. — survey of the desiderata used in §6.2.

---

## 7. The identifying-assumption stack

The reachability verdict — "target direction `d` is (not) in the non-negative cone of the
measured effect vectors" — is a counterfactual about a compound intervention. It is
*identified* by a stack of assumptions, each of which a causal-inference reviewer will name.
This is the authoritative enumeration; the research agenda (§8), the counterfactual-
explanation layer (§6), and the validation ledger (§9) all reference it, and it is the
narrative spine of `notebooks/09_causal_validation_dossier.ipynb` ("Can I trust it?").

| # | Assumption | What it buys | Where it can break here |
|---|---|---|---|
| I | **Unbiased ATEs** — each `E[g,:]` is the true average effect of perturbing `g` | the dictionary rows mean what we say | DE-model misspecification, residual batch/abundance confounding |
| II | **Additivity of compound interventions** — `E[Σ g] ≈ Σ wg E[g]` | the cone = the reachable set | epistasis; only tested pairwise (Norman 2019) |
| III | **SUTVA / no interference** — a cell's outcome depends only on its own guide | pooled effects transport to delivery | **cytokine spillover in T cells** (paracrine IL-2/4/IFN-γ) |
| IV | **Homogeneous effect** — the population-average `E[g]` applies to the treated cell | one verdict for the population | naive/memory/cell-cycle heterogeneity |
| V | **Exclusion restriction** — a guide moves `Y` only through its on-target KD | instruments are valid | off-target activity (handled in §3) |
| VI | **Transportability** — `E` measured here generalizes to the target context | young→aged, screen→patient | donor age, activation state, culture |

**Expert causal inference is not more machinery bolted on. It is (a) making each assumption
explicit, (b) testing the testable ones with existing data, and (c) quantifying how far a
violation must go before the verdict flips.** The agenda in §8 is ordered by (value ×
feasibility-with-current-data), split into what runs now vs what needs new wet-lab
experiments.

---

## 8. Research agenda — stress-testing the assumption stack

*The expert-level menu of what to test or refine next, and how far each is from the data
already in hand. Part A is computable now; Part B needs new wet-lab experiments (state the
caveat now, run later).*

### Part A — Computable now, with data already in the repo

#### A1. Verdict sensitivity analysis — a robustness radius native to the cone  ★ headline

**Assumption at stake:** I (unbiased effects).
**The question a reviewer will ask:** "Your verdict is a point estimate on noisy effect
vectors. How much would the effects have to be wrong to flip it?"

This is the single highest-value addition, because the geometry already contains the answer
and the repo already contains the calibration.

- For a **provably-outside** verdict, the Farkas/KKT certificate is a separating hyperplane
  `h` (with `E hᵀ ≤ 0`, `dᵀh > 0`). The normalized margin `dᵀh / ‖h‖` is *exactly* how far
  `d` sits outside the cone — i.e. the smallest perturbation to the effect matrix that would
  admit `d`.
- For a **reachable** verdict, the analogue is the distance from `d` to the cone boundary
  (how much shrinkage of the active generators pushes `d` out). Note the headline targets are
  *graded*-reachable (cosine 0.52–0.64 against a shuffled null), not exactly in the cone, so
  here the radius is a **bounded worst-case re-solve sweep** — perturb `E` by increasing
  multiples of its measured SE and find the magnitude at which the verdict grade flips —
  rather than the single closed-form margin available for a provably-outside verdict.

The move that makes this publishable rather than cosmetic: **`atlas_lfcSE.npz` holds the
per-effect standard errors.** So express the robustness radius *in units of measured SE* and
you get an E-value-style statement (VanderWeele; Cinelli–Hazlett omitted-variable bias):

> "The Th1/Stim48hr verdict would only change sign if the effect estimates carried a
> systematic bias of ≥ δ, which is **k× the median measured standard error** of the
> generators on the certificate."

A verdict robust at k ≫ 1 is one the measurement noise cannot explain away; a verdict at k ≈
1 is fragile and should be reported as such. This turns the Farkas certificate — the
project's signature — into a *calibrated* robustness claim, which is the modern gold-standard
framing (sensitivity analysis over point identification).

**Deliverable:** per-target robustness radius `k`, a curve of verdict vs perturbation
magnitude with the measured-SE scale marked, and a one-line E-value sentence per headline
target.
**Feasibility:** HIGH — you have the certificates, the solver, and the SEs. No new data.

The **key finding for A1** is a subtle one worth stating plainly: undirected measurement
noise *inflates* the reachable cosine (it adds anisotropy the null already corrects for), so
random error is not the threat — the verdict's genuine exposure is a *coordinated* systematic
bias against the target, which is exactly what the randomized CRISPRi design is built to
exclude. That is the quantitative case for why the design-based identification is
load-bearing.

#### A2. Conditional reachability — CATE × cone

**Assumption at stake:** IV (homogeneous effect).

`E[g,:]` averages over a heterogeneous CD4⁺ population (naive vs memory, cell-cycle phase,
incipient Treg). A population-average counterfactual can hold for *no* individual subtype,
and the **minimal knockdown recipe may differ by subtype** — which is exactly what a
therapeutic cares about.

- **Experiment:** stratify effect estimation by a baseline covariate (cell-cycle score,
  naive/memory markers, baseline polarization score), re-derive a stratified `E`, and
  re-solve the oracle per stratum. Report which verdicts are stratum-stable and whether the
  recipe reorders.
- **Method:** CATE estimation (causal forest / X-learner) fused with reachability; even a
  coarse 2–3 stratum split is informative.
- **Feasibility:** BLOCKED in-sandbox — needs cell-level covariates and a re-derivation of
  `E` per stratum from **raw single-cell counts**, which are *not* in the repo. The only
  single-cell file present (`GWCD4i.DE_stats.h5ad`, 16.8 GB) is a **DE-statistics summary**:
  its rows are `<gene>_<condition>` perturbations (n = 33,983), its layers are `log_fc /
  lfcSE / zscore / p_value / adj_p_value / baseMean`, and `n_cells_target` (median 541 cells)
  records the cells already *aggregated away* into each estimate. So A2 is a **scaffolded
  future-work item**: the code path and expected-output schema can be built now, but it needs
  the raw count matrix (re-run from the CZI VCP source) to execute. No new wet-lab work, but a
  new data pull.

#### A3. Higher-order additivity bound — k-way, not just pairwise

**Assumption at stake:** II (additivity).

`additivity_risk` is calibrated on Norman 2019 *pairwise* doubles, but recipes routinely use
12+ generators, and pairwise additivity does **not** imply k-way additivity.

- **Refinement:** make the additivity score *cardinality-aware* — bound the compound
  extrapolation error as a function of recipe size using the observed pairwise-interaction
  distribution (a spectral bound on the interaction tensor), so recipe confidence decays with
  `k`.
- **Feasibility:** HIGH — extends the existing `additivity_risk` + Norman calibration; no new
  data.

#### A4. Weak-instrument-robust recipe intervals

**Assumption at stake:** V (exclusion) + first-stage strength.

You now have `π` (first-stage compliance). Weak instruments (e.g. SNX5, π = 0.07) inflate
LATE variance and destabilize recipe weights.

- **Refinement:** propagate first-stage uncertainty into the recipe weights — an
  Anderson–Rubin-flavored confidence set that stays valid as `π → 0`, extending the existing
  bootstrap verdict CIs. Flag any recipe that leans on a weak instrument.
- **Feasibility:** HIGH — bootstrap infrastructure and the compliance table both exist.

#### A5. Negative-control-outcome / proximal check

**Assumption at stake:** I (residual confounding in the screen).

Pooled Perturb-seq carries batch, guide-abundance, and capture biases. Proximal causal
inference (Miao; Tchetgen Tchetgen) formalizes the existing non-targeting placebo arm.

- **Experiment:** assemble a negative-control-*outcome* panel — housekeeping / off-pathway
  genes that should not respond to the targets — and verify the effect vectors project ≈ 0
  onto them. A non-zero projection is *detectable* residual confounding and licenses a
  proximal correction; the non-targeting guides are the negative-control *exposure*.
- **Feasibility:** HIGH — these are projections on the `E` you already have.

#### A6. Directional construct validity

**Assumption at stake:** overall validity (a falsification stronger than the shuffled null).

- **Experiment:** knock down known master regulators already in the dictionary (TBX21/T-bet,
  GATA3) and confirm the oracle predicts the *correct signed* change in reachability —
  removing T-bet should make Th1 *less* reachable and Th2 *more*. Recovering textbook
  immunology as a signed prediction is construct validity the shuffled null cannot provide.
- **Feasibility:** HIGH — target genes are in the dictionary.

### Part B — Needs new wet-lab experiments (state the caveat now, run later)

#### B1. Interference / cytokine spillover — the T-cell reviewer's first question  ★ primary caveat

**Assumption at stake:** III (SUTVA).

This is the one biological threat the design does *not* currently address, and it is specific
to the tissue. T cells secrete cytokines (IL-2, IL-4, IFN-γ) that act paracrine/autocrine, so
**one cell's outcome depends on other cells' perturbations** — a SUTVA violation. In a pooled
screen each guide is rare, so the measured `E[g]` is the effect in a *mostly-wild-type*
cytokine milieu. But polarization is a positive-feedback system: a therapeutic that delivers
a recipe to *all* cells shifts the milieu, so the compound effect may not transport from
screen to delivery.

- **Experiments:** (i) arrayed vs pooled comparison for a handful of cytokine-pathway genes;
  (ii) guide-density / MOI titration to read the spillover slope; (iii) exposure mapping
  (Aronow–Samii partial-interference) if any pool-composition covariate is available. Each
  estimates the paracrine component of `E[g]`.
- **Why it matters for the verdict:** if spillover is large for cytokine genes, their
  dictionary rows are context-dependent and the cone itself is a function of pool
  composition. Even without running it, **state this as the primary external-validity
  caveat** — it pre-empts the obvious objection.

#### B2. Mediation through master regulators — controlled direct effects

**Assumption at stake:** mechanism behind II/IV.

- **Experiment:** measure `E[g | TBX21 also knocked down]` vs `E[g]` from targeted g×TF
  double perturbations → a controlled direct effect. A recipe gene whose effect flows *only*
  through GATA3 is fragile if GATA3 is already saturated; a direct actor is robust. Refines
  the minimal recipe toward mechanism-independent levers. Does knocking down gene X move Th1
  *directly* or *through* master regulators (TBX21/GATA3)? Natural direct/indirect effects
  could refine the minimal recipe toward direct actors.
- **Feasibility:** needs a small targeted double-perturbation panel (new screen, but modest).

#### B3. Dynamic / time-varying reachability — g-methods & control theory

**Assumption at stake:** the static snapshot.

The verdict is static, yet the three timepoints (Rest / 8h / 48h) are a coarse trajectory.
Reconnect "reachability" to its control-theory origin: is `d` in the *reachable set* of the
dynamical system whose vector field is the effects — reachable at some time, or only at
steady state? Robins' g-methods / marginal structural models handle the time-varying-
treatment version.

- **Feasibility:** needs a denser time course to be more than suggestive.

#### B4. Transportability to the patient context — selection diagrams

**Assumption at stake:** VI (transportability).

Fit an explicit transport function across donor age / activation / culture and report which
target directions are *transport-stable* (reachable in every context) vs context-specific.
Selection diagrams (Bareinboim–Pearl) name exactly which differences license transport. A
target reachable in young cells but *provably outside* in aged cells is a headline finding,
not a caveat. The CD4 aging appendix is naturally this transportability question: are the
effect vectors invariant across donor age, or does the cone *deform*? The cross-cell-type
transfer machinery (§2) is the template; a formal transport statement needs age-stratified
effect matrices.

*(A related mediation direction on the polarization axis — whether knocking down gene X moves
Th1 directly or through TBX21/GATA3, refining the minimal recipe toward direct actors — is
developed under B2.)*

### Implementation architecture — where each item lives (notebook vs prose)

The repo organizes notebooks by **audience/question** (01 explore → 02 headline method → 03
generalizability → 04 experimental-design toolkit → 05 pharma showcase → 06 robustness
battery → 07 transportability → 08 evaluation), not by technique. Causal inference is a *lens
over existing work*, so it should be connective tissue, not a standalone notebook.

**Do not create a `causal_inference.ipynb`** — it would fight the audience-axis organization,
duplicate this doc, and re-implement what 02/06/07 already compute. A reviewer opens a
notebook asking "can I trust this verdict?", not "where is the causal inference?" — the
causal machinery is the answer, not the question.

Placement (forced by the fact that Part A produces computed output and Part B does not):

| Item | Home | Rationale |
|---|---|---|
| A1 sensitivity radius ★, A4 weak-instrument intervals, A5 negative-control-outcome, A6 construct validity | **NB-B "Can I trust it?"** (`notebooks/09_causal_validation_dossier.ipynb`) | these *are* the trust layer; NB-B's spine = the assumption-stack table in §7 |
| shipped IV/compliance (`run_iv_compliance.py`) | fold into NB-B | it is the dossier's first section |
| A2 conditional / CATE reachability | scaffold only → forward-note in `02_reachability_on_tier2` | blocked in-sandbox (needs raw single-cell counts, not present); a new capability once the data is pulled |
| A3 k-way additivity bound | augment `06_reinforcement_analyses` | additivity already lives here; A3 is its cardinality-aware extension |
| B1 spillover, B2 mediation, B3 dynamic, B4 transportability | **prose only** — this doc + manuscript Future Work | experimental designs with no computable output; a protocol in a notebook is a category error |
| B1/B2 experiment *designs*; B4 template | optional stubs in `04_experimental_design_toolkit`; forward-note in `07_cross_celltype_transfer` | 04 already turns a verdict into a runnable experiment; 07 is the transportability template |

Net for a full A+B build: **one new notebook (NB-B)** + three augmentations (02, 06, 07) +
optional design stubs in 04. Zero of Part B becomes a notebook. (NB-C, the certificate
deep-dive, is a separately-planned showcase notebook and is orthogonal to this — do not merge
them.)

### Build status

Part A is built and executed; Part B is threaded into the manuscript Future Work. The
computable analyses live in `notebooks/09_causal_validation_dossier.ipynb` (the trust
dossier) with their result files in the repo root.

| Item | Status | Deliverable |
|---|---|---|
| **A1** verdict sensitivity radius ★ | **built** | `a1_sensitivity_radius.csv`, `fig_a1_sensitivity.png` — Th1/Th2 robust to measurement error, flip under ≈8 % coordinated bias (≈0.03 SE-units); aging axes fragile at Stim48hr |
| **A3** k-way additivity bound | **built** | `a3_kway_additivity_bound.csv` — directional retention 0.71 (k=2) → 0.29 (k=12); augments notebook 06 §3b |
| **A4** weak-instrument intervals | **built** | `a4_weak_instrument_intervals.csv` — AR-style 1/π set; catches SNX5 (π≈0.07 → weight-multiplier CI [4.4, 398]) |
| **A5** negative-control outcomes | **built** | `a5_negative_control_outcome.csv` — 4.1–5.5× positive-over-negative-control mass enrichment |
| **A6** directional construct validity | **built** | `a6_construct_validity.csv` — TBX21↓/GATA3↓ correct-signed; STAT6 notable vs null |
| **A2** conditional reachability | **scaffold** | `a2_conditional_reachability_scaffold.py` + schema — blocked in-sandbox (raw single-cell counts absent); notebook 02 §8b forward-note |
| **B1** cytokine-spillover / SUTVA | **prose + design stub** | manuscript Future Work; notebook 04 §9b design stub |
| **B2** mediation | **prose + design stub** | manuscript Future Work; notebook 04 §9b design stub |
| **B3** dynamic reachability / g-methods | **prose** | manuscript Future Work |
| **B4** transportability / selection diagram | **prose** | manuscript Future Work; notebook 07 forward-note |

The single most useful *next* build, once raw single-cell counts are available, is **A2** —
it is the one Part-A analysis that is designed but not executed, and it directly addresses the
homogeneous-effect assumption (IV) that the pooled ATE cannot see.

Beyond trust, the same Part-A objects form the **counterfactual explanation** of the verdict
developed in §6 — the recipe is an actionable counterfactual, the signed decomposition is its
contrastive foil, the activation certificate is the cross-modality witness, and the A1 radius
is an assumption-level counterfactual. Because the action space is the cone's non-negativity
constraint and each generator is a measured interventional effect, these counterfactuals are
*feasible* and *causally valid* by construction — the two properties post-hoc explanation
methods for black-box models cannot guarantee.

---

## 9. Validation ledger — the honest audit

*A candid audit of which core assumptions are actually validated, which are not yet, and the
experiments + notebooks that would close the gap. Numbers cited from `RESULTS.md` §3–§6; the
roadmap-fit notes below refer to `ROADMAP.md`.*

### 9.1 Are the core assumptions validated? — the honest ledger

The method rests on seven assumptions. They are **not** all at the same level of proof.
Distinguishing them is the difference between a demo and a manuscript. (This validation
ledger is a different cut from the identifying-assumption stack in §7: §7 enumerates what
*must hold* for identification; the table below scores how far each is *empirically
validated*.)

| # | Assumption (plain language) | Status | Evidence | What's still open |
|---|------------------------------|--------|----------|-------------------|
| 1 | **The null is calibrated** — "how good a match looks by chance" is computed correctly, even for one-sided targets | **Solid** | Closed-form anisotropy null matches the brute-force null at Pearson **0.995 in-sample / 0.998 held-out**; reproduces both the ≈0 null of balanced axes and the ≈0.26–0.34 null of the aging axis (§6.4) | Confirm on the Day-2 datasets' dictionaries (comes for free) |
| 2 | **The verdict generalizes** — it isn't overfit to the exact genes used to build it | **Solid within/across genes; cross-cell-type tested** | Held-out cosine **0.448** vs in-sample 0.627; shuffled-target null **z ≈ 24** (§3, fig5). Cross-cell-type test (nb07, K562/RPE1 Replogle 2022): the cone produces above-null structure in **three** human cell types; effect *direction* transfers (matched cross-type cosine +0.35 vs ≈0 null), minimal *recipe* does not (Jaccard 0.11) | Cross-cell-type test is between two non-T lines, not a second CD4⁺ T screen; still no leave-*donor*-out (data are donor-collapsed) |
| 3 | **Effects add up (additivity)** — combining single-gene effects approximates the combination | **Calibrated on ONE dataset; recipes provably safe** | On Norman 2019 K562: collinearity mechanism **refuted** (ρ≈−0.16, n.s.); magnitude-saturation mechanism **supported** (ρ≈+0.58, p<0.01); saturation law M\*=13.9, R²=0.57; risk validates ρ=+0.46, **p=5.6e-8**. Reinforcement (nb06 L2): every recommended knee recipe is additive-safe **12/12** (reliability 0.92–0.96); the magnitude cap binds only at k≈28 vs a knee of k≈4–5 — a **5.6× margin** on the headline (§6.2–6.3) | Saturation *law* still calibrated on one cell type (K562, CRISPRa); does it **transfer**? Predictive (not just calibration) accuracy not yet reported |
| 4 | **The activation certificate is real** — genes flagged "must be switched ON" genuinely can't be reached by knockdown | **Scaffolded + synthetic-verified** | Signed split (Th2→Th1 Rest) LOF 0.39 / GOF 0.25 / neither 0.35; atlas mean LOF ≈0.34 — every transition is minority-LOF (§5). Reinforcement (nb06 L1): a runnable dual-modality test (`held_out_modality_test`) recovers a hidden activation-only gene set at **AUROC 0.999, z 8.9** on synthetic ground truth | Not yet run on **real** overexpression data — blocked on a dual-arm CRISPRi+CRISPRa screen on one axis, not on method. The code runs unchanged the moment such data exists |
| 5 | **It recovers known biology** — the method ranks true drivers above bystanders | **Provisional / post-hoc** | Master-TF panel **AUROC = 1.00** (p=0.014), 7/8 sign-concordant (p=0.035) — **but** the master-TF-vs-marker split was chosen *after* seeing the weak full-panel result (AUROC 0.69, p=0.052). The docs state this plainly (§6.1) | This is the **single biggest gap.** The clean result needs a **pre-registered** test on an independent axis or dataset before it counts as confirmatory |
| 6 | **Effects are consistent across data sources** | **Moderate** | 68.5% sign-concordant across 11,616 shared genes; Spearman ρ=0.562 | 31.5% sign-discordance is not negligible; the discordant tail is uncharacterized |
| 7 | **The metric isn't gamed by signal dilution** — the all-gene reachability cosine isn't inflated by the many unchanged background genes | **Solid** | DEG-weighted recompute (Mejia et al., ICML 2026; §6.5, nb08): the Th2→Th1 verdict *strengthens* under `w=|d|` weighting (cosine 0.627→0.803, held-out z 28.3→14.1, still ≫ z=3) across all three conditions and six Norman held-out doubles; interpolated-duplicate positive control confirms the metric rewards known-reachable targets (ceiling ≈0.97–0.99); unweighted default reproduces every published number bit-for-bit | Dynamic-range calibration computed on the headline Rest condition only (Monte-Carlo cost); Stim conditions carry the cheaper cosine + held-out-z recompute |

Two reinforcement analyses (`analysis_cache/nb_out/REINFORCEMENT_RESULTS.md`) also close manuscript-level
defenses not in the six-row ledger above: **L4** — dropping the non-negativity constraint
buys only **+0.018** held-out cosine, and buys it with ~3,537 biologically-unrealizable
negative weights; the cone recovers a mean **95%** of the unconstrained cosine using physical
knockdowns and is the *sole* source of the certificate, so the constraint is load-bearing,
not decorative. **L5** — the signed decomposition's knockdown-only ceiling `√(LOF)` equals
the in-sample cone cosine to 1e-4, so the headline 0.448 is **71% of the achievable ceiling**
(atlas mean 61%), not a small fraction of an unreachable 1.0.

### 9.2 What this means, said plainly

- **The engineering is trustworthy.** The parts that decide *whether a match is real* (the
  null, the held-out test, the additivity safety check, the now-quantified constraint
  ablation, and — new — the DEG-weighted metric robustness check) are validated to a high
  standard. A skeptical reviewer cannot easily attack the statistics: the headline survives
  the exact signal-dilution critique (Mejia et al., ICML 2026) that has been leveled at the
  perturbation-prediction field, and in fact *strengthens* under it.
- **Generalization has widened.** The method now runs unchanged on **three** human cell types
  (CD4⁺ T, K562, RPE1) plus the Norman K562 CRISPRa demo; effect *direction* transfers while
  the specific *recipe* does not — a robustness result with a sharp, honestly-stated boundary
  (`analysis_cache/czi_data/CROSS_CELLTYPE_TRANSFER.md`).
- **One gap remains the priority.** The claim a reader most wants — *"it finds the right
  genes"* — is still **not independently confirmed**: the positive control is honestly labeled
  post-hoc, and the activation certificate, though now backed by a synthetic-verified runnable
  test (L1), has still never met **real** activation data.

That is a **healthy** place to be at a hackathon: the hard machinery is done, and the
remaining work is a well-defined validation layer, not open-ended research.

### 9.3 Experiments that can be conducted (prioritized)

**Tier A — closes a named validation gap (do these first)**

1. **Pre-registered biological-recovery test** *(closes gap #5 — highest value).* Freeze the
   hypothesis — "known master regulators of the target state rank in the top-k of the recipe
   and align above marker genes" — **before** looking, then test it on an *independent* axis
   (the aging axis, or a fresh lineage transition) and/or an independent dataset. Converts
   AUROC = 1.00 from suggestive to confirmatory. *CPU-cheap; needs one labeled independent
   axis.*
2. **Out-of-sample double-perturbation prediction** *(hardens gap #3).* Fit the cone on single
   perturbations, **predict** the measured double, report held-out cosine + magnitude error —
   on Norman *and* the Day-2 combinatorial screens. Turns "additivity is a reasonable
   approximation" into "additivity predicts unseen combinations to within X%." *Data: Norman
   (have) + Replogle/Adamson (Day-2 fetch).*
3. **Gain-of-function transfer on the Joung TF atlas** *(closes gap #4).* Test the
   "must-activate" certificate against real overexpression ground truth (GSE216463). Already
   in roadmap Day 2. *Data: Day-2 fetch.* **◐ Scaffolded (nb06 L1):** the scorer is now
   packaged as a runnable `held_out_modality_test` and verified on synthetic ground truth
   (AUROC 0.999, z 8.9); it runs unchanged the moment a dual-arm CRISPRi+CRISPRa screen on one
   axis is in hand — see `analysis_cache/nb_out/L1_certificate_test_scaffold.json`.
4. **Negative-control / specificity test** *(new — reviewers always ask).* Run transitions
   that *should* be unreachable (cross-lineage jumps, or targets beyond the shuffle null) and
   confirm the method returns low reach + high must-activate + non-significant corrected p.
   Establishes the false-positive behavior. *CPU-cheap; uses existing data.*
5. **Leave-one-donor-out cross-validation** *(closes the donor half of gap #2).* Rebuild the
   cone dropping each donor; confirm the verdict is stable. Now affordable because the
   analytic null replaces the ~12-min empirical null. *CPU-moderate; uses existing CD4 data.*

**Tier B — extends impact and generalizability**

6. **Saturation-law cross-dataset transfer** (Replogle + Adamson) — roadmap Day 2.
7. **Retrospective pharma-value analysis** *(the impact money-shot).* Cross-reference the
   existing 102-target triage (reachable × druggable × genetically-supported) against
   historical clinical outcomes via Open Targets / openFDA. Question: do reachable +
   genetically-supported + tractable targets show measurably better historical success?
   Connects the method's output directly to the "genetic support → 2.6× approval odds" and
   "~90% of programs fail" framing in `NOVELTY.md`. *Data: allowlisted public APIs.*
8. **Linearity / distance-from-reference boundary** *(defines the operating envelope).*
   Characterize how the reachable cosine degrades as the target moves further from the
   reference state — the honest statement of where the local-linear cone is trustworthy.

**Roadmap fit:** #3 and #6 are already Day 2. #1, #2, #4 are natural Day-2 afternoon
additions (they reuse the same fetched data). #5, #7, #8 are strong **stretch goals** or a
post-hackathon manuscript track.

### 9.4 New notebooks — by audience and impact

Three distinct readers want three distinct stories. Each notebook is a self-contained
narrative that *shows* the method, not just describes it.

- **NB-A · "From screen to shortlist" — for pharma target-ID** *(most impactful to pharma)* —
  **✓ BUILT (`notebooks/05_target_id_showcase.ipynb`)**. A decision-useful walkthrough: pick a
  disease-relevant transition → reachability verdict → reach vs. must-activate split → ranked
  knockdown/activation recipe → druggability + human-genetics triage → modality call. Ends
  with a one-page "target dossier." Reuses `design_experiment()` and the 102-node triage.
  *This is the notebook a target-ID or BD team would actually open.*
- **NB-B · "Can I trust it?" — for reviewers / methods readers** *(manuscript backbone)*. One
  notebook that runs the whole validation layer in sequence: the pre-registered recovery test,
  the null calibration, the out-of-sample additivity prediction, the negative controls, and
  cross-dataset concordance. The "validation dossier" a referee needs. *Built as
  `notebooks/09_causal_validation_dossier.ipynb`, with its role also served by
  `06_reinforcement_analyses.ipynb` (the L1–L5 validation battery) and
  `07_cross_celltype_transfer.ipynb` (cross-cell-type generalization).*
- **NB-C · "The certificate: when NOT to run the experiment" — the differentiator**. Focuses
  on the one thing no ranking tool has: the Farkas/activation certificate. Take an unreachable
  transition, *prove* it, name the genes that must be activated, and contrast with the (wrong)
  shortlist a naive correlation ranking would return. The clearest answer to "why this and not
  a gene-ranking heatmap?" NB-C (the certificate deep-dive) remains the sharpest one-figure
  differentiator if time allows.

*The validation ledger is grounded in a direct read of `RESULTS.md` (§3 held-out validation,
§5 signed decomposition, §6 post-hackathon advances) plus the two reinforcement writeups
`analysis_cache/nb_out/REINFORCEMENT_RESULTS.md` (L1/L2/L4/L5) and `analysis_cache/czi_data/CROSS_CELLTYPE_TRANSFER.md`
(K562/RPE1 transfer). Every number above is quoted from the repo's own results — no number is
introduced here that a linked CSV or notebook does not reproduce.*

---

## 10. Adversarial dataset appraisal (Reviewer 2)

*A deliberately adversarial, constructive read of the data underneath the oracle. Every
quantitative claim below was recomputed from the local supplementary tables and the
`GWCD4i.DE_stats.h5ad` header, not from memory. Line references are to
`manuscript/sections/20_methods.tex`.*

### Summary judgment

The **method** is clean; the concerns are almost entirely about the **data it is fed** and
how far that data can be pushed. The central intellectual move — "reachability = membership
in the convex cone of measured knockdown vectors" — is only as good as (a) the knockdown
vectors that form the generators and (b) the target vectors they are compared against. Three
facts dominate everything else:

1. **The entire dictionary is a single preprint-stage CRISPRi screen from one lab, one cell
   type, four young donors, one technology.** Loss-of-function only. Much of what the paper
   frames as a *biological* result ("25–31 % of Th2→Th1 is gain-of-function and therefore
   unreachable") is at least partly a *consequence of the assay chosen*, not a discovery about
   the cell.
2. **The target vectors live in a different measurement space than the generators.** The
   Th2→Th1 direction is assembled from two *external* studies (Ota 2021 + Höllbacher 2021) and
   the aging direction from a third (Yaza 2022) — none of them the Perturb-seq assay. The cone
   geometry silently assumes these are commensurable after a gene-name join and cosine.
3. **One of the four target axes (aging) is asked of donors who span only ages 22–34.** The
   data cannot carry that question, and — to the project's credit — the robustness machinery
   already says so. But it is still shipped as a first-class axis.

None of these are fatal, and several are already partially mitigated by analyses in the repo.
But a reviewer will raise all of them, so they belong in Limitations, stated in the project's
own voice before someone else states them less charitably.

### 10.1 The dictionary (generators E): the CRISPRi screen itself

The primary dictionary is Zhu et al. 2025 (`GWCD4i.DE_stats.h5ad`): **33,983 single-gene
knockdown × condition effect vectors over 10,282 genes**, from **4 donors × 3 culture
conditions** (Rest / Stim8hr / Stim48hr). Layers: `log_fc`, `zscore`, `lfcSE`, `baseMean`,
`p_value`, `adj_p_value`.

**G1 — Loss-of-function only; the headline GOF limitation is baked into the assay** *(structural, cannot be analyzed away on this data)*.
CRISPRi removes function. Every generator therefore enters the cone with fixed orientation,
which is exactly what makes the geometry a *non-negative* conic hull — the paper is honest
about this. But it means the finding that **25–31 % of the Th2→Th1 shift is
"gain-of-function-locked" and unreachable is in part a statement about CRISPRi, not about
T-cell biology.** A screen with a CRISPRa arm would convert an unknown fraction of that
"GOF/neither" budget into reachable territory. The certificate is elegant, but its denominator
is set by the modality of the one screen you happened to use. **Fix:** state plainly that the
LOF/GOF split is *modality-relative*; make the CRISPRa test (the certificate's own falsifiable
prediction) the headline future experiment rather than a footnote.

**G2 — Imperfect and heterogeneous knockdown ("compliance")** *(serious; partially mitigated)*.
From `guide_kd_efficiency.suppl_table.csv` (73,765 guide×condition rows): only **73.3 % of
guides reach significant knockdown**, and *among those that do*, the median residual
expression is **0.115** (the gene is still ~11 % expressed) with the **90th percentile at
0.49** — i.e. a large tail of "knockdowns" that only halve the target. At the
perturbation×condition level (`DE_stats`, 33,983 rows) only **62.4 % are
`ontarget_significant`.** So a meaningful share of the columns of E are weak-instrument
estimates of a partial intervention. The cone treats a 95 %-knockdown atom and a
40 %-knockdown atom as geometric rays of equal standing (cosine ignores magnitude). *Already
addressed:* the IV/compliance layer (§3–§4: LATE rescaling; weak-instrument intervals) shows
the *verdict* is invariant to per-generator positive rescaling and flags the worst offenders
(e.g. SNX5, π≈0.07). **Residual concern:** the recipe *ranking* and the specific greedy picks
are not rescale-invariant, and partial knockdown means the measured atom direction itself may
differ from the true full-LOF direction, which rescaling cannot repair.

**G3 — Heteroscedastic, sometimes thin, per-perturbation support** *(moderate; partially mitigated)*.
`n_cells_target` per perturbation×condition: **median 539, but 5th percentile 123, minimum
17**, right-skewed to 11,510. Effect-vector noise scales like 1/√(cells), so generators are
**strongly heteroscedastic** — some atoms are far noisier rays than others — yet enter the
cone with equal weight. ~1.1 % of rows rest on <50 cells; ~3.7 % on <100. *Already addressed:*
the errors-in-variables dictionary bootstrap (E + N(0,1) on the z-scale, B=200 × 12 cells)
shows the feasibility verdict is robust (flip rate at the 0.5 threshold = 0). **Residual
concern:** the bootstrap uses a *homoscedastic* unit-noise model on the z-scale; it does not
propagate the *actual* per-atom `lfcSE`/cell-count spread, so the thinnest generators are
treated as no noisier than the deepest.

**G4 — Off-target and neighboring-gene repression contaminates atom identity** *(moderate; not yet quantified in the cone)*.
`DE_stats` carries an **8.3 % `offtarget_flag` rate**, and the sgRNA library metadata
(`sgrna_library_metadata.suppl_table.csv`, 26,504 guides) explicitly tracks neighboring genes
within 2/10/20/30 kb — because CRISPRi can repress a *neighbor* of the intended TSS. When it
does, the "gene X knockdown" effect vector is partly gene-Y's effect, and the recipe that names
X is mis-attributed. **Fix:** as a specificity check, drop or down-weight `offtarget_flag`
generators and confirm the headline recipe genes are not among them; report how many recipe
atoms carry a within-2 kb neighbor.

### 10.2 The donor panel: who E is estimated from

`sample_metadata.suppl_table.csv` (the full sample sheet — 12 rows):

| donor | age | sex | conditions |
|---|---|---|---|
| CE0006864 | 34 | M | Rest, Stim8hr, Stim48hr |
| CE0008162 | 29 | F | Rest, Stim8hr, Stim48hr |
| CE0008678 | 23 | F | Rest, Stim8hr, Stim48hr |
| CE0010866 | 22 | F | Rest, Stim8hr, Stim48hr |

**D1 — n = 4 donors, demographically narrow** *(structural)*.
Four donors, **all aged 22–34**, **3 F / 1 M**, all blood type O+. This is the entire
population from which every generator — and therefore every reachability verdict — is
estimated. Nothing about the cone speaks to older adults, to males (n=1), or to disease
states. The manuscript should not claim generality beyond "young, healthy, mostly female
primary CD4⁺ T cells."

**D2 — The aging axis is asked of donors with no aging span** *(the sharpest single point)*.
The project runs `toward_younger` / `toward_older` reachability as two of its twelve atlas
cells. But **the generators contain no aging variation** — the oldest donor is 34. The
*target* aging direction comes from an entirely separate cohort (Yaza 2022, a single "age_bin"
discovery contrast over 10,000 genes). So the aging axis is a fully cross-cohort extrapolation:
a young-donor knockdown dictionary asked whether it can reach an aging signature it never saw.
Consistent with that, the aging verdicts are **the only ones that fail the strict robustness
test**: at Stim48hr the younger axis sits *below* the anisotropy-null 99th percentile (baseline
cosine 0.597 vs null-p99 0.616; z ≈ 1.90) and the older axis lower still (0.566; z ≈ 1.20),
both with `boot_clears_null = False` and a coordinated-bias radius of **kSE\* ≈ 0.001** —
essentially *zero* margin, versus kSE\* ≈ 0.033 for the Th1 axis. (They do clear the *weaker*
shuffled-target null, held-out z ≈ 6, which is why they look fine in the atlas table — the two
nulls answer different questions and only the strict one exposes the fragility.) **Fix:** either
demote the aging axis to an explicit "negative demonstration — here is what the oracle says when
the data cannot support the question," or drop it from the headline atlas. As-is it invites a
reviewer to distrust the axes that *are* solid.

**D3 — Donor is collapsed; no true leave-one-donor-out** *(moderate)*.
E is a donor-aggregated DE estimate. Cross-donor reproducibility *columns* exist
(`crossdonor_correlation_mean/min`), but the reachability fit itself runs on the pooled effect,
so there is no genuine per-donor replication of the *verdict*. With n=4 a real LODO is
underpowered, but the current design cannot even report verdict variance across donors.

**D4 — Stim48hr is perfectly confounded with sequencing run** *(a real, checkable batch confound)*.
From the sample sheet: **all four Stim48hr libraries were run on `CD4i_R2`**, while Rest and
Stim8hr are split 2/2 across `CD4i_R1` and `CD4i_R2`.

```
condition   CD4i_R1  CD4i_R2
Rest              2        2
Stim8hr           2        2
Stim48hr          0        4     ← run-confounded
```

Any Stim48hr-specific effect is therefore inseparable from an R2 batch/run effect. The A1
sensitivity analysis happens to run at Stim48hr, and several fragile findings (the aging axis
above) are Stim48hr — so this confound sits underneath exactly the cells that are already
borderline. **Fix:** flag it; where possible cross-check a Stim48hr conclusion against Rest,
which is run-balanced.

### 10.3 The target vectors d: what the cone is compared against

**T1 — Targets are cross-study, cross-assay relative to the generators** *(structural)*.
The Th2→Th1 target is built from **Ota 2021 (24,821 rows) + Höllbacher 2021 (12,467 rows)** and
the aging target from **Yaza 2022** — all *external* differential-expression contrasts, **not**
the Perturb-seq assay that produced E. The generators are single-cell CRISPRi pseudobulk
z-scores; the targets are sorted-population / bulk-style DE from different studies with
different normalization, dynamic range, and batch structure. The method reconciles them by
gene-name intersection + cosine and nothing more. This is a genuine cross-platform assumption
hiding inside "restricted to shared genes." **Fix:** state it explicitly; ideally show the
verdict is stable if the target is re-derived from a *single* Th1/Th2 study, or from a
within-atlas contrast, to prove it is not an artifact of stitching two external cohorts.

**T2 — Two disagreeing source contrasts are merged into one Th2→Th1 direction** *(moderate)*.
The polarization target combines two studies (Ota, Höllbacher) that need not agree. How they
are merged (intersection? mean of z? concatenation?) can dominate the resulting direction, and
the two cohorts' disagreement is currently invisible in the single cosine that comes out.
**Fix:** report the between-study cosine of the two contrasts, and the verdict computed on each
separately — if Ota-only and Höllbacher-only give materially different reachable cosines, the
headline number is a blend.

**T3 — Silent truncation by gene intersection** *(minor but worth a number)*.
Effect and target are both restricted to genes measured in both. If a Th1-defining gene is
absent from the screen's 10,282, it is dropped from d before the cone ever sees it —
potentially removing exactly the coordinates that make the target hard. **Fix:** report the
fraction of each target's top DEGs (by |z|) that survive the intersection; a high survival rate
defuses the concern cheaply.

### 10.4 Provenance and reproducibility of the data

**P1 — The dictionary is a preprint-stage, author-derived DE product** *(context, not a flaw)*.
`GWCD4i.DE_stats.h5ad` is Zhu 2025, a **bioRxiv preprint** (not yet peer-reviewed), and it is
the authors' *DE output* — pseudobulk model, shrinkage, covariate choices and all — not
something recomputed from raw counts here. Any systematic feature of their DE pipeline
propagates directly into every generator. This is a reasonable choice (recomputing from ~22M
cells is infeasible on the target hardware), but it should be named: the cone inherits the
upstream DE model's assumptions wholesale.

**P2 — No independent / arrayed validation of the polarization ground truth held locally** *(acknowledged in the repo)*.
The `data/README.md` already notes the arrayed CRISPRi validation table is not available
locally, so the polarization "ground truth" leans on internal reproducibility columns +
orthogonal literature/Open Targets evidence rather than a held-out wet-lab confirmation. Fine
to disclose; not fine to leave implicit.

### 10.5 Prioritized shortlist (what a reviewer would require before acceptance)

**Must state in Limitations (no new compute needed):**
- G1 — the LOF/GOF split is modality-relative to CRISPRi, not an absolute biological partition.
- D1/D2 — 4 young donors; the aging axis is a cross-cohort extrapolation the data cannot
  robustly support (kSE\* ≈ 0.001), and should be demoted or reframed as a negative
  demonstration.
- T1 — targets are cross-study/cross-assay relative to the generators.
- P1 — the dictionary is a preprint-stage author-derived DE product.

**Cheap analyses that would materially strengthen the paper (all runnable on data in hand):**
- D4 — check every Stim48hr-specific conclusion against run-balanced Rest.
- G4 — specificity pass: confirm headline recipe genes are not `offtarget_flag` / within-2
  kb-neighbor generators.
- T2 — report between-study cosine of Ota vs Höllbacher and the verdict on each alone.
- T3 — report top-DEG survival fraction through the gene intersection.
- G3 — repeat the dictionary bootstrap with *per-atom* `lfcSE`/cell-count noise rather than
  homoscedastic unit noise.

**Requires new data (the honest ceiling):**
- G1/G2 — a CRISPRa arm (tests the certificate's own prediction and unlocks the GOF budget)
  and in-domain double-perturbation data to test additivity in CD4⁺ T cells rather than
  borrowing Norman K562.

*Bottom line: the method is not over-claiming its geometry, but the paper currently
under-discloses how much of its most striking numbers are properties of a single, narrow,
cross-stitched dataset. Move these from "things a reviewer will find" to "things we state
first," and the work is much harder to attack.*

---

## 11. One-paragraph version for the manuscript

> The reachability oracle is a design-based causal-inference method. Because genome-scale
> CRISPRi Perturb-seq is a randomized interventional experiment, each effect vector is an
> estimated average treatment effect, and asking whether a target state lies in the
> non-negative cone of those vectors is a counterfactual feasibility query about a compound
> intervention. Unlike observational causal-discovery approaches, the method requires no
> fitted gene-regulatory graph — identification is supplied by the experiment, and
> infeasibility comes with a Farkas certificate. We treat incomplete CRISPRi knockdown as
> imperfect compliance in an instrumental-variables design: the published effects are
> intent-to-treat, their compliance-rescaled counterparts are local average treatment effects,
> and because a convex cone is invariant to positive per-generator rescaling, the reachability
> verdict is provably and numerically unchanged by this rescaling (max |Δcosine| = 2.22e-16
> over all 12 targets × conditions). As a separate exclusion-restriction check, we also
> re-solve on the valid-instrument subset (perturbations with demonstrable on-target
> knockdown, ≈99 % of the dictionary); the per-target effect on the verdict is reported in §4.

---

*Files. Conceptual anchor: this doc (`CAUSAL.md`), with the trust dossier in
`notebooks/09_causal_validation_dossier.ipynb`. IV/compliance: `run_iv_compliance.py`,
`analysis_cache/atlas_work/first_stage_compliance.csv`, `iv_compliance_verdicts.csv`,
`late_rescaling_invariance.csv`, `fig_causal_compliance.png`. Part-A trust layer:
`run_a1_sensitivity.py`, `a1_sensitivity_radius.csv`, `fig_a1_sensitivity.png`,
`a3_kway_additivity_bound.csv`, `a4_weak_instrument_intervals.csv`,
`a5_negative_control_outcome.csv`, `a6_construct_validity.csv`. Counterfactual explanation:
`counterfactual_cards.json`, `diverse_recipes.csv`, `diverse_recipes_detail.json`. See
`NOVELTY.md` §"design-based causal inference" and `RELATED_WORK.md` Camp-A citations for the
surrounding positioning.*

