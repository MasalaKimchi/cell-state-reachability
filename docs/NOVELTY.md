# Novelty, impact, and field positioning — what is new, why it matters, and where it sits

*This is the positioning dossier for the reachability oracle: what is scientifically new
(stated as a precise delta against the right prior art), how to sharpen it, where it has
real-world impact, and how it sits relative to the field's current activity. It absorbs and
supersedes the former `IMPACT.md` (the industry-scientist case — attrition economics, the two
odds-moving levers, concrete users and decisions) and `TIDE_VS_WAVE.md` (the field-positioning
argument — the forward-prediction tide vs the interventional-decision wave). Companion to
`README.md` (framing), `RELATED_WORK.md` (the 91-method survey), `RESULTS.md` (the built method,
atlas, and validation), and `CAUSAL.md` (the design-based trust layer). Numbers, citations, and
figure captions are carried verbatim from the source docs; the provenance note in §5.7 states
which impact figures were read at sentence level versus confirmed at citation level.*

---

## 1. The one sentence — a verdict, not a ranking

> **We turn a genome-scale CRISPRi screen into a feasibility oracle for cell-state
> engineering: given a desired transcriptional shift, we return a *verdict* — reachable
> by some combination of knockdowns, or provably not — together with the smallest
> knockdown set that reaches it and a certificate for what it cannot.**

The word doing the work is **verdict**. Everything else in the field returns a *ranking*
(which gene is most associated) or a *prediction* (what will this perturbation do). We
return a **decision with a proof attached**: reachable, or here is the direction no
knockdown can produce.

![Reachability as convex geometry: a feasibility verdict, not a ranking]({{artifact:cc2d1434-0603-4ae7-9872-32d852e484f8}})

*Left: the target shift `d` lies inside the cone spanned by the (non-negative) knockdown
effect vectors — reachable, and the mixing weights are the minimal recipe. Right: `d`
lies outside the cone — no non-negative combination reaches it; the residual is a
direction that requires gene **activation**, and the separating hyperplane is the
infeasibility certificate. That right-hand panel is the output no existing method
produces.*

---

## 2. Question → geometry → verdict (the framing that earns the cone)

The machinery (convex cone, non-negative combination, residual) only makes sense once the
question it answers is on the table. Stated the wrong way round — machinery first — a reader
hits "convex cone" with no anchor for *why a cone*, and the geometry reads as decoration on top
of a ranking; it also buries the single most novel output, the *provably-outside* verdict, in a
subordinate clause. The framing that works runs **question → why the geometry is forced → what
falls out** (this is the framing to use, in preference to the machinery-first paragraph in
`CLAUDE.md §0`):

> A genome-scale CRISPRi screen measures, for every gene, the causal effect of knocking
> it down on the whole transcriptome. Read each of those measured effects as a vector,
> and a natural engineering question appears: **to move a cell from state A toward a
> target state B, which knockdowns — and in what combination — get closest, and is it
> even possible?**
>
> Two facts about CRISPRi make this question geometric. Knockdown only *removes* function
> (loss-of-function), and you can apply several knockdowns together but you cannot apply a
> "negative knockdown." So the reachable set of transcriptome shifts is exactly the
> **non-negative combinations of the measured effect vectors** — a convex cone. The
> target either lies inside that cone (**reachable** — and the combination weights are the
> recipe) or outside it (**provably unreachable by knockdown alone** — and the leftover
> direction tells you what would have to be *activated* instead).
>
> That verdict is the product. It is not a differential-expression list (what differs
> between A and B) and not a similarity score (which single gene looks most like the
> target). It is a **reachability decision with a residual and a confidence score**, and
> either answer — reachable or not — is a result you can take to the bench.

This version leads with the *question* and earns "convex cone" by first stating the two facts
that force it, rather than opening on the machinery. ("Dictionary" is fine as a technical term —
it is standard sparse-coding vocabulary for a matrix of atoms — but it belongs in the method
section, not the opening sentence a non-specialist reads first.)

---

## 3. The novelty, stated as a precise delta against the RIGHT prior art

Novelty is a claim about a neighbourhood. It is tempting to position against **Mogrify** and
**CellOracle** — both about *reprogramming factor discovery* — and stop there. That misses the
field this method actually lives in, and against which the delta is cleanest: the **network-control
theory of cell fate**, a decade-old, active line of work.

### 3.1 The comparison class that matters: control of cell fate

This field asks *the exact question* — what is the minimal set of interventions that
drives a cell to a target fate — and has a mature toolbox for it:

| Prior method | What it does | What it controls | The basis it needs |
|---|---|---|---|
| **Zañudo & Albert 2015** (*PLoS Comput Biol*) — stable-motif control | Minimal node interventions that switch a Boolean network's attractor | discrete attractors | a hand-built **Boolean model** |
| **NETISCE 2022** (*npj Syst Biol Appl*) — feedback-vertex-set control | Ranks intervention sets that steer a signalling network to a desired steady state | steady states | an **inferred/curated network** + simulation |
| **CellOracle 2023** (*Nature*) | Simulates KO/overexpression through an inferred GRN into a "vector map of transitions" | cell-identity transitions | a **GRN inferred from WT scRNA-seq + ATAC** |
| **"Reverse control of cell-fate" 2024** (*npj Syst Biol Appl*) | Derives control targets from cell-fate landscape structure | fate attractors | a **model of the fate landscape** |
| **Mogrify 2016** (*Nat Genet*) | Predicts TF set for a transdifferentiation | cell-type conversion | expression + **regulatory-network** heuristic |

**Every one of these runs on an *inferred or hand-built model* of the regulatory
network.** That is the field's structural bottleneck, and it is stated openly in the
literature: the control set you get out is only as trustworthy as the network you put in,
and GRN inference from observational data is famously underdetermined. These methods
answer "what should we perturb" by *simulating a model* of what perturbations would do.

### 3.2 The delta — four things, in priority order

**(1) The effect vectors are measured, not modelled — this collapses the field's
bottleneck.** The dictionary rows are the *empirically observed* genome-scale CRISPRi
effects on the real transcriptome. We are not inferring what a knockdown does and then
controlling the inference; we are composing measurements. This is the difference between
"control the model and hope it transfers" and "compose interventions whose effects were
directly observed in the target cell type." *This is the headline — lead with it.*

**(2) A convex-cone reachability *verdict with a certificate*, not a ranked set.** The
control-theory methods return ranked/heuristic intervention sets. None of them return a
*provable infeasibility result*. Because the reachable set is literally a convex cone,
"unreachable" is not "we didn't find a good set" — it is a theorem: the target has a
component in the orthogonal complement of the non-negative span, and the separating
hyperplane is a checkable certificate. **"Provably outside the knockdown cone → this
transition requires activation" is an output no method in the table above produces.** A
negative result becomes a positive, falsifiable claim.

**(3) Held-out-gene validity + a screen-native confidence decomposition.** We validate
by fitting weights on a random half of the target signature's genes and scoring alignment
on the held-out half — a generalisation test that cannot be inflated by overfitting, and
that the association/ranking framings structurally cannot run. Confidence decomposes into
the screen's own reproducibility metrics (cross-guide, cross-donor), stability-selection
frequency, and orthogonal literature/genetics evidence.

**(4) A modality-resolved intervention verdict — the geometry answers "by which
modality," and that answer is cross-verified against the druggable genome.** Because the
cone is signed (knockdown = LOF, sign-flipped = GOF), every per-node requirement carries a
*modality*: a node reachable in the LOF cone is a small-molecule / degrader / siRNA program;
a node that only enters after allowing activation is an agonist / cytokine / cell-engineering
program; a node reachable by neither is not a drug target at any modality. Crossing that
requirement with Open Targets tractability turns the verdict into a **triage** — and the
striking, decision-useful finding is the collision: **44% of the required knockdown nodes
across the atlas are hard-to-drug, and the single strongest genetically-supported nomination
among them (IRF1, 17 immune-disease genetic associations, top Th2 knockdown node) has no
conventional handle at all.** No ranking, association, or model-control method returns
"reachable + genetically supported + undruggable-at-this-modality," because none of them
represent the modality requirement geometrically in the first place. (The full triage over
102 real nominations — three actionable tiers with named genes — is in §5.4.)

### 3.3 What must NOT be claimed as novel (credibility depends on this)

- **Not the Th1/Th2 or aging regulators.** Zhu et al. already report these. Recovering
  GATA3 is a *positive control*, not a discovery. (`RESULTS.md` is already honest
  about this — keep it that way.)
- **Not "combining perturbation effects linearly."** Additive composition of effect
  vectors is standard; the 2025 *Nat Methods* baseline paper (Ahlmann-Eltze et al.) is
  a friend here — it licenses the linear model as the right primary model, but linear
  composition itself is not the novel part.
- **Not "convex optimisation on omics."** NNLS/LASSO on expression data is old. The novel
  composite is *measured LOF effect vectors + non-negativity-as-biology + a feasibility
  verdict with a certificate*, applied to cell-state engineering.

### 3.4 The defensible novelty sentence (drop-in)

> Prior cell-fate control methods (stable-motif control, NETISCE, CellOracle, reverse
> control) compute intervention sets by simulating an **inferred model** of the regulatory
> network. We instead compose the **measured** genome-scale CRISPRi effect vectors
> directly, and — because knockdown is loss-of-function and combinations are non-negative
> — cast reachability as membership in a convex cone. This yields what model-simulation
> approaches do not: a **provable reachable / outside-the-cone verdict with an
> infeasibility certificate**, validated by held-out-gene generalisation. We claim novelty
> for this measured-effect, certificate-carrying reachability decision — not for the
> regulator lists, which the source screen already reports.

That sentence, as a figure — screen → convex cone → verdict + activation certificate →
cross-dataset transfer:

![Figure 1 — the reachability oracle end to end. (A) A genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells gives a measured effect vector per knockdown (no inferred network). (B) Reachability = membership in the convex cone of non-negative combinations of those measured vectors; the Th2→Th1 target splits into an in-cone reachable component (teal) and an out-of-cone residual (amber), with "outside" certified by a separating hyperplane (Farkas/KKT). (C) The Th2→Th1 verdict: partly reachable (held-out cosine 0.448, null z ≈ 24, KKT residual 1.1×10⁻¹¹), signed decomposition LOF 39% / GOF 25% / neither 35%, and an activation certificate naming falsifiable CRISPRa hypotheses. (D) The operator transfers unchanged to the Norman 2019 K562 CRISPRa screen (held-out CEBPA cosine 0.878). This is the measured-effect, certificate-carrying reachability decision the sentence above claims as novel.]({{artifact:724c5ece-a685-440a-8850-b3617cb1adcb}})

### 3.5 Positioning vs prior work — the full survey

The table in §3.1 is the argument in miniature. The full, citation-grounded version — **91
prior methods across three research communities, each verified against a primary database record** —
is in [`RELATED_WORK.md`](RELATED_WORK.md), with a machine-readable capability matrix in
`method_comparison_matrix.csv` and every DOI/PMID in `references.csv`. The one figure to carry
from it is the landscape map: it shows that the two capabilities our method combines —
*grounded in measured effects* and *answers the feasibility question with a certificate* — are
**disjoint across the entire prior literature.**

![Where a measured-effect reachability oracle sits in the cell-fate-control landscape. Perturbation-ML methods (amber) use measured data but only predict responses; network-control methods (blue) reach toward feasibility but run on inferred or hand-built models; GRN-inference methods (green) sit lower-left. The measured-and-certified corner (top-right) is empty before this work.]({{artifact:891d8a45-8d5c-44fe-91e3-abab7ff9e434}})

The counts that back the claim, over the **91 prior methods** (2011–2026, ours excluded): **14**
operate on measured effects, but all 14 only predict or rank; **0** return a feasibility verdict on
measured data, **0** emit a data-grounded certificate, and **0** can return a provable-infeasibility
result. The 13–14 "partial" verdict/certificate entries are model-theoretic (network-control and a
few GRN/perturbation-prediction methods), whose guarantees hold only inside an inferred or hand-built
model, not on measured data. Our row is the only entry at the measured-and-certified intersection —
which is the precise, defensible meaning of "new" here.

---

## 4. How to make the novelty *stronger* (ranked by value-per-effort)

These are ways to widen the delta, not just defend it. ★ = highest-value, do-this-first;
the rest are the "future work" slide.

> **Status update (expansion round).** Items (a), (b), (c), and (e) are **done** and
> are now the backbone of the result: the activation certificate ships as a
> first-class constructive output (a); both nulls (shuffled-target + held-out-gene,
> plus a random-perturbation null) gate every one of the 12 atlas verdicts (b); the
> directional decomposition became the full **signed LOF/GOF/neither split**
> (`signed_reachability()`, c); and condition-resolved cones are the atlas's second
> axis — Rest is measurably more reachable than Stim8hr/Stim48hr for every transition
> (e). The expansion also added a delta not on this list — **modality triage** (see
> §3.2 item 4). Items (f), (g) remain future work; (g)'s low-dimensional-cone
> caveat is now quantified and stated in `RESULTS.md`.
>
> **Status update (post-hackathon round).** Item **(d) is done** — with a course
> correction worth stating: the hypothesised mechanism (collinearity → sub-additivity)
> was **refuted** on the Norman doubles, and the data-supported mechanism, **magnitude
> saturation**, was calibrated instead and shipped as `additivity_risk()` +
> `epistasis_penalty`. The positive control was also broadened from 2 genes to a
> regulator panel (master TFs give AUROC = 1.00), and the empirical shuffled-target null
> now has a **closed-form anisotropy correction** (`analytic_anisotropy_null()`) that
> ties directly to the resolution caveat in (g): it quantifies *exactly* how much of a
> "partially outside" verdict is the free anisotropy floor rather than signal. See
> `RESULTS.md` §6 for all three. Item (f) — cross-cone transfer to a second atlas —
> remains the standing future-work headline (the Norman K562 CRISPRa demo is a
> single-target cross-dataset step toward it, not the full atlas-level claim).

**★ (a) Ship the infeasibility certificate as a first-class output, not a residual
number.** Rather than leaving "outside the cone" as a scalar residual, make it
*constructive*: when the target is unreachable, return the **specific genes in the
residual direction that would need activation** — i.e. the dual variables / the
separating hyperplane's heavy coordinates. That converts "not reachable" into "not
reachable by knockdown; here are the 3 genes a CRISPRa arm would need to test," which is
a concrete, wet-lab-actionable hypothesis. **This is the single highest-value addition**
— it is the output that most sharply distinguishes the method from every ranking method, and it
directly motivates the Schmidt 2022 CRISPRa validation already cited.

**★ (b) Report the reachability *spectrum* as the headline figure, with both nulls.**
Alignment-vs-k against a shuffled-target null *and* a random-perturbation null. This is
what makes the verdict honest — the null band is the achievable-by-chance floor at every
sparsity level — and it is the figure that survives expert scrutiny. Without it, any single
cosine is dismissable as overfitting.

**★ (c) Make non-negativity earn its keep with a directional decomposition.** Split the
target `d` into its *in-cone* (knockdown-reachable) and *out-of-cone* (activation-
requiring) components and report the fraction of the target norm in each. "72% of the
Th2→Th1 shift is reachable by knockdown; 28% requires activation, concentrated in
{TBX21, STAT4, ...}" is a far richer, more falsifiable claim than a single verdict, and
it is a direct consequence of the geometry already in hand.

**(d) Turn the additivity caveat into a model feature. ✅ done (post-hackathon).**
The original plan was to penalise co-selecting genes with *collinear* effect vectors
(shared downstream program → predicted sub-additivity). Calibrating against the 126
measured Norman doubles **refuted that specific mechanism** (Spearman ρ ≈ −0.16 with
non-additivity, n.s.; collinear pairs are if anything more additive) and identified the
real one: **magnitude saturation** — the combined effect magnitude falls below the
additive sum, growing with combined magnitude (ρ ≈ +0.58). This ships as a per-recipe
`additivity_risk()` score and an `epistasis_penalty` on `reachability_spectrum` (risk
validates against measured deficits at ρ = +0.46, p = 5.6e-8). Re-annotating the atlas
showed its recipes sit safely in the additive regime (risk 0.04–0.08 ≪ 0.5), so the
additivity assumption is now *calibrated and validated* rather than merely asserted. The
designed-experiment queue (rank the next-most-informative double) remains the natural
next step on top of this calibration.

**(e) Condition-resolved cones.** Rest / Stim8hr / Stim48hr are available. Ask whether a target
is reachable in one condition but not another — "reachable only in activated T cells" is a
biologically meaningful and novel kind of statement, and it uses data already in hand.

**(f) Cross-cone transfer as the generality claim.** Nothing in the method is
T-cell-specific. Demonstrating the *same* machinery on a second public Perturb-seq atlas
(even a small one) upgrades the contribution from "a T-cell result" to "a general operator
for perturbation atlases." This is the strongest *scientific* claim and the natural
future-work headline.

**(g) Bound the cone's resolution honestly.** With ~2,000 perturbations spanning a
10,000-gene space, the cone is low-dimensional relative to the ambient space, so *most*
random targets will look "partially outside." State this — it is why the shuffled-target
null (b) is essential, and acknowledging it pre-empts the sharpest reviewer objection.

**(h) Name the method as design-based causal inference — and ship its trust layer. `[built]`**
The strongest *framing* upgrade, and it is done. The oracle is not an observational model:
the effect vectors are average treatment effects from a randomized CRISPRi experiment, and a
reachability verdict is a counterfactual query on a compound intervention. Saying so places the
method in the **design-based** camp (Fisher/Neyman/Rubin, LATE/IV) rather than the
observational-discovery camp (PC/GES, NOTEARS, causal representation learning) — the camp whose
identification does *not* rest on a fitted graph, which is precisely the field's weak point that
the headline delta (§3.2, point 1) already exploits. The novelty this unlocks is not the reframe
alone but the **trust layer** it makes natural: a verdict sensitivity radius in measured-SE units
(an E-value for a geometric feasibility claim — A1), an instrumental-variables treatment of
incomplete knockdown as noncompliance (ITT vs LATE, with the cone's rescaling-invariance making the
verdict provably compliance-robust), negative-control-outcome and signed-construct-validity
falsifications (A5/A6), and weak-instrument-robust recipe intervals (A4). No prior
cell-fate-control method ships a *calibrated* feasibility verdict with its assumptions enumerated
and stress-tested. See [`CAUSAL.md`](CAUSAL.md) and
`notebooks/09_causal_validation_dossier.ipynb`.

---

## 5. Real-world impact — where a feasibility oracle is actually useful, and to whom

The impact is not "we found a drug target." It is that **a reachable/unreachable verdict
changes what a lab does next**, and it does so before any wet-lab spend. Making that case to
an industry scientist needs three things established first: the problem is large and expensive,
the specific failure mode this method attacks is the dominant one, and the "measured, not
inferred" design is aligned with the only two things shown to move drug-development odds.

### 5.1 The problem is failure, and failure is the cost

Drug development is not primarily expensive because success is expensive — it is expensive
because **failure is the norm and it is paid for late.** The canonical figures:

- **~90% of clinical programs never reach approval.** Only about 10% of clinical
  programmes eventually receive approval (Minikel et al., *Nature* 2024; the BIO/Thomas
  2006–2015 analysis puts the cumulative rate at 9.6%). The cost of drug discovery and
  development is *driven primarily by failure*.
- **The failure concentrates at Phase II, and its cause is efficacy.** Phase II has the
  lowest transition rate of any phase — **30.7%**, versus 63.2% (Phase I) and 58.1%
  (Phase III) in the BIO dataset (Thomas et al. 2016). And **lack of efficacy is the
  reason for approximately half of all Phase II and Phase III failures** (Harrison, *Nat
  Rev Drug Discov* 2016, analysing 174 attrition events with a stated reason over
  2013–2015) — i.e. wrong-target / wrong-mechanism, not safety or operations, is the
  dominant late-stage failure mode.
- **Autoimmune indications — this dataset's home turf — sit near the bottom.** Wong, Siah
  & Lo (*Biostatistics* 2019, 406,038 trial records) put autoimmune candidates at roughly
  **15% probability of success**, well below the cross-indication average.

**The implication for a computational method:** the highest-leverage place to intervene is
*target selection* — the decision, made years before Phase II, of *which* gene to go
after. A wrong target is an efficacy failure waiting to happen, and it is only discovered
after the most expensive experiments have run.

### 5.2 Two things move the odds — and both are "measured, not inferred"

The literature is unusually clear about what *reduces* efficacy failure, and both levers
point the same way as this project's core design choice.

**Human genetic support ~2.6× the odds of approval.** Drug mechanisms with human genetic
support are **2.6 times more likely** to succeed from Phase I to approval than those without
(Minikel et al., *Nature* 2024); the earlier estimate (Nelson et al., *Nat Genet* 2015) was
~2×. Critically, the authors frame genetic support not as a gate but as **an enrichment signal
and a probabilistic tool for portfolio prioritisation.** That is exactly the register a
reachability oracle operates in: it does not promise a target works — it *re-weights the
portfolio* toward mechanisms with independent, measured support. *This project's disease layer
already carries this signal:* the local autoimmune-enrichment table links perturbation clusters
to **17 immune-mediated diseases with 185 significant GWAS-gene enrichments (FDR < 0.05)**. A
nomination that is *both* knockdown-reachable *and* genetically supported for the target disease
is precisely the "high relative-success" quadrant Minikel describes.

**Inferred preclinical claims don't reproduce — measured effects are the antidote.** The
reason "wrong target" is so common upstream is that the published preclinical literature it
draws on is not reliable. The two landmark industry audits:

- **Amgen: 6 of 53 (11%)** landmark preclinical cancer papers could be reproduced (Begley
  & Ellis, *Nature* 2012).
- **Bayer: ~25%** of target-validation projects reproduced cleanly (Prinz et al., *Nat Rev
  Drug Discov* 2011).

This is the crux of the "measured, not inferred" argument in §3.2, stated in business terms:
**the field's target ideas are largely built on inferred networks and individually-published
effects that fail to replicate.** A method whose dictionary rows are *directly measured
genome-scale CRISPRi effects in the relevant primary human cell type* is drawing on exactly the
kind of systematic, internally-reproducible measurement (cross-guide, cross-donor) that the
inferred literature lacks. That is not a nice-to-have — it is the specific defect that makes 90%
of programs fail.

**The one-sentence pharma pitch:** *the two interventions shown to cut efficacy failure —
human-genetic support and measured (not inferred) target effects — are the two things this
method is built on. It scores knockdown nominations by measured effect and cross-references
them to disease genetics, in the primary cell type where the biology actually happens.*

### 5.3 Pain point vs current practice — what the existing toolbox cannot do

The attrition numbers above describe *the cost of a wrong target*. The computational
toolbox a discovery team currently reaches for to avoid that cost does not, in fact, answer the
question that would avoid it. The **91 verified prior methods** (2011–2026) across
cell-fate-control theory, GRN inference, and single-cell perturbation modelling — full comparison
in [`RELATED_WORK.md`](RELATED_WORK.md) and the machine-readable `method_comparison_matrix.csv` —
yield this structural finding, in business terms:

| What a target-selection team needs | What current methods deliver | The pain point |
|---|---|---|
| An answer grounded in **measured** effects in the relevant cell type | 14 of 91 methods use measured data — **but all 14 only predict or rank responses** | measured-data methods never say *whether a goal is reachable* |
| A **yes/no achievability** answer, not another ranking | **0 of 91** return a feasibility verdict on measured data | teams get optimistic rankings; every tool only ever says *go* |
| A credible **"stop — not reachable by this modality"** | **0 of 91** can return a provable-infeasibility result | the expensive *stop* decision (the STOP verdict, §5.4) has no computational support |
| Achievability reasoning that isn't hostage to an **inferred model** | the 13 methods that reason about controllability run on inferred/hand-built networks | the control set is only as good as the network, which is underdetermined |

The pain point in one sentence: **current practice can predict what a perturbation will
do, or reason about control inside a model it cannot fully trust — but it cannot tell a team,
from measured data, that a desired cell-state change is achievable, and it certainly cannot
tell them to stop before the modality is proven wrong.** Because late efficacy failure is the
dominant cost (§5.1), the single missing capability is also the most expensive one to lack. The
reachability oracle is built to fill exactly that hole: a measured-grounded reachable /
provably-outside verdict, with a constructive certificate for what is missing and a
modality-resolved triage over the druggable genome (§5.4). *The methods-level evidence for each
"0 of 91" above is in [`RELATED_WORK.md`](RELATED_WORK.md) §2 and the capability matrix.*

### 5.4 Where it is useful — concrete users and decisions

The use cases below reduce to one decision diagram and one funnel: *what the oracle tells a
team to do* (GO / STOP / REDIRECT, then a modality triage), and *where each value lever acts*
on the attrition funnel.

![How the reachability oracle is useful to a drug-development team. (A) The decision it makes — a desired cell-state change is routed to GO (reachable → minimal ranked knockdown set → arrayed screen), STOP (provably outside the cone → don't spend a CRISPRi arm), or REDIRECT (activation-required → switch to CRISPRa; the certificate names the genes to test). The GO branch feeds a modality triage over 102 real nominations (44% hard-to-drug; 10 with a clinical-grade drug): green-light (JAK2, ICOS, MAPK14, CD3D), tractable-but-untried (IL7R, ZAP70, TET2), and required-but-undruggable (IRF1, 17 immune-disease genetic associations, no conventional handle). (B) Four value levers mapped onto the drug-development attrition funnel — measured-not-inferred effects (attacks 11–25% preclinical reproducibility), disease-genetics cross-reference (2.6× approval odds), the provable "unreachable" STOP (redirects before the expensive phase; ~half of Phase II/III failures are efficacy), and the minimal ranked recipe (allocates scarce screen slots). ≈10% of programmes reach approval; Phase II is the lowest transition at 30.7%. Nominations are wet-lab hypotheses, not validated targets.]({{artifact:3a67feb4-18d0-4fbd-b88d-69317752dfdd}})

**Target triage for T-cell engineering & cell therapy (the nearest application).**
The featured dataset is primary human CD4+ T cells, and **T-cell state engineering is a
market in the middle of a step-change.** CAR-T, proven in oncology, is moving into
autoimmunity fast: as of December 2024 there were **116 clinical trials** evaluating CAR-T
against autoimmune conditions (Clinical & Translational Science, 2025), and trials in
autoimmune rheumatic disease alone **peaked at 25 in 2024** (Frontiers in Immunology,
2025), with early SLE cohorts reaching durable, steroid-free remission.
A programme deciding *which knockdowns* to engineer into a T-cell product — for
persistence, for Treg stability, for Th-bias correction — faces a combinatorial design
space and a fixed screening budget. **This is the decision the reachability oracle is
built for:** "this state shift is reachable by this ranked 3-gene knockdown set — put it in
the arrayed screen" and, just as valuably, "this shift is *not* reachable by knockdown —
don't spend a CRISPRi arm on it." The users are cell-therapy target-selection teams and
academic immuno-engineering labs; the decision is which constructs enter a screen that
costs weeks and animal cohorts per arm.

**Autoimmune / allergic disease — the Gladstone-prize hook.**
The local data already links perturbation clusters to the **17 immune-mediated diseases with
185 significant GWAS-gene enrichments** of §5.2 (Crohn's, IBD, asthma, atopic eczema, psoriasis,
MS, T1D, lupus, RA, ...). Th1/Th2 balance is a validated therapeutic axis — Th2-skew
drives allergy/atopy, Th1/Th17 drives autoimmunity. A method that nominates the *minimal,
ranked, confidence-scored* knockdowns to rebalance that axis — **and states when a
rebalancing is not achievable by knockdown at all** — is a target-selection instrument for
immunotherapy, with the disease link already in the data rather than asserted.

**The "don't run this experiment" verdict — the underrated value.**
Every ranking or association tool only ever says *go*. The distinctive output here is a
**provable stop**: when the target lies outside the knockdown cone, the tool returns not
just "no" but the **specific genes that would have to be *activated* instead** (the
constructive certificate in `reachability.py`). For a lab, that converts a dead end into a
redirected experiment: *stop the loss-of-function screen, switch to CRISPRa / an agonist,
and here are the 3 genes to test first.* Given that late efficacy failure is the dominant
cost (§5.1), **a credible early "this modality cannot reach the goal" is worth more than
another optimistic ranking** — it moves resources before they are sunk. Negative results
that redirect resources are where a decision tool earns its cost, and this is the use case
no ranking or association method can serve.

**Portfolio prioritisation across a Perturb-seq atlas.**
Framed the way Minikel frames genetic support — a probabilistic prioritisation signal, not
a gate — the reachability verdict + confidence decomposition is a **portfolio-ranking
instrument.** For any (target state, disease) pair, it returns a reachable/outside verdict,
a minimal ranked set, a screen-native confidence score, and a disease-genetics
cross-reference. A discovery team can rank a whole slate of desired cell-state changes by
*how reachable they are and how well-supported the drivers are* — before committing a
single wet-lab FTE.

**Modality triage against the druggable genome (delivered, not promised).**
The expansion turned the "which modality" verdict into a **cross-verified triage
over 102 real nominations**, and the empirical result is the sharpest impact claim in
this document. Every knockdown node the atlas nominates was crossed against Open
Targets tractability and immune-disease human genetics (see `RESULTS.md`). The finding:
**44% (45/102) of the required knockdown nodes are hard-to-drug**, only 10 have a
clinical-grade drug, and — the decision-relevant collision — the strongest
genetically-supported nomination among the hard-to-drug set, **IRF1** (17 immune-disease
genetic associations, a top Th2-axis knockdown node), has *no conventional handle at
all*. For a portfolio team this is exactly the signal that saves spend: it separates three
actionable tiers: **green-light** nominations already carrying a clinical-grade drug
(JAK2, ICOS, MAPK14, CD3D — strong immune genetics *and* approved/candidate drugs);
**tractable-but-untried** nominations with a plausible modality and strong genetics but
no drug yet (IL7R, antibody-tractable, the single highest genetic support; ZAP70 and
TET2, SM-tractable) — the highest-value *new* leads; and the
**required-but-undruggable** nominations (IRF1 and the other degrader-only/undruggable
nodes) that should route to a degrader-discovery or cell-engineering effort *before* a
doomed small-molecule campaign is funded. This is the genetic-support odds argument of §5.2 and the modality-stop
argument (STOP verdict, above), fused into one triage table.

**A reusable operator for the Perturb-seq era (the general claim).**
Nothing in the method is T-cell-specific. Perturb-seq atlases are proliferating across
tissues and disease models; any atlas + any target state (disease→healthy, aged→young,
cell-type A→B) plugs into the same convex-cone machinery. The durable contribution is a
**general reachability operator for perturbation atlases** — the transferability, not the
one T-cell result, is the scientific and commercial headline.

### 5.5 Quantified value proposition (order-of-magnitude, honestly framed)

The value is not a claimed hit rate — it is **shifting where a program spends and where it
stops.** In the terms of §5.1–§5.2:

| Lever this method pulls | The number it speaks to | Direction |
|---|---|---|
| Prioritise measured over inferred target effects | 11–25% preclinical reproducibility (Amgen/Bayer) | attacks the root cause of wrong-target efficacy failure |
| Cross-reference nominations to disease genetics | 2.6× approval odds with genetic support (Minikel 2024) | moves the portfolio toward the high-success quadrant |
| Provable "unreachable by knockdown" verdicts | ~48% of Phase II failures are efficacy; autoimmune PoS ~15% | redirects modality *before* the expensive phase |
| Rank minimal knockdown sets for a screen | screen slots / animal cohorts per arm | allocates the scarcest preclinical resource |

**The framing that survives scrutiny:** this is a **hypothesis-prioritisation and
experiment-triage instrument**, not a target-validation engine. It tells a team where to
look first and where not to look at all, using measured effects and disease genetics — the
two signals with published evidence of improving odds. It does not replace the wet lab; it
decides what the wet lab does next.

### 5.6 Honesty guardrails (state these plainly — they are what make the impact credible)

- **Nominations are hypotheses for wet-lab testing, not validated targets.** A
  transcriptomic-signature match is not proof of functional rescue, let alone clinical
  efficacy.
- **CRISPRi is loss-of-function only.** The tool can nominate knockdowns and *flag* what
  needs activation; it cannot test the activation hypotheses itself (that is the CRISPRa
  validation arm, e.g. Schmidt et al. 2022).
- **Multi-gene sets assume additivity (no epistasis).** Every combination is an
  extrapolation to be tested, and should be flagged as such.
- **One primary-cell system, four donors.** External validity across genotypes, tissues,
  and disease contexts is a claim to be earned, not assumed.

An oracle that over-promises is worthless to a drug-development team — the honest framing
above is precisely what makes the (real, large) upside believable to an expert.

### 5.7 References (impact figures)

1. Minikel, Painter, Dong & Nelson. *Refining the impact of genetic evidence on clinical
   success.* **Nature** (2024). doi:10.1038/s41586-024-07316-0. — verified quotes: "only
   about 10% of clinical programmes eventually receiving approval"; "probability of
   success for drug mechanisms with genetic support is 2.6 times greater than those
   without."
2. Nelson et al. *The support of human genetic evidence for approved drug indications.*
   **Nat Genet** 47, 856–860 (2015). — original "genetic evidence doubles the success
   rate" (~2×) estimate, as restated in ref [1].
3. Wong, Siah & Lo. *Estimation of clinical trial success rates and related parameters.*
   **Biostatistics** 20(2):273–286 (2019). doi:10.1093/biostatistics/kxx069. — 406,038
   records over 21,143 compounds; autoimmune PoS ~15%; oncology 3.4%.
4. Thomas et al. *Clinical Development Success Rates 2006–2015.* BIO / Biomedtracker /
   Amplion (2016). — 9.6% cumulative success; Phase II transition 30.7% (vs Phase I 63.2%,
   Phase III 58.1%).
5. Harrison. *Phase II and phase III failures: 2013–2015.* **Nat Rev Drug Discov** 15,
   817–818 (2016). doi:10.1038/nrd.2016.184. — verified quote: "lack of efficacy is the
   reason for approximately half of all phase II and phase III failures" (174 failures
   with a stated reason).
6. Begley & Ellis. *Drug development: raise standards for preclinical cancer research.*
   **Nature** 483, 531–533 (2012). doi:10.1038/483531a. — Amgen reproduced 6 of 53 (11%)
   landmark papers.
7. Prinz, Schlange & Asadullah. *Believe it or not: how much can we rely on published data
   on potential drug targets?* **Nat Rev Drug Discov** 10, 712 (2011). — Bayer ~25%
   reproducibility.
8. CAR-T in autoimmunity: 116 trials as of Dec 2024 (Clinical & Translational Science,
   2025); autoimmune-rheumatic-disease trials peaked at 25 in 2024 (Frontiers in
   Immunology, 2025).

*Provenance note (honest about verification level). Read at sentence level from the
primary source this session: the 2.6× / ~10%-approval figures (ref 1), the ~50%-of-
failures-are-efficacy figure (ref 5, Harrison 2016), Amgen 6/53 (ref 6), and the CAR-T
counts of 116 and 25-in-2024 (ref 8). Consistent with well-known published figures but
confirmed here only at title/citation level, not re-read from the source text this
session: the 9.6% cumulative rate and the 30.7% / 63.2% / 58.1% phase-transition rates
(ref 4, BIO/Thomas 2016), the ~15% autoimmune PoS (ref 3), the ~2× Nelson 2015 estimate
(ref 2), and Bayer ~25% (ref 7). Confirm the latter set against the publisher record
before any external submission.*

---

## 6. Field positioning — tide vs wave

This section answers a harder question a senior judge or reviewer will ask: **what is the
overarching problem in this field, which current activity is a passing tide and which is the
durable wave, and what would make this work relevant to the bottleneck rather than to a
leaderboard.** (Strategic companion to `RELATED_WORK.md` — the 91-method survey — and
[`CAUSAL.md`](CAUSAL.md) — the trust layer.)

### 6.1 The one-paragraph answer

The field is drowning in interventional data and has spent most of its effort on the wrong
question. The loud, well-funded activity — **AI Virtual Cell foundation models that predict a
cell's transcriptome after a perturbation** — is a *tide*: it rises on scale and attention,
and as of 2025 it does not beat deliberately simple linear baselines at its own benchmark. The
*wave* underneath it is quieter and structural: the shift from **observational association to
design-based, interventional, counterfactual reasoning** — using AI as an *instrument that
interrogates measured interventions and returns a decision with a proof and its assumptions
stated*, rather than as an *oracle that emits a predicted expression vector*. This repo is a
specimen of the wave. The highest-value work now is not to make the verdict more accurate; it
is to **sharpen the wave-claim against the tide, and to close the falsification loop** so the
instrument earns its trust at the bench. Five concrete moves do that (§6.5).

### 6.2 The overarching problem

Perturb-seq gives biology something it almost never has: **randomized interventions at genome
scale**, in the native cell type, read out transcriptome-wide. That is an extraordinary
substrate for causal reasoning. The field's response has been to pour its energy into a single
task shape:

> **Forward prediction.** Given a perturbation `g`, predict the resulting expression profile
> `Y(g)` — ideally for a gene, cell type, or dose never seen in training.

This is the task behind scGPT, scFoundation, Geneformer, STATE, Tahoe-x1, the Virtual Cell
Challenge, and the diffusion/flow-matching wave of 2025–26. It is a legitimate task. But it is
**not the question a drug or cell-engineering program actually faces**, which is inverse,
combinatorial, and decision-shaped:

> **Inverse feasibility.** Given a *target cell state* I want to reach, what is the minimal set
> of interventions that gets there, **is it even possible**, and **how much can I trust that
> answer**?

The bottleneck is not model capacity. It is that the field has been optimizing forward
prediction accuracy — a quantity that (a) simple baselines already match, and (b) does not
convert into a decision — while the inverse, feasibility-shaped, trust-calibrated question that
a lab would actually spend money on remains under-served. **Prediction is not
understanding, and neither is a decision.**

### 6.3 Tide vs wave — the call, with the evidence

**The tide: the AI Virtual Cell as a prediction oracle.**
The dominant, best-capitalized activity is the "AI Virtual Cell" (AIVC) — foundation models
that promise to predict perturbation responses and serve as in-silico screens (Bunne et al.,
*Cell* 2024; Rood et al., *Cell* 2024; the Arc Institute Virtual Cell Challenge, Roohani et
al., *Cell* 2025). The 2025–26 literature is a flood of these: STATE (Adduri et al. 2025),
Tahoe-x1, xVERSE, and generative variants (Squidiff, CellFlow, PerturbDiff, scDFM).

Why it reads as a **tide, not (yet) a wave**:

1. **It loses to linear baselines at its own game.** Ahlmann-Eltze, Huber & Anders
   (*Nat. Methods* 22:1657–1661, 2025; doi:10.1038/s41592-025-02772-6) benchmarked five
   foundation models + two other DL models against deliberately simple baselines for single-
   and double-perturbation prediction. **None outperformed the baselines.** For unseen genes,
   DL did not beat "predict the training mean," and a simple linear model reliably won. Their
   stated hypothesis: the pre-training data is *observational*, so the models learn
   associations, not intervention effects.
2. **The community's own competition confirms it.** The Arc Institute Virtual Cell Challenge
   2025 wrap-up: "Purely AI-based approaches did not consistently outperform statistical
   baselines"; the winning entries *fused* DL with classical statistical features, and pure
   end-to-end learning "is yet to solve this problem."
3. **The evaluations don't transport.** A Dec-2025 AIVC review ("Insights into Artificial
   Intelligence Virtual Cells," arXiv:2510.12498) notes evaluations remain predominantly
   within single datasets, transport across labs/platforms is limited, some splits are
   vulnerable to leakage, and dose/time/combination effects are not systematically handled.
4. **The models capture associations, not mechanism** ("Virtual Cells as Causal World Models,"
   AI4D3 @ NeurIPS 2025). Observational pre-training cannot, by construction, identify an
   intervention effect (DoFormer, bioRxiv 2026.05.02.722054, makes the same point:
   observational RNA-seq alone is insufficient for causal modeling; perturbational data is
   essential).

None of this means foundation models are worthless — they may become the wave once
intervention-aware and transport-tested. It means **the current leaderboard is a tide**:
scaling a forward-prediction objective that baselines already saturate and that does not answer
a decision.

**The wave: observational → interventional, prediction → decision, oracle → instrument.**
The durable shift — visible across the causal-ML perspective in *Nature Genetics*
(57:797–808, 2025; doi:10.1038/s41588-025-02124-2), the design-based experimental-design work
(IterPert, Genentech), and the intervention-aware modeling turn (DoFormer, CINEMA-OT,
the Wisconsin Perturb-seq mediation thesis 2025) — is this:

| | Tide (oracle) | Wave (instrument) |
|---|---|---|
| **Object returned** | predicted expression vector | a **decision** (reachable / not) with a **proof** and a **residual** |
| **Data used as** | observational corpus to pre-train | **randomized interventions** to compose |
| **Question** | forward: `g → Y` | inverse: `target state → intervention set + verdict` |
| **Trust** | held-out accuracy on one dataset | stated **identifying assumptions** + calibrated **robustness radius** |
| **Failure mode** | confident wrong prediction | an **explicit, falsifiable** infeasibility claim |
| **What ships** | the trained model | **the understanding**: a map + its falsification machinery |

This is precisely the **AI-as-instrument** thesis: computation used to interrogate the
structure of measured interventional data, where the *understanding* is the product — not a
supervised model whose emitted labels are the product. AI derives real insight only inside a
falsification loop; outside one it derives confident nonsense faster than any prior tool.
**This repo is already on the wave.** The remaining work is to make the wave-claim unmistakable
and to close the loop.

### 6.4 What the repo already gets right (so §6.5 builds, not repeats)

- **Feasibility verdict + certificate**, not a ranking or a prediction — the one output no
  method in the 91-method survey produces on measured data.
- **Measured, not modelled, effect vectors** — sidesteps the exact weakness (observational
  pre-training) that sinks the foundation models.
- **Design-based causal framing + a built trust dossier** (A1 sensitivity radius in
  measured-SE units, A4 weak-instrument intervals, A5 negative-control outcomes, A6 signed
  construct validity, the IV/compliance layer; see §4h and [`CAUSAL.md`](CAUSAL.md)) — the
  assumption-explicit, robustness-calibrated posture the field is only now converging on.
- **Modality triage** (signed LOF/GOF/neither) crossed with tractability.
- **Honest, graded transfer** (direction ports; recipe does not, Jaccard 0.11) — stated as a
  bound, not oversold.

That is a strong position. The five moves below are what turn "a strong T-cell result with a
good causal story" into "a statement about the field's bottleneck."

### 6.5 Five moves that make this relevant to the bottleneck (ranked)

Each is tied to a current field signal, uses data/machinery already in hand, and is a *wave*
move (sharpens the idea or closes the loop) rather than a *tide* move (chases a metric). These
extend the "make it stronger" list of §4 with field-relevance framing; where a move builds on an
item already delivered there, it is cross-referenced rather than restated.

**Move 1 — Draw the line against the tide, explicitly. `[framing, do first]`**
`RELATED_WORK.md` positions against network-control methods. The positioning that matters most
in 2025–26 is **against the AIVC/foundation-model prediction paradigm.** The argument writes
itself and is devastating in the best way:

> The field's flagship models optimize forward prediction — and lose to linear baselines at it
> (Ahlmann-Eltze 2025; VCC 2025), because observational pre-training identifies association,
> not intervention. We take the *one thing that does work* — linear composition of measured
> effects — and redirect it at the question that actually carries a decision: **inverse
> feasibility with a certificate.** Same linear insight the benchmarks vindicated; a question
> worth answering.

This reframes Ahlmann-Eltze from "a citation that licenses our linear model" (its role in §3.3
above) into **the field-level indictment that motivates the whole reframe.** It is the single
highest-leverage paragraph available, and it costs nothing to run.

**Move 2 — Close the loop: turn the certificate into an experiment-design acquisition function. `[build]`**
The field's design-based frontier (IterPert; "active learning for optimal intervention design")
is converging on the right question: **not "predict the atlas" but "which intervention should I
measure next?"** This method is uniquely positioned to answer it *principledly*, because the
geometry already names what is missing: the **residual direction of an infeasibility certificate
is exactly the most informative next generator to add to the screen.** Concretely: given the
current cone and a target, rank candidate next perturbations by how much each would (i) shrink
the infeasibility residual, or (ii) tighten the verdict's robustness radius (A1). That converts
the static oracle (building on the constructive certificate of §4a) into a **closed-loop screen
designer** — the inverse of the existing "don't run this experiment" value. It is the difference
between a verdict and an *instrument a screening program plans around*: the tool decides what to
ask next.

**Move 3 — Elevate the modality gap into a critique of how the field generates data. `[framing + one analysis]`**
Across the 12-cell atlas, knockdown is *never the majority modality* (mean LOF fraction 0.34;
headline Th2→Th1 only 39% in the LOF cone). That is easy to treat as a limitation of the method.
**It is actually a finding about the field's instrument.** As a field-level claim:

> Most cell-state-engineering goals are not reachable by loss-of-function at all. A field whose
> perturbation atlases are overwhelmingly CRISPRi is therefore structurally blind to the
> majority of the reachability question. The reachability map tells you, per goal, *which
> modality* (LOF / GOF / neither) is required — i.e. **which screen the community should run
> next**, not just which gene.

This is a "right question" contribution that no prediction model can make, because none of them
represent the modality requirement geometrically. It reframes the method's honest ceiling
(√f_LOF) as a *prescription for data generation*. One atlas-wide LOF/GOF/neither figure, already
computable, carries it.

**Move 4 — Make the falsification loop the deliverable, not the outlook. `[design artifact now, wet-lab later]`**
Running one CRISPRa arm against a certificate's named genes would provide the first direct test
of an infeasibility prediction — the framework's most distinctive and least-validated claim.
**Package that as a pre-registered, falsifiable prediction now**, as a first-class artifact: the
named activation gene set, the exact contrast, the predicted signed effect, the null it must
clear, and the decision rule that would *refute* the oracle. A pre-registered prediction that can
be proven wrong is a stronger scientific object than any additional in-silico number — and it is
the physical embodiment of the AI-as-instrument-inside-a-falsification-loop thesis. It also turns
the certificate's current honest caveat (top genes are state-markers/negative regulators,
"activate" may be sign-wrong for some) into the *content* of a sharp, testable hypothesis rather
than a hedge.

**Move 5 — Make "portable vs local" a first-class verdict axis. `[build, medium]`**
The field's *stated* #1 open problem is generalization to novel conditions (*Nat. Genet.*
2025 causal-ML perspective; the AIVC review's transport critique). The transfer result already
shows direction ports but recipe does not (§4f). Promote this from a caveat to an **output
axis**: tag every verdict `transport-stable` (reachable across cell types/donors/conditions) vs
`context-specific` (reachable here only). A target that is *provably outside* in aged cells but
reachable in young ones is a headline, not a footnote. This directly answers the problem the
field says it cares most about, using the transfer machinery already built (07) and the
condition-resolved cones of §4e.

### 6.6 What is explicitly *not* the priority (so effort goes to the wave)

- **Not a higher held-out cosine.** Chasing the reachability number up is tide-work; the number
  is already ~71% of its modality-imposed ceiling. Accuracy is not the bottleneck.
- **Not a bigger model / a neural cone.** The linear, measured-effect composition is the part
  the field's own benchmarks vindicate. Replacing it with a learned model trades away the one
  defensible edge for the one thing that loses to baselines.
- **Not more datasets for their own sake.** A second atlas matters only insofar as it turns
  transfer (Move 5) from an existence proof into a generalization curve — i.e. in service of the
  field's stated open problem, not as a leaderboard entry.

*The discipline: every addition should either sharpen the inverse/feasibility/trust idea or
close the falsification loop. If it only moves a metric, it is a tide.*

---

## 7. One-line summary for the pitch

> Everyone else ranks genes or predicts what a perturbation will do. **We answer a
> yes/no engineering question with a proof: can you get this cell to that state by
> knockdowns, what's the smallest set, and — when you can't — what must you activate
> instead.** Measured effects, convex geometry, a certificate either way.
