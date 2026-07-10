# Related work — how cell-state reachability differs from prior cell-fate-control methods

*A citation-grounded literature review and methods comparison. Companion to `NOVELTY.md`
(what is new, why it matters to a drug-development team, and the field positioning). This file does
the middle job an expert reviewer asks for: it maps the neighbourhood of prior methods,
states the pain point precisely, and shows — against the actual literature — what our
approach adds and which limitations it addresses.*

*Every method discussed below was retrieved and verified against a primary database record
(OpenAlex / PubMed); DOIs, venues, years, and citation counts are listed in
[§8 References](#8-references) and in the machine-readable `references.csv`. The capability
classifications behind the two comparison figures are in `method_comparison_matrix.csv`.*

---

## 1. Why this matters — the problem the field is actually trying to solve

Cell-state engineering asks an inverse-design question: **given a starting cell state and a
desired target state, what set of interventions moves the cell there — and is it even
possible?** This is the question behind directed differentiation, transdifferentiation,
CAR-T and Treg engineering, and disease-state reversal. It is also, framed slightly
differently, the central question of pharmaceutical target selection: *which gene, knocked
down or activated, will move a diseased cell toward a healthy phenotype?*

Our answer, in one figure: compose the *measured* effect vectors from a genome-scale screen,
ask whether the target lies in their convex cone, and return a verdict with a checkable
certificate — then show it transfers to a second screen unchanged.

![Figure 1 — the reachability oracle end to end. (A) A genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells yields a measured effect vector per knockdown — the empirical "dictionary", with no inferred network. (B) Reachability becomes membership in a convex cone of non-negative combinations of those measured vectors; the Th2→Th1 target decomposes into an in-cone reachable component (teal) and an out-of-cone residual (amber), with "outside" certified by a separating hyperplane (Farkas/KKT). (C) The verdict for Th2→Th1: partly reachable (held-out cosine 0.448, null z ≈ 24, KKT residual 1.1×10⁻¹¹), a signed decomposition (LOF 39% / GOF 25% / neither 35%), and an activation certificate naming falsifiable CRISPRa hypotheses. (D) The same operator transfers unchanged to the Norman 2019 K562 CRISPRa screen (held-out CEBPA cosine 0.878).]({{artifact:724c5ece-a685-440a-8850-b3617cb1adcb}})

That question is expensive to get wrong. Only about 10% of clinical programmes reach
approval, and the cost of drug development is driven primarily by failure that is paid for
late (Minikel et al., *Nature* 2024). The failure concentrates at Phase II — the lowest
transition rate of any phase — and roughly half of Phase II/III failures are for lack of
efficacy, i.e. the wrong target or wrong mechanism rather than safety (Harrison, *Nat Rev
Drug Discov* 2016). The two interventions with published evidence of improving those odds
are **human genetic support** (≈2.6× the odds of approval; Minikel et al. 2024) and
**measured rather than inferred target effects** — the direct antidote to the 11–25%
preclinical reproducibility documented by the Amgen and Bayer audits (Begley & Ellis,
*Nature* 2012; Prinz et al., *Nat Rev Drug Discov* 2011). *(Full citations for these
pharma-attrition figures are in `NOVELTY.md` §References.)*

So the field needs a method that (a) works from measured effects in the relevant cell type,
(b) answers the *achievability* question rather than merely ranking candidates, and (c) says
so honestly when a goal is not achievable. The rest of this document shows that the large,
active body of prior work does (a) **or** approaches (b), but no prior method does both — and
none does (c).

---

## 2. The landscape at a glance

We surveyed **91 verified methods (2011–2026)** spanning three research communities —
control-theoretic network control, gene-regulatory-network (GRN) inference with in-silico
perturbation, and deep perturbation-response prediction / foundation models — plus the
combinatorial-screen, optimal-transport, experimental-design, and (for the 2025–2026 cohort)
virtual-cell and benchmark literature that borders them. Of these, 34 form the core survey
(2011–2024) and 57 are verified 2025–2026 additions (see [§9](#9-recent-work-20252026-the-wave-arrived-the-gap-did-not-close)).
Positioning them on two axes that matter for this problem makes the gap visible:

- **What question does it answer?** Forward prediction (*what will this perturbation do?*)
  versus inverse design / feasibility (*can this target state be reached, and how?*).
- **What is it grounded in?** An inferred or hand-built model of the regulatory network
  versus directly measured perturbation effects, with or without a checkable certificate.

![Method-landscape positioning map over 91 methods (2011–2026). Perturbation-response, foundation and virtual-cell models (amber/purple) are grounded in measured data but answer only the forward-prediction question; network-control methods (blue) reach toward feasibility but operate on inferred or hand-built models; GRN-inference/TF predictors (green) sit in the lower-left; benchmarks (brown) cluster at the measured-but-forward left edge. Larger dots are the 2011–2024 core; smaller faded dots are the 2025–2026 additions. The top-right region — a data-grounded feasibility verdict with a certificate — is empty across all 91 prior methods.]({{artifact:891d8a45-8d5c-44fe-91e3-abab7ff9e434}})

The two capabilities are **disjoint across the entire field**: the methods that use measured
effects answer only forward questions, and the methods that reach toward feasibility do so by
simulating a model. Our convex-cone reachability oracle is the only entry at their
intersection.

The same picture, resolved to individual capabilities and aggregated by family, is a heatmap
(Panel A) plus a per-method disjointness scatter (Panel B):

![Capability comparison over 91 prior methods (2011–2026). Panel A: the six method families (rows) against seven capabilities (columns); each cell is the fraction of that family's methods with the capability (● = all). The measured-effect / forward families (perturbation-response, foundation/virtual-cell, combinatorial, benchmarks) populate the left columns but are empty on feasibility verdict, certificate, and provable infeasibility; only the network-control family shows verdict/certificate — and only at ~35–38%, and only because those guarantees hold inside an assumed model (its measured-effects score is ~0.09). Panel B: every one of the 91 methods placed by grounding (y: measured effects) against feasibility (x: verdict/certificate); no prior method occupies the top-right, which holds only this work. The single genuinely all-empty column for every prior method is provable-infeasibility.]({{artifact:827b6ea2-f6c7-4ea5-b03f-efe938f49703}})

The headline counts, computed over the **91 prior methods** (our own row excluded):

| Capability | prior methods with it |
|---|---|
| Operates directly on **measured** perturbation effects | 14 full, 35 partial (of 91) |
| Quantitative **held-out** validation | 51 full, 30 partial |
| Handles **combinatorial / epistatic** perturbations | 9 full, 53 partial |
| Returns a **feasibility verdict** (reachable / not) | **0 full** (13 partial, all model-theoretic) |
| Emits a data-grounded **feasibility certificate** | **0 full** (14 partial, all model-theoretic) |
| Can return **provable infeasibility** + what is missing | **0 full of 91** (1 partial) |
| Integrates target **druggability / modality** | 1 full, 2 partial |

The precise reading of the two feasibility columns matters, because the heatmap does show
faint marks there. **No prior method — none of the 91 — scores a full *yes* on feasibility
verdict, certificate, or provable infeasibility.** The 13–14 *partial* marks on
verdict/certificate are all *model-theoretic*: their guarantee holds only inside a hand-built,
fitted, or inferred model, never over measured effects. Most are the network-control family
(§3, structural/Boolean control and control-from-data such as Ronquist [9]), but the set also
includes an inferred single-cell GRN — CEFCON [22], which imports feedback-vertex-set
driver-control onto a *reconstructed* network — and one causally-inspired perturbation
predictor (PAIRING [49]); the common thread is the assumed model, not a single family. Not one
of the 15 is grounded in measured effects: their measured-effects scores are "no" (13) or
"partial" (2), never "yes". And the disjunction is exact in the other direction: **all 14
methods that operate on measured effects score "no" on both verdict and certificate.** So the
two capabilities this work combines — measured grounding and a feasibility
verdict-with-certificate — are held by disjoint sets of prior methods, and the
provable-infeasibility column is empty (no full *yes*) for the entire field. That disjunction
is the gap.

---

## 3. Control-theoretic network control of cell fate

This is the field that asks *our exact question* — the minimal intervention set that drives a
system to a target state — and has the most mature theory for it.

**Structural and target controllability.** Liu, Slotine & Barabási [1] cast controllability
of a network as a maximum-matching problem and compute the minimum set of *driver nodes* that
render a linear time-invariant (LTI) system steerable in principle; Gao et al. [2] extend
this to controlling a chosen *target subset* of nodes. These are foundational network-science
results (Liu et al. is cited >3,000 times) and they *do* carry a controllability guarantee.
But the guarantee is about an **assumed linear model on the topology**, not about whether a
specific measured expression state is reachable from real perturbations; there is no
infeasibility certificate for a biological target, and the work is validated on network
topologies, not cell-fate experiments. Mochizuki & Fiedler [3] and Kim et al. [4] similarly
identify determining node sets (feedback vertex sets, "control kernels") from network
*structure*.

**Boolean attractor control.** Zañudo & Albert's stable-motif control [5] and the
feedback-vertex-set/source control of Zañudo, Yang & Albert [6] compute node interventions
that switch a Boolean network between attractors; Joo et al. [7] assess the relative
stability of Boolean cell states, and NETISCE [8] steers a signalling network toward a
desired steady state. These methods return exactly the kind of *intervention set* an engineer
wants — but every one requires a **hand-built or curated Boolean/logical model** as input, and
the control set is only as trustworthy as that model.

**Control models fit to data.** A newer strand fits a *dynamical* model to measured omics and
then applies control theory: Ronquist et al. [9] fit an LTI system to Hi-C / time-series data;
cSTAR (Rukhlenko et al. [10]) builds a state-transition vector plus a mechanistic network from
omics to identify perturbation targets; Wytock & Motter [11] transfer-learn functional
transcriptional networks for reprogramming design. cSTAR is the closest conceptual neighbour
to our work in *vocabulary* — it also speaks of "cell-state transitions" and which nodes to
perturb — but it does so by **fitting and then controlling a model**, and it returns target
nodes, not a reachable/not-reachable verdict with a certificate. Where methods in this strand
do carry a guarantee — Ronquist et al.'s LTI formulation admits a controllability argument, so
it scores a *partial* on our certificate axis — that guarantee is a property of the fitted
linear model, not a data-grounded infeasibility proof about the measured effects; the
distinction our verdict turns on is exactly this one.

**TF-reprogramming predictors.** Mogrify (Rackham et al. [12]) and CellNet (Cahan et al. [13])
predict the transcription-factor set for a cell-type conversion using an expression +
regulatory-network heuristic. They are influential (CellNet cited >500 times) and
experimentally motivated, but they output a *ranked TF set* from network influence, not a
feasibility decision, and they run on an inferred regulatory network.

> **The shared limitation.** Every method in this family computes its answer by *reasoning
> over a model* — a Boolean network, an LTI system, an inferred GRN, a network-influence
> heuristic. This is the field's structural bottleneck, stated openly in its own literature:
> the control set you get out is only as trustworthy as the network you put in, and network
> inference from observational data is underdetermined. CEFCON (Wang et al. [22]) even imports
> feedback-vertex-set driver-control ideas directly onto an *inferred* single-cell GRN,
> inheriting exactly this dependency.

**Our delta against this family.** We do not model the network. We compose the **measured**
genome-scale CRISPRi/CRISPRa effect vectors directly, so the "dictionary" of what each
perturbation does is empirical, not inferred. And because knockdown is loss-of-function and
combinations are non-negative, reachability becomes membership in a **convex cone** — which
means "unreachable" is a theorem (a separating hyperplane / Farkas certificate), not "we
didn't find a good set." No method above returns a provable infeasibility result grounded in
measured effects.

---

## 4. GRN inference + in-silico perturbation

The single-cell community's answer to the same question runs an inferred network *forward*.

CellOracle (Kamimoto et al. [14]) is the exemplar: it infers a GRN from scRNA-seq + ATAC, then
simulates a transcription-factor knockout or overexpression and propagates it into a "vector
field" of predicted cell-identity transitions. It is widely used (>600 citations) and its
in-silico perturbation is genuinely useful for hypothesis generation. scTenifoldKnk (Osorio et
al. [18]) performs a "virtual knockout" on a single-cell GRN; ANANSE (Xu et al. [17]) ranks
key TFs from enhancer networks for a fate transition; SCENIC and SCENIC+ (Aibar et al. [15];
Bravo González-Blas et al. [20]) and Dictys (Wang et al. [21]) infer regulons and dynamic
networks. BEELINE (Pratapa et al. [16]) is the benchmark that quantifies how hard GRN
inference is, and dyngen (Cannoodt et al. [19]) is the simulator used to generate
ground-truth for such benchmarks.

**What they share.** The perturbation effect is *predicted by propagating it through an
inferred network*. The output is a predicted response or a driver ranking, and — critically —
its reliability is bounded by GRN-inference accuracy, which BEELINE [16] shows is modest and
method-dependent. None returns a feasibility verdict, and none returns a certificate: an
in-silico knockout that fails to reach the target is not distinguishable from an inference
error.

**Our delta.** Our reachability verdict does not depend on reconstructing the network at all.
It asks a strictly weaker, answerable question — is the target in the non-negative span of the
*observed* effects? — and answers it with a numerically checkable optimality certificate
(KKT/Farkas residual ≈1e-11 in `reachability.py`). Where CellOracle says "here is the
predicted transition if my inferred network is right," we say "this target is reachable by
this measured combination — or provably is not, and here is what must be activated instead."

---

## 5. Deep perturbation-response prediction & foundation models

This family is *measured-grounded* — it learns from real Perturb-seq data — but it is
resolutely **forward**: it predicts what a perturbation will do, not whether a target is
reachable.

GEARS (Roohani et al. [26]) predicts the transcriptional outcome of novel *multi-gene*
perturbations using a gene-similarity graph, and is validated on held-out combinations —
genuinely strong on the combinatorial axis. CPA and chemCPA (Lotfollahi et al. [25]; Hetzel et
al. [24]) learn a compositional latent space of perturbations (chemCPA extends this to unseen
drug structures, the one method in our survey that touches druggability). scGen (Lotfollahi et
al. [23]) predicts responses by latent vector arithmetic; Biolord (Piran et al. [27])
disentangles known attributes for prediction. The single-cell foundation models — Geneformer
(Theodoris et al. [28]), scGPT (Cui et al. [29]), scFoundation (Hao et al. [30]) — are
pretrained on massive expression corpora and fine-tuned to, among other tasks, predict
perturbation responses or perform in-silico deletion.

**What they share.** Every one maps *perturbation → predicted response*. This is the inverse
of our question. A response predictor can tell you what GATA3 knockdown will do; it cannot tell
you whether *any* combination of available knockdowns reaches Th1, and it structurally cannot
tell you that a target is unreachable — a predictor always returns *some* prediction. None
emits a feasibility certificate; held-out validation here means held-out-perturbation response
accuracy, a different quantity from held-out-gene *reachability* generalisation.

**Adjacent: combinatorial data, optimal transport, experimental design.** Norman et al. [31]
map genetic-interaction manifolds from combinatorial Perturb-seq — the measured epistasis data
that a reachability model's additivity assumption should ultimately be tested against (and the
dataset our method transfers to; see §6). Waddington-OT (Schiebinger et al. [32]) and CellRank
(Lange et al. [33]) infer developmental *trajectories* and fate probabilities by optimal
transport / Markov chains on observed data — they describe how cells *do* move, not which
interventions *would* move them. IterPert (Huang et al. [34]) is the nearest thing to a
decision tool: it sequentially selects which perturbation experiment to run next. But its
output is an acquisition choice for screen design, not a reachability verdict for a target
state.

**Adjacent: the evaluation-methodology critique.** A parallel line of work asks not what these
models predict but whether the *metrics* used to judge them are sound. The linear-baseline
benchmark [83] and Systema [88] showed that deep perturbation models often fail to beat simple
baselines; Mejia et al. [92] trace much of that result to **signal dilution** — when a
perturbation moves only a few genes, an unweighted transcriptome-wide error is dominated by the
thousands of unchanged background genes and rewards a predictor that merely reproduces the
dataset mean. Their remedy is DEG-aware metrics (weighted MSE and weighted ΔR²) reported
against *both* a negative and a positive control (their "interpolated duplicate"), summarised as
a scale-free **dynamic-range fraction**. That critique targets *forward-prediction* evaluation,
but its core concern — that a transcriptome-wide similarity can be inflated by the quiet
background — applies directly to our reachability *cosine*, which is likewise scored over all
genes. We therefore adopted its machinery to stress-test our own verdict (RESULTS §6.5 and
`notebooks/08_deg_weighted_evaluation.ipynb`): under DEG weighting the Th2→Th1 conclusion is
preserved and in fact strengthens (reachable cosine 0.63 → 0.80, held-out z 28 → 14, both far
above the significance floor), confirming it is not a dilution artifact. The unweighted default
reproduces every previously published number bit-for-bit, so the addition is non-breaking.

**Our delta.** We turn the measured effects these models predict *into an inverse-design
oracle*. Prediction is necessary but not sufficient: knowing every perturbation's response
does not tell you whether a target is in their reachable set. Convex geometry closes that gap —
and, uniquely, closes it in a way that can return "no, and here is the missing direction."

---

## 6. The pain point, and how our methodology addresses it

**The pain point, stated once.** Existing methods either (i) rank or predict on *measured*
data but never certify achievability, or (ii) reason about achievability but only inside an
*inferred or hand-built model* whose error they cannot bound. A target-selection team is
therefore left with optimistic rankings and model-dependent simulations — and, crucially, with
no principled way to be told *"stop: this goal is not reachable by this modality."* Since late
efficacy failure is the dominant cost of drug development (§1), the missing capability — a
credible, early, measured-grounded *stop* — is exactly the expensive one.

**How the convex-cone reachability oracle addresses it.** Four design choices map onto the
four gaps:

1. **Measured, not modelled (addresses the §3–§4 model-dependence).** The cone's generators
   are the empirically observed genome-scale CRISPRi effect vectors in the relevant primary
   human cell type. There is no network to infer and no inference error to propagate.
2. **A verdict with a certificate (addresses the §3–§5 "ranking/prediction only").**
   Non-negative least squares over the effect cone is a convex program; its KKT/Farkas
   conditions give a numerically checkable optimality certificate. "Reachable" comes with the
   mixing weights (the recipe); "unreachable" comes with a separating hyperplane and the
   specific genes in the residual direction that would need to be *activated* — a constructive
   output no surveyed method produces.
3. **Held-out-gene validation + a target-specific null (addresses inflated in-sample fit).**
   We fit weights on a random half of the target signature's genes and score alignment on the
   held-out half, and we compare every verdict against a shuffled-target null. On Th2→Th1 in
   resting CD4⁺ T cells the held-out cosine is 0.448 at a null z of 24.4; the target decomposes
   into 39% reachable by knockdown / 31% activation-required / 30% orthogonal. The method
   transfers *unchanged* to the Norman K562 CRISPRa screen [31], where held-out CEBPA is
   partially reachable (cosine 0.878, null z ≈37, held-out-gene z ≈23.5) and measured
   double-perturbations let us confirm additivity is a bounded approximation (median cosine
   0.71, sum-of-singles vs. measured) — a check the singles-only T-cell screen could not do.
4. **A modality-resolved, druggability-cross-referenced triage (addresses the §5 single
   druggability entry).** Because the cone is signed (LOF vs. GOF), every required node carries
   a modality; crossing that with Open Targets tractability turns the verdict into a triage.
   The decision-relevant finding: 44% of the required knockdown nodes across the atlas are
   hard-to-drug, and the strongest genetically-supported node among them (IRF1, 17 immune-
   disease genetic associations) has no conventional drug handle at all.

The same four design choices, read as *decisions a drug-development team makes* rather than
method properties, are what turn the oracle into a triage instrument:

![How the oracle is useful to a drug-development team. (A) The decision it makes — a desired cell-state change is routed to GO (reachable → put the minimal ranked knockdown set in the screen), STOP (provably outside the cone → don't spend a CRISPRi arm), or REDIRECT (activation-required → switch to CRISPRa; the certificate names the genes). The GO branch feeds a modality triage over 102 real nominations (44% hard-to-drug; 10 with a clinical-grade drug), separating green-light (JAK2, ICOS, MAPK14, CD3D), tractable-but-untried (IL7R, ZAP70, TET2), and required-but-undruggable (IRF1) tiers. (B) Four value levers mapped onto the drug-development attrition funnel (≈10% reach approval; Phase II is the lowest transition at 30.7%; ~half of Phase II/III failures are lack of efficacy). Nominations are wet-lab hypotheses, not validated targets; this is a prioritisation-and-triage instrument, not a claimed hit rate.]({{artifact:3a67feb4-18d0-4fbd-b88d-69317752dfdd}})

### The two results that show the machinery working

The reachability *spectrum* is the honesty figure — the observed alignment rises far above a
shuffled-target null band at every recipe size, then plateaus well below 1, which is the
correct behaviour for a partially-reachable target:

![Reachability spectrum for Th2→Th1 by knockdown. The observed alignment (blue) rises to ~0.46 and plateaus, with a knee at k=7; the shuffled-target null band (grey) sits far below (~0.05–0.15) at every recipe size. The gap between the two is the signal; the plateau below 1 is the honest ceiling of loss-of-function alone.]({{artifact:5cd1b5a4-4e3a-4e70-b23f-d475f05e837a}})

And the decomposition + **activation certificate** is the output no prior method produces —
the target splits into what knockdown can reach and what it provably cannot, and the
unreachable direction is resolved into named genes that a CRISPRa arm would have to test:

![Left: the Th2→Th1 target splits into 39% reachable by knockdown, 31% activation-required, 30% orthogonal residual. Right: the activation certificate names the genes carrying the unmet upward demand (LYAR, IKZF3, CRTAM, LAG3, …) — concrete, wet-lab-testable CRISPRa hypotheses, with immune-function genes marked.]({{artifact:a63d1f60-b868-4f2c-a3d0-0668aae4d0e3}})

---

## 7. Positioning summary — the one-paragraph delta

Prior cell-fate-control methods either simulate an **inferred or hand-built model** to reason
about achievability (structural/Boolean control [1–8], control-from-data [9–11], TF predictors
[12–13], GRN in-silico perturbation [14–22]) or **predict responses** from measured data
without ever asking achievability (perturbation-ML [23–27], foundation models [28–30],
combinatorial/OT/design [31–34]). We compose the **measured** effect vectors directly and, via
the convex geometry that loss-of-function non-negativity forces, return a **reachable /
provably-outside-the-cone verdict with a constructive certificate**, validated by held-out-gene
generalisation and resolved to intervention modality. We claim novelty for this
measured-effect, certificate-carrying reachability *decision* — not for the Th1/Th2 or CEBPA
regulator lists, which the source screens already report and which we use only as positive
controls.

*See `NOVELTY.md` §2 for the full novelty argument and its impact section for the translational
use cases.*

---

## 9. Recent work (2025–2026) — the wave arrived, the gap did not close

The 2024 survey above could have aged badly: 2025–2026 has been the most active two-year
period this field has seen, dominated by the "virtual cell" programme and a flood of large
perturbation models. We therefore re-ran the search over 2025–2026 only and verified **57
additional methods** against primary records (real DOIs/PMIDs; machine-readable in
`references.csv`, cohort `2025-2026`). The headline is that the wave populated the *same*
regions of the landscape and left the decisive corner empty. Across all 91 methods, the count
of prior methods that return a feasibility verdict, emit a data-grounded certificate, or prove
infeasibility is still **zero**. The prose below names representative exemplars per family; the
complete set of 57 (entries [35]–[91]) is in [§8](#8-references) and `references.csv`.

**The virtual-cell / large-model surge (forward predictors, at scale).** Arc Institute's
**State** [47] is an explicit perturbation model trained across cell contexts; **TranscriptFormer**
[64], **scPRINT** [63], **CellFM** [69], **Nicheformer** [67] and the Cell2Sentence-scale LLM
approach **C2S-Scale** [66] extend the foundation-model line; **GET** [62] predicts expression
from regulatory sequence; **TxPert** [60] targets transcriptional perturbation response. The
community also formalised the goal itself — the **Virtual Cell Challenge** [74] and the
"AI virtual cell" roadmap pieces [73] — and released the data substrate for it, notably the
**Tahoe-100M** perturbation atlas [91] and the **moscot** scaling of optimal transport to
atlas-scale [81]. Every one of these is, by construction, a *forward* model or a dataset: it
predicts a response, embeds a cell, or provides training data. None returns a reachable /
provably-outside verdict, a minimal knockdown recipe, or a Farkas/KKT certificate — the gap our
method fills is orthogonal to what scaling a predictor buys.

**The critical-benchmark literature (which strengthens our design choice).** The most
decision-relevant 2025 development is a cluster of careful benchmarks reporting that current
perturbation-prediction models — including the foundation models — often **fail to beat simple
mean or linear baselines** on genuinely unseen perturbations: Ahlmann-Eltze et al. [83],
Wong et al. [54], Csendes et al. [84], Kedzierska et al. [85], Wu et al. [90], Systema
(Viñas et al. [88]), Kernfeld et al. [86], Wei et al. [89], and the Open Problems effort
[87]. This is not a criticism we need to make from the outside; the field is making it of
itself. It also sharpens our positioning: if forward prediction on unseen perturbations is
this fragile, then a method whose output is a *feasibility verdict computed directly from the
measured effects* — one that does not extrapolate to unseen perturbations at all, and that
comes with a checkable certificate and a shuffled-target null — is attacking exactly the
reliability gap these benchmarks expose. Our claim is deliberately weaker and therefore more
defensible than "predict the response of an unseen perturbation."

**Network control did not stand still either.** New Boolean/structural-control work — biobalm
succession diagrams (Trinh et al. [38]), modular Boolean control (Murrugarra et al. [37]),
trap-space control (Cifuentes-Fontanals et al. [35]), stable-motif control in multicellular
models (Metzcar et al. [36]), and biomolecular feedback control (Vecchio et al. [39]) —
continues to produce genuine control strategies, and these are the 2025–2026 methods that
score *partial* on our verdict/certificate axes. As with their pre-2024 predecessors, the
guarantee is a theorem about a hand-built or inferred model, not a certificate over measured
effects; they sit in the same lower/right model-theoretic region of the landscape, not in the
measured-and-certified corner.

**Net effect on the positioning.** Adding 57 verified 2025–2026 methods did not move a single
entry into the top-right quadrant. It did the opposite: it enlarged the measured-but-forward
cloud (more predictors), enlarged the model-theoretic-feasibility cloud (more Boolean control),
and — through the benchmark papers — supplied independent evidence that the forward-prediction
strategy those clouds represent is unreliable on unseen perturbations. The measured-effect,
certificate-carrying reachability verdict remains, as of 2026, unoccupied by any other method
we could verify.

---

## 8. References

*All records verified against OpenAlex / PubMed; DOIs, PMIDs, venues, and citation counts as
retrieved. Machine-readable copy in `references.csv`. Citation counts are a snapshot and will
drift.*

**Structural & target network controllability**

1. **Structural controllability (minimum driver-node set)** — Liu, Slotine & Barabási *Controllability of complex networks.* Nature (2011). doi:10.1038/nature10011. PMID:21562557. (cited by 3,298)
2. **Target control of complex networks** — Gao, Liu, Slotine & Barabási *Target control of complex networks.* Nature Communications (2014). doi:10.1038/ncomms6415. PMID:25388503. (cited by 422)
3. **Feedback vertex set / determining nodes** — Mochizuki, Fiedler et al. *Dynamics and control at feedback vertex sets. II: A faithful monitor to determine the diversity of molecular activities in regulatory networks.* Journal of Theoretical Biology (2013). doi:10.1016/j.jtbi.2013.06.009. PMID:23774067. (cited by 193)
4. **Control kernel of biomolecular regulatory networks** — Kim, Park & Cho *Discovery of a kernel for controlling biomolecular regulatory networks.* Scientific Reports (2013). doi:10.1038/srep02223. PMID:23860463. (cited by 117)

**Boolean / logical attractor control**

5. **Stable-motif control of Boolean networks** — Zañudo & Albert *Cell Fate Reprogramming by Control of Intracellular Network Dynamics.* PLoS Computational Biology (2015). doi:10.1371/journal.pcbi.1004193. PMID:25849586. (cited by 248)
6. **Structure-based (feedback-vertex-set + source) control** — Zañudo, Yang & Albert *Structure-based control of complex networks with nonlinear dynamics.* Proceedings of the National Academy of Sciences (2017). doi:10.1073/pnas.1617387114. PMID:28655847. (cited by 282)
7. **Boolean relative dynamic stability of cell states** — Joo, Zhou, Huang & Cho *Determining Relative Dynamic Stability of Cell States Using Boolean Network Model.* Scientific Reports (2018). doi:10.1038/s41598-018-30544-0. PMID:30104572. (cited by 56)
8. **NETISCE (attractor/steady-state steering)** — Marazzi, Shah et al. *NETISCE: a network-based tool for cell fate reprogramming.* npj Systems Biology and Applications (2022). doi:10.1038/s41540-022-00231-y. PMID:35725577. (cited by 17)

**Control-theoretic reprogramming fit to data**

9. **Control-theoretic cellular reprogramming (LTI from Hi-C/time-series)** — Ronquist, Patterson, Rajapakse et al. *Algorithm for cellular reprogramming.* Proceedings of the National Academy of Sciences (2017). doi:10.1073/pnas.1712350114. PMID:29078370. (cited by 41)
10. **cSTAR (cell-State Transition Assessment & Regulation)** — Rukhlenko, Halász, Kholodenko et al. *Control of cell state transitions.* Nature (2022). doi:10.1038/s41586-022-05194-y. PMID:36104561. (cited by 89)
11. **Transfer-learned functional-network control** — Wytock & Motter *Cell reprogramming design by transfer learning of functional transcriptional networks.* Proceedings of the National Academy of Sciences (2024). doi:10.1073/pnas.2312942121. PMID:38437548. (cited by 8)

**Network-influence TF-reprogramming predictors**

12. **Mogrify (network-influence TF predictor)** — Rackham et al. (Mogrify) *A predictive computational framework for direct reprogramming between human cell types.* Nature Genetics (2016). doi:10.1038/ng.3487. PMID:26780608. (cited by 306)
13. **CellNet (GRN cell-identity assessment)** — Cahan et al. (CellNet) *CellNet: Network Biology Applied to Stem Cell Engineering.* Cell (2014). doi:10.1016/j.cell.2014.07.020. PMID:25126793. (cited by 545)

**GRN inference + in-silico perturbation**

14. **CellOracle** — Kamimoto et al. *Dissecting cell identity via network inference and in silico gene perturbation.* Nature (2023). doi:10.1038/s41586-022-05688-9. PMID:36755098. (cited by 656)
15. **SCENIC** — Aibar et al. *SCENIC: single-cell regulatory network inference and clustering.* Nature Methods (2017). doi:10.1038/nmeth.4463. PMID:28991892. (cited by 6,993)
16. **BEELINE** — Pratapa et al. *Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data.* Nature Methods (2020). doi:10.1038/s41592-019-0690-6. PMID:31907445. (cited by 825)
17. **ANANSE** — Xu et al. *ANANSE: an enhancer network-based computational approach for predicting key transcription factors in cell fate determination.* Nucleic Acids Research (2021). doi:10.1093/nar/gkab598. PMID:34244796. (cited by 98)
18. **scTenifoldKnk** — Osorio et al. *scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation.* Patterns (2022). doi:10.1016/j.patter.2022.100434. PMID:35510185. (cited by 157)
19. **dyngen** — Cannoodt et al. *Spearheading future omics analyses using dyngen, a multi-modal simulator of single cells.* Nature Communications (2021). doi:10.1038/s41467-021-24152-2. PMID:34168133. (cited by 129)
20. **SCENIC+** — Bravo Gonzalez-Blas et al. *SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks.* Nature Methods (2023). doi:10.1038/s41592-023-01938-4. PMID:37443338. (cited by 695)
21. **Dictys** — Wang et al. *Dictys: dynamic gene regulatory network dissects developmental continuum with single-cell multiomics.* Nature Methods (2023). doi:10.1038/s41592-023-01971-3. PMID:37537351. (cited by 107)
22. **CEFCON** — Wang et al. *Deciphering driver regulators of cell fate decisions from single-cell transcriptomics data with CEFCON.* Nature Communications (2023). doi:10.1038/s41467-023-44103-3. PMID:38123534. (cited by 25)

**Deep perturbation-response prediction**

23. **scGen** — Lotfollahi et al. *scGen predicts single-cell perturbation responses.* Nature Methods (2019). doi:10.1038/s41592-019-0494-8. PMID:31363220. (cited by 696)
24. **chemCPA** — Hetzel et al. *Predicting Cellular Responses to Novel Drug Perturbations at a Single-Cell Resolution.* NeurIPS 2022 (arXiv:2204.13545) (2022). doi:10.48550/arXiv.2204.13545. (cited by 42)
25. **CPA (Compositional Perturbation Autoencoder)** — Lotfollahi et al. *Predicting cellular responses to complex perturbations in high-throughput screens.* Molecular Systems Biology (2023). doi:10.15252/msb.202211517. PMID:37154091. (cited by 272)
26. **GEARS** — Roohani et al. *Predicting transcriptional outcomes of novel multigene perturbations with GEARS.* Nature Biotechnology (2023). doi:10.1038/s41587-023-01905-6. PMID:37592036. (cited by 309)
27. **Biolord** — Piran et al. *Disentanglement of single-cell data with biolord.* Nature Biotechnology (2024). doi:10.1038/s41587-023-02079-x. PMID:38225466. (cited by 63)

**Single-cell foundation models**

28. **Geneformer** — Theodoris et al. *Transfer learning enables predictions in network biology.* Nature (2023). doi:10.1038/s41586-023-06139-9. PMID:37258680. (cited by 1,051)
29. **scGPT** — Cui et al. *scGPT: toward building a foundation model for single-cell multi-omics using generative AI.* Nature Methods (2024). doi:10.1038/s41592-024-02201-0. PMID:38409223. (cited by 1,073)
30. **scFoundation** — Hao et al. *Large-scale foundation model on single-cell transcriptomics.* Nature Methods (2024). doi:10.1038/s41592-024-02305-7. PMID:38844628. (cited by 473)

**Combinatorial data, optimal transport & experimental design**

31. **Genetic-interaction manifolds (Norman combinatorial Perturb-seq)** — Norman et al. *Exploring genetic interaction manifolds constructed from rich single-cell phenotypes.* Science (2019). doi:10.1126/science.aax4438. PMID:31395745. (cited by 376)
32. **Waddington-OT** — Schiebinger et al. *Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming.* Cell (2019). doi:10.1016/j.cell.2019.01.006. PMID:30712874. (cited by 855)
33. **CellRank** — Lange et al. *CellRank for directed single-cell fate mapping.* Nature Methods (2022). doi:10.1038/s41592-021-01346-6. PMID:35027767. (cited by 783)
34. **IterPert (sequential optimal experimental design of perturbation screens)** — Huang et al. *Sequential Optimal Experimental Design of Perturbation Screens Guided by Multi-modal Priors.* RECOMB 2024 (Lecture Notes in Computer Science) (2024). doi:10.1007/978-1-0716-3989-4_2. (cited by 7)


### Recent additions (2025–2026)

*57 additional methods, verified against primary records this round (real DOIs/PMIDs). Citation counts are a snapshot; brand-new preprints may show none.*


**Network control (2025–2026)**

35. **Trap-space-based control strategy identification (ASP)** — Laura Cifuentes-Fontanals, Elisa Tonello et al. *Node and edge control strategy identification via trap spaces in Boolean networks.* BMC Bioinformatics (2025). doi:10.1186/s12859-025-06135-y. PMID:41057758. (cited by 1)
36. **Boolean network stable-motif/exhaustive control + multiscale ABM translation** — John Metzcar, Katie Pletz et al. *Translating and evaluating single-cell Boolean network interventions in the multiscale setting.* arXiv (Cornell University) (2025). doi:10.48550/arxiv.2501.16052. (new/preprint)
37. **Modular Boolean Network Control** — David Murrugarra, Alan Veliz‐Cuba et al. *Modular Control of Boolean Network Models.* Bulletin of Mathematical Biology (2025). doi:10.1007/s11538-025-01471-9. PMID:40461704. (cited by 1)
38. **biobalm (succession diagram mapper)** — Van-Giang Trinh, Kyu Hyong Park et al. *Mapping the attractor landscape of Boolean networks with biobalm.* Bioinformatics (2025). doi:10.1093/bioinformatics/btaf280. PMID:40327535. (cited by 5)
39. **Biomolecular feedback/feedforward control design** — Domitilla Del Vecchio *Control systems for synthetic biology and a case-study in cell fate reprogramming.* Open MIND (2026). doi:10.48550/arxiv.2601.20135. (new/preprint)
40. **Computational blueprints for cell fate programming (review/framework)** — Pengyi Yang *Computational blueprints for cell fate programming.* Stem Cell Reports (2026). doi:10.1016/j.stemcr.2026.102929. PMID:42208532. (new/preprint)

**GRN inference + in-silico perturbation (2025–2026)**

41. **Airqtl** — Matthew W. Funk, Yuhe Wang et al. *Airqtl dissects cell state-specific causal gene regulatory networks with efficient single-cell eQTL mapping.* bioRxiv (Cold Spring Harbor Laboratory) (2025). doi:10.1101/2025.01.15.633041. (cited by 1)
42. **SMOGT (heterogeneous graph transformer for HRNet)** — Yuhong Huang, Chao Liu et al. *Deciphering hierarchical regulatory network of cell fate via an epigenetics-informed heterogeneous graph transformer on single-cell multi-omics data.* Briefings in Bioinformatics (2025). doi:10.1093/bib/bbaf664. PMID:41643202. (new/preprint)
43. **Regulator-Program-Trait causal graph** — Mineto Ota, Jeffrey P. Spence et al. *Causal modelling of gene effects from regulators to programs to traits.* Nature (2025). doi:10.1038/s41586-025-09866-3. PMID:41372418. (cited by 10)
44. **IGNITE (kinetic Ising GRN inference)** — Clelia Corridori, Merrit Romeike et al. *Unveiling gene perturbation effects through gene regulatory networks inference from single-cell transcriptomic data..* PLoS computational biology (2026). doi:10.1371/journal.pcbi.1014067. PMID:41984780. (new/preprint)
45. **CellPolaris** — Guihai Feng, Xin Qin et al. *CellPolaris: Transfer Learning for Gene Regulatory Network Construction to Guide Cell State Transitions.* Advanced Science (2026). doi:10.1002/advs.202508697. PMID:41498638. (cited by 1)
46. **scMagnifier** — Zhenhui He, Kangning Dong *scMagnifier: Resolving fine-grained cell subtypes via GRN-informed perturbations and consensus clustering.* PLoS Computational Biology (2026). doi:10.1371/journal.pcbi.1014167. PMID:42313807. (new/preprint)

**Perturbation-response prediction (2025–2026)**

47. **State (transformer perturbation model)** — Abhinav Adduri, Dhruv Gautam et al. *Predicting cellular responses to perturbation across diverse contexts with State.* bioRxiv (Cold Spring Harbor Laboratory) (2025). doi:10.1101/2025.06.26.661135. (cited by 63)
48. **PertSpectra** — Seowon Chang, Anna Shcherbina et al. *PertSpectra: Interpretable Matrix Factorization for Predicting Functional Impact of Genetic Perturbation Experiments.* Pacific Symposium on Biocomputing (2025). doi:10.1142/9789819824755_0033. PMID:41758161. (new/preprint)
49. **PAIRING** — Young-Hyun Han, Hyunjin Kim et al. *Identifying an optimal perturbation to induce a desired cell state by generative deep learning.* Cell Systems (2025). doi:10.1016/j.cels.2025.101405. PMID:40997798. (new/preprint)
50. **Squidiff (diffusion)** — Siyu He, Yuefei Zhu et al. *Squidiff: predicting cellular development and responses to perturbations using a diffusion model.* Nature Methods (2025). doi:10.1038/s41592-025-02877-y. PMID:41184550. (cited by 21)
51. **CellFlow (flow matching)** — Dominik Klein, Jonas Simon Fleck et al. *CellFlow enables generative single-cell phenotype modeling with flow matching.* bioRxiv (Cold Spring Harbor Laboratory) [preprint] (2025). doi:10.1101/2025.04.11.648220. (cited by 34)
52. **SCCVAE** — Emily Liu, Jiaqi Zhang et al. *Learning Genetic Perturbation Effects with Variational Causal Inference.* bioRxiv (Cold Spring Harbor Laboratory) (2025). doi:10.1101/2025.06.05.657988. PMID:40501829. (cited by 1)
53. **Large Perturbation Model (LPM)** — Djordje Miladinovic, Tobias Höppe et al. *In silico biological discovery with large perturbation models.* Nature Computational Science (2025). doi:10.1038/s43588-025-00870-1. PMID:41094040. (cited by 5)
54. **Baseline mean/linear control for perturbation prediction** — Daniel R. Wong, Abby S. Hill et al. *Simple controls exceed best deep learning algorithms and reveal foundation model effectiveness for predicting genetic perturbations.* Bioinformatics (2025). doi:10.1093/bioinformatics/btaf317. PMID:40407144. (cited by 12)
55. **GPerturb** — Hanwen Xing, Christopher Yau *GPerturb: Gaussian process modelling of single-cell perturbation data.* Nature Communications (2025). doi:10.1038/s41467-025-61165-7. PMID:40593897. (cited by 5)
56. **PerturbNet** — Hengshi Yu, Weizhou Qian et al. *PerturbNet predicts single-cell responses to unseen chemical and genetic perturbations.* Molecular Systems Biology (2025). doi:10.1038/s44320-025-00131-3. PMID:40640612. (cited by 18)
57. **PertAdapt** — Ding Bai, Le Song et al. *PertAdapt: unlocking single-cell foundation models for genetic perturbation prediction via condition-sensitive adaptation..* Bioinformatics (Oxford, England) (2026). doi:10.1093/bioinformatics/btag307. PMID:42412811. (new/preprint)
58. **PRIM (Prior-guided Response Inference Model)** — Xiuhao Fu, Chao Yang et al. *Incorporating valuable prior knowledge to improve deep learning prediction of genetic perturbation responses.* Genome Research (2026). doi:10.1101/gr.281523.125. PMID:41887798. (new/preprint)
59. **TF Perturb-seq rejuvenation screen** — Janine Sengstack, Jiashun Zheng et al. *Systematic identification of single transcription factor perturbations that drive cellular and tissue rejuvenation.* Proceedings of the National Academy of Sciences (2026). doi:10.1073/pnas.2515183123. PMID:41512022. (cited by 9)
60. **TxPert** — Frederik Wenkel, Wilson Tu et al. *TxPert: using multiple knowledge graphs for prediction of transcriptomic perturbation effects.* Nature Biotechnology (2026). doi:10.1038/s41587-026-03113-4. PMID:42067667. (cited by 1)

**Single-cell foundation models (2025–2026)**

61. **Foundation Model Benchmarking Study (scGPT/scGPT-like, Geneformer, Seurat v5)** — Srijan Atti, Shankar Subramaniam *Fundamental Limitations of Foundation Models in Single-Cell Transcriptomics.* bioRxiv (Cold Spring Harbor Laboratory) (2025). doi:10.1101/2025.06.26.661767. (cited by 5)
62. **GET (general expression transformer)** — Xi Fu, Shentong Mo et al. *A foundation model of transcription across human cell types.* Nature (2025). doi:10.1038/s41586-024-08391-z. PMID:39779852. (cited by 98)
63. **scPRINT** — Jérémie Kalfon, Jules Samaran et al. *scPRINT: pre-training on 50 million cells allows robust gene network predictions.* Nature Communications (2025). doi:10.1038/s41467-025-58699-1. PMID:40240364. (cited by 36)
64. **TranscriptFormer** — James D Pearce, Sara E. Simmonds et al. *A Cross-Species Generative Cell Atlas Across 1.5 Billion Years of Evolution: The TranscriptFormer Single-cell Model.* bioRxiv (Cold Spring Harbor Laboratory) [preprint] (2025). doi:10.1101/2025.04.25.650731. (cited by 27)
65. **GRNFormer** — Mufan Qiu, Xi-Jiang Hu et al. *GRNFormer: A Biologically-Guided Framework for Integrating Gene Regulatory Networks into RNA Foundation Models.* arXiv (Cornell University) (2025). doi:10.48550/arxiv.2503.01682. (new/preprint)
66. **C2S-Scale (Cell2Sentence, scaling LLMs)** — Syed Asad Rizvi, Daniel Lévine et al. *Scaling Large Language Models for Next-Generation Single-Cell Analysis.* bioRxiv (Cold Spring Harbor Laboratory) [preprint] (2025). doi:10.1101/2025.04.14.648850. PMID:41279114. (cited by 32)
67. **Nicheformer** — Alejandro Tejada-Lapuerta, Anna C. Schaar et al. *Nicheformer: a foundation model for single-cell and spatial omics.* Nature Methods (2025). doi:10.1038/s41592-025-02814-z. PMID:41168487. (cited by 51)
68. **scGPT-spatial** — Chloe Xueqi Wang, Haotian Cui et al. *scGPT-spatial: Continual Pretraining of Single-Cell Foundation Model for Spatial Transcriptomics.* bioRxiv (Cold Spring Harbor Laboratory) [preprint] (2025). doi:10.1101/2025.02.05.636714. (cited by 37)
69. **CellFM** — Yuansong Zeng, Jiancong Xie et al. *CellFM: a large-scale foundation model pre-trained on transcriptomics of 100 million human cells.* Nature Communications (2025). doi:10.1038/s41467-025-59926-5. PMID:40393991. (cited by 59)
70. **scLong** — Ding Bai, Shentong Mo et al. *scLong: a billion-parameter foundation model for capturing long-range gene context in single-cell transcriptomics.* Nature Communications (2026). doi:10.1038/s41467-026-69102-y. PMID:41639087. (cited by 5)
71. **CellxPert** — Andaç Demir, Erik W. Anderson et al. *CellxPert: Inference-Time MCMC Steering of a Multi-Omics Single-Cell Foundation Model for In-Silico Perturbation.* arXiv (Cornell University) (2026). doi:10.48550/arxiv.2605.00930. (new/preprint)

**Virtual-cell models & data (2025–2026)**

72. **ARTEMIS** — Sayali Alatkar, Daifeng Wang *ARTEMIS integrates autoencoders and Schrödinger Bridges to predict continuous dynamics of gene expression, cell population, and perturbation from time-series single-cell data.* Bioinformatics (2025). doi:10.1093/bioinformatics/btaf218. PMID:40662824. (cited by 3)
73. **Grow AI virtual cells (position/roadmap)** — Liujia Qian, Zhen Dong et al. *Grow AI virtual cells: three data pillars and closed-loop learning.* Cell Research (2025). doi:10.1038/s41422-025-01101-y. PMID:40128605. (cited by 40)
74. **Virtual Cell Challenge** — Yusuf Roohani, Tony J. Hua et al. *Virtual Cell Challenge: Toward a Turing test for the virtual cell.* Cell (2025). doi:10.1016/j.cell.2025.06.008. PMID:40578317. (cited by 59)
75. **SCALE (set-aware flow, conditional transport)** — Shuizhou Chen, Lang Yu et al. *SCALE:Scalable Conditional Atlas-Level Endpoint transport for virtual cell perturbation prediction.* arXiv (Cornell University) (2026). doi:10.48550/arxiv.2603.17380. (new/preprint)
76. **VCBench (VC benchmarking framework)** — Xinjie Mao, Songming Zhang et al. *Benchmarking virtual cell models for in-the-wild perturbation response.* arXiv (Cornell University) (2026). doi:10.48550/arxiv.2604.27646. (cited by 1)
77. **SBB Principles (Signal, Bounds, Baselines)** — Michael Vollenweider, Peter Bühlmann *Signal, Bounds, and Baselines: Principles for Evaluating Virtual Cell Perturbation Models.* bioRxiv (Cold Spring Harbor Laboratory) (2026). doi:10.64898/2026.04.20.719650. (cited by 1)
78. **AROMA** — Zhenyu Wang, Geyan Ye et al. *AROMA: Augmented Reasoning Over a Multimodal Architecture for Virtual Cell Genetic Perturbation Modeling.* arXiv (Cornell University) (2026). doi:10.48550/arxiv.2604.20263. (new/preprint)
79. **Lingshu-Cell (masked discrete diffusion world model)** — Han Zhang, Guo-Hua Yuan et al. *Lingshu-Cell: A generative cellular world model for transcriptome modeling toward virtual cells.* arXiv (Cornell University) (2026). doi:10.48550/arxiv.2603.25240. (new/preprint)

**Combinatorial / optimal transport / design (2025–2026)**

80. **Causally-inspired combinatorial perturbation prediction** — Guadalupe Gonzalez, Xiang Lin et al. *Combinatorial prediction of therapeutic perturbations using causally inspired neural networks.* Nature Biomedical Engineering (2025). doi:10.1038/s41551-025-01481-x. PMID:40925962. (cited by 18)
81. **moscot (multi-omics single-cell optimal transport)** — Dominik Klein, Giovanni Palla et al. *Mapping cells through time and space with moscot.* Nature (2025). doi:10.1038/s41586-024-08453-2. PMID:39843746. (cited by 127)
82. **CAR-OT (conditional optimal transport for CAR response)** — Alice Driessen, Jannis Born et al. *Modeling chimeric antigen receptor response at the single-cell level with conditional optimal transport.* Cell Systems (2026). doi:10.1016/j.cels.2026.101591. PMID:42025164. (cited by 1)

**Benchmarks & critical evaluations (2025–2026)**

83. **Linear-baseline benchmark of perturbation models** — Constantin Ahlmann-Eltze, Wolfgang Huber et al. *Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines.* Nature Methods (2025). doi:10.1038/s41592-025-02772-6. PMID:40759747. (cited by 96)
84. **Benchmark of foundation models for post-perturbation RNA-seq** — Gerold Csendes, Gema Sanz et al. *Benchmarking foundation cell models for post-perturbation RNA-seq prediction.* BMC Genomics (2025). doi:10.1186/s12864-025-11600-2. PMID:40269681. (cited by 34)
85. **Zero-shot evaluation of single-cell foundation models** — Kasia Z. Kedzierska, Lorin Crawford et al. *Zero-shot evaluation reveals limitations of single-cell foundation models.* Genome biology (2025). doi:10.1186/s13059-025-03574-x. PMID:40251685. (cited by 65)
86. **Comparison of expression-forecasting (perturbation) methods** — Eric Kernfeld, Yunxiao Yang et al. *A comparison of computational methods for expression forecasting.* Genome biology (2025). doi:10.1186/s13059-025-03840-y. PMID:41250104. (cited by 6)
87. **Open Problems in single-cell analysis (benchmark platform)** — Malte D. Luecken, Scott Gigante et al. *Defining and benchmarking open problems in single-cell analysis.* Nature Biotechnology (2025). doi:10.1038/s41587-025-02694-w. PMID:40595413. (cited by 31)
88. **Systema** — Ramón Viñas, Maciej Wiatrak et al. *Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation.* Nature Biotechnology (2025). doi:10.1038/s41587-025-02777-8. PMID:40854979. (cited by 31)
89. **Benchmark of generalizable single-cell perturbation prediction** — Zhiting Wei, Yiheng Wang et al. *Benchmarking algorithms for generalizable single-cell perturbation response prediction.* Nature Methods (2025). doi:10.1038/s41592-025-02980-0. PMID:41381899. (cited by 15)
90. **Biology-driven insights into single-cell foundation models** — Jialu Wu, Qing Ye et al. *Biology-driven insights into the power of single-cell foundation models.* Genome biology (2025). doi:10.1186/s13059-025-03781-6. PMID:41044630. (cited by 15)
91. **Tahoe-100M** — Jesse Zhang, Airol A Ubas et al. *Tahoe-100M : A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling.* bioRxiv (Cold Spring Harbor Laboratory) [preprint] (2025). doi:10.1101/2025.02.20.639398. (cited by 85)

**Evaluation methodology incorporated into this work** *(not part of the 91-method
cell-fate-control landscape survey; adopted to stress-test our own reachability metric)*

92. **Needles in the Haystack (DEG-aware metrics + control calibration)** — Gabriel Mateo Mejia, Henry E. Miller, Francis J. A. Leblanc, Bo Wang, Brendan Swain, Lucas Paulo de Lima Camillo. *Needles in the Haystack: Addressing Signal Dilution Improves scRNA-seq Perturbation Response Modeling and Evaluation.* ICML 2026 (Poster). Preprint: *Deep Learning-Based Genetic Perturbation Models Do Outperform Uninformative Baselines on Well-Calibrated Metrics*, bioRxiv (Cold Spring Harbor Laboratory) [preprint] (2025). doi:10.1101/2025.10.20.683304. (cited by 3)