# Five-day push — submission-readiness summary

*Cell-state reachability: a feasibility oracle with an infeasibility certificate for
cell-state engineering. This document maps every deliverable built in the five-day push
to the five questions that motivated it, and states honestly what is new work, what is
consolidation, and what remains a proposal.*

*Generated 2026-07-10. Companion to `README.md`, `NOVELTY.md`, `RESULTS.md`.*

---

## The five questions, answered

The push was organized around five questions: **(1) usefulness, (2) feasibility,
(3) impact, (4) more data / more testing, (5) SOTA comparison — and can it even be done.**
Each is answered below with the concrete artifacts that answer it.

---

## Q5 — Does it need a SOTA comparison, and can one be achieved?

**Yes, and yes — but the honest comparison is a *task-placement* argument, not a leaderboard win.**

The decisive finding from the benchmark-era literature is that deep and foundation models for
perturbation prediction (GEARS, scGPT, scFoundation, CPA, STATE) **do not yet beat simple
linear / additive / mean baselines** — and for double perturbations they do *worse* than an
additive baseline (Ahlmann-Eltze, Huber & Anders, *Nature Methods* 2025). That result is what
makes a measured-effect, non-negative-additive method well-founded rather than naive.

We ran the CPU-cheap comparison the field's own protocol calls for, on the Norman 2019 K562
double-perturbation dataset, as **two explicitly distinct tasks**:

- **Task 1 — predict an unseen double (no method sees the target).** The additive baseline
  (MSE 0.200) is the one Ahlmann-Eltze show deep models don't beat; our operator does **not**
  claim to win here — it does not attempt point prediction.
- **Task 2 — fit on 80% of a target's genes, predict the rest (the constraint's value).** The
  non-negative reachability model **wins outright**: MSE 0.125 / PCC 0.724 / cosine 0.742,
  beating the additive-rescaled baseline (MSE 0.144) on 98% of doubles *and* beating the
  unconstrained linear model (MSE 0.129) — which achieves its number using **34% negative
  weight mass** (unrealizable "anti-perturbations"). Non-negativity is a regularizer that
  *improves* generalization, not a handicap.

![Benchmark on Norman K562 doubles, two distinct tasks: additive baseline on Task 1 (blind prediction) vs. the non-negative reachability model winning Task 2 (constraint value).]({{artifact:434da92e-9498-49cf-8869-80707040ca8b}})

The positioning is made structurally by the **task taxonomy**: methods placed by whether they
are grounded in *measured* effects (x) and by output type (y: point prediction → ranking →
feasibility verdict + certificate). The top-right cell — measured-grounded *and* emitting a
feasibility verdict with a machine-checkable certificate — is occupied by this method alone.

![Task-taxonomy positioning of ~92 methods by measured-effect grounding (x) and output type (y); the measured-grounded feasibility-verdict-with-certificate cell is occupied by this method alone.]({{artifact:6ec1bbe3-d788-4c44-ac03-6f6cbfba4eb8}})

**Citations are real and verified.** All **12** benchmark-era citations were resolved live
(8 via PubMed PMID, 2 arXiv, 2 bioRxiv) with **zero fabrications**; each carries a `verified_via`
evidence trail. Deliverables: [benchmark_comparison.csv]({{artifact:d43a3360-75ef-45d3-94e9-267752bc170a}}), [benchmark_citations.csv]({{artifact:5a63e69f-d2ac-4f78-8b57-52fd56c67c37}}),
[task_taxonomy_data.csv]({{artifact:c705908c-d71c-4625-8830-dd15d1ed25ca}}), and a drop-in [manuscript_benchmark_section.tex]({{artifact:916560f6-43e4-4662-940f-7269cd058570}})
(~740 words, cites all 12 keys, explicit that we do **not** claim to beat prediction models at
their task).

**On the survey counts.** The benchmark pass surfaced exactly **2 net-new citations**
(PerturBench 2024, scGenePT 2024). Both score `no` on feasibility-verdict and certificate, so
the "only entry at the intersection" claim is **intact whether or not they are counted**. Rather
than re-thread the "91 prior / 92 total" convention through ~10 documents (a count-drift risk
with no scientific change), the two net-new citations are held as a documented, ready-to-apply
proposal in [reconciliation_note.md]({{artifact:8a55b27d-0c83-4060-8e5b-b5a0d8908aed}}) with exact deltas for both scenarios.

---

## Q4 — More data and more testing

**A second, independent public dataset confirms the central generality claim.**

On Replogle et al. 2022 genome-wide Perturb-seq (CRISPRi) — a different lab, different cell
types (K562 leukemia and RPE1 epithelium), 843 perturbations over 2,832 shared genes — the
result is **"reachability transfers, recipes don't":**

- **The verdict transfers.** Per-perturbation reach cosine correlates across cell types
  (Spearman ρ = 0.57, p = 1.8 × 10⁻⁷³); cross-basis reach (fitting a K562 target with only the
  RPE1 cone, and vice versa) stays far above the shuffled-gene null (medians 0.50 / 0.73 vs null
  0.058); binary verdict agreement is 99.3% / 100%.
- **The recipe does not.** Same-target recipe overlap (Jaccard) is 0.11, barely above the 0.053
  null — the *direction* of an achievable shift is portable biology; the *specific minimal recipe*
  is basis-specific.

![Second-dataset generality on Replogle 2022 K562/RPE1: the reachability verdict transfers across cell types (rho=0.57) while the specific recipe does not (Jaccard 0.11, near the 0.053 null).]({{artifact:0d8628e3-6176-4ffe-8aa8-4743baf9026d}})

*Honest scope:* this K562/RPE1 result was **first established in earlier project work**
(notebook `07_cross_celltype_transfer`); the push **independently reproduces it from raw effect
vectors to max|Δ| = 0.0000** and packages it as a figure-backed generality claim. A genuinely new
genome-wide download was infeasible under this session's compute budget (no GPU, ~0.5 GB free
RAM) — the plan anticipated this fallback. Deliverables:
[generality_second_dataset_summary.csv]({{artifact:5dc0f8a7-e70e-46df-9ff5-f76cb134fe30}}), [generality_second_dataset_per_perturbation.csv]({{artifact:a9728e31-3302-43e9-b751-326ef8f54011}}),
[GENERALITY_SECOND_DATASET.md]({{artifact:e80c82b2-4536-4804-b3a6-44e09e82632d}}) (with the provenance note stated up front).

---

## Q1 — Usefulness

**Four artifacts turn the method from a result into something a researcher can *use*.**

- **Per-transition design cards** — one decision-ready page per transition: verdict, confidence,
  the LOF/GOF/neither split, the ranked knockdown *and* activation recipes, the optimal
  next-screen library, druggability counts, and a modality message. Rendered for all 12 atlas
  cells. [design_cards.pdf]({{artifact:953a7cfc-e58b-4070-9b67-6f7b4132ef32}}) (12 pages) and the 4 resting-condition cards below.

![Per-transition design cards for the four resting-condition transitions: verdict, confidence, LOF/GOF/neither split, ranked knockdown and activation recipes, next-screen library, and druggability.]({{artifact:c46ad5fa-fa53-4758-9d53-2c3e66020ced}})

- **Bring-your-own-target notebook** — [bring_your_own_target.ipynb]({{artifact:ea440e5b-85cf-4f75-9447-b98861c1c6bd}}), executed clean
  end-to-end (6 cells, 0 errors, ~2.6 s). A researcher supplies a target signature (or picks a
  measured perturbation) and gets a live verdict + recipe + certificate. The worked example
  (the measured AHR+FEV double) returns **partially reachable, cosine 0.853**, the right-triangle
  identity checks to 1.000, and the top-2 knockdown recipe **recovers exactly AHR and FEV** — a
  built-in proof the method works.
- **Interactive explorer** — [reachability_explorer.html]({{artifact:dee30f53-69ea-4bca-8100-8380ef191b12}}), a single self-contained file
  (verified: zero external references, data embedded inline) that lets anyone interrogate the
  cone geometry for any of the 12 atlas cells. Staged for GitHub Pages deployment.
- **Reader-question notebook order** — [notebooks_README_reading_order.md]({{artifact:6bc125d6-8c0b-4f4a-b87a-260abc044607}}) reorders the
  9 notebooks so a cold reader hits the payoff (what can it tell me?) before the plumbing.

---

## Q3 — Impact and implication

**Eight per-target disease dossiers tie the method to the autoimmune / cell-therapy thesis.**

[target_dossiers.pdf]({{artifact:357ba5b0-98a8-4aef-a44a-914dbbfd4c4c}}) (11 pages) and [target_dossiers.md]({{artifact:ad86a6de-f997-4dcb-b229-577268204596}}) turn the strongest
nominations into decision-ready one-pagers, chosen to span the full triage surface:

- **Green-light (druggable + genetic):** JAK2, ICOS, CD3D, and — deliberately — MAPK14 as the
  *cautionary* case (chemically tractable, repeated clinical-efficacy failure).
- **Tractable-but-untried (the highest-value new leads):** IL7R (genetic score 0.959),
  ZAP70 (0.882), TET2 (0.988 — flagged as an **oncology safety signal, not a green light**).
- **Required-but-undruggable:** IRF1 (asthma, 17 associations, degrader-only) — the collision the
  method exists to surface: *don't fund the small-molecule campaign; route to CRISPRa/degrader.*

Every number traces to an input table; the thesis is grounded in the published evidence that
human genetic support (~2.6× approval odds; Minikel 2024) and *measured* target effects both
improve target-selection odds, against a backdrop of 116 CAR-T autoimmune trials by Dec 2024.
Deliverables include [dossier_targets.csv]({{artifact:ffe6d4da-a824-4b3e-8df7-b04bc25f44b7}}) with the selection rationale. The honest
caveats are stated throughout: several nodes flip modality by destination, and several activation-
certificate genes (IKZF3/Aiolos, CBLB, LAG3) are *negative* regulators — a certificate hit is a
direction *hypothesis*, not a guaranteed "activate this."

---

## Q2 — Feasibility (that it works, and reproducibly)

**The method reproduces from a clean build with one command.**

- **`reachability.py` depends only on numpy + scipy.** Verified: in a clean environment with
  only `numpy==2.4.6 scipy==1.17.1 pytest`, the module imports and its full 38-assert self-test
  passes.
- **Test suite** — [test_reachability.py]({{artifact:c1317297-3747-43a4-85a3-e961ff68de49}}), 11 tests (the packaged self-test + 10
  independent property tests re-deriving the invariants from scratch: the right-triangle identity
  `reachable_cosine² + residual_norm² = 1`, KKT/Farkas certification < 1e-5, the signed
  decomposition summing to 1, non-negative weights, spectrum monotonicity, the one-call design
  API). **All 11 pass in 16 s** in the minimal environment.
- **One-command reproduction** — [reproduce.sh]({{artifact:f782c370-0cdf-4b06-bd85-b1012293e71d}}) (`--plain` / `--venv` / `--conda`) builds
  the environment, runs the tests, and regenerates the headline verdict banner; tested end-to-end.
- **Environment files** — [environment.yml]({{artifact:2eb9eaf8-fff4-4770-bb10-a9abe327c098}}) and [requirements.txt]({{artifact:d06ecba7-14ff-4bdf-807d-f1983552ba27}}), tiered
  (core method = numpy+scipy; notebooks add pandas/matplotlib/anndata/scanpy), with the validated
  exact version pins.

---

## Submission-readiness status

| Question | Status | Primary deliverable |
|---|---|---|
| Q1 Usefulness | ✅ complete | design cards, BYO-target notebook, interactive explorer |
| Q2 Feasibility | ✅ complete | pytest suite (11/11 pass) + one-command repro |
| Q3 Impact | ✅ complete | 8 per-target disease dossiers (11-page PDF) |
| Q4 More data / testing | ✅ complete | second-dataset generality (Replogle K562/RPE1) |
| Q5 SOTA comparison | ✅ complete | benchmark (2 tasks) + task-taxonomy + verified citations |

**What is genuinely new in this push:** the two-task benchmark against field baselines with
bootstrap CIs; the task-taxonomy positioning figure; 12 independently-verified benchmark-era
citations; the per-target dossiers; the bring-your-own-target notebook; the design cards; and the
full reproducibility bundle (tests + repro script + environment files), verified in a clean env.

**What is consolidation, stated honestly:** the K562/RPE1 generality result reproduces prior
in-repo work rather than adding a new download (compute-constrained; disclosed in the artifact).

**What remains a proposal (not applied):** the 2 net-new citations and the 91→93 count re-thread,
held ready-to-apply in `reconciliation_note.md` because the disjointness claim is unaffected.
