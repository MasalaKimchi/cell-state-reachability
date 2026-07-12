# Cell-State Reachability

**Built with Claude — Life Sciences Hackathon · Research / Lab Track**

## Can knockdown point a cell where you want it to go?

A GPS does more than suggest a route: it can tell you when the road you need is not on the
map. This project asks the analogous question for cell engineering. Given a library of
effects measured in a real perturbation screen, can non-negative combinations of those
effects point the transcriptome toward a desired state—and where does the measured
intervention class fall short?

This is a **directional feasibility test relative to a measured screen**, not an oracle of
biological truth. It returns a held-out, null-calibrated alignment; a ranked perturbation
mixture and greedy sparse candidate panel; and, when the target lies outside the measured
cone, a Farkas/KKT dual certificate for the full mismatched direction.

> **New here?** Start with the [2-minute pitch](./SUBMISSION.md), then open the
> [interactive walkthrough](./app/index.html). The [validation report](./docs/VALIDATION_REPORT.md)
> separates what is verified, what is model-dependent, and what still needs a wet-lab test.

---

## The result in 30 seconds

For the flagship **Th2 → Th1** direction in resting primary human CD4⁺ T cells:

- Fixed-split held-out cosine: **0.448** (in-sample 0.627).
- Shuffled-target control: larger than all **60** shuffles (plus-one empirical
  **p = 1/61**; descriptive z ≈ 24).
- Staged modality proxy: **39% measured LOF / 25% sign-flipped GOF proxy / 35% neither**.
- KKT violation: **1.1 × 10⁻¹¹**, certifying the numerical cone projection—not biological
  efficacy.

The result is **partly directionally reachable by knockdown**, not fully reachable as a
cell phenotype. Positive residual coordinates rank unmet readouts for CRISPRa or
de-repression follow-up; an individual gene on that list is not thereby proven causal.

Across 12 transition–condition cases, knockdown is never the majority component (mean LOF
fraction 0.34). Those atlas cells each exceed eight shuffled screening controls, but the
eight-shuffle z range is descriptive rather than a formal permutation tail test.

## The method in plain language

Each CRISPRi perturbation has a measured transcriptional fingerprint. A fingerprint can
contain both increases and decreases in gene expression. The non-negativity constraint does
**not** mean “knockdown only lowers expression”; it means the complete measured effect may
be combined with non-negative weight but may not be reversed into an unmeasured
“anti-knockdown.”

Under a first-order additivity approximation, these combinations form a convex cone:

1. Build a dictionary `E` of measured perturbation effects.
2. Define a target direction `d = destination − source`.
3. Project `d` onto `{Eᵀw : w ≥ 0}` with non-negative least squares.
4. Evaluate the fitted weights on held-out genes and against shuffled targets.
5. Use a greedy spectrum to propose compact panels. The knee is a practical sparse panel,
   not a proof of globally minimum cardinality.
6. If the residual is non-zero, use it as a separating direction. The complete residual is
   the mathematical certificate; its largest positive coordinates are follow-up hypotheses.

The signed decomposition is also explicitly model-based: after the measured LOF fit, the
remaining residual is fitted with `-E`, treating the negative of each knockdown effect as a
hypothetical activation effect. Real CRISPRa effects need not be exact mirrors of CRISPRi.

## What the data actually contain

The primary source is Zhu et al. 2025, *Genome-scale Perturb-seq in primary human CD4⁺ T
cells* (Marson and Pritchard labs; CZI Virtual Cells Platform, MIT license).

- Dataset card: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1

The derived Tier-2 H5AD is **33,983 perturbation–condition profiles × 10,282 readout
genes**, representing 11,526 distinct targeted genes over Rest, Stim8hr, and Stim48hr. It is
not 33,983 unique knockdown genes. After source-screen quality filters, the working
dictionaries contain 6,871, 7,155, and 7,195 generators respectively.

- **Tier 1** (~38 MB): seven supplementary CSVs for target signatures, quality control,
  donor metadata, guide design, and immune-disease enrichment.
- **Tier 2** (16.8 GB): `data/GWCD4i.DE_stats.h5ad`, read selectively and cached; the raw
  ~22-million-cell dataset is not required.

See [`data/README.md`](./data/README.md) for acquisition and exact grain. Raw data remain
gitignored and must never be committed.

## What is novel—and what is not

Minimal regulator sets, network control, and cell-state conversion are established fields.
The source paper already reports Th1/Th2 and aging regulators; their recovery here is a
positive control, not a discovery.

The claimed delta is the combination of:

- directly measured perturbation-effect vectors rather than an inferred regulatory network;
- a target-specific convex-cone directional verdict rather than only a similarity ranking;
- a dual separating certificate for an outside-the-cone target direction;
- held-out-gene and shuffled-target calibration; and
- explicit modality and additivity caveats attached to every candidate panel.

A structured survey of 91 prior methods found no prior entry combining measured-effect
grounding with this target-specific certificate. That is a bounded literature-survey claim,
not proof that no adjacent method exists.

## Scope, stated before the repository map

- Transcriptomic alignment is not functional rescue.
- The cone is unbounded and the main score is directional; it does not establish an
  achievable biological dose or magnitude.
- Multi-gene combinations assume approximate additivity. In 126 Norman K562 CRISPRa
  doubles, measured-versus-additive direction has median cosine 0.71; CD4⁺ combinations
  remain unvalidated.
- The primary effects are donor-collapsed across four donors, so true leave-one-donor-out
  validation is unavailable.
- The CEBPA transfer is one held-out state in a different screen; full cross-atlas transfer
  remains future work.
- Open Targets and literature annotations are saved snapshots and prioritization aids.
- Every nomination is a wet-lab hypothesis, never a validated target.

## Start here

| If you want… | Open |
|---|---|
| **2-minute pitch** | [`SUBMISSION.md`](./SUBMISSION.md) |
| **Guided interactive story** | [`app/index.html`](./app/index.html) |
| **Paper** | [`manuscript/main.pdf`](./manuscript/main.pdf) |
| **Validation report** | [`docs/VALIDATION_REPORT.md`](./docs/VALIDATION_REPORT.md) |
| **Full technical record** | [`docs/Technical_Dossier.pdf`](./docs/Technical_Dossier.pdf) · [source](./docs/Technical_Dossier.md) |
| **Verify core software** | [`reproduce.sh`](./reproduce.sh) |
| **Re-run analyses** | [`notebooks/`](./notebooks/) · [`scripts/`](./scripts/) |

## Repository map

```text
cell-state-reachability/
├── SUBMISSION.md            # accessible 2-minute pitch
├── README.md                # orientation, scope, and reproduction routes
├── AGENTS.md / CLAUDE.md    # enduring facts and scientific guardrails
├── DEMO_VIDEO_SCRIPT.md     # 3-minute demo plan
├── reachability.py          # cone fit, dual certificate, nulls, spectrum, additivity risk
├── reproduce.sh             # software verification: 11 tests + synthetic self-test
├── requirements.txt         # pinned core and analysis dependencies
│
├── docs/
│   ├── VALIDATION_REPORT.md # audit findings, verified claims, and residual risks
│   ├── Technical_Dossier.pdf
│   ├── Technical_Dossier.md
│   └── figures/
│
├── manuscript/
│   ├── main.pdf
│   ├── main.tex + sections/
│   ├── figures/
│   └── manuscript_facts.json
│
├── app/                     # seven self-contained interactive explorers + guided index
├── notebooks/               # 01–09 + bring_your_own_target
├── scripts/                 # batch analysis and validation drivers
├── results/                 # committed result tables and fact ledgers
├── tests/                   # 11 method/property tests
├── data/                    # local, gitignored Tier-1/Tier-2 data
└── analysis_cache/          # large local intermediates and small tracked summaries
```

## Verify the software

```bash
bash reproduce.sh
```

This runs the 11-test suite and the module's synthetic invariant battery. It verifies the
implementation, KKT conditions, staged decomposition, null machinery, and design API. It
does **not** regenerate the real-data headline from a fresh clone because the 16.8 GB Tier-2
matrix is intentionally not committed.

The synthetic test suite currently takes several minutes on the reference laptop; the core
full-data cone solve itself is much faster once the reduced matrix is cached.

## Reproduce the analysis

After obtaining the data and building `analysis_cache/atlas_work/inputs.npz` through
notebook 02:

```bash
python scripts/run_atlas.py
python scripts/run_nulls.py
python scripts/run_bootstrap.py
python scripts/run_a1_sensitivity.py
python scripts/run_iv_compliance.py
python scripts/run_deg_weighted_eval.py
python scripts/run_heldout_split_stability.py
```

The notebook route is documented in [`notebooks/README.md`](./notebooks/README.md): `01`
EDA → `02` primary atlas → `03` K562 generalization → `04` design toolkit → `05` target-ID
showcase → `06` reinforcement → `07` cross-cell-type transfer → `08` DEG weighting → `09`
causal validation, plus `bring_your_own_target`.

## License

MIT. Data: MIT via the CZI Virtual Cells Platform. Please cite Zhu et al. 2025.
