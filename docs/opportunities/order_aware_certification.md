# Order-Aware Certification: Monotonicity Theory (I3) + Order-3 Harness & Wet Spec (E1)

**Scope label (read first).** Every quantitative result in §1–§3 is **COMPUTED
from the repo's SYNTHETIC substrate** `synth_triple_screen.npz` (planted 3-way
epistasis; 14 singles, 91 doubles, 120 labelled triples). It validates the
*code path and its discrimination*, and it makes an order-aware **property**
(monotonicity) explicit. It is **never** evidence about real biological 3-way
epistasis — the Norman screen measures zero triples, which is exactly why the
substrate is synthetic. §4 (the wet-lab acquisition spec) is **DESIGN-ONLY**: a
protocol to *replace* this synthetic validation with an in-domain one, not a
result.

The monotonicity **statement** in §1 is a mathematical property of the
projection geometry (a theorem, holding for any cone and any target); the
**numbers** that instantiate it (80→36, etc.) are computed on the synthetic
substrate.

---

## 1. Monotonicity statement (I3) — order-aware certified sets nest

### 1.1 The property (holds for any order, any data)

CombiCone certifies a measured combination *t* as emergent relative to a
reference cone `cone(A)` = the non-negative conic hull of a set of effect atoms
`A`. `reachability.project_cone` solves the non-negative least-squares
projection of *t* onto `cone(A)` and returns the **residual fraction**
`ρ(A) = ‖t − fit‖ / ‖t‖`, under a KKT-certified optimum
(`kkt_violation < 1e-8`). The solver enforces the Pythagorean projection
identity as an invariant:

> `‖t‖² = ‖fit‖² + ‖residual‖²`  (residual ⟂ fitted; `project_cone` rejects any
> solve whose `projection_identity` term exceeds the certification tolerance).

**Claim (residual monotonicity).** If `A ⊆ B` (the richer cone `B` contains
every atom of `A`), then for every target *t*:

> `ρ(B) ≤ ρ(A)`.

*Proof sketch.* `cone(A) ⊆ cone(B)`, so the best non-negative approximation of
*t* available in `B` is at least as good as the best in `A`; the NNLS objective
`‖t − fit‖` can only decrease. The denominator `‖t‖` is a property of the target
alone and is identical under both cones. Hence `ρ = ‖residual‖/‖t‖` is
non-increasing. ∎

**Claim (certified-set nesting).** The two-bar certified-emergent set is

> `C(A) = { t : p_noise(t; A) < α  AND  floor_ratio(t; A) ≥ 1.9 }`.

The **geometry** bar (raw residual, bar-a's substrate) is monotone by the above.
The effect-size bar-b (`floor_ratio = ρ / noise_floor`) shares the same monotone
numerator against a target-and-noise-defined floor. Empirically on this
substrate the enrichment produces **no new certifications** and the set nests:

> `C(B) ⊆ C(A)`  when `A ⊆ B`  —  *"certified emergence is a property of the
> cone you certify against, and enriching the cone can only shrink (never grow)
> the certified set."*

This is the order-aware generalization of the repo's stated order-2 caveat.

### 1.2 The numbers (COMPUTED on the synthetic substrate)

Enriching the reference cone from **singles (14 atoms)** to **singles + doubles
(105 atoms)** and re-certifying all 120 triples with the *unchanged* two-bar
verdict (analytic null, α=0.05, floor≥1.9):

| quantity | cone = singles | cone = singles + doubles |
|---|---:|---:|
| triples certified emergent (two-bar) | **80** | **36** |
| shrink | — | **55%** |
| strict nesting `C(rich) ⊆ C(sparse)` | — | **True** |
| newly certified under richer cone | — | **0** |
| flipped emergent → reachable | — | **44** |
| residual strictly non-increasing (per triple) | — | **True** (max Δ = −5.6×10⁻⁴) |

**Why it shrinks — resolved by planted class** (the mechanism, not a black box):

| planted class | n | certified vs singles | certified vs singles+doubles |
|---|---:|---:|---:|
| additive | 40 | 0 | 0 |
| **reducible** (2-way structure only) | 40 | 40 | **0** |
| **emergent** (genuine 3-way) | 40 | 40 | **36** |

The 80→36 shrink is driven **entirely** by the **40 reducible triples** (whose
departure from the singles cone is fully explained by their constituent pairwise
effects — once the doubles enter the cone their residual drops from ~0.20–0.26
to ~0.05–0.08, below the noise floor, so they fail bar-b) **plus 4 emergent
triples** declined at the noise floor (the four smallest planted terms,
‖Tt‖ = 1.50–1.92; the certified emergent span ‖Tt‖ = 2.06–7.00). All 120
triples remain **geometrically outside both cones** (0 land inside tolerance);
the shrink is **effect-size driven, not geometric** — precisely the pattern the
repo reports for order-2 (40→16).

**This directly generalizes the repo's order-2 result** (real Norman doubles,
`docs/FINDINGS.md §8`): enriching 105→235 atoms shrinks certified doubles
40→16, 0 newly certified, residual strictly non-increasing. Order-3 reproduces
the *same qualitative law* (strict nesting, 0 new certs, monotone residual) on
an independent substrate and higher order.

### 1.3 Monte-Carlo cross-check

The analytic null is conservative by construction (it can only withhold a
certificate, never inflate one). Re-running with `method="montecarlo"`
(n_boot=200, seed=0) gives 80→**37** (one borderline emergent triple clears),
strict nesting still holds, 0 newly certified. The analytic 36 is the frozen
headline count in `docs/FINDINGS.md`; the two paths bracket the truth as
expected.

---

## 2. Order-3 triage harness (E1) — recovery of planted emergence

`order3_harness.py` runs order-3 certification with the two-bar verdict (per-
triple `noise_sd`) and measures recovery of `triple_label_emergent` against
random. **All numbers below use the order-aware enriched cone (singles+doubles)
— the cone that isolates genuine 3-way structure from pairwise-reducible
structure.**

### 2.1 Discrimination (COMPUTED, synthetic; analytic null)

| scorer | AUROC | top-40 precision |
|---|---:|---:|
| raw residual (unreachable fraction) | 0.926 | 0.850 |
| **noise-aware z** | **1.000** | **1.000** |
| floor_ratio | 1.000 | 1.000 |

The raw residual leaks the magnitude confound (AUROC 0.926); the **noise-aware
z removes it (AUROC 1.000)** — the order-3 echo of the repo's central
size-vs-confidence result. Random top-40 precision = base rate = 0.333, so the
noise-aware z-ranked triage is a **3.0× enrichment** over random at both top-20
and top-40.

### 2.2 Two-bar verdict vs planted label (COMPUTED, synthetic)

| null | TP | FP | TN | FN | sensitivity | specificity | precision |
|---|---:|---:|---:|---:|---:|---:|---:|
| **analytic** (frozen headline) | **36** | **0** | 80 | 4 | 0.900 | 1.000 | 1.000 |
| montecarlo (n_boot=200) | 37 | 0 | 80 | 3 | 0.925 | 1.000 | 1.000 |

**36/40 planted-emergent certified, 0/80 additive-or-reducible false
positives** — reproducing the repo's frozen k-way headline exactly. The 4
analytic misses are the smallest planted terms (‖Tt‖ 1.50–1.92), correctly
declined at the noise floor.

### 2.3 Honest boundary — *prospective* triage does NOT transfer to this substrate

The repo's **training-free `−cos` triage** (singles-only, no measured triple)
enriches 2.4× on **real Norman doubles**. On this **synthetic** substrate it is
**at chance**:

- AUROC(−cos, planted label) = **0.53**;
  Spearman(−cos score, planted label) = 0.048 (p=0.60, n.s.);
  top-20 precision 0.35 vs 0.33 base rate = **1.05×**.

Reason (verified): the synthetic generator plants 3-way emergence as a term
**orthogonal to the pairwise-cosine axis**, so the pairwise-cosine heuristic has
no signal to find here. **This is an honest limitation of the synthetic
substrate, not of the method** — it is precisely why prospective triage must be
recalibrated per screen (the same per-screen-recalibration boundary the repo
reports for ALMANAC synergy and CaRPool). The **post-measurement certification**
(the harness's actual product) is unaffected: it recovers planted emergence
perfectly (AUROC 1.000) because it *measures* the triple rather than predicting
it from singles.

`order3_harness.py` reports both the certification metrics (its product) and,
as a documented negative control, flags that prospective `−cos` on this
substrate is uninformative.

---

## 3. Reusable harness — `order3_harness.py`

numpy/scipy-only, fail-closed, matches the `combicone` house style. Imports
`combicone` + `reachability`; adds no new dependency (dependency-free AUROC /
rank utilities).

```python
import order3_harness as h
data = h.load_synth_triple_screen("synth_triple_screen.npz")

mono = h.monotonicity_report(data, method="analytic")
print(mono.n_certified_sparse, "->", mono.n_certified_rich)   # 80 -> 36
print(mono.is_subset, mono.n_newly_certified)                 # True 0
print(mono.residual_non_increasing)                           # True

tri = h.triage_report(data, method="analytic", cone="rich")
print(tri.auroc["z"])                                         # 1.0
print(tri.two_bar_confusion)                                  # TP=36 FP=0 ...
```

CLI:

```
python order3_harness.py --npz synth_triple_screen.npz --method analytic \
    --cone rich --json order_aware_certification.json
```

Public surface: `load_synth_triple_screen`, `certify_triples`, `two_bar_mask`,
`monotonicity_report`, `triage_report`, plus `auroc` / `top_k_precision` /
`random_top_k_precision`. Every entry point carries the synthetic-scope
disclaimer; the CLI prints it as a banner.

---

## 4. Wet-lab acquisition spec for MEASURED triples (DESIGN-ONLY)

**Purpose.** The order-3 result above is synthetic. This spec is the concrete
experiment that would replace it with an in-domain one — measuring a modest
panel of real triple perturbations and testing whether the certified-set
shrinkage seen from singles→doubles (§1) continues to order 3, and whether real
3-way emergence exists above the noise floor. **Nothing here is a result; it is a
protocol.** This is exactly the "measure a modest panel of triples" experiment
the manuscript Outlook (`sections/50_discussion.tex`) names as the highest-value
next step.

### 4.1 Platform

Two established multiplexed CRISPR single-cell platforms both deliver ≥3
perturbations per cell with a transcriptome readout on the additive
(log-normalized) scale CombiCone consumes:

| platform | mechanism | why it fits | reference |
|---|---|---|---|
| **CROP-seq / multiplexed CROP-seq** (CRISPRa or CRISPRi, e.g. dCas9-based à la Norman) | polymerized guide array, guide captured in the transcript | matches the Norman CRISPRa modality the singles/doubles cone is built on — keeps the dictionary on one mechanism | CROP-seq / Norman combinatorial CRISPRa method papers (verify citations before external use) |
| **CaRPool-seq** (Cas13d RNA-targeting knockdown) | one array encodes up to 4 guides; direct capture; high multiplicity | the repo already ingests CaRPool as its transfer screen; native to ≥3-way arrays; loss-of-function complement | CaRPool-seq method paper, as cited in the repo (`docs/FINDINGS.md`: "Wessels 2023, GSE213957") — confirm the accession/citation against the source before external use |

**Recommendation:** CROP-seq-multi (CRISPRa) to keep the whole singles→doubles→
triples ladder on one mechanism (cleanest test of the monotonicity law on real
data); CaRPool-seq if a loss-of-function 3-way panel is the scientific target.

### 4.2 Design and scale

Order-3 certification requires, for every measured triple, that **all three
constituent singles and (for the order-aware / monotonicity test) the three
constituent doubles are also measured in the same screen** (this is what lets
the cone be enriched singles→singles+doubles, and what `screen_ingest.py`
resolves via `has_all_singles`).

**Combinatorial budget** (why the panel must be focused):

| N singles | C(N,2) doubles | C(N,3) triples |
|---:|---:|---:|
| 10 | 45 | 120 |
| 12 | 66 | 220 |
| 20 | 190 | 1,140 |
| 50 | 1,225 | 19,600 |

**Recommended focused panel** — a self-contained ladder that measures *every*
lower-order term (so the monotonicity test is exact) at a tractable size:

- **N = 10–12 single-gene targets**, chosen from genes with the largest,
  most orthogonal single-gene effects (the regime where 3-way emergence is most
  likely and most certifiable).
- **Full ladder:** measure all singles + all C(N,2) doubles + all C(N,3)
  triples. At **N=10: 10 + 45 + 120 = 175 conditions**; at **N=12: 12 + 66 +
  220 = 298 conditions**; **+1 combinatorial control** (`NT_NT_NT` /
  triple-safe-harbor).
- If C(N,3) is still too large, measure a **triage-selected subset** of triples
  (≥120, to match the synthetic substrate's statistical power) plus **all**
  their constituent singles and doubles.

**Guides / replicates.**

- **2 guides per gene** minimum (independent guides guard against
  guide-specific off-target artifacts; the repo's guide-pair transfer analysis
  is built on exactly this redundancy). A triple then has 2³ = 8 guide-
  combinations; measure ≥2–3 distinct guide arrays per gene-triple as biological
  replicates of the *condition*.
- **No separate bulk replicate needed for the noise model** — the certificate's
  per-gene `noise_sd` is estimated by a **cell split-half within each
  condition** (`screen_ingest.split_half_noise_sd`, `|e₁−e₂|/2`). This is the
  single most important sizing constraint below.

**Cells per condition.**

- **Hard floor: ≥ 12 QC-passing cells per condition** (`min_cells_per_half = 6`
  per half; conditions below `2×min` return `NaN` noise and cannot be
  certified).
- **Target: 150–300 QC-passing cells per condition** for a stable pseudobulk
  effect vector and a low-variance split-half SE (the Norman/CaRPool substrates
  sit in this range). The control needs the same or more.
- At N=12 (298 conditions + control) × 200 cells ≈ **~60,000 QC cells**; budget
  ~2× for guide-assignment/QC dropout → **~120,000 cells targeted** (2–3 10x
  lanes). N=10 (175 conditions) ≈ ~35k QC / ~70k targeted.

### 4.3 Exact arrays `screen_ingest.py` needs

`ingest_screen` (arrays path) turns the raw screen into the CombiCone inputs.
Supply:

```python
import screen_ingest as si
sub = si.ingest_screen(
    expression = X,            # (n_cells, n_genes) float, ALREADY log-normalized
                               #   (CP10k/log1p or equivalent additive scale)
    conditions = labels,       # (n_cells,) str, one label per cell (grammar below)
    gene_names = genes,        # (n_genes,) str, shared gene axis
    control_label = "NT_NT_NT",# combinatorial control token (or "ctrl")
    separator = "+",           # gene-join token: "A+B+C"
    arm_handling = "merge",    # pool GENE+ctrl / ctrl+GENE spellings of a single
    compute_noise = True,      # split-half per-gene SD (the certificate's noise_sd)
    min_cells_per_half = 6,    # >=12 cells/condition required to certify
)
```

**Condition-label grammar** (`screen_ingest.parse_condition`):

| kind | label examples |
|---|---|
| control | `"NT_NT_NT"`, `"ctrl"`, or `"NT_NT"` (set `control_label=`) |
| single | `"AHR"`, or `"AHR+ctrl"` / `"ctrl+AHR"` (arm with control collapses to the single) |
| double | `"AHR+KLF1"` |
| triple | `"AHR+KLF1+SET"` (three genes joined by `separator`, none the control token) |

`ingest_screen` then yields a `ScreenSubstrate` with exactly:

- `atoms` **(n_atoms, n_genes)** — single-gene effect vectors = the cone;
- `atom_names` **(n_atoms,)**;
- `genes` **(n_genes,)** — shared gene axis;
- `combos` — a `ComboRecord` per measured combination, each with `effect`
  (n_genes,), `noise_sd` (n_genes,) from the split-half, `constituent_singles`,
  and `has_all_singles` (True required for order-aware certification).

Then, per measured triple:

```python
kw = sub.certify_ready("AHR+KLF1+SET")          # dict ready to splat
cc.certify_emergence(cone_atoms=kw["cone_atoms"],   # singles cone
                     measured_combo=kw["measured_combo"],
                     noise_sd=kw["noise_sd"])
```

For the **monotonicity test on real data**, build the enriched cone yourself
from the ingested doubles and re-certify — exactly what `order3_harness.py`
does on the synthetic substrate:

```python
singles = sub.atoms
doubles = np.vstack([c.effect for c in sub.combos if len(c.genes) == 2])
cone_rich = np.vstack([singles, doubles])
cc.certify_emergence(cone_atoms=cone_rich, measured_combo=triple_effect,
                     noise_sd=triple_noise_sd)
```

### 4.4 Pre-registered analysis (what the measured screen would test)

1. **Monotonicity on real data (the I3 test):** certify each measured triple
   against cone(singles) then cone(singles+doubles); confirm strict nesting and
   residual non-increase (the synthetic law of §1 predicts this must hold
   geometrically — the *biological* question is how many real triples survive
   bar-b under the enriched cone, i.e. whether genuine 3-way emergence exists).
2. **Certified 3-way emergence rate:** two-bar certified fraction under the
   enriched cone = the real-data analog of the synthetic 36/40.
3. **Prospective triage recalibration:** fit the triage model on the
   singles+doubles of *this* screen and test top-k triple emergence rate vs
   random (recalibrated per screen, per the repo's stated boundary — do **not**
   assume the Norman `−cos` rule transfers).

---

## 5. Files

| file | kind | contents |
|---|---|---|
| `order_aware_certification.png` | figure (300 dpi) | (a) certified-set shrink 80→36 by class; (b) per-triple residual monotonicity scatter; (c) order-3 triage ROC (z AUROC 1.000) |
| `order3_harness.py` | code | reusable order-3 certification + monotonicity + triage harness (numpy/scipy only) |
| `order_aware_certification.json` | results | full computed report (monotonicity + triage), analytic null |
| `order_aware_certification.md` | this doc | monotonicity statement + wet spec |

**Honesty ledger.** §1.2/§1.3/§2 are **COMPUTED on SYNTHETIC data**. §1.1 is a
**theorem** (holds for any data). §4 is **DESIGN-ONLY**. No real-data 3-way
emergence claim is made anywhere; the only real-data order-k evidence in the
repo remains the order-2 Norman shrink (40→16), which this work generalizes to
order 3 *in property* and *on a synthetic substrate*.
