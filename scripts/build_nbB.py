"""Builder for 09_causal_validation_dossier.ipynb (NB-B, "Can I trust it?").

Narrative spine = the six-assumption identifying stack. Each Part-A analysis is mapped to
the assumption it stress-tests. Numbers are loaded from the frozen result files so the
notebook's tables are bit-identical to the saved artifacts.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

C = []
md = lambda s: C.append(new_markdown_cell(s))
co = lambda s: C.append(new_code_cell(s))

# ---------------------------------------------------------------- title
md(r"""# Can I trust the verdict?  — a causal-validation dossier

**Notebook 09 — the reviewer's notebook.**

The reachability oracle returns a falsifiable verdict: *is this cell-state shift reachable
by CRISPRi knockdowns, and if so, by which minimal recipe?* Every other notebook answers
**"what does the oracle say?"** This one answers the question a referee asks next:
**"under what assumptions is that answer causal, and how hard would each assumption have to
fail to overturn the verdict?"**

The oracle is not an observational model — it never fits a gene-regulatory graph. Its inputs
are **average treatment effects** measured in a genome-scale randomized CRISPRi Perturb-seq
experiment (Zhu et al. 2025), and a reachability verdict is a **counterfactual query on a
compound intervention**: *does there exist a non-negative combination of measured single-gene
interventions whose effect points at the target state?* That places the method in the
**design-based** camp of causal inference (Fisher/Neyman randomization, Rubin potential
outcomes, LATE/IV), not the **observational-discovery** camp (PC/GES, NOTEARS, causal
representation learning). Identification comes from the experiment's randomized guide
assignment, and infeasibility comes with a **Farkas certificate** — an exact proof, not a
p-value.

This notebook makes the identifying assumptions **explicit**, tests the testable ones, and —
where a test is impossible in silico — quantifies **how far the assumption must fail** to flip
the verdict. That last move (sensitivity analysis in place of point identification) is the
modern gold standard, and it is what makes the trust layer publishable rather than cosmetic.

> Companion prose: `../CAUSAL.md` — the causal positioning and the full A1–A6 / B1–B4
> agenda. This notebook executes the computable half (Part A).""")

# ---------------------------------------------------------------- setup
md(r"""## 0 · Setup

We load the frozen effect matrices and target directions (`atlas_work/inputs.npz`, the same
`E` and `d` the headline notebooks build), the per-effect standard errors
(`atlas_lfcSE.npz`), and the Part-A result files. The oracle itself is the unchanged
`reachability.py`.""")

co(r"""import os, sys, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Image, display

REPO = os.environ.get("REACH_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
import reachability as R

WORK    = os.path.join(REPO, "analysis_cache", "atlas_work")
RESULTS = os.path.join(REPO, "results")
def rp(name):                       # result-file path: prefer repo root, else atlas_work
    for base in (REPO, WORK, RESULTS):
        cand = os.path.join(base, name)
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(name)

inp = np.load(os.path.join(WORK, "inputs.npz"), allow_pickle=True)
COND = "Stim48hr"                   # the differentiated state used across Part A
E = np.asarray(inp[f"E_{COND}"]); gene = np.asarray(inp[f"gene_{COND}"]).astype(str)
targets = {n: np.asarray(inp[f"t_{n}"]) for n in
           ["toward_Th1", "toward_Th2", "toward_younger", "toward_older"]}
print(f"E[{COND}] = {E.shape}  |  {len(targets)} target axes  |  reachability.py loaded")""")

# ---------------------------------------------------------------- assumption spine
md(r"""## 1 · The identifying stack — six assumptions the verdict rests on

A reachability verdict is causal **only if** a chain of assumptions holds. Naming them is the
contribution of this notebook; the rest of it tests them one by one. Each row lists the
assumption, whether it is testable with the data in hand, and the analysis that does the
testing.""")

co(r"""stack = pd.DataFrame([
 ["I",   "Unbiased ATEs",           "guide assignment is randomized; E[g,:] is an unconfounded average treatment effect",
        "partially — negative controls (A5), and the IV/compliance rescaling", "A5, IV layer"],
 ["II",  "Additivity / no epistasis","a compound recipe's effect is the non-negative sum of its parts",
        "yes — calibrated on Norman 2019 doubles", "A3 (+ existing additivity_risk)"],
 ["III", "SUTVA / no interference",  "one cell's perturbation does not change another's outcome (no cytokine spillover)",
        "no in silico — needs arrayed-vs-pooled wet lab", "B1 (Future Work)"],
 ["IV",  "Homogeneous effect",       "the population-average effect holds within each cell subtype",
        "no in-sandbox — needs raw single-cell counts", "A2 (scaffold)"],
 ["V",   "Exclusion restriction",    "guides act on the outcome only through their on-target knockdown",
        "yes — off-target back-route dropped; weak-instrument intervals", "IV layer, A4"],
 ["VI",  "Transportability",         "an effect measured in this assay carries to the target context",
        "partially — cross-cell-type transfer; aging is the open case", "B4 (Future Work), NB07"],
], columns=["#","assumption","statement","testable in silico?","analysis"])
stack.style.hide(axis="index") if hasattr(stack.style,"hide") else stack""")

md(r"""**Reading the stack.** Assumptions **I, II, V** are testable now and are the subject of
§2–§6 below. **IV** (homogeneous effect) needs raw single-cell counts that are not in the
repo, so A2 ships as a scaffold. **III** (SUTVA) and **VI** (transportability) need wet-lab
or new context data and are carried in the manuscript's Future Work (B1, B4). The honest
headline is that **no assumption is assumed away** — each is either tested or explicitly
booked as a caveat.""")

# ---------------------------------------------------------------- IV/compliance
md(r"""## 2 · Assumptions I & V — the instrument, compliance, and the exclusion restriction

Each CRISPRi guide is an **instrument** `Z`; the realized on-target knockdown is the
**treatment** `T`; the transcriptome shift is the **outcome** `Y`. The published effect
`E[g,:]` is an **intent-to-treat** (ITT) contrast (assigned-guide vs control). Dividing by the
first-stage compliance `π` (realized fractional knockdown) would give the **LATE** —
the effect among compliers.

The key structural fact for this method: the reachable cone is **invariant to positive
per-generator rescaling**. Rescaling `E[g,:] → E[g,:]/π_g` cannot move the verdict — only the
recipe weights. We verified this to machine precision across all 12 target×condition cells
(`run_iv_compliance.py`).""")

co(r"""iv = pd.read_csv(rp("iv_compliance_verdicts.csv"))
show = iv.rename(columns={"cond":"condition","cos_itt":"cos_ITT","cos_valid":"cos_valid",
                          "delta":"Δcosine (drop invalid)","kkt":"kkt_max","pct_drop":"% invalid dropped"})
cols = ["target","condition","cos_ITT","cos_valid","Δcosine (drop invalid)","kkt_max","% invalid dropped"]
print("IV / compliance verdict table (ITT vs valid-instrument), 12 cells:")
display(show[cols].round(6))
print(f"\\nExclusion-restriction drop arm: max |Δcosine| = {iv['delta'].abs().max():.1e}, "
      f"max KKT = {iv['kkt'].abs().max():.1e}")
inv = pd.read_csv(rp("late_rescaling_invariance.csv"))
mx = float(np.max(np.abs(inv["abs_delta"])))
print(f"LATE 1/π rescaling invariance: max |Δcosine| over 12 cells = {mx:.2e}  "
      f"(verdict unchanged by compliance rescaling; scale factors up to {inv['scale_max'].max():.1f}×)")""")

md(r"""**Exclusion restriction (V).** Identification assumes a guide reaches the outcome *only*
through its on-target knockdown — the off-target back-route `Z→Y` is excluded. We enforce this
operationally by dropping generators with no significant on-target knockdown (the
**valid-instrument** arm). Removing them moves the reachable cosine by at most **|Δ| = 4×10⁻⁴**
across all 12 cells (table above), because such guides carry negligible weight. The
compliance-annotated recipe (§6) then flags any generator whose LATE weight leans on a
**weak** instrument.""")

# ---------------------------------------------------------------- A1 sensitivity
md(r"""## 3 · Assumption I, quantified — the verdict sensitivity radius (headline)

The verdict is a point estimate on noisy effects. **A1** asks: *how far would the effects have
to be wrong to flip it?* — answered in units of the **measured standard error**, an
E-value-style statement (VanderWeele; Cinelli–Hazlett).

Two arms:
- **(M) measurement error** — resample `E ~ N(Ê, SE²)`, re-solve, report the bootstrap CI of
  the reachable cosine against the null threshold.
- **(C) coordinated bias** — the worst-case adversary that attenuates the target-aligned
  component of every generator by a fraction β (re-solve-proof); we bisect for the β\* that
  drops the cosine to the null p99 and express it in median-SE units (`kSE*`).""")

co(r"""a1 = pd.read_csv(rp("a1_sensitivity_radius.csv"))
disp = a1[["target","baseline_cosine","null_p99","null_z","boot_lo","boot_hi",
           "boot_clears_null","beta_star","kse_star"]].round(4)
print("A1 — verdict sensitivity radius (condition = Stim48hr):")
display(disp)
th1 = a1[a1.target=="Th1"].iloc[0]
print(f"\\nHeadline (Th1): the boot 95% CI [{th1.boot_lo:.3f},{th1.boot_hi:.3f}] CLEARS the null "
      f"p99={th1.null_p99:.3f} (robust to measurement error),")
print(f"but a COORDINATED bias of only β*={th1.beta_star:.1%} of the target-aligned signal "
      f"= {th1.kse_star:.3f} median-SE units flips it.")""")

co(r"""# render the A1 figure if present
try:
    display(Image(filename=rp("fig_a1_sensitivity.png")))
except FileNotFoundError:
    print("fig_a1_sensitivity.png not yet built (see Phase-5 figure step)")""")

md(r"""**Why this is the right probe.** Undirected measurement noise *adds anisotropy*, which
makes any target easier to reach by chance — so random error tends to **inflate** the reachable
cosine, and the null already corrects for it. The verdict's genuine exposure is therefore not
random error but **coordinated systematic bias** against the target direction — exactly the
failure mode a randomized design rules out. That the polarization axes (Th1/Th2) clear the null
under bootstrap but flip under an ~8% coordinated attenuation is the quantitative case for why
the design-based identification is load-bearing.

The aging axes (**younger/older**) at Stim48hr sit *below* the null p99 (β\*≈0.002): the oracle
honestly reports them as **not separated from chance** in the differentiated state — a caveat
the raw cosine (0.57–0.60, superficially higher than Th1's 0.53) would have hidden.""")

# ---------------------------------------------------------------- A5 negative controls
md(r"""## 4 · Assumption I — negative-control outcomes

If the effect estimates carried broad confounding (e.g. a technical axis correlated with guide
assignment), the reachable vector would deposit mass on genes with **no business** responding to
a polarization push — housekeeping, proteasome, structural genes. **A5** is the
negative-control-outcome test: project the fitted reachable vector onto a housekeeping/off-pathway
panel and compare per-gene mass against a positive-control panel of canonical polarization genes
(TBX21, GATA3, IFNG, IL4, STAT1/4/6, …).""")

co(r"""a5 = pd.read_csv(rp("a5_negative_control_outcome.csv"))
disp = a5[["target","pos_per_gene","neg_per_gene","enrich_pos_over_neg"]].round(5)
disp.columns = ["target","pos-ctrl mass/gene","neg-ctrl mass/gene","enrichment (pos/neg)"]
print("A5 — negative-control-outcome projection:")
display(disp)
print(f"\\nAcross all four targets the reachable vector puts "
      f"{a5.enrich_pos_over_neg.min():.1f}×–{a5.enrich_pos_over_neg.max():.1f}× more mass per gene "
      f"on polarization positive controls than on housekeeping/proteasome negative controls.")""")

md(r"""**Interpretation.** A valid design puts near-zero mass on negative-control outcomes. The
4–5.5× enrichment of positive over negative controls is what we want: the reachable vector
concentrates on biology, not on housekeeping machinery. (Ribosomal and mitochondrial genes are
absent from the readout — they were filtered during signature construction — so the panel is
housekeeping + proteasome.)""")

# ---------------------------------------------------------------- A6 construct validity
md(r"""## 5 · Assumption I — directional construct validity

A stronger falsification than the shuffled-target null: does the oracle recover **known signed
immunology**? T-bet (**TBX21**) is the master Th1 driver; **GATA3** the master Th2 driver.
Knocking each down should point the cell *away* from its own lineage. **A6** measures the signed
cosine of each master regulator's own knockdown-effect vector with the Th1 axis, against the
generator-wide null (mean ≈ 0.004, sd ≈ 0.038).""")

co(r"""a6 = pd.read_csv(rp("a6_construct_validity.csv"))
print("A6 — directional construct validity (signed cosine vs Th1 axis; targets are antiparallel):")
display(a6.round(4))
canon = a6[a6.regulator.isin(["TBX21","GATA3"])]
print("\\nCanonical check:")
for _,r in canon.iterrows():
    print(f"  {r.regulator}: KD cos→Th1 = {r.kd_cos_toward_Th1:+.4f}  (null z = {r.null_z:+.2f})")""")

md(r"""**Interpretation — honest.** Both canonical regulators carry the **correct sign** (TBX21↓
points away from Th1; GATA3↓ points away from Th2, hence toward Th1), and 4 of 5 lineage drivers
do. But only **STAT6** exceeds the null band (z ≈ +3.0); the others are correct-signed yet within
noise. This is *expected*, not disappointing: the method's own premise is that **no single
generator dominates** the target direction (that is why a whole recipe is needed), so a
single-knockdown alignment *should* be small. Construct validity here means the signs line up
with immunology — they do — not that any one gene reaches the target alone. (The single-generator
*drop*-and-re-solve test returns ≈0 for every regulator: with ~7,195 collinear generators, removing
one changes nothing — a documented redundancy, and the reason we score the effect vector directly.)""")

# ---------------------------------------------------------------- A4 weak instruments
md(r"""## 6 · Assumption V, quantified — weak-instrument-robust recipe intervals

The exclusion restriction is only as good as the first stage. If a generator's on-target
knockdown `π` is small and noisily estimated, its LATE recipe weight (`∝ 1/π`) is poorly
identified — the classic **weak-instrument** problem. **A4** propagates first-stage `π`
uncertainty into an Anderson–Rubin-flavored confidence set for each recipe weight, valid as
`π → 0`.""")

co(r"""a4 = pd.read_csv(rp("a4_weak_instrument_intervals.csv"))
print("A4 — headline Th1/Stim48hr recipe, weight-multiplier (1/π) intervals:")
display(a4.round(3))
demo = pd.read_csv(rp("a4_weak_instrument_demo.csv"))
print("\\nStress test — the method on a genuinely weak instrument (SNX5, π≈0.07):")
display(demo)
print("SNX5's weight-multiplier CI runs to the hundreds — the AR signature of a set that is "
      "effectively unbounded as π→0. The headline recipe does not rely on such an instrument.")""")

md(r"""**Interpretation.** The headline recipe is clean — no generator's first-stage lower CI
dips into weak territory — though RARA (π≈0.46) and KIF5B (π≈0.44) have wide (~4.5×) weight
multipliers and should be read as *directionally* right but quantitatively soft. The stress test
confirms the interval **does** blow up for a real weak instrument (SNX5), so the check has teeth.""")

# ---------------------------------------------------------------- A2 scaffold + summary
md(r"""## 7 · Assumption IV — conditional reachability (scaffold)

The population-average counterfactual (assumption IV) can hold for no individual subtype. **A2**
would stratify effect estimation by a baseline covariate (naive/memory, cell-cycle phase) and
re-solve the oracle per stratum, reporting stratum-stability and recipe reordering. This needs
**raw single-cell counts**, which are not in the repo — the in-sandbox `GWCD4i.DE_stats.h5ad` is
a DE-summary with the cells already aggregated away. A2 therefore ships as a **scaffold** (code
path + expected-output schema) and is booked as future work pending a CZI VCP raw-count pull.""")

co(r"""schema = json.load(open(rp("a2_expected_output_schema.json")))
print("A2 expected-output schema (scaffold — not executed):")
print(json.dumps(schema, indent=1))""")

md(r"""## 8 · Verdict on trust — the stack, scored

| # | assumption | status in this notebook |
|---|------------|--------------------------|
| I | unbiased ATEs | **tested** — negative-control outcomes clean (A5), construct validity signs correct (A6), verdict robust to measurement error (A1-M) |
| II | additivity | **bounded** — pairwise-calibrated saturation score + k-way directional decay (A3); recipe confidence declines with cardinality |
| III | SUTVA | **booked** — cytokine-spillover screen specified (B1, Future Work) |
| IV | homogeneous effect | **scaffolded** — A2 code path ready; needs raw counts |
| V | exclusion restriction | **tested + quantified** — off-target route dropped (|Δ|≤4e-4), weak-instrument intervals (A4) |
| VI | transportability | **partial** — cross-cell-type transfer holds direction (NB07); aging is the open case (B4) |

The verdict for the **polarization axes (Th1/Th2)** survives every in-silico test: it is robust
to measurement error, its recipe is compliance-clean, its reachable vector is biologically
concentrated, and it recovers signed immunology. Its one honest vulnerability is a **small
coordinated bias** (~0.03 SE-units) — which is precisely the failure mode the randomized CRISPRi
design is built to exclude. The **aging axes** are flagged fragile at Stim48hr and should be
reported as such.

> **Bottom line for a referee:** the reachability verdict is a design-based causal claim whose
> assumptions are enumerated, tested where testable, and quantified (in SE-units) where not. The
> Farkas certificate makes infeasibility a proof; this notebook makes feasibility a *calibrated*
> claim.""")

nb = new_notebook(cells=C, metadata={"kernelspec":
     {"display_name":"Python 3","language":"python","name":"python3"},
     "language_info":{"name":"python"}})
out = os.path.join(REPO, "notebooks", "09_causal_validation_dossier.ipynb")
nbf.write(nb, out)
print(f"wrote {out} with {len(C)} cells")
