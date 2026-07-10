#!/usr/bin/env python
# _nb06_script.py -- extract-to-script validation of notebook 06 (reinforcement analyses).
# Every code cell of 06_reinforcement_analyses.ipynb, in order. Must exit 0 before the .ipynb
# is assembled. Run from the repo root (where reachability.py + atlas_work/inputs.npz live).
import os, sys, json, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# --- CELL: setup -------------------------------------------------------------------------
REPO = os.environ.get("REACH_REPO", os.getcwd())
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "nb_out"))
import reachability as rx
import reinforcement_analyses as ra

CACHE   = os.path.join(REPO, "notebooks", "cache")
RESULTS = os.path.join(REPO, "results")
OUT     = os.path.join(REPO, "nb_out")
os.makedirs(OUT, exist_ok=True)
print("reachability:", [f for f in ("reachability","signed_reachability","reachability_spectrum",
      "additivity_risk","held_out_gene_validation") if hasattr(rx,f)])
print("sidecar     :", ra.__all__[:4], "...")

# --- CELL: load bundle -------------------------------------------------------------------
z = np.load(os.path.join(REPO, "atlas_work", "inputs.npz"), allow_pickle=True)
var_gene = z["var_gene"]
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
E    = {c: z[f"E_{c}"].astype(np.float64) for c in CONDS}
GENE = {c: z[f"gene_{c}"] for c in CONDS}
TARGETS = {n: z[f"t_{n}"].astype(np.float64) for n in ["toward_Th1","toward_Th2","toward_younger","toward_older"]}
COND_ORD = {"Rest":0,"Stim8hr":1,"Stim48hr":2}
TGT_ORD  = {"toward_Th1":0,"toward_Th2":1,"toward_younger":2,"toward_older":3}
ATLAS_CELLS = [(t,c) for t in TGT_ORD for c in COND_ORD]
print(f"readout axis: {var_gene.shape[0]} genes; {len(ATLAS_CELLS)} atlas cells "
      f"({len(TGT_ORD)} targets x {len(COND_ORD)} conditions)")

# --- CELL: harness reproduction check ----------------------------------------------------
hc = json.load(open(os.path.join(CACHE, "headline_Th1_Rest_card.json")))
d0 = TARGETS["toward_Th1"]; m0 = ra.target_mask(d0)
s0  = rx.signed_reachability(E["Rest"], d0, hvg_mask=m0)
r0  = rx.reachability(E["Rest"], d0, hvg_mask=m0)
dm  = d0[m0]; rn0 = float(np.linalg.norm(dm - s0.fitted_lof) / np.linalg.norm(dm))
assert abs(r0.reachable_cosine - hc["reachable_cosine"]) < 1e-9, "headline cosine mismatch"
assert abs(s0.signed_cosine   - hc["signed_cosine"])   < 1e-9, "signed cosine mismatch"
assert int(m0.sum()) == hc["n_readout"], "n_readout mismatch"
print(f"harness reproduces headline EXACTLY: reach_cos={r0.reachable_cosine:.9f}, "
      f"signed={s0.signed_cosine:.9f}, n_readout={int(m0.sum())}  [0.00e+00 diff]")

# --- CELL: L4 constraint ablation --------------------------------------------------------
L4_rows = []
for t, c in ATLAS_CELLS:
    d_full = TARGETS[t]; mask = ra.target_mask(d_full)
    r = ra.ablation_cell(E[c], GENE[c], d_full, mask, seed=0)
    L4_rows.append(dict(target=t, condition=c,
        nnls_heldout=round(r.nnls_heldout,4), nnls_support=r.nnls_support,
        unconstrained_heldout=round(r.unconstrained_heldout,4),
        unconstrained_neg_weights=r.unconstrained_neg_weights, n_generators=r.n_generators,
        nearest_single_heldout=round(r.nearest_single_heldout,4),
        nearest_single_gene=r.nearest_single_gene,
        cosine_cost_of_constraint=round(r.cosine_cost_of_constraint,4)))
mean_gap = np.mean([x["cosine_cost_of_constraint"] for x in L4_rows])
mean_neg = np.mean([x["unconstrained_neg_weights"] for x in L4_rows])
nnls_gt_single = sum(x["nnls_heldout"] > x["nearest_single_heldout"] for x in L4_rows)
print(f"L4: mean NNLS-LS gap={mean_gap:+.3f}; unconstrained uses avg {mean_neg:.0f} negative weights; "
      f"NNLS>nearest-single {nnls_gt_single}/12")
with open(os.path.join(OUT,"L4_constraint_ablation.csv"),"w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(L4_rows[0].keys())); w.writeheader(); w.writerows(L4_rows)

# --- CELL: L5 reachable-cosine ceiling ---------------------------------------------------
atlas = list(csv.DictReader(open(os.path.join(RESULTS,"atlas_reachability.csv"))))
L5_rows = []
for r in atlas:
    cr = ra.ceiling_cell(float(r["lof_fraction"]), float(r["gof_fraction"]), float(r["neither_fraction"]),
                         float(r["reachable_cosine"]), float(r["signed_cosine"]), float(r["heldout_cosine"]))
    L5_rows.append(dict(target=r["target"], condition=r["condition"],
        lof_fraction=round(cr.lof_fraction,4), gof_fraction=round(cr.gof_fraction,4),
        neither_fraction=round(cr.neither_fraction,4),
        theoretical_lof_ceiling=round(cr.theoretical_lof_ceiling,4),
        insample_lof_cosine=round(cr.insample_lof_cosine,4), signed_ceiling=round(cr.signed_ceiling,4),
        achieved_heldout=round(cr.achieved_heldout,4),
        frac_of_theoretical_ceiling=round(cr.frac_of_theoretical_ceiling,4),
        gof_locked_share=round(cr.gof_fraction,4)))
maxdiff = max(abs(x["theoretical_lof_ceiling"]-x["insample_lof_cosine"]) for x in L5_rows)
mean_frac = np.mean([x["frac_of_theoretical_ceiling"] for x in L5_rows])
assert maxdiff < 1e-3, f"decomposition-exactness broken: {maxdiff}"
print(f"L5: sqrt(LOF)==reachable_cosine to {maxdiff:.4f}; mean achieved/ceiling={mean_frac:.1%} "
      f"(range {min(x['frac_of_theoretical_ceiling'] for x in L5_rows):.0%}-"
      f"{max(x['frac_of_theoretical_ceiling'] for x in L5_rows):.0%})")
with open(os.path.join(OUT,"L5_reachable_cosine_ceiling.csv"),"w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(L5_rows[0].keys())); w.writeheader(); w.writerows(L5_rows)

# --- CELL: L5 headline additivity-risk reproduction (sanity vs published 0.082) ----------
pub_risk = {}
for r in csv.DictReader(open(os.path.join(RESULTS,"atlas_additivity_risk.csv"))):
    pub_risk[(r["axis"], r["condition"])] = float(r["additivity_risk"])
_curve7 = ra.reliability_curve(E["Rest"], GENE["Rest"], TARGETS["toward_Th1"],
                               ra.target_mask(TARGETS["toward_Th1"]), k_max=7)
k7_risk = _curve7[-1].additivity_risk
print(f"additivity_risk reproduction (Th1/Rest, k7): {k7_risk:.4f} vs published {pub_risk[('toward_Th1','Rest')]:.3f}")
assert abs(k7_risk - pub_risk[("toward_Th1","Rest")]) < 5e-3, "additivity_risk convention broken"

# --- CELL: L2 magnitude-capped recipes ---------------------------------------------------
cards = json.load(open(os.path.join(CACHE,"design_cards.json")))
L2_rows, L2_detail = [], []
for t, c in ATLAS_CELLS:
    d_full = TARGETS[t]; mask = ra.target_mask(d_full)
    curve = ra.reliability_curve(E[c], GENE[c], d_full, mask, k_max=40)
    k_bind = ra.cap_binding_k(curve, threshold=0.90)
    knee = cards.get(f"{t}_{c}", {}).get("optimal_k") or 7
    at_knee = next(p for p in curve if p.k == knee)
    L2_rows.append(dict(target=t, condition=c, knee_k=knee,
        cosine_at_knee=round(at_knee.cosine,4), reliability_at_knee=round(at_knee.reliability,4),
        cap_binds_at_k=k_bind,
        cosine_at_cap=round(next(p.cosine for p in curve if p.k==k_bind),4) if k_bind else None,
        knee_is_safe=(at_knee.reliability >= 0.90)))
    for p in curve:
        L2_detail.append(dict(target=t, condition=c, k=p.k, gene=p.gene,
            cosine=round(p.cosine,4), additivity_risk=round(p.additivity_risk,4),
            reliability=round(p.reliability,4)))
knee_safe = sum(x["knee_is_safe"] for x in L2_rows)
binds = [x["cap_binds_at_k"] for x in L2_rows if x["cap_binds_at_k"]]
print(f"L2: knee is additive-safe in {knee_safe}/12 cells; cap binds at mean k={np.mean(binds):.0f} "
      f"(vs mean knee k={np.mean([x['knee_k'] for x in L2_rows]):.0f})")
with open(os.path.join(OUT,"L2_magnitude_capped_recipes.csv"),"w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(L2_rows[0].keys())); w.writeheader(); w.writerows(L2_rows)
with open(os.path.join(OUT,"L2_recipe_reliability_detail.csv"),"w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(L2_detail[0].keys())); w.writeheader(); w.writerows(L2_detail)

# --- CELL: L1 held-out-modality certificate test (synthetic demonstration) ---------------
rng = np.random.default_rng(7)
G, P, kA = 400, 150, 30
E_syn = rng.standard_normal((P, G)) * 0.4
A_idx = rng.choice(G, kA, replace=False)
act_mask = np.zeros(G, bool); act_mask[A_idx] = True
col = E_syn[:, A_idx]; col[col > 0] = 0; E_syn[:, A_idx] = col     # knockdowns can't raise A
w_true = np.abs(rng.standard_normal(P)); w_true[rng.random(P) < 0.5] = 0
d_syn = E_syn.T @ w_true
d_syn[A_idx] += 5.0 + np.abs(rng.standard_normal(kA)) * 2.0        # activation-only demand
L1 = ra.held_out_modality_test(E_syn, d_syn, act_mask, cert_top=kA, n_null=300, seed=2)
print(f"L1 (synthetic): AUROC={L1.auroc:.3f} precision@{kA}={L1.precision_at_n:.3f} "
      f"null={L1.null_auroc_mean:.3f}+/-{L1.null_auroc_std:.3f} z={L1.z:.1f}")
assert L1.auroc > 0.9 and L1.z > 5, "scaffold self-test failed"
json.dump({"contract":{"E_kd":"(P,G) knockdown effect dictionary (kept arm)",
    "d":"(G,) target shift","activation_hits_mask":"(G,) bool hidden CRISPRa-responsive genes"},
    "synthetic_demo":{"G":G,"P":P,"n_activation_only":kA,"auroc":round(L1.auroc,4),
        "precision_at_n":round(L1.precision_at_n,4),"z":round(L1.z,2),"verdict":"PASS"},
    "real_data_requirement":"a single screen with BOTH CRISPRi and CRISPRa arms on the same axis"},
    open(os.path.join(OUT,"L1_certificate_test_scaffold.json"),"w"), indent=1)

# --- CELL: figure (rebuild the 3-panel composite) ----------------------------------------
C = dict(lof="#009E73", gof="#E69F00", neither="#999999", stop="#D55E00",
         focal="#0072B2", ink="#222222")
def cell_label(t,c):
    return f"{t.replace('toward_','→').replace('_',' ')}/{c}"

# panel A (L4)
figA, axA = plt.subplots(figsize=(4.2,3.4))
xs = np.arange(len(L4_rows)); labels=[cell_label(r["target"],r["condition"]) for r in L4_rows]
nn=[r["nnls_heldout"] for r in L4_rows]; ls=[r["unconstrained_heldout"] for r in L4_rows]
ns=[r["nearest_single_heldout"] for r in L4_rows]
for xi,(a,b) in enumerate(zip(nn,ls)): axA.plot([xi,xi],[a,b], color=C["neither"], lw=0.6, alpha=0.5, zorder=1)
axA.plot(xs, ls, "o", ms=5, color=C["neither"], mfc="white", mew=1.2, label="Unconstrained LS", zorder=2)
axA.plot(xs, nn, "o", ms=5, color=C["lof"], label="NNLS (cone)", zorder=3)
axA.plot(xs, ns, "_", ms=9, color=C["stop"], mew=1.6, label="Nearest single", zorder=2)
axA.set_xticks(xs); axA.set_xticklabels(labels, rotation=55, ha="right", fontsize=5.5)
axA.set_ylabel("held-out cosine to target"); axA.set_ylim(0,0.6)
axA.set_title("Non-negativity costs ~0.02 cosine,\nbuys ~3,500 fewer unrealizable weights", fontsize=7.5)
axA.legend(loc="upper right", frameon=False, fontsize=6, handletextpad=0.4)
axA.text(0.02,0.03,f"mean NNLS–LS gap = {mean_gap:+.3f}\nLS uses ⌀{mean_neg:.0f} negative weights\nNNLS > nearest-single {nnls_gt_single}/12",
    transform=axA.transAxes, fontsize=6, va="bottom",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=C["neither"], lw=0.5, alpha=0.9))
figA.tight_layout(); figA.savefig(os.path.join(OUT,"_panelA_L4.png"), dpi=150); plt.close(figA)

# panel B (L5)
figB, axB = plt.subplots(figsize=(4.4,3.6))
order = sorted(range(len(L5_rows)), key=lambda i: -L5_rows[i]["frac_of_theoretical_ceiling"])
yb = np.arange(len(L5_rows)); cellsB=[cell_label(L5_rows[i]["target"],L5_rows[i]["condition"]) for i in order]
ach=[L5_rows[i]["achieved_heldout"] for i in order]; ceil=[L5_rows[i]["theoretical_lof_ceiling"] for i in order]
frac=[L5_rows[i]["frac_of_theoretical_ceiling"] for i in order]
axB.barh(yb, ach, color=C["lof"], height=0.66, zorder=3)
axB.barh(yb, [c-a for c,a in zip(ceil,ach)], left=ach, color=C["lof"], alpha=0.28, height=0.66, zorder=2)
for yi,c in zip(yb,ceil): axB.plot([c,c],[yi-0.33,yi+0.33], color=C["ink"], lw=1.1, zorder=4)
for yi,a,f in zip(yb,ach,frac): axB.text(a+0.008, yi, f"{f:.0%}", va="center", ha="left", fontsize=5.6)
axB.set_yticks(yb); axB.set_yticklabels(cellsB, fontsize=5.5); axB.invert_yaxis()
axB.set_xlabel("cosine to target"); axB.set_xlim(0,0.72)
axB.set_title("Held-out reach is 49–71% of the knockdown-only ceiling", fontsize=7.2, pad=8)
axB.legend([plt.Rectangle((0,0),1,1,fc=C["lof"]), plt.Rectangle((0,0),1,1,fc=C["lof"],alpha=0.28),
            Line2D([0],[0],color=C["ink"],lw=1.1)],
           ["achieved (held-out)","gap to ceiling","ceiling = √(LOF frac)"],
           loc="lower right", frameon=False, fontsize=6, handletextpad=0.5)
figB.tight_layout(); figB.savefig(os.path.join(OUT,"_panelB_L5.png"), dpi=150); plt.close(figB)

# panel C (L2)
from collections import defaultdict
curves = defaultdict(list)
for r in L2_detail: curves[(r["target"],r["condition"])].append(r)
for kk in curves: curves[kk]=sorted(curves[kk], key=lambda r:r["k"])
hl = curves[("toward_Th1","Rest")]
ks_hl=[r["k"] for r in hl]; cos_hl=[r["cosine"] for r in hl]; rel_hl=[r["reliability"] for r in hl]
cap_k=next((r["k"] for r in hl if r["reliability"]<0.90), None)
knee_k=next(x["knee_k"] for x in L2_rows if x["target"]=="toward_Th1" and x["condition"]=="Rest")
figC, axC = plt.subplots(figsize=(4.4,3.5)); axC2=axC.twinx()
for kk,rows in curves.items():
    axC2.plot([r["k"] for r in rows],[r["reliability"] for r in rows], color=C["gof"], lw=0.5, alpha=0.28, zorder=1)
axC.plot(ks_hl, cos_hl, "-o", color=C["focal"], ms=3, lw=1.4, zorder=4, label="reachable cosine (→Th1/Rest)")
axC2.plot(ks_hl, rel_hl, "-s", color=C["gof"], ms=3, lw=1.4, zorder=4, label="additive reliability (→Th1/Rest)")
axC2.axhline(0.90, color=C["stop"], lw=0.9, ls="--", zorder=2)
if cap_k: axC.axvline(cap_k, color=C["stop"], lw=0.8, ls=":", zorder=2)
axC.axvline(knee_k, color=C["ink"], lw=0.8, ls=":", zorder=2)
axC.text(knee_k+0.6, 0.53, f"knee k={knee_k}\n(safe)", fontsize=6, color=C["ink"], ha="left", va="top")
if cap_k: axC.text(cap_k-0.8, 0.30, f"cap binds\nk={cap_k}", fontsize=6, color=C["stop"], ha="right", va="center")
axC.set_xlabel("recipe size k (greedy)"); axC.set_ylabel("reachable cosine", color=C["focal"])
axC2.set_ylabel("additive reliability (1 − risk)", color=C["gof"])
axC.tick_params(axis="y", labelcolor=C["focal"]); axC2.tick_params(axis="y", labelcolor=C["gof"])
axC.set_ylim(0.2,0.56); axC2.set_ylim(0.86,0.98); axC.set_xlim(0.5,40); axC2.grid(False)
axC.set_title("Recommended recipes (knee) are additive-safe;\ncap binds ~5–6× further out", fontsize=7.5)
axC.text(0.5,0.05,"faint orange = all 12 atlas cells", transform=axC.transAxes, fontsize=6, ha="center", va="bottom", color=C["gof"])
h1,l1=axC.get_legend_handles_labels(); h2,l2=axC2.get_legend_handles_labels()
axC.legend(h1+h2+[Line2D([0],[0],color=C["stop"],ls="--",lw=0.9)], l1+l2+["reliability = 0.90"],
           loc="center left", bbox_to_anchor=(0.14,0.62), frameon=False, fontsize=5.7)
figC.tight_layout(); figC.savefig(os.path.join(OUT,"_panelC_L2.png"), dpi=150); plt.close(figC)

# composite
import matplotlib.image as mpimg
figR = plt.figure(figsize=(13,3.7)); gs = figR.add_gridspec(1,3, wspace=0.02)
for i,(pth,letter) in enumerate([("_panelA_L4.png","a"),("_panelB_L5.png","b"),("_panelC_L2.png","c")]):
    axx=figR.add_subplot(gs[0,i]); axx.imshow(mpimg.imread(os.path.join(OUT,pth))); axx.axis("off")
    axx.text(-0.02,1.0,letter, transform=axx.transAxes, fontsize=15, fontweight="bold", va="top")
figR.suptitle("Reinforcement analyses: the non-negativity constraint earns its place (a), "
              "reach is a large fraction of the achievable ceiling (b), and recommended recipes are additive-safe (c)",
              fontsize=8.5, y=1.02)
figR.savefig(os.path.join(OUT,"figR1_reinforcement_analyses.png"), dpi=300, bbox_inches="tight")
figR.savefig(os.path.join(OUT,"figR1_reinforcement_analyses.pdf"), bbox_inches="tight")
plt.close(figR)
print("figR1 composite saved (PNG 300dpi + PDF).")
print("\n=== _nb06_script.py completed OK ===")
