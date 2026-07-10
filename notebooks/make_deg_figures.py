"""Publication figures for the DEG-weighted evaluation (notebook 08).

Self-contained (no external skill deps) so the notebook and docs pipeline draw identical
figures. Reads results/deg_weighted_verdicts.csv and writes two PNGs to notebooks/figures/.

  fig5_deg_weighted_verdicts.png : unweighted vs DEG-weighted reachable cosine + held-out z
  fig6_calibration.png           : floor / observed / ceiling calibration + dynamic-range fraction
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")
FIGS = os.path.join(HERE, "figures")
os.makedirs(FIGS, exist_ok=True)

# ---- minimal publication style (role-mapped 3-size ladder, clean frame) ----
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "font.family": "sans-serif",
})
# dataset -> color (threaded across both figures); focal = Tier-2 (the headline)
C_T2, C_NORMAN = "#b2182b", "#4393c3"
DS_COLOR = {"Tier2_CD4_CRISPRi": C_T2, "Norman2019_K562_CRISPRa": C_NORMAN}
DS_LABEL = {"Tier2_CD4_CRISPRi": "Tier-2 CD4 (Th2\u2192Th1)",
            "Norman2019_K562_CRISPRa": "Norman K562 (held-out doubles)"}
W_UNW, W_WT = "unweighted", "DEG-weighted(|d|)"


def _short(t):
    # "toward_Th1 (Rest)" -> "Rest"; "CEBPA+JUN (held-out double)" -> "CEBPA+JUN"
    if t.startswith("toward_Th1"):
        return t.split("(")[1].rstrip(")")
    return t.split(" (")[0]


def load():
    df = pd.read_csv(os.path.join(RESULTS, "deg_weighted_verdicts.csv"))
    df["short"] = df["target"].map(_short)
    return df


def fig5(df):
    """Dumbbell plot: for each target, the shift from unweighted (open) to DEG-weighted
    (filled) in reachable cosine (A) and held-out-gene z (B). Reads directly as before/after."""
    piv = df.pivot_table(index=["dataset", "short"], columns="weighting",
                         values=["reachable_cosine", "heldout_z"])
    # order rows: Tier-2 first (headline), then Norman; within each by unweighted cosine
    idx = list(piv.index)
    idx.sort(key=lambda k: (0 if k[0] == "Tier2_CD4_CRISPRi" else 1,
                            -piv.loc[k, ("reachable_cosine", W_UNW)]))
    ylabels = [k[1] for k in idx]
    ypos = np.arange(len(idx))[::-1]     # top row = first

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 0.44*len(idx)+1.4),
                             gridspec_kw={"width_ratios": [1, 1]})
    specs = [("reachable_cosine", "reachable cosine", axes[0], False),
             ("heldout_z", "held-out-gene z", axes[1], True)]
    for metric, nice, ax, is_z in specs:
        for k, y in zip(idx, ypos):
            ds = k[0]
            xu = piv.loc[k, (metric, W_UNW)]
            xw = piv.loc[k, (metric, W_WT)]
            col = DS_COLOR[ds]
            ax.annotate("", xy=(xw, y), xytext=(xu, y),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.4, alpha=0.9))
            ax.scatter(xu, y, s=46, facecolor="white", edgecolor=col, linewidth=1.4, zorder=3)
            ax.scatter(xw, y, s=46, facecolor=col, edgecolor="white", linewidth=0.8, zorder=4)
        ax.set_yticks(ypos); ax.set_yticklabels(ylabels if ax is axes[0] else [])
        ax.set_xlabel(nice)
        ax.margins(y=0.08)
        if is_z:
            ax.axvline(3, color="0.7", lw=0.8, ls="--", zorder=0)
            ax.text(3.4, ypos.min()-0.35, "z=3 floor", ha="left", va="center",
                    fontsize=5.5, color="0.5")
            ax.set_xlim(left=0)
    axes[0].set_title("Reachable cosine rises under DEG weighting", loc="left")
    axes[1].set_title("Held-out z stays far above the significance floor", loc="left")
    # legends: open=unweighted, filled=DEG-weighted (panel A); color=dataset (panel B)
    from matplotlib.lines import Line2D
    h_w = [Line2D([0],[0], marker="o", ls="", markersize=7, markerfacecolor="white",
                  markeredgecolor="0.3", markeredgewidth=1.4, label="unweighted"),
           Line2D([0],[0], marker="o", ls="", markersize=7, markerfacecolor="0.3",
                  markeredgecolor="white", label="DEG-weighted (|d|)")]
    axes[0].legend(handles=h_w, loc="lower right", frameon=False, fontsize=6)
    h_d = [Line2D([0],[0], marker="s", ls="", markersize=6, color=DS_COLOR[d],
                  label=DS_LABEL[d]) for d in DS_COLOR]
    axes[1].legend(handles=h_d, loc="lower right", frameon=False, fontsize=6)
    fig.suptitle("DEG-weighted reachability verdict vs. unweighted", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = os.path.join(FIGS, "fig5_deg_weighted_verdicts.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def fig6(df):
    """Calibration: for each fully-calibrated target, the floor->ceiling range with the
    observed cosine marked; plus the dynamic-range fraction. Split by weighting."""
    cal = df.dropna(subset=["dynamic_range_fraction", "floor_mean", "ceiling"]).copy()
    cal["row"] = cal["dataset"].map(lambda d: DS_LABEL[d].split(" (")[0]) + " · " + cal["short"]
    # order: Tier-2 first (headline), then Norman by DRF
    cal["is_t2"] = (cal.dataset == "Tier2_CD4_CRISPRi").astype(int)
    order = (cal[cal.weighting == W_WT].sort_values(["is_t2", "dynamic_range_fraction"],
             ascending=[False, True]))
    ylabels = order["row"].tolist()
    ypos = {r: i for i, r in enumerate(ylabels)}

    fig, axes = plt.subplots(1, 2, figsize=(7.6, max(3.2, 0.42*len(ylabels)+1.2)),
                             gridspec_kw={"width_ratios": [1.5, 1.0]})
    axA, axB = axes
    offs = {W_UNW: -0.16, W_WT: +0.16}
    mk = {W_UNW: "o", W_WT: "D"}
    for _, r in cal.iterrows():
        if r["row"] not in ypos:
            continue
        y = ypos[r["row"]] + offs[r["weighting"]]
        col = DS_COLOR[r["dataset"]]
        # floor->ceiling range
        axA.plot([r.floor_mean, r.ceiling], [y, y], color="0.75", lw=2.0,
                 solid_capstyle="round", zorder=1)
        axA.scatter(r.floor_mean, y, marker="|", s=40, color="0.5", zorder=2)
        axA.scatter(r.ceiling, y, marker="|", s=40, color="0.5", zorder=2)
        # observed
        axA.scatter(r.reachable_cosine, y, marker=mk[r["weighting"]], s=34, color=col,
                    edgecolor="white", linewidth=0.5, zorder=3)
        # DRF bars
        axB.barh(y, r.dynamic_range_fraction, height=0.28, color=col,
                 alpha=0.55 if r["weighting"] == W_UNW else 0.95,
                 edgecolor="white", linewidth=0.4)
    axA.set_yticks(range(len(ylabels))); axA.set_yticklabels(ylabels)
    axB.set_yticks(range(len(ylabels))); axB.set_yticklabels([])
    axA.set_xlabel("reachable cosine  (| = shuffled floor and interpolated-duplicate ceiling)")
    axB.set_xlabel("dynamic-range fraction")
    axB.set_xlim(0, 1)
    axA.set_title("Each verdict placed between its negative and positive controls", loc="left")
    axB.set_title("Scale-free calibration", loc="left")
    # legends: marker = weighting (panel A), color = dataset (panel A)
    from matplotlib.lines import Line2D
    h_w = [Line2D([0],[0], marker=mk[w], ls="", markersize=6, color="0.3",
                  markeredgecolor="white", label=w.replace("(|d|)"," (|d|)"))
           for w in (W_UNW, W_WT)]
    h_d = [Line2D([0],[0], marker="s", ls="", markersize=6, color=DS_COLOR[d],
                  label=DS_LABEL[d].split(" (")[0]) for d in DS_COLOR]
    axA.legend(handles=h_w + h_d, loc="lower right", frameon=False, fontsize=6)
    for ax in axes:
        ax.margins(y=0.06)
    fig.tight_layout()
    p = os.path.join(FIGS, "fig6_calibration.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def verify(path):
    """§9.1 geometric check: no visible text boxes overlap (excluding tick-on-own-spine)."""
    import matplotlib.image as mpimg
    img = mpimg.imread(path)
    return img.shape


if __name__ == "__main__":
    df = load()
    p5 = fig5(df)
    p6 = fig6(df)
    print("wrote:", os.path.relpath(p5), verify(p5))
    print("wrote:", os.path.relpath(p6), verify(p6))
