"""Build the single repo-facing Cell-State Reachability storyboard.

The figure deliberately separates what was measured, what the mathematical model
does, what the retrospective case study showed, and what remains to be tested. All
displayed numbers are loaded from the canonical findings ledger.

Outputs: fig_at_a_glance.png and fig_at_a_glance.pdf
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("SOURCE_DATE_EPOCH", "0")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results"
findings = json.loads((RESULTS / "findings.json").read_text(encoding="utf-8"))
if findings.get("schema_version") != "1.0.0":
    raise ValueError("unsupported findings schema")
headline = findings["headline"]
observed = float(headline["historical_fixed_split_cosine"])
n_null = int(headline["diagnostic_target_shuffles"])
split_values = [float(value) for value in headline["split_values"]]
split_mean = float(headline["held_out_cosine_mean"])
split_sd = float(headline["held_out_cosine_sd"])
finding_by_id = {entry["id"]: entry for entry in findings["updated_findings"]}
target_scope = finding_by_id["target_observation_scope"]["values"]
target_total = int(target_scope["target_genes_total"])
target_measured = int(target_scope["target_genes_in_screen"])
top_50_measured = int(target_scope["top_50_surviving"])
if len(split_values) != 12 or n_null != 60:
    raise ValueError("unexpected registered split/shuffle counts")


# restrained, print-safe palette
NAVY = "#17324D"
INK = "#31506A"
MUTE = "#526675"
TEAL = "#1F6A5C"
TEAL_LIGHT = "#E7F2EF"
BLUE = "#3D6697"
BLUE_LIGHT = "#EAF0F8"
GOLD = "#805600"
GOLD_LIGHT = "#F6F0E2"
GRAY = "#A9B1B7"
GRAY_LIGHT = "#F0F2F3"
BORDER = "#CCD8DF"
WHITE = "#FFFFFF"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "axes.unicode_minus": False,
    }
)

fig = plt.figure(figsize=(14, 7), dpi=200, facecolor=WHITE)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")


def box(x, y, w, h, *, fc=WHITE, ec=BORDER, lw=1.3, radius=0.8, z=1):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        mutation_aspect=0.5,
        zorder=z,
    )
    ax.add_patch(patch)
    return patch


def label(
    x,
    y,
    text,
    *,
    size=10,
    color=INK,
    weight="normal",
    ha="left",
    va="center",
    style="normal",
    z=5,
    linespacing=1.2,
):
    artist = ax.text(
        x,
        y,
        text,
        fontsize=size,
        color=color,
        fontweight=weight,
        ha=ha,
        va=va,
        style=style,
        zorder=z,
    )
    artist.set_linespacing(linespacing)
    return artist


def step_header(x, number, title, subtitle):
    label(x + 1.7, 78.2, str(number), size=11, color=WHITE, weight="bold", ha="center")
    ax.add_patch(plt.Circle((x + 1.7, 78.2), 1.45, facecolor=NAVY, edgecolor=NAVY, zorder=3))
    label(x + 4.0, 79.4, title, size=11.5, color=NAVY, weight="bold")
    label(x + 4.0, 76.3, subtitle, size=8.6, color=MUTE)


# Title
label(4, 94.0, "Cell-State Reachability", size=25, color=NAVY, weight="bold")
label(
    4,
    89.1,
    "A screen-relative test of transcriptomic direction — not a state-conversion claim",
    size=13.2,
    color=TEAL,
    weight="bold",
)
label(
    96,
    93.5,
    "MEASURED  →  MODELLED  →  TEST NEXT",
    size=9.0,
    color=MUTE,
    weight="bold",
    ha="right",
)
ax.add_line(Line2D([4, 96], [85.2, 85.2], color=BORDER, lw=1.2))


# Four peer panels
xs = [4.0, 27.5, 51.0, 74.5]
pw, py, ph = 21.5, 23.7, 59.0
for x in xs:
    box(x, py, pw, ph, fc=WHITE, ec=BORDER, lw=1.3, radius=0.8)


# 1. Inputs
x = xs[0]
step_header(x, 1, "Inputs", "two separately sourced objects")

box(x + 1.7, 60.0, pw - 3.4, 11.0, fc=BLUE_LIGHT, ec="#C9D7E8", lw=1.0, radius=0.6)
label(x + 3.0, 68.0, "TARGET DIRECTION  d", size=8.4, color=BLUE, weight="bold")
label(x + 3.0, 64.5, "Constructed external Th1-vs-Th2\ndirection (two source contrasts)", size=9.2, color=NAVY)

box(x + 1.7, 44.7, pw - 3.4, 11.0, fc=TEAL_LIGHT, ec="#BDD9D2", lw=1.0, radius=0.6)
label(x + 3.0, 52.7, "EFFECT DICTIONARY  E", size=8.4, color=TEAL, weight="bold")
label(x + 3.0, 49.2, "Screen-derived CRISPRi DE z-scores\nin resting primary CD4 cells", size=9.2, color=NAVY)

label(
    x + 1.8,
    36.8,
    f"Only {target_measured:,} / {target_total:,} target genes\n"
    f"intersect the screen ({top_50_measured} / top 50)\n"
    "Dictionary is donor-collapsed; not\nmeasured in polarized Th2 cells",
    size=8.7,
    color=MUTE,
    style="italic",
    va="top",
    linespacing=1.35,
)


# 2. Projection geometry
x = xs[1]
step_header(x, 2, "Projection", "what the model computes")

origin = np.array([x + 4.2, 34.0])
cone = np.array(
    [
        origin,
        [x + 18.8, 46.5],
        [x + 15.2, 65.0],
    ]
)
ax.add_patch(Polygon(cone, closed=True, facecolor=BLUE_LIGHT, edgecolor="none", zorder=1))

for end in [(x + 17.5, 45.4), (x + 15.9, 53.0), (x + 14.4, 62.0)]:
    ax.add_patch(
        FancyArrowPatch(
            origin,
            end,
            arrowstyle="-|>",
            mutation_scale=10,
            color=BLUE,
            lw=1.6,
            zorder=3,
        )
    )

target = np.array([x + 10.0, 69.5])
projection = np.array([x + 13.3, 57.0])
ax.add_patch(
    FancyArrowPatch(origin, target, arrowstyle="-|>", mutation_scale=12, color=NAVY, lw=2.3, zorder=4)
)
ax.add_patch(
    FancyArrowPatch(origin, projection, arrowstyle="-|>", mutation_scale=12, color=TEAL, lw=2.8, zorder=5)
)
ax.add_patch(
    FancyArrowPatch(projection, target, arrowstyle="-|>", mutation_scale=10, color=GOLD, lw=2.0, zorder=5)
)
label(target[0] - 0.2, target[1] + 2.1, "target", size=8.2, color=NAVY, weight="bold", ha="center")
label(x + 13.5, 54.5, "projected\nmodel fit", size=8.1, color=TEAL, weight="bold")
label(x + 11.6, 64.1, "unmatched\nresidual", size=8.1, color=GOLD, weight="bold")
label(
    x + 1.7,
    28.0,
    "Non-negative linear combinations\nin screen z-score space\nResidual is model-relative",
    size=8.6,
    color=MUTE,
    va="bottom",
    linespacing=1.3,
)


# 3. Retrospective challenge
x = xs[2]
step_header(x, 3, "Challenge", "frozen retrospective diagnostics")

label(x + 1.8, 66.5, f"{split_mean:.3f} ± {split_sd:.3f}", size=23, color=TEAL, weight="bold")
label(x + 1.9, 61.7, "mean ± SD, 12 fixed random-gene splits", size=8.8, color=NAVY, weight="bold")
label(x + 1.9, 58.5, "external Th1-like direction", size=8.5, color=MUTE)

# compact split-stability strip
sx0, sx1, sy = x + 2.0, x + pw - 2.0, 48.5
lo, hi = 0.425, 0.465
ax.add_line(Line2D([sx0, sx1], [sy, sy], color=BORDER, lw=1.5, zorder=2))
for value in split_values:
    px = sx0 + (value - lo) / (hi - lo) * (sx1 - sx0)
    ax.add_line(Line2D([px, px], [sy - 1.5, sy + 1.5], color=BLUE, lw=1.2, zorder=3))
mean_x = sx0 + (split_mean - lo) / (hi - lo) * (sx1 - sx0)
ax.add_line(Line2D([mean_x, mean_x], [sy - 2.4, sy + 2.4], color=NAVY, lw=2.5, zorder=4))
label(x + 1.9, 44.3, f"historical fixed split: {observed:.3f}", size=8.7, color=NAVY, weight="bold")

box(x + 1.7, 31.1, pw - 3.4, 8.5, fc=GRAY_LIGHT, ec="#D8DDE0", lw=1.0, radius=0.6)
label(x + 3.0, 37.5, f"{n_null}/{n_null} diagnostic target shuffles\nbelow observed", size=8.5, color=NAVY, weight="bold", linespacing=1.1)
label(x + 3.0, 32.8, f"plus-one empirical p = 1/{n_null + 1}", size=8.4, color=MUTE)

label(
    x + 1.8,
    29.0,
    "Retrospective directional evidence;\nnot donor or functional validation",
    size=8.4,
    color=MUTE,
    style="italic",
    va="top",
)


# 4. Decision boundary
x = xs[3]
step_header(x, 4, "Decide", "what this result changes")

rows = [
    (TEAL_LIGHT, TEAL, "REPLICATION NEEDED", "Donor-resolved, polarized\nTh2 starting-state study"),
    (GOLD_LIGHT, GOLD, "MEASURE NEXT", "Matched modalities, combinations\n+ functional readouts"),
    (GRAY_LIGHT, MUTE, "NOT ESTABLISHED", "State conversion, rescue,\nor target validation"),
]
for i, (fill, accent, heading, body) in enumerate(rows):
    ry = 59.5 - i * 14.4
    box(x + 1.7, ry, pw - 3.4, 11.2, fc=fill, ec=accent, lw=1.0, radius=0.6)
    ax.add_patch(Rectangle((x + 1.7, ry), 0.8, 11.2, facecolor=accent, edgecolor="none", zorder=3))
    label(x + 3.2, ry + 8.2, heading, size=8.0, color=accent, weight="bold")
    label(x + 3.2, ry + 4.0, body, size=9.2, color=NAVY, weight="bold", linespacing=1.25)


# Footer: the durable boundary
box(4, 7.0, 92, 11.8, fc=GRAY_LIGHT, ec=BORDER, lw=1.1, radius=0.7)
label(6.0, 15.1, "CLAIM BOUNDARY", size=8.8, color=GOLD, weight="bold")
label(
    6.0,
    11.1,
    "Screen-relative directional support  •  measured subset only  •  model-relative geometry ≠ biological efficacy",
    size=10.7,
    color=NAVY,
    weight="bold",
)
label(
    96,
    3.4,
    "Values loaded from results/findings.json; source hashes recorded in results/manifest.json",
    size=8.2,
    color=MUTE,
    ha="right",
)

fig.savefig(HERE / "fig_at_a_glance.png", dpi=200, facecolor=WHITE, bbox_inches="tight", pad_inches=0.12)
fig.savefig(HERE / "fig_at_a_glance.pdf", facecolor=WHITE, bbox_inches="tight", pad_inches=0.12)
print("wrote fig_at_a_glance.png / .pdf")
