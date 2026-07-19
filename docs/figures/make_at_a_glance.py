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
if findings.get("schema_version") != "3.2.0":
    raise ValueError("unsupported findings schema")
alignment = findings["source_bound_alignment"]
split_values = [float(value) for value in alignment["split_values"]]
split_mean = float(alignment["held_out_cosine_mean"])
split_sd = float(alignment["held_out_cosine_sd"])
target_scope = findings["target_lineage"]
target_total = int(target_scope["union_genes"])
target_measured = int(target_scope["shared_screen_genes"])
target_analyzed = int(target_scope["registered_genes"])
transfer = findings["cross_source_directional_transfer"]
donor = findings["donor_pair_transfer"]["run_balanced"]
arce = findings["arce_external_validation"]["spearman"]
schmidt = findings["schmidt_external_validation"]
schmidt_top = schmidt["source_selected_top_k_transfer"]
guide = findings["guide_position_transfer"]
guide_same_source = guide["same_source_guide_position_transfer"]
arrayed = findings["zhu_arrayed_followup"]
goudy = findings["goudy_cross_experiment_stress"]
goudy_primary = goudy["primary_transcriptome"]
goudy_reliability = goudy["reliability"]
if len(split_values) != 12:
    raise ValueError("unexpected source-bound split count")
if goudy["status"]["geometric_model"] != "FAILS_DECLARED_GEOMETRIC_MODEL":
    raise ValueError("unexpected Goudy geometric status")
if (
    schmidt["tier"] != "STRESS"
    or schmidt["primary_contract"]["common_complete_genes"] != 18568
    or schmidt["primary_contract"]["top_k"] != 200
):
    raise ValueError("unexpected Schmidt stress contract")
if (
    guide["tier"] != "STRESS"
    or guide["benchmark_role"]
    != "negative released positional-guide reciprocal-transfer stress"
    or guide["guide_identity_status"]
    != "POSITIONAL_MODALITIES_ONLY_IDENTITIES_NOT_RELEASED"
    or guide["positional_modalities"] != 2
    or guide["common_rest_atoms"] != 8323
    or guide["fit_count"] != 12
    or guide["challenge_rows"] != 24
    or guide_same_source["challenge_rows"] != 12
    or guide_same_source["unique_fits"] != 12
):
    raise ValueError("unexpected released positional-guide stress contract")
guide_within_cosine = float(guide_same_source["median_within_training_guide_cone_cosine"])
guide_reciprocal_cosine = float(
    guide_same_source["median_reciprocal_held_guide_cone_cosine"]
)
guide_reciprocal_delta = float(guide_same_source["median_reciprocal_minus_within_cone_cosine"])
guide_positive_fraction = float(
    guide_same_source[
        "fraction_reciprocal_held_guide_cosine_improvement_positive_over_training_best_single"
    ]
)
if not np.isfinite(
    [
        guide_within_cosine,
        guide_reciprocal_cosine,
        guide_reciprocal_delta,
        guide_positive_fraction,
    ]
).all():
    raise ValueError("non-finite released positional-guide stress values")
guide_positive_count_float = guide_positive_fraction * guide_same_source["challenge_rows"]
guide_positive_count = round(guide_positive_count_float)
if not np.isclose(guide_positive_count_float, guide_positive_count, atol=1e-12, rtol=0.0):
    raise ValueError("released positional-guide positive fraction is not an exact row count")
guide_reciprocal_label = f"{guide_reciprocal_cosine:.3f}"
guide_delta_label = f"{guide_reciprocal_delta:+.3f}"


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
    "Point-estimate cone geometry — not a calibrated state-conversion claim",
    size=13.2,
    color=TEAL,
    weight="bold",
)
label(
    96,
    93.5,
    "MEASURED  →  MODELLED  →  VALIDATE",
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
step_header(x, 1, "Inputs", "source objects reused from Zhu workflow")

box(x + 1.7, 60.0, pw - 3.4, 11.0, fc=BLUE_LIGHT, ec="#C9D7E8", lw=1.0, radius=0.6)
label(x + 3.0, 68.0, "TARGET DIRECTION  d", size=8.4, color=BLUE, weight="bold")
label(x + 3.0, 64.5, "Source-study-reused population contrast\n(not independent or a trajectory)", size=9.2, color=NAVY)

box(x + 1.7, 44.7, pw - 3.4, 11.0, fc=TEAL_LIGHT, ec="#BDD9D2", lw=1.0, radius=0.6)
label(x + 3.0, 52.7, "EFFECT DICTIONARY  E", size=8.4, color=TEAL, weight="bold")
label(x + 3.0, 49.2, "Screen-derived CRISPRi DE z-scores\nin post-expansion Rest CD4 cells", size=9.2, color=NAVY)

label(
    x + 1.8,
    36.8,
    f"Target union: {target_total:,}; shared screen: {target_measured:,}\n"
    f"Final merged analysis: {target_analyzed:,} genes\n"
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
    "Non-negative linear combinations\nin declared effect space\nResidual is model-relative",
    size=8.6,
    color=MUTE,
    va="bottom",
    linespacing=1.3,
)


# 3. Retrospective challenge
x = xs[2]
step_header(x, 3, "Challenge", "source-bound + measured transfer tests")

label(x + 1.8, 66.5, f"{split_mean:.3f} ± {split_sd:.3f}", size=23, color=TEAL, weight="bold")
label(x + 1.9, 61.7, "mean ± SD, 12 hash-frozen gene splits", size=8.8, color=NAVY, weight="bold")
label(x + 1.9, 58.5, "selected population contrast; not trajectory", size=8.5, color=MUTE)

# compact split-stability strip
sx0, sx1, sy = x + 2.0, x + pw - 2.0, 48.5
lo, hi = min(split_values) - 0.005, max(split_values) + 0.005
ax.add_line(Line2D([sx0, sx1], [sy, sy], color=BORDER, lw=1.5, zorder=2))
for value in split_values:
    px = sx0 + (value - lo) / (hi - lo) * (sx1 - sx0)
    ax.add_line(Line2D([px, px], [sy - 1.5, sy + 1.5], color=BLUE, lw=1.2, zorder=3))
mean_x = sx0 + (split_mean - lo) / (hi - lo) * (sx1 - sx0)
ax.add_line(Line2D([mean_x, mean_x], [sy - 2.4, sy + 2.4], color=NAVY, lw=2.5, zorder=4))
label(x + 1.9, 44.3, f"range: {min(split_values):.3f}–{max(split_values):.3f}", size=8.7, color=NAVY, weight="bold")

box(x + 1.7, 30.6, pw - 3.4, 9.0, fc=GRAY_LIGHT, ec="#D8DDE0", lw=1.0, radius=0.6)
label(x + 3.0, 37.2, "Frozen transfer challenges", size=8.5, color=NAVY, weight="bold")
label(
    x + 3.0,
    33.4,
    f"Source gain +{transfer['ota_to_hollbacker']['mean_cosine_improvement_over_test_selected_better_baseline']:.3f} / "
    f"+{transfer['hollbacker_to_ota']['mean_cosine_improvement_over_test_selected_better_baseline']:.3f}\n"
    f"Donor-pair gain +{donor['median_cosine_improvement_over_training_best_single']:.3f}; "
    f"nRMSE {donor['median_cone_normalized_rmse']:.3f} vs "
    f"{donor['median_training_best_single_normalized_rmse']:.3f}\n"
    f"Released positional: cone {guide_within_cosine:.3f} within → "
    f"{guide_reciprocal_label} reciprocal\n"
    f"Δ {guide_delta_label}; {guide_positive_count}/{guide_same_source['challenge_rows']} beat best single",
    size=6.5,
    color=MUTE,
    linespacing=1.0,
)

label(
    x + 1.8,
    29.6,
    f"Arce ranks {arce['Resting_Teff']:.3f}/{arce['Stimulated_Teff']:.3f}/{arce['Resting_Treg']:.3f} | "
    f"Zhu 9/9; cosine {arrayed['profile_replication']['panel_centered_median_cosine']:.3f}\n"
    f"Schmidt ρ(top-200): donor {schmidt_top['same_screen_held_donor']['median_signed_spearman']:.2f} | "
    f"donor+M/L {schmidt_top['donor_plus_modality_library_same_context']['median_signed_spearman']:.2f} | "
    f"donor+ctx {schmidt_top['donor_plus_cross_context_cytokine_cell_type_same_modality']['median_signed_spearman']:.2f}\n"
    "RNA→flow 0.72–0.85; descriptive, not donor/guide/state validation",
    size=6.5,
    color=MUTE,
    style="italic",
    va="top",
    linespacing=1.0,
)


# 4. Decision boundary
x = xs[3]
step_header(x, 4, "Limits", "what remains unresolved")

rows = [
    (TEAL_LIGHT, TEAL, "REPLICATION GATE", "Donor-resolved, polarized\nTh2 starting-state study"),
    (GOLD_LIGHT, GOLD, "CALIBRATION GATE", "Structured holdouts, uncertainty\n+ matched baselines"),
    (GRAY_LIGHT, MUTE, "NOT ESTABLISHED", "State conversion, rescue,\nor target validation"),
]
for i, (fill, accent, heading, body) in enumerate(rows):
    ry = 59.5 - i * 14.4
    box(x + 1.7, ry, pw - 3.4, 11.2, fc=fill, ec=accent, lw=1.0, radius=0.6)
    ax.add_patch(Rectangle((x + 1.7, ry), 0.8, 11.2, facecolor=accent, edgecolor="none", zorder=3))
    label(x + 3.2, ry + 8.2, heading, size=8.0, color=accent, weight="bold")
    label(x + 3.2, ry + 4.0, body, size=9.2, color=NAVY, weight="bold", linespacing=1.25)

label(
    x + 1.8,
    29.0,
    "Goudy cross-experiment stress: cone "
    f"{goudy_primary['component_cone_median_cosine']:.3f};\n"
    "triple donor cosine "
    f"{goudy_reliability['triple_pairwise_donor_median_cosine']:.3f} — model fails",
    size=7.1,
    color=GOLD,
    weight="bold",
    va="top",
    linespacing=1.25,
)


# Footer: the durable boundary
box(4, 7.0, 92, 11.8, fc=GRAY_LIGHT, ec=BORDER, lw=1.1, radius=0.7)
label(6.0, 15.1, "CLAIM BOUNDARY", size=8.8, color=GOLD, weight="bold")
label(
    6.0,
    11.1,
    "Point-estimate directional geometry  •  correlated holdout diagnostic  •  model cone ≠ biological reachability",
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
