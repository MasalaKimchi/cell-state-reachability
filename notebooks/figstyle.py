"""figstyle.py — the single source of truth for the manuscript figure look.

Every reader-facing figure in this repo (manuscript/, docs/, README) must import
from this module so the palette, role→colour mapping, typography, and axes style
are identical across the whole figure set. Before this module existed the
manuscript figures were polished by hand while the doc/notebook figures used
ad-hoc rcParams, and the two drifted (wrong palette, red-as-data, a stale
decomposition number). This file exists so that never happens again.

Usage
-----
    import sys; sys.path.insert(0, "<repo>/notebooks")
    from figstyle import apply_figure_style, OI, ROLE, set_frame, panel_letter, verify_bbox
    apply_figure_style()                 # call once, before plotting
    ax.bar(..., color=ROLE["LOF"])       # colour by SEMANTIC ROLE, never a raw hue

The locked palette is Okabe-Ito (colour-vision-deficiency safe). Do NOT hand-pick
hues in figure code — reference ``ROLE`` so the same entity is the same colour in
every panel of every figure.

ROLE semantics (the cross-reference the whole figure set relies on)
-------------------------------------------------------------------
    LOF       green      loss-of-function / reachable by knockdown (CRISPRi)
    GOF       orange     gain-of-function / activation-required (CRISPRa/ORF)
    neither   grey       orthogonal residual / neither direction
    observed  blue       the observed reachability curve / in-sample series
    null      grey       shuffled / random-perturbation null band
    control   purple     positive-control marker (e.g. GATA3)
    transfer  sky-blue   cross-dataset transfer series (Norman / Replogle)
    stop      vermillion ALARM ONLY — infeasible / off-target / anomaly. Never a
                         routine data series (reserve red so it stays meaningful).
"""

import os
import sys
import glob

import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Locked Okabe-Ito palette (colour-vision-deficiency safe). Reference by name.
# ---------------------------------------------------------------------------
OI = dict(
    black="#000000",
    orange="#E69F00",
    skyblue="#56B4E9",
    green="#009E73",
    yellow="#F0E442",
    blue="#0072B2",
    vermillion="#D55E00",
    purple="#CC79A7",
    grey="#999999",
)

# Semantic role → colour. Figure code colours by ROLE, never by raw hue, so the
# same entity is the same colour across every figure in the repo.
ROLE = dict(
    LOF=OI["green"],
    GOF=OI["orange"],
    neither=OI["grey"],
    stop=OI["vermillion"],
    observed=OI["blue"],
    null=OI["grey"],
    control=OI["purple"],
    transfer=OI["skyblue"],
)

# Confidence tiers (used by the design-card / reliability figures). Ordinal, so a
# single-hue ramp on the observed-blue family — NOT green/gold/red, which would
# collide with the LOF/GOF role colours.
CONFIDENCE = dict(
    high=OI["blue"],
    medium=OI["skyblue"],
    low="#BFD9EC",   # light tint of skyblue for the lowest tier
)

# Neutral grey for meta text (goodness cues, "y = x" guides, faded annotations).
META_GREY = "#888888"


# ---------------------------------------------------------------------------
# House style
# ---------------------------------------------------------------------------
def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    """Set global rcParams to the manuscript house style.

    Call once before plotting. ``sizes=(base, secondary, tick)`` is the
    role-mapped three-step font ladder (title/label = base, legend/annotation =
    secondary, tick labels = tick). ``frame='open'`` shows only left+bottom
    spines; 'boxed' shows all four; 'none' shows none.
    """
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    # Register any bundled fonts shipped with the conda env, if present.
    try:
        import matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,   # editable text in vector PDF
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def set_frame(ax, style="open"):
    """Show a subset of spines on one axis. 'open' = left+bottom only."""
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    """Bold panel letter in axes coordinates (top-left, outside the box)."""
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


def panel_letter_margin(ax, letter, extra_left=0.0, fontsize=11):
    """Bold panel letter in the FIGURE's left margin (clears a left-aligned title).

    Preferred over ``panel_letter`` when the panel has a long sentence-title that
    would otherwise collide with an in-axes letter.
    """
    fig = ax.figure
    bb = ax.get_position()
    x = max(0.004, bb.x0 - 0.050 - extra_left)
    y = min(0.990, bb.y1 + 0.006)
    fig.text(x, y, letter, fontsize=fontsize, fontweight="bold", ha="left", va="bottom")


def goodness_arrow(ax, text="higher = better", loc="upper left", axis="y", fontsize=None):
    """Small upright direction-of-goodness cue in a corner of the axes."""
    if fontsize is None:
        fontsize = plt.rcParams["legend.fontsize"]
    pos = {"upper left": (0.02, 0.98), "upper right": (0.98, 0.98),
           "lower left": (0.02, 0.02), "lower right": (0.98, 0.02)}[loc]
    ha = "left" if "left" in loc else "right"
    va = "top" if "upper" in loc else "bottom"
    arrow = "\u2191 " if axis == "y" else "\u2192 "
    ax.text(pos[0], pos[1], arrow + text, transform=ax.transAxes,
            fontsize=fontsize, color=META_GREY, ha=ha, va=va)


def verify_bbox(fig, min_area=6.0):
    """Geometric render-then-verify check (figure-style §9.1).

    Returns ``(overlaps, out_of_bounds)``: a list of colliding visible text/spine
    box pairs and a list of text that falls outside the figure. Empty lists mean
    the figure passes the geometric check. A tick label sitting on its own spine
    is not counted.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    tickmap = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
               for ax in fig.axes}

    def inter(a, b):
        dx = min(a.x1, b.x1) - max(a.x0, b.x0)
        dy = min(a.y1, b.y1) - max(a.y0, b.y0)
        return dx * dy if (dx > 0 and dy > 0) else 0

    real = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if inter(ba, bb) > min_area:
                real.append((a.get_text()[:20], b.get_text()[:20], round(inter(ba, bb))))
    for t, bt in texts:
        for s, bs in spines:
            if t in tickmap.get(s.axes, set()):
                continue
            if inter(bt, bs) > min_area:
                real.append((t.get_text()[:20], "<spine>", round(inter(bt, bs))))
    fb = fig.bbox
    oob = [t.get_text()[:20] for t, bt in texts
           if not (bt.x0 >= fb.x0 - 1 and bt.y0 >= fb.y0 - 1
                   and bt.x1 <= fb.x1 + 1 and bt.y1 <= fb.y1 + 1)]
    return real, oob


def save_fig(fig, stem, outdir="."):
    """Save a figure as both 300-dpi PNG and vector PDF at ``outdir/stem.{png,pdf}``."""
    png = os.path.join(outdir, f"{stem}.png")
    pdf = os.path.join(outdir, f"{stem}.pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    return png, pdf
