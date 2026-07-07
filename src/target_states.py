"""Build target-state direction vectors from the dataset's own signature tables.

A "target direction" d in gene space encodes the desired transcriptomic shift, e.g.
"become more Th1-like" or "look transcriptionally younger". We align these vectors to
the same gene ordering as the perturbation dictionary so the solver can operate.

SIGN CONVENTIONS (read before trusting any direction):
  * The polarization table's contrast is **Th2_vs_Th1**: a positive z-score means the
    gene is HIGHER in Th2. Therefore moving *toward Th1* is the NEGATIVE of the table,
    and moving *toward Th2* is the table as-is. (This was previously inverted.)
  * The aging table's contrast is **aged_vs_young**: positive means higher in aged
    cells, so moving *toward young* is the negative of the table.

ROBUSTNESS (polarization only):
  The polarization signature ships TWO independent source contrasts (Ota 2021 and
  Höllbacher 2021) stacked in one file. A naive dict(zip(gene, value)) silently keeps
  only the last contrast. We instead pivot to gene x contrast, then build the target
  from the **sign-concordant core** (genes both sources agree on), which is the
  cross-source robustness check the project claims. `concordance_report()` exposes the
  agreement so it can be reported, not just assumed.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

POLARIZATION_CSV = DATA_DIR / "Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv"
AGING_CSV = DATA_DIR / "CD4T_aging_signature_DE_results_full.suppl_table.csv"


def _gene_col(df: pd.DataFrame) -> str:
    """Return the gene-identifier column, tolerant of schema drift.

    The aging table has `gene_name`; the polarization table only has `variable`.
    """
    for c in ("gene_name", "variable", "gene", "symbol"):
        if c in df.columns:
            return c
    raise KeyError(
        f"No gene-identifier column found. Looked for gene_name/variable/gene/symbol; "
        f"got {list(df.columns)}"
    )


def _align_to_genes(values: dict[str, float], genes: np.ndarray) -> np.ndarray:
    """Project a {gene: value} map onto the solver's gene ordering (missing -> 0)."""
    return np.array([values.get(str(g), 0.0) for g in genes], dtype=np.float32)


# --------------------------------------------------------------------------- #
# Polarization (Th1 <-> Th2)
# --------------------------------------------------------------------------- #
def _polarization_wide(value_col: str = "zscore") -> pd.DataFrame:
    """gene x contrast matrix of the requested value column (one row per gene)."""
    sig = pd.read_csv(POLARIZATION_CSV)
    gcol = _gene_col(sig)
    # median collapses any within-contrast duplicate gene rows deterministically.
    wide = sig.pivot_table(
        index=gcol, columns="contrast", values=value_col, aggfunc="median"
    )
    return wide


@dataclass
class ConcordanceReport:
    n_shared: int
    n_concordant: int
    frac_concordant: float
    spearman: float
    pearson: float
    contrasts: list[str]


def concordance_report(value_col: str = "zscore") -> ConcordanceReport:
    """Quantify Ota-vs-Höllbacher agreement — the cross-source robustness check."""
    wide = _polarization_wide(value_col).dropna(how="any")
    cols = list(wide.columns)
    if len(cols) < 2:
        return ConcordanceReport(len(wide), len(wide), 1.0, float("nan"), float("nan"), cols)
    a, b = wide[cols[0]].to_numpy(), wide[cols[1]].to_numpy()
    concordant = np.sign(a) == np.sign(b)
    return ConcordanceReport(
        n_shared=int(len(wide)),
        n_concordant=int(concordant.sum()),
        frac_concordant=float(concordant.mean()),
        spearman=_corr(_rankdata(a), _rankdata(b)),  # numpy-only Spearman (no scipy)
        pearson=_corr(a, b),
        contrasts=cols,
    )


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / denom) if denom > 0 else float("nan")


def _rankdata(x: np.ndarray) -> np.ndarray:
    """Average-rank of x (ties averaged), numpy-only."""
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1)
    # average ties
    _, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    return (sums / counts)[inv]


def polarization_target(
    genes: np.ndarray,
    direction: str = "toward_Th1",
    mode: str = "concordant",
    value_col: str = "zscore",
) -> np.ndarray:
    """Th1<->Th2 polarization target direction, aligned to `genes`.

    direction : "toward_Th1" (default) = negative of the Th2_vs_Th1 signature;
                "toward_Th2" = signature as-is.
    mode      : "concordant" (default) keeps only genes both source contrasts agree on
                (sign-concordant core; discordant genes set to 0);
                "mean" averages across available contrasts;
                a literal contrast name uses that single source only.
    """
    wide = _polarization_wide(value_col)
    cols = list(wide.columns)

    if mode == "concordant":
        shared = wide.dropna(how="any")
        concordant = np.sign(shared[cols[0]]) == np.sign(shared[cols[1]])
        core = shared[concordant]
        series = core.mean(axis=1)
    elif mode == "mean":
        series = wide.mean(axis=1, skipna=True)
    elif mode in cols:
        series = wide[mode].dropna()
    else:
        raise ValueError(f"mode must be 'concordant', 'mean', or one of {cols}; got {mode!r}")

    values = {str(g): float(v) for g, v in series.items()}
    d = _align_to_genes(values, genes)  # in Th2_vs_Th1 orientation (positive = Th2)

    if direction == "toward_Th1":
        return -d
    if direction == "toward_Th2":
        return d
    raise ValueError(f"direction must be 'toward_Th1' or 'toward_Th2'; got {direction!r}")


# --------------------------------------------------------------------------- #
# Aging (aged -> young-like)
# --------------------------------------------------------------------------- #
def aging_target(genes: np.ndarray, direction: str = "toward_young") -> np.ndarray:
    """Reverse-aging direction from the CD4+ T-cell aging signature table.

    The aging signature encodes aged-vs-young change; to move *toward young* we invert.
    NOTE: the local donors are all young (ages 22-34), so this axis has NO in-sample
    aged reference and is exploratory only — see ROADMAP limitations.
    """
    sig = pd.read_csv(AGING_CSV)
    gcol = _gene_col(sig)
    values = {str(g): float(v) for g, v in zip(sig[gcol], sig["zscore"])}
    d_aging = _align_to_genes(values, genes)
    return -d_aging if direction == "toward_young" else d_aging


def normalize(d: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(d)
    return d / n if n > 0 else d
