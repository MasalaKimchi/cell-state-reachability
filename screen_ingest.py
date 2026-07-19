"""Ingest a real combinatorial perturbation screen into CombiCone inputs.

This adapter turns a raw, *already-preprocessed* single-cell screen — an AnnData
``.h5ad`` or a pair of CSV/TSV tables — into the exact arrays the
:mod:`combicone` triage/certify API consumes:

    atoms          (n_atoms, n_genes)   single-perturbation effect vectors (the cone)
    atom_names     (n_atoms,)           gene name per atom
    combos         list[ComboRecord]    each measured k-way combination, with its
                                        effect vector, per-gene split-half noise SD,
                                        and the constituent single names it resolves to
    genes          (n_genes,)           shared gene axis

It is a thin, honest layer over :func:`effect_dictionary.build_effect_dictionary`:
it computes pooled ``mean(condition) - mean(control)`` effects, parses the
condition grammar to separate singles from combinations, matches every
combination to its constituent measured singles, and (optionally) estimates the
per-gene measurement noise from a random split-half of the cells — the standard
error the noise-aware certificate needs.

What it does NOT do (by design; these are upstream scientific choices that must
match the intended estimand): normalize raw counts, build replicate-aware
pseudobulks, correct donor/batch effects, or decide that cells are independent
replicates. Feed it data already on an additive (log-normalized) scale.

Condition grammar
-----------------
Each cell carries one condition label. The default grammar matches the public
Perturb-seq convention used by the Norman and CaRPool-seq screens:

  * a **control** label (default ``"ctrl"``; a combinatorial control such as
    ``"NT_NT"`` is also supported via ``control_label=``);
  * a **single** perturbation is one gene, optionally paired with the control
    token on one arm (``"AHR+ctrl"``, ``"NT_AHR"``);
  * a **combination** is ``order`` genes joined by ``separator`` (``"AHR+KLF1"``,
    ``"A+B+C"``), none of them the control token.

Both ``"GENE+ctrl"`` / ``"ctrl+GENE"`` and a bare ``"GENE"`` are read as the
single for that gene; a combination arm equal to the control token collapses the
label to the remaining single(s).

Public API
----------
``ingest_screen``          top-level: file/AnnData/arrays -> ScreenSubstrate
``ScreenSubstrate``        dataclass holding atoms, combos, genes, and provenance
``parse_condition``        the label parser (control / single / combination)
``triage_ready``           (atoms, atom_names) for :func:`combicone.triage_combinations`
``certify_ready``          per-combo (cone_atoms, measured_combo, noise_sd, names)
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

try:  # SciPy sparse is optional at import time; AnnData X is often sparse.
    from scipy import sparse
except Exception:  # pragma: no cover - scipy is a hard dep of the project
    sparse = None  # type: ignore

import effect_dictionary as _ed


__all__ = [
    "ScreenSubstrate",
    "ComboRecord",
    "ingest_screen",
    "parse_condition",
    "split_half_noise_sd",
]


# --------------------------------------------------------------------------- #
# Label grammar
# --------------------------------------------------------------------------- #
def parse_condition(
    label: str,
    *,
    control_label: str = "ctrl",
    separator: str = "+",
) -> tuple[str, ...]:
    """Return the tuple of perturbed genes encoded by one condition label.

    The control token is stripped from every arm, so ``"AHR+ctrl"`` -> ``("AHR",)``
    and ``"AHR+KLF1"`` -> ``("AHR", "KLF1")``. The pure control returns ``()``.
    Genes are returned in the order they appear (not sorted), so a caller can see
    the raw arm order; downstream matching is order-insensitive.

    Raises
    ------
    ValueError
        If ``label`` is empty/whitespace or an arm is empty (e.g. ``"AHR+"``).
    """
    if not isinstance(label, str) or not label.strip():
        raise ValueError("condition label must be a non-empty string")
    raw_arms = label.split(separator)
    arms: list[str] = []
    for arm in raw_arms:
        arm = arm.strip()
        if arm == "" and len(raw_arms) > 1:
            raise ValueError(f"empty arm in condition label {label!r}")
        arms.append(arm)
    genes = tuple(a for a in arms if a and a != control_label)
    return genes


def _classify(
    conditions: np.ndarray,
    *,
    control_label: str,
    separator: str,
) -> tuple[dict[str, str], dict[tuple[str, ...], str]]:
    """Map gene -> single-condition-label and gene-tuple -> combo-condition-label.

    Singles are keyed by their one gene; combinations by the sorted tuple of their
    genes, so ``A+B`` and ``B+A`` collapse to one canonical key.
    """
    singles: dict[str, str] = {}
    combos: dict[tuple[str, ...], str] = {}
    for label in conditions:
        genes = parse_condition(
            label, control_label=control_label, separator=separator
        )
        if len(genes) == 0:
            continue  # control
        if len(genes) == 1:
            singles.setdefault(genes[0], label)
        else:
            combos.setdefault(tuple(sorted(genes)), label)
    return singles, combos


# --------------------------------------------------------------------------- #
# Split-half noise
# --------------------------------------------------------------------------- #
def split_half_noise_sd(
    matrix: Any,
    condition_labels: np.ndarray,
    target_labels: Sequence[str],
    *,
    control_label: str,
    seed: int = 7,
    min_cells_per_half: int = 6,
) -> dict[str, np.ndarray]:
    """Per-gene measurement-noise SD for each target condition, via a cell split-half.

    For each condition we randomly partition its cells into two halves, form the
    effect (half-mean minus the matching control half-mean) on each half, and take
    ``|e1 - e2| / 2`` as the per-gene standard-error estimate — the same construction
    the emergence certificate's noise-injection null expects. Conditions with fewer
    than ``2 * min_cells_per_half`` cells return ``NaN`` (too few to split honestly).
    """
    rng = np.random.default_rng(seed)
    is_sparse = sparse is not None and sparse.issparse(matrix)
    ctrl_rows = np.where(condition_labels == control_label)[0]
    if ctrl_rows.size < 2 * min_cells_per_half:
        raise ValueError(
            "too few control cells to build a split-half noise estimate"
        )
    rng.shuffle(ctrl_rows)
    ctrl_a = ctrl_rows[: ctrl_rows.size // 2]
    ctrl_b = ctrl_rows[ctrl_rows.size // 2 :]

    def _mean(rows: np.ndarray) -> np.ndarray:
        if is_sparse:
            return np.asarray(matrix[rows].mean(axis=0)).ravel()
        return matrix[rows].mean(axis=0)

    ctrl_mean_a = _mean(ctrl_a)
    ctrl_mean_b = _mean(ctrl_b)

    out: dict[str, np.ndarray] = {}
    n_genes = matrix.shape[1]
    for label in target_labels:
        rows = np.where(condition_labels == label)[0]
        if rows.size < 2 * min_cells_per_half:
            out[label] = np.full(n_genes, np.nan)
            continue
        rows = rows.copy()
        rng.shuffle(rows)
        half_a = rows[: rows.size // 2]
        half_b = rows[rows.size // 2 :]
        e1 = _mean(half_a) - ctrl_mean_a
        e2 = _mean(half_b) - ctrl_mean_b
        out[label] = np.abs(e1 - e2) / 2.0
    return out


# --------------------------------------------------------------------------- #
# Substrate container
# --------------------------------------------------------------------------- #
@dataclass
class ComboRecord:
    """One measured combination and everything needed to certify it."""

    name: str
    genes: tuple[str, ...]
    effect: np.ndarray
    noise_sd: np.ndarray | None
    constituent_singles: tuple[str, ...]
    has_all_singles: bool


@dataclass
class ScreenSubstrate:
    """CombiCone-ready arrays parsed from a real screen.

    Attributes
    ----------
    atoms : (n_atoms, n_genes) float array
        Single-perturbation effect vectors — the cone for triage and certification.
    atom_names : (n_atoms,) str array
        Gene name per atom row.
    genes : (n_genes,) str array
        Shared gene axis.
    combos : list[ComboRecord]
        Every measured combination (order >= 2), each resolved to its constituent
        singles and (if noise was requested) its per-gene split-half SD.
    control_label, separator : str
        The grammar used to parse condition labels.
    provenance : dict
        Free-form counts and settings for the run (cell type, n cells, etc.).
    """

    atoms: np.ndarray
    atom_names: np.ndarray
    genes: np.ndarray
    combos: list[ComboRecord]
    control_label: str = "ctrl"
    separator: str = "+"
    provenance: dict[str, Any] = field(default_factory=dict)

    # -- summaries ---------------------------------------------------------- #
    @property
    def n_atoms(self) -> int:
        return int(self.atoms.shape[0])

    @property
    def n_combos(self) -> int:
        return len(self.combos)

    @property
    def n_genes(self) -> int:
        return int(self.genes.shape[0])

    def coverage(self) -> float:
        """Fraction of combinations whose every constituent single was measured."""
        if not self.combos:
            return float("nan")
        return float(np.mean([c.has_all_singles for c in self.combos]))

    def triage_ready(self) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(atoms, atom_names)`` for :func:`combicone.triage_combinations`."""
        return self.atoms, self.atom_names

    def certify_ready(self, combo_name: str) -> dict[str, Any]:
        """Return the kwargs to certify one combination against the singles cone.

        The returned dict is ready to splat into
        ``combicone.certify_emergence(cone_atoms=..., measured_combo=...,
        noise_sd=...)``.
        """
        rec = next((c for c in self.combos if c.name == combo_name), None)
        if rec is None:
            raise KeyError(f"no measured combination named {combo_name!r}")
        return {
            "cone_atoms": self.atoms,
            "measured_combo": rec.effect,
            "noise_sd": rec.noise_sd,
            "genes": self.genes,
            "constituent_singles": rec.constituent_singles,
        }

    def summary(self) -> dict[str, Any]:
        d = {
            "n_atoms": self.n_atoms,
            "n_combos": self.n_combos,
            "n_genes": self.n_genes,
            "singles_coverage": self.coverage(),
            "control_label": self.control_label,
            "separator": self.separator,
        }
        d.update(self.provenance)
        return d


# --------------------------------------------------------------------------- #
# Readers
# --------------------------------------------------------------------------- #
def _read_anndata(
    path_or_adata: Any,
    *,
    condition_key: str,
    gene_name_key: str | None,
    layer: str | None,
) -> tuple[Any, np.ndarray, np.ndarray]:
    """Return ``(matrix, condition_labels, gene_names)`` from an AnnData or .h5ad."""
    try:
        import anndata as ad
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "reading .h5ad requires the 'anndata' package; install it or pass "
            "arrays directly via ingest_screen(expression=..., conditions=...)"
        ) from exc
    adata = (
        path_or_adata
        if not isinstance(path_or_adata, (str, Path))
        else ad.read_h5ad(str(path_or_adata))
    )
    if condition_key not in adata.obs:
        raise ValueError(
            f"condition_key {condition_key!r} not in adata.obs "
            f"(available: {list(adata.obs.columns)[:12]})"
        )
    conditions = np.asarray(adata.obs[condition_key].astype(str).values)
    matrix = adata.layers[layer] if layer is not None else adata.X
    if gene_name_key is not None and gene_name_key in adata.var:
        genes = np.asarray(adata.var[gene_name_key].astype(str).values)
    else:
        genes = np.asarray(adata.var_names.astype(str).values)
    return matrix, conditions, genes


def _read_csv_pair(
    expression_csv: str | Path,
    conditions_csv: str | Path,
    *,
    condition_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a dense expression CSV (cells x genes, header = gene names) plus a
    per-cell conditions CSV (one column named ``condition_key``)."""
    expression_csv = Path(expression_csv)
    with expression_csv.open() as fh:
        reader = csv.reader(fh)
        header = next(reader)
    # First column may be a cell-id index; detect by name.
    skip_first = header[0].strip().lower() in {"", "cell", "cell_id", "index", "barcode"}
    gene_start = 1 if skip_first else 0
    genes = np.asarray([h.strip() for h in header[gene_start:]], dtype=str)
    matrix = np.loadtxt(
        expression_csv,
        delimiter=",",
        skiprows=1,
        usecols=range(gene_start, len(header)),
    )
    if matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    conds: list[str] = []
    with Path(conditions_csv).open() as fh:
        dreader = csv.DictReader(fh)
        if condition_key not in (dreader.fieldnames or []):
            raise ValueError(
                f"condition_key {condition_key!r} not a column in "
                f"{conditions_csv} (found: {dreader.fieldnames})"
            )
        for row in dreader:
            conds.append(str(row[condition_key]))
    conditions = np.asarray(conds, dtype=str)
    if conditions.size != matrix.shape[0]:
        raise ValueError(
            f"cell count mismatch: expression has {matrix.shape[0]} rows, "
            f"conditions has {conditions.size}"
        )
    return matrix, conditions, genes


# --------------------------------------------------------------------------- #
# Top-level entry point
# --------------------------------------------------------------------------- #
def ingest_screen(
    source: Any = None,
    *,
    # AnnData / arrays
    expression: Any = None,
    conditions: Sequence[str] | np.ndarray | None = None,
    gene_names: Sequence[str] | np.ndarray | None = None,
    condition_key: str = "condition",
    gene_name_key: str | None = None,
    layer: str | None = None,
    # CSV pair
    conditions_csv: str | Path | None = None,
    # grammar
    control_label: str = "ctrl",
    separator: str = "+",
    arm_handling: str = "merge",
    # noise
    compute_noise: bool = True,
    noise_seed: int = 7,
    min_cells_per_half: int = 6,
    dtype: np.dtype | type = np.float64,
) -> ScreenSubstrate:
    """Ingest a combinatorial screen into a :class:`ScreenSubstrate`.

    Dispatch by input:
      * ``ingest_screen("screen.h5ad", condition_key="guide_identity")`` — AnnData file
      * ``ingest_screen(adata, condition_key=...)`` — an in-memory AnnData
      * ``ingest_screen("expr.csv", conditions_csv="meta.csv", condition_key=...)`` — CSVs
      * ``ingest_screen(expression=X, conditions=labels, gene_names=g)`` — raw arrays

    Parameters mirror :func:`effect_dictionary.build_effect_dictionary` for the
    pooled-mean core, plus the combinatorial grammar (``control_label``,
    ``separator``) and the split-half noise controls (``compute_noise``,
    ``noise_seed``, ``min_cells_per_half``).

    Returns
    -------
    ScreenSubstrate
    """
    # ---- resolve the (matrix, conditions, genes) triple -------------------- #
    if expression is not None:
        if conditions is None:
            raise ValueError("pass conditions alongside expression")
        matrix = expression
        cond = np.asarray([str(c) for c in conditions])
        genes = (
            np.asarray([str(g) for g in gene_names])
            if gene_names is not None
            else np.asarray([f"gene_{i}" for i in range(matrix.shape[1])], dtype=str)
        )
    elif conditions_csv is not None:
        if source is None:
            raise ValueError("pass the expression CSV path as the first argument")
        matrix, cond, genes = _read_csv_pair(
            source, conditions_csv, condition_key=condition_key
        )
    elif source is not None:
        matrix, cond, genes = _read_anndata(
            source,
            condition_key=condition_key,
            gene_name_key=gene_name_key,
            layer=layer,
        )
    else:
        raise ValueError(
            "provide one of: source .h5ad/AnnData, source CSV + conditions_csv, "
            "or expression + conditions arrays"
        )

    # ---- canonicalize every cell's label to its parsed gene-tuple ---------- #
    # This is what makes pooling correct: a screen may spell the same biological
    # condition several ways (``BAK1+ctrl`` and ``ctrl+BAK1``; ``A+B`` and
    # ``B+A``). We collapse each to one canonical label BEFORE pooling. Canonical
    # forms:
    #   control      -> control_label
    #   single gene  -> "<gene>" (bare)
    #   combination  -> separator.join(sorted(genes))
    #
    # ``arm_handling`` decides how a SINGLE spelled with the control on one arm is
    # pooled when both arms are present (common in Perturb-seq; 47/105 Norman
    # genes have both). It never affects controls or true combinations:
    #   "merge"         (default) pool cells from every arm — most cells, the most
    #                   statistically efficient single-gene estimate.
    #   "control_left"  keep only ``ctrl+GENE`` cells (control on the left arm).
    #                   Reproduces the shipped Norman reference substrate.
    #   "control_right" keep only ``GENE+ctrl`` cells (control on the right arm).
    if arm_handling not in {"merge", "control_left", "control_right"}:
        raise ValueError(
            "arm_handling must be 'merge', 'control_left', or 'control_right', "
            f"got {arm_handling!r}"
        )

    def _arm_of(raw: str) -> str:
        """Classify a single-gene raw label's control arm: 'left', 'right', or 'bare'."""
        arms = [a.strip() for a in raw.split(separator)]
        if len(arms) == 1:
            return "bare"
        if arms[0] == control_label:
            return "left"
        if arms[-1] == control_label:
            return "right"
        return "bare"  # unusual: control token elsewhere; treat as no-arm-info

    # First pass: for each single gene, discover which arms exist. The arm
    # PREFERENCE only bites when the preferred arm is actually present for that
    # gene; a gene measured on only one arm is always kept (no choice to make).
    gene_arms: dict[str, set[str]] = {}
    for raw in np.unique(cond):
        g = parse_condition(raw, control_label=control_label, separator=separator)
        if len(g) == 1:
            gene_arms.setdefault(g[0], set()).add(_arm_of(str(raw)))

    def _preferred_single_arm(raw: str, gene: str) -> bool:
        """Whether this raw single-label should contribute under arm_handling."""
        if arm_handling == "merge":
            return True
        arm = _arm_of(str(raw))
        if arm == "bare":
            return True
        want = "left" if arm_handling == "control_left" else "right"
        arms_present = gene_arms.get(gene, set())
        # Fall back to whatever arm exists if the preferred one is absent.
        if want not in arms_present:
            return True
        return arm == want

    canon_cond = np.empty(cond.shape, dtype=object)
    canon_to_genes: dict[str, tuple[str, ...]] = {}
    raw_variants: dict[str, set[str]] = {}
    _DROP = "\x00drop"  # cells excluded by arm_handling (not control, not a pert)
    for i, raw in enumerate(cond):
        g = parse_condition(raw, control_label=control_label, separator=separator)
        if len(g) == 0:
            canon = control_label
        elif len(g) == 1:
            if not _preferred_single_arm(str(raw), g[0]):
                canon_cond[i] = _DROP
                continue
            canon = g[0]
        else:
            g = tuple(sorted(g))
            canon = separator.join(g)
        canon_cond[i] = canon
        if canon != control_label:
            canon_to_genes[canon] = g
            raw_variants.setdefault(canon, set()).add(str(raw))
    canon_cond = canon_cond.astype(str)

    # Drop cells excluded by arm_handling before pooling.
    if np.any(canon_cond == _DROP):
        keep = canon_cond != _DROP
        matrix = matrix[keep]
        canon_cond = canon_cond[keep]

    # ---- pooled effects via the audited low-level adapter ------------------ #
    eff = _ed.build_effect_dictionary(
        matrix,
        canon_cond,
        control_label=control_label,
        gene_names=genes,
        dtype=dtype,
    )
    E, perts, genes = eff["E"], eff["perts"], eff["genes"]
    pert_index = {p: i for i, p in enumerate(perts)}

    # ---- split canonical perturbations into singles vs combinations -------- #
    singles_map = {c: c for c, g in canon_to_genes.items() if len(g) == 1 and c in pert_index}
    combos_map = {canon_to_genes[c]: c for c, g in canon_to_genes.items() if len(g) >= 2 and c in pert_index}

    # atoms = the single-gene effects, ordered by gene name for determinism
    atom_genes = sorted(singles_map)
    atom_rows = [pert_index[g] for g in atom_genes]
    if not atom_rows:
        raise ValueError("no single-gene atoms parsed from conditions")
    atoms = E[atom_rows]
    atom_names = np.asarray(atom_genes, dtype=str)
    atom_set = set(atom_genes)

    # ---- optional split-half noise for the combos we will certify ---------- #
    noise_targets = list(combos_map.values())
    noise_by_label: dict[str, np.ndarray] = {}
    if compute_noise and noise_targets:
        noise_by_label = split_half_noise_sd(
            matrix,
            canon_cond,
            noise_targets,
            control_label=control_label,
            seed=noise_seed,
            min_cells_per_half=min_cells_per_half,
        )

    # ---- assemble combo records ------------------------------------------- #
    combos: list[ComboRecord] = []
    for genes_tuple in sorted(combos_map):
        label = combos_map[genes_tuple]
        constituents = tuple(g for g in genes_tuple if g in atom_set)
        combos.append(
            ComboRecord(
                name=label,
                genes=genes_tuple,
                effect=E[pert_index[label]],
                noise_sd=noise_by_label.get(label),
                constituent_singles=constituents,
                has_all_singles=(len(constituents) == len(genes_tuple)),
            )
        )

    provenance = {
        "n_cells": int(matrix.shape[0]),
        "n_conditions": int(perts.size),
        "n_singles_parsed": len(singles_map),
        "n_combos_parsed": len(combos_map),
        "n_multi_variant_labels": int(sum(1 for v in raw_variants.values() if len(v) > 1)),
        "arm_handling": arm_handling,
        "noise_computed": bool(compute_noise and noise_targets),
        "noise_seed": noise_seed,
    }
    return ScreenSubstrate(
        atoms=np.asarray(atoms, dtype=float),
        atom_names=atom_names,
        genes=np.asarray(genes, dtype=str),
        combos=combos,
        control_label=control_label,
        separator=separator,
        provenance=provenance,
    )
