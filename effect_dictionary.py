"""Build and validate portable perturbation effect dictionaries.

The reachability geometry expects one dense matrix with orientation
``(perturbations, genes)`` plus labels for both axes::

    E       floating array, shape (n_perturbations, n_genes)
    perts   string array, shape (n_perturbations,)
    genes   string array, shape (n_genes,)
    ncells  optional integer array, shape (n_perturbations,)

``build_effect_dictionary`` is a deliberately small adapter for already
preprocessed cell-by-gene data.  It computes a pooled arithmetic mean for each
condition and subtracts the pooled control mean.  It does *not* normalize raw
counts, construct replicate-aware pseudobulks, correct donors/batches, estimate
uncertainty, or establish that cells are independent replicates.  Those choices
belong upstream and must match the intended scientific estimand.

Dense NumPy arrays and SciPy sparse cell matrices are supported.  The saved
``.npz`` is intentionally limited to non-object NumPy arrays and is always read
with ``allow_pickle=False``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from io import BytesIO
from pathlib import Path
import tempfile
from typing import Any
from zipfile import ZIP_STORED, ZipFile, ZipInfo

import numpy as np
from scipy import sparse


REQUIRED_KEYS = ("E", "perts", "genes")
OPTIONAL_KEYS = ("ncells",)
ALLOWED_KEYS = frozenset((*REQUIRED_KEYS, *OPTIONAL_KEYS))


def _coerce_input_labels(
    values: Sequence[str] | np.ndarray,
    *,
    name: str,
    expected_size: int,
) -> np.ndarray:
    """Convert user labels to a safe Unicode vector without stringifying nulls."""

    raw = np.asarray(values)
    if raw.ndim != 1 or raw.size != expected_size:
        raise ValueError(f"{name} must be 1D with {expected_size} entries")
    items = raw.tolist()
    if any(not isinstance(item, (str, np.str_)) for item in items):
        raise ValueError(f"{name} must contain strings (missing/numeric labels are invalid)")
    labels = np.asarray(items, dtype=str)
    if any(not label.strip() or label != label.strip() for label in labels.tolist()):
        raise ValueError(f"{name} must contain non-empty, whitespace-trimmed strings")
    return labels


def _dense_row_mean(matrix: Any, row_mask: np.ndarray) -> np.ndarray:
    """Return a 1D mean while keeping a sparse input sparse until reduction."""

    if sparse.issparse(matrix):
        return np.asarray(matrix[row_mask].mean(axis=0), dtype=float).reshape(-1)
    return np.asarray(matrix[row_mask].mean(axis=0), dtype=float).reshape(-1)


def build_effect_dictionary(
    expression: np.ndarray | sparse.spmatrix,
    conditions: Sequence[str] | np.ndarray,
    *,
    control_label: str = "ctrl",
    gene_names: Sequence[str] | np.ndarray | None = None,
    dtype: np.dtype | type = np.float32,
) -> dict[str, np.ndarray]:
    """Build condition-minus-control effects from preprocessed cell data.

    Parameters
    ----------
    expression
        Dense or SciPy sparse ``(cells, genes)`` matrix on an additive scale
        for which arithmetic means and subtraction are scientifically valid.
    conditions
        One non-empty string condition label per cell.
    control_label
        Label whose pooled cell mean defines the baseline.
    gene_names
        One unique string per expression column.  Deterministic ``gene_N``
        placeholders are used only when labels truly are unavailable.
    dtype
        Floating dtype for the returned dense effect matrix.

    Notes
    -----
    Rows are sorted lexicographically by perturbation label and the control is
    excluded.  ``ncells`` records pooled cell counts for provenance only; it is
    not a replicate count and is never used as a statistical weight.
    """

    if not isinstance(control_label, str) or not control_label.strip():
        raise ValueError("control_label must be a non-empty string")
    if control_label != control_label.strip():
        raise ValueError("control_label must be whitespace-trimmed")

    output_dtype = np.dtype(dtype)
    if not np.issubdtype(output_dtype, np.floating):
        raise ValueError("dtype must be a floating dtype")

    if sparse.issparse(expression):
        if expression.ndim != 2:
            raise ValueError("expression must be 2D (cells, genes)")
        matrix = expression.astype(float, copy=False).tocsr()
        if not np.all(np.isfinite(matrix.data)):
            raise ValueError("expression contains non-finite values")
    else:
        matrix = np.asarray(expression, dtype=float)
        if matrix.ndim != 2:
            raise ValueError("expression must be 2D (cells, genes)")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("expression contains non-finite values")

    n_cells, n_genes = matrix.shape
    if n_cells == 0 or n_genes == 0:
        raise ValueError("expression must contain at least one cell and one gene")

    condition_labels = _coerce_input_labels(
        conditions, name="conditions", expected_size=n_cells
    )
    if gene_names is None:
        genes = np.asarray([f"gene_{index}" for index in range(n_genes)], dtype=str)
    else:
        genes = _coerce_input_labels(
            gene_names, name="gene_names", expected_size=n_genes
        )
    if np.unique(genes).size != genes.size:
        raise ValueError("gene_names must be unique")

    control_rows = condition_labels == control_label
    if not np.any(control_rows):
        raise ValueError(
            f"no control cells found for control_label={control_label!r}; "
            "a measured baseline is required"
        )
    control_mean = _dense_row_mean(matrix, control_rows)

    perturbations = np.asarray(
        sorted(label for label in np.unique(condition_labels) if label != control_label),
        dtype=str,
    )
    if perturbations.size == 0:
        raise ValueError("no non-control conditions found")

    effects = np.empty((perturbations.size, n_genes), dtype=float)
    cell_counts = np.empty(perturbations.size, dtype=np.int64)
    for index, perturbation in enumerate(perturbations):
        rows = condition_labels == perturbation
        cell_counts[index] = int(np.count_nonzero(rows))
        effects[index] = _dense_row_mean(matrix, rows) - control_mean

    result = {
        "E": effects.astype(output_dtype),
        "perts": perturbations,
        "genes": genes,
        "ncells": cell_counts,
    }
    problems = validate_effect_dictionary(result)
    if problems:
        raise ValueError(
            "constructed effect dictionary is not representable:\n  - "
            + "\n  - ".join(problems)
        )
    return result


def _label_problems(
    values: Any,
    *,
    name: str,
    expected_size: int,
) -> list[str]:
    """Validate labels without coercing object arrays into apparently safe strings."""

    labels = np.asarray(values)
    problems: list[str] = []
    if labels.ndim != 1 or labels.size != expected_size:
        return [f"{name} must be 1D with {expected_size} entries; got shape {labels.shape}"]
    if labels.dtype.kind not in {"U", "S"}:
        return [f"{name} must use a non-object string dtype; got {labels.dtype}"]
    as_text = labels.astype(str)
    items = as_text.tolist()
    if any(not label.strip() or label != label.strip() for label in items):
        problems.append(f"{name} must contain non-empty, whitespace-trimmed strings")
    if np.unique(as_text).size != as_text.size:
        problems.append(f"{name} contains duplicate labels")
    return problems


def validate_effect_dictionary(dictionary: Mapping[str, Any]) -> list[str]:
    """Return format problems; an empty list means the dictionary is well formed.

    This is a transport and alignment gate, not a scientific validity check.
    It cannot determine whether upstream normalization, controls, replicates,
    batches, or perturbation estimates are appropriate.
    """

    if not isinstance(dictionary, Mapping):
        return ["effect dictionary must be a mapping"]

    problems: list[str] = []
    missing = [key for key in REQUIRED_KEYS if key not in dictionary]
    problems.extend(f"missing required key: {key!r}" for key in missing)
    unexpected = sorted(set(dictionary) - ALLOWED_KEYS, key=str)
    problems.extend(f"unexpected key: {key!r}" for key in unexpected)
    if missing:
        return problems

    raw_effects = dictionary["E"]
    if sparse.issparse(raw_effects):
        problems.append("E must be a dense array in the portable NPZ format")
        return problems
    effects = np.asarray(raw_effects)
    if effects.ndim != 2:
        problems.append(
            f"E must be 2D (perturbations, genes); got shape {effects.shape}"
        )
        return problems
    n_perturbations, n_genes = effects.shape
    if n_perturbations == 0 or n_genes == 0:
        problems.append("E must contain at least one perturbation and one gene")
    if not np.issubdtype(effects.dtype, np.floating):
        problems.append(f"E must use a floating dtype; got {effects.dtype}")
    else:
        try:
            finite = bool(np.all(np.isfinite(effects)))
        except TypeError:
            finite = False
        if not finite:
            problems.append("E contains non-finite values (NaN or Inf)")

    problems.extend(
        _label_problems(
            dictionary["perts"], name="perts", expected_size=n_perturbations
        )
    )
    problems.extend(
        _label_problems(dictionary["genes"], name="genes", expected_size=n_genes)
    )

    if "ncells" in dictionary:
        counts = np.asarray(dictionary["ncells"])
        if counts.ndim != 1 or counts.size != n_perturbations:
            problems.append(
                f"ncells must be 1D with {n_perturbations} entries; got shape {counts.shape}"
            )
        elif (
            np.issubdtype(counts.dtype, np.bool_)
            or not np.issubdtype(counts.dtype, np.integer)
        ):
            problems.append(f"ncells must use an integer dtype; got {counts.dtype}")
        elif np.any(counts <= 0):
            problems.append("ncells entries must be positive")
    return problems


def _portable_payload(dictionary: Mapping[str, Any]) -> dict[str, np.ndarray]:
    """Return canonical little-endian, C-contiguous arrays in fixed field order."""

    effects = np.asarray(dictionary["E"])
    effects_dtype = effects.dtype.newbyteorder("<")
    perturbations = np.asarray(dictionary["perts"]).astype(str)
    genes = np.asarray(dictionary["genes"]).astype(str)
    payload = {
        "E": np.ascontiguousarray(effects, dtype=effects_dtype),
        "perts": np.ascontiguousarray(
            perturbations, dtype=f"<U{max(len(value) for value in perturbations.tolist())}"
        ),
        "genes": np.ascontiguousarray(
            genes, dtype=f"<U{max(len(value) for value in genes.tolist())}"
        ),
    }
    if "ncells" in dictionary:
        payload["ncells"] = np.ascontiguousarray(dictionary["ncells"], dtype="<i8")
    return payload


def _npy_bytes(array: np.ndarray) -> bytes:
    buffer = BytesIO()
    np.lib.format.write_array(buffer, array, allow_pickle=False)
    return buffer.getvalue()


def save_effect_dictionary(path: str | Path, dictionary: Mapping[str, Any]) -> Path:
    """Atomically save a pickle-free, byte-stable ``.npz`` effect dictionary.

    Members are uncompressed ``.npy`` payloads with fixed ZIP timestamps, permissions,
    ordering, endianness, and C layout. Avoiding implementation-dependent compression makes
    identical dictionaries reproduce the same full-file SHA-256 across runs and platforms.
    """

    destination = Path(path)
    if destination.suffix.lower() != ".npz":
        raise ValueError("effect dictionary path must end in .npz")
    problems = validate_effect_dictionary(dictionary)
    if problems:
        raise ValueError(
            "refusing to save malformed effect dictionary:\n  - "
            + "\n  - ".join(problems)
        )
    payload = _portable_payload(dictionary)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=destination.parent, prefix=f".{destination.name}.", delete=False
        ) as temporary:
            temporary_name = temporary.name
            with ZipFile(temporary, mode="w", compression=ZIP_STORED, allowZip64=True) as archive:
                for key in ("E", "perts", "genes", "ncells"):
                    if key not in payload:
                        continue
                    info = ZipInfo(f"{key}.npy", date_time=(1980, 1, 1, 0, 0, 0))
                    info.compress_type = ZIP_STORED
                    info.create_system = 3
                    info.external_attr = 0o600 << 16
                    archive.writestr(info, _npy_bytes(payload[key]))
        Path(temporary_name).replace(destination)
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
    return destination


def load_effect_dictionary(path: str | Path) -> dict[str, np.ndarray]:
    """Load and validate an effect dictionary without enabling pickle.

    Object arrays are refused rather than deserialized.  This makes the format
    suitable for ordinary numeric/string NPZ exchange; it does not make an
    untrusted archive immune to resource-exhaustion attacks.
    """

    source = Path(path)
    try:
        loaded = np.load(source, allow_pickle=False)
    except (OSError, ValueError) as error:
        raise ValueError(f"could not open {source} as a safe NPZ archive") from error
    if not isinstance(loaded, np.lib.npyio.NpzFile):
        raise ValueError(f"{source} is not an NPZ archive")
    try:
        try:
            dictionary = {key: loaded[key].copy() for key in loaded.files}
        except ValueError as error:
            raise ValueError(
                f"{source} contains an object/pickled array; only numeric and string arrays are allowed"
            ) from error
    finally:
        loaded.close()

    problems = validate_effect_dictionary(dictionary)
    if problems:
        raise ValueError(
            f"{source} is not a valid effect dictionary:\n  - "
            + "\n  - ".join(problems)
        )
    return dictionary
