#!/usr/bin/env python3
"""Rebuild the registered retrospective effect dictionaries from frozen sources.

The Norman and Replogle inputs contain already-processed cell-by-gene expression.
Their builder computes pooled condition-minus-control means on that processed scale; it
does not normalize raw counts, estimate replicate-aware effects, or provide biological
inference.  The Zhu builder extracts the registered donor-collapsed ``log_fc`` DE layer.

Every source is checked against its registered byte length and full-file SHA-256 before
it is opened.  Outputs use :func:`effect_dictionary.save_effect_dictionary`, whose
pickle-free NPZ encoding is byte stable.  ``--check`` rebuilds into temporary storage and
requires the complete output bytes, SHA-256, and array schema to match the config without
modifying the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Mapping, Sequence
from pathlib import Path
import re
import sys
import tempfile
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from effect_dictionary import (  # noqa: E402
    build_effect_dictionary,
    load_effect_dictionary,
    save_effect_dictionary,
)
from reachability import InputError  # noqa: E402
from scripts.run_source_reconstruction import (  # noqa: E402
    _load_effect_subset,
    build_profile,
)


DEFAULT_CONFIG = ROOT / "configs" / "library_coverage_crossdataset.json"
DATASETS = ("zhu", "norman", "replogle")
PORTABLE_FORMAT = "portable_effect_dictionary_v1"
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def _progress(message: str) -> None:
    print(f"[cache-reconstruction] {message}", file=sys.stderr, flush=True)


def sha256_file(path: str | Path) -> str:
    """Return a streaming full-file SHA-256 digest."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_exact_keys(
    value: Any, *, required: set[str], context: str
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise InputError(f"{context} must be an object")
    observed = set(value)
    missing = sorted(required - observed)
    unexpected = sorted(observed - required)
    if missing or unexpected:
        raise InputError(
            f"{context} keys differ from the contract: "
            f"missing={missing}, unexpected={unexpected}"
        )
    return value


def _registered_path(project_root: Path, relative_path: Any, *, context: str) -> Path:
    if not isinstance(relative_path, str) or not relative_path.strip():
        raise InputError(f"{context}.path must be a non-empty repository-relative string")
    relative = Path(relative_path)
    if relative.is_absolute():
        raise InputError(f"{context}.path must be repository-relative")
    root = project_root.resolve()
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root):
        raise InputError(f"{context}.path escapes the repository root")
    return resolved


def _validate_identity_fields(spec: Mapping[str, Any], *, context: str) -> None:
    expected_bytes = spec["bytes"]
    if type(expected_bytes) is not int or expected_bytes <= 0:
        raise InputError(f"{context}.bytes must be a positive integer")
    expected_hash = spec["sha256"]
    if not isinstance(expected_hash, str) or _SHA256_PATTERN.fullmatch(expected_hash) is None:
        raise InputError(f"{context}.sha256 must be 64 lowercase hexadecimal characters")


def verify_registered_source(
    source_spec: Mapping[str, Any], *, project_root: Path = ROOT
) -> Path:
    """Verify a raw source's complete registered identity before any data loading."""

    spec = _require_exact_keys(
        source_spec,
        required={"path", "bytes", "sha256"},
        context="cache_build source",
    )
    _validate_identity_fields(spec, context="cache_build source")
    path = _registered_path(project_root, spec["path"], context="cache_build source")
    if not path.is_file():
        raise FileNotFoundError(f"registered source is absent: {path}")
    actual_bytes = path.stat().st_size
    if actual_bytes != spec["bytes"]:
        raise InputError(
            f"registered source byte length differs: expected {spec['bytes']}, "
            f"found {actual_bytes}"
        )
    actual_hash = sha256_file(path)
    if actual_hash != spec["sha256"]:
        raise InputError("registered source SHA-256 differs from the frozen identity")
    return path


def validate_cache_build_section(
    config: Mapping[str, Any], *, project_root: Path = ROOT
) -> Mapping[str, Any]:
    """Fail closed unless the cache-builder declarations have the exact v1 schema."""

    if not isinstance(config, Mapping) or "cache_build" not in config:
        raise InputError("config is missing cache_build")
    section = _require_exact_keys(
        config["cache_build"],
        required={"format", "writer_contract", *DATASETS},
        context="cache_build",
    )
    if section["format"] != PORTABLE_FORMAT:
        raise InputError(f"cache_build.format must be {PORTABLE_FORMAT!r}")
    if not isinstance(section["writer_contract"], str) or not section[
        "writer_contract"
    ].strip():
        raise InputError("cache_build.writer_contract must be non-empty")

    zhu = _require_exact_keys(
        section["zhu"],
        required={
            "builder",
            "source_config",
            "layer",
            "gene_universe",
            "effect_rows",
        },
        context="cache_build.zhu",
    )
    if zhu["builder"] != "zhu_de_layer":
        raise InputError("cache_build.zhu.builder must be 'zhu_de_layer'")
    if zhu["layer"] != "log_fc" or zhu["gene_universe"] != "shared_screen_genes":
        raise InputError(
            "the registered Zhu builder must use log_fc on shared_screen_genes"
        )
    _registered_path(
        project_root,
        zhu["source_config"],
        context="cache_build.zhu.source_config",
    )

    cell_required = {
        "builder",
        "source",
        "condition_field",
        "control_field",
        "control_label",
        "gene_field",
        "gene_axis_note",
        "source_matrix_semantics",
    }
    for dataset in ("norman", "replogle"):
        required = set(cell_required)
        if dataset == "replogle":
            required.add("feature_name_note")
        spec = _require_exact_keys(
            section[dataset],
            required=required,
            context=f"cache_build.{dataset}",
        )
        if spec["builder"] != "condition_minus_control_mean":
            raise InputError(
                f"cache_build.{dataset}.builder must be "
                "'condition_minus_control_mean'"
            )
        for field in ("condition_field", "control_field", "control_label", "gene_field"):
            if not isinstance(spec[field], str) or not spec[field].strip():
                raise InputError(f"cache_build.{dataset}.{field} must be non-empty")
            if spec[field] != spec[field].strip():
                raise InputError(f"cache_build.{dataset}.{field} must be whitespace-trimmed")
        for field in (
            "source_matrix_semantics",
            "gene_axis_note",
            *(('feature_name_note',) if dataset == "replogle" else ()),
        ):
            if not isinstance(spec[field], str) or not spec[field].strip():
                raise InputError(f"cache_build.{dataset}.{field} must be non-empty")
        source = _require_exact_keys(
            spec["source"],
            required={
                "path",
                "bytes",
                "sha256",
                "url",
                "etag_observed",
                "last_modified_observed",
                "version_id_observed",
                "retrieved_on",
            },
            context=f"cache_build.{dataset}.source",
        )
        _validate_identity_fields(source, context=f"cache_build.{dataset}.source")
        for field in (
            "url",
            "etag_observed",
            "last_modified_observed",
            "version_id_observed",
            "retrieved_on",
        ):
            if not isinstance(source[field], str) or not source[field].strip():
                raise InputError(
                    f"cache_build.{dataset}.source.{field} must be non-empty"
                )
    return section


def _string_vector(values: Sequence[Any], *, name: str) -> np.ndarray:
    items = list(values)
    if any(not isinstance(value, (str, np.str_)) for value in items):
        raise InputError(f"{name} must contain strings without missing values")
    result = np.asarray(items, dtype=str)
    if result.ndim != 1:
        raise InputError(f"{name} must be one-dimensional")
    if any(not value.strip() or value != value.strip() for value in result.tolist()):
        raise InputError(f"{name} must contain non-empty, whitespace-trimmed strings")
    return result


def validate_control_annotation(
    conditions: Sequence[Any],
    controls: Sequence[Any],
    *,
    control_label: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate independent control annotations and return a semantic boolean axis.

    The registered benchmark H5ADs encode this logical field as integer 0/1.  Native
    booleans and that exact integer encoding are accepted; floats, strings, missing
    values, and any other integers fail closed.
    """

    condition_array = _string_vector(conditions, name="condition annotation")
    raw_controls = np.asarray(controls)
    if raw_controls.ndim != 1 or raw_controls.shape != condition_array.shape:
        raise InputError("condition and control annotations must share one cell axis")
    if raw_controls.dtype.kind == "b":
        control_array = raw_controls.astype(bool, copy=False)
    elif raw_controls.dtype.kind in {"i", "u"} and np.all(
        (raw_controls == 0) | (raw_controls == 1)
    ):
        control_array = raw_controls.astype(bool)
    else:
        raise InputError("control annotation must be boolean or integer encoded as 0/1")
    expected = condition_array == control_label
    if not np.array_equal(control_array, expected):
        mismatches = int(np.count_nonzero(control_array != expected))
        raise InputError(
            "boolean control annotation disagrees with condition == "
            f"{control_label!r} for {mismatches} cells"
        )
    return condition_array, control_array


def build_condition_minus_control_dictionary(
    build_spec: Mapping[str, Any], *, project_root: Path = ROOT
) -> dict[str, np.ndarray]:
    """Build pooled effects from a registered processed-expression H5AD.

    Source verification deliberately precedes importing AnnData or opening the H5AD.
    """

    source_identity = {
        key: build_spec["source"][key] for key in ("path", "bytes", "sha256")
    }
    source_path = verify_registered_source(source_identity, project_root=project_root)
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - exercised in the external environment
        raise RuntimeError(
            "cache reconstruction requires requirements-external.txt"
        ) from exc

    adata = ad.read_h5ad(source_path)
    condition_field = build_spec["condition_field"]
    control_field = build_spec["control_field"]
    gene_field = build_spec["gene_field"]
    missing_obs = [
        field for field in (condition_field, control_field) if field not in adata.obs.columns
    ]
    if missing_obs:
        raise InputError(f"source H5AD is missing obs fields: {missing_obs}")
    if adata.X is None:
        raise InputError("source H5AD has no processed expression matrix in X")

    condition_series = adata.obs[condition_field]
    control_series = adata.obs[control_field]
    if bool(condition_series.isna().any()):
        raise InputError("condition annotation contains missing values")
    if bool(control_series.isna().any()):
        raise InputError("control annotation contains missing values")
    raw_controls = control_series.to_numpy()
    conditions, _ = validate_control_annotation(
        condition_series.tolist(),
        raw_controls,
        control_label=build_spec["control_label"],
    )

    if gene_field in adata.var.columns:
        gene_series = adata.var[gene_field]
        if bool(gene_series.isna().any()):
            raise InputError(f"source H5AD var.{gene_field} contains missing values")
        gene_values = gene_series.tolist()
    elif adata.var.index.name == gene_field:
        gene_values = adata.var.index.tolist()
    else:
        raise InputError(
            f"source H5AD is missing configured var field/index {gene_field!r}"
        )
    genes = _string_vector(gene_values, name=f"var.{gene_field}")
    if np.unique(genes).size != genes.size:
        raise InputError(f"source H5AD var.{gene_field} labels must be unique")

    return build_effect_dictionary(
        adata.X,
        conditions,
        control_label=build_spec["control_label"],
        gene_names=genes,
        dtype=np.float32,
    )


def build_zhu_dictionary(
    build_spec: Mapping[str, Any], *, project_root: Path = ROOT
) -> dict[str, np.ndarray]:
    """Extract the registered Zhu log-fold-change dictionary and its exact axes."""

    source_config_path = _registered_path(
        project_root,
        build_spec["source_config"],
        context="cache_build.zhu.source_config",
    )
    if not source_config_path.is_file():
        raise FileNotFoundError(f"Zhu source config is absent: {source_config_path}")
    source_config = json.loads(source_config_path.read_text(encoding="utf-8"))
    if source_config.get("schema_version") != "1.0.0":
        raise InputError("unsupported source-reconstruction config schema")

    for name, source_spec in source_config.get("inputs", {}).items():
        if not isinstance(source_spec, Mapping) or "path" not in source_spec:
            raise InputError(f"source_reconstruction.inputs.{name} is malformed")
        _registered_path(
            project_root,
            source_spec["path"],
            context=f"source_reconstruction.inputs.{name}",
        )
    # build_profile verifies all raw file identities before profile_h5ad loads data.
    _, internals = build_profile(
        source_config, verify_hash=True, project_root=project_root
    )
    axes = internals["axes"]
    lineage = internals["lineages"][build_spec["layer"]]
    genes = tuple(lineage[build_spec["gene_universe"]])
    source_h5ad = _registered_path(
        project_root,
        source_config["inputs"]["de_stats"]["path"],
        context="source_reconstruction.inputs.de_stats",
    )
    effects, _ = _load_effect_subset(
        source_h5ad,
        build_spec["layer"],
        axes["selected_rows"],
        axes["genes"],
        genes,
    )
    return {
        "E": np.asarray(effects, dtype=np.float32),
        "perts": np.asarray(axes["selected_perturbations"], dtype=str),
        "genes": np.asarray(genes, dtype=str),
    }


def build_dataset_dictionary(
    dataset: str,
    cache_build: Mapping[str, Any],
    *,
    project_root: Path = ROOT,
) -> dict[str, np.ndarray]:
    """Dispatch one validated cache declaration to its deterministic builder."""

    if dataset not in DATASETS:
        raise InputError(f"unknown cache dataset: {dataset!r}")
    build_spec = cache_build[dataset]
    if dataset == "zhu":
        return build_zhu_dictionary(build_spec, project_root=project_root)
    return build_condition_minus_control_dictionary(build_spec, project_root=project_root)


def assert_cache_contract(path: str | Path, input_spec: Mapping[str, Any]) -> dict[str, Any]:
    """Require exact artifact identity and complete portable-array schema."""

    cache_path = Path(path)
    spec = _require_exact_keys(
        input_spec,
        required={
            "path",
            "bytes",
            "sha256",
            "effect_array",
            "label_array",
            "required_arrays",
        },
        context="registered output",
    )
    if spec["effect_array"] != "E" or spec["label_array"] != "perts":
        raise InputError("portable output must declare effect_array='E' and label_array='perts'")
    _validate_identity_fields(spec, context="registered output")
    required_arrays = spec["required_arrays"]
    if not isinstance(required_arrays, Mapping) or not required_arrays:
        raise InputError("registered output required_arrays must be a non-empty object")

    actual_bytes = cache_path.stat().st_size
    if actual_bytes != spec["bytes"]:
        raise InputError(
            f"rebuilt cache byte length differs: expected {spec['bytes']}, "
            f"found {actual_bytes}"
        )
    actual_hash = sha256_file(cache_path)
    if actual_hash != spec["sha256"]:
        raise InputError("rebuilt cache SHA-256 differs from the registered output")

    dictionary = load_effect_dictionary(cache_path)
    if set(dictionary) != set(required_arrays):
        raise InputError(
            "rebuilt cache arrays differ from the complete registered schema: "
            f"expected={sorted(required_arrays)}, found={sorted(dictionary)}"
        )
    for key, expected in required_arrays.items():
        expected_spec = _require_exact_keys(
            expected,
            required={"shape", "dtype"},
            context=f"registered output required_arrays.{key}",
        )
        array = dictionary[key]
        if list(array.shape) != expected_spec["shape"]:
            raise InputError(
                f"rebuilt cache {key} shape differs: "
                f"expected {expected_spec['shape']}, found {list(array.shape)}"
            )
        if str(array.dtype) != expected_spec["dtype"]:
            raise InputError(
                f"rebuilt cache {key} dtype differs: "
                f"expected {expected_spec['dtype']}, found {array.dtype}"
            )
    return {
        "path": str(cache_path),
        "bytes": actual_bytes,
        "sha256": actual_hash,
        "arrays": sorted(dictionary),
    }


def reconstruct_dataset(
    dataset: str,
    config: Mapping[str, Any],
    *,
    project_root: Path = ROOT,
    output_dir: Path | None = None,
    check: bool = False,
) -> dict[str, Any]:
    """Build and validate one cache, optionally installing it after validation."""

    cache_build = validate_cache_build_section(config, project_root=project_root)
    if "inputs" not in config or dataset not in config["inputs"]:
        raise InputError(f"config is missing inputs.{dataset}")
    input_spec = config["inputs"][dataset]
    registered_output = _registered_path(
        project_root, input_spec.get("path"), context=f"inputs.{dataset}"
    )

    if check:
        with tempfile.TemporaryDirectory(prefix=f"{dataset}-cache-check-") as temporary:
            staged = Path(temporary) / registered_output.name
            _progress(f"{dataset}: verifying sources and rebuilding in temporary storage")
            dictionary = build_dataset_dictionary(
                dataset, cache_build, project_root=project_root
            )
            save_effect_dictionary(staged, dictionary)
            status = assert_cache_contract(staged, input_spec)
            status["path"] = str(registered_output)
        return {"dataset": dataset, "status": "PASS", **status}

    destination = (
        Path(output_dir).resolve() / registered_output.name
        if output_dir is not None
        else registered_output
    )
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite existing cache: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{destination.name}.build-", dir=destination.parent
    ) as temporary:
        staged = Path(temporary) / destination.name
        _progress(f"{dataset}: verifying sources and rebuilding {destination.name}")
        dictionary = build_dataset_dictionary(
            dataset, cache_build, project_root=project_root
        )
        save_effect_dictionary(staged, dictionary)
        status = assert_cache_contract(staged, input_spec)
        try:
            os.link(staged, destination)
        except FileExistsError as exc:
            raise FileExistsError(
                f"refusing to overwrite existing cache: {destination}"
            ) from exc
    return {"dataset": dataset, "status": "BUILT", **status, "path": str(destination)}


def run(
    config: Mapping[str, Any],
    *,
    dataset: str = "all",
    project_root: Path = ROOT,
    output_dir: Path | None = None,
    check: bool = False,
) -> list[dict[str, Any]]:
    """Reconstruct one or all configured effect dictionaries."""

    validate_cache_build_section(config, project_root=project_root)
    selected = DATASETS if dataset == "all" else (dataset,)
    return [
        reconstruct_dataset(
            name,
            config,
            project_root=project_root,
            output_dir=output_dir,
            check=check,
        )
        for name in selected
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--dataset", choices=("all", *DATASETS), default="all")
    destinations = parser.add_mutually_exclusive_group()
    destinations.add_argument(
        "--output-dir",
        type=Path,
        help="write the registered cache filename under this directory",
    )
    destinations.add_argument(
        "--check",
        action="store_true",
        help="rebuild temporarily and require exact registered identity/schema",
    )
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    results = run(
        config,
        dataset=args.dataset,
        output_dir=args.output_dir,
        check=args.check,
    )
    for result in results:
        print(
            f"{result['dataset']}: {result['status']} "
            f"bytes={result['bytes']} sha256={result['sha256']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
