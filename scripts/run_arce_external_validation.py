#!/usr/bin/env python3
"""Run the frozen Arce external IL2RA and supplied-score benchmarks.

The screen outcomes are used only after guide-count eligibility and Zhu-generator
availability/admission have been fixed. The benchmark is deliberately narrow: it
tests cross-study, cross-modality rank alignment and makes no state-reachability,
donor-generalization, or causal-validation claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import struct
import sys
from typing import Any
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
from scipy.stats import kendalltau, rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reachability import InputError  # noqa: E402


DEFAULT_CONFIG = ROOT / "configs" / "arce_external_validation.json"
_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_DOC_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"

S14_COLUMNS = [
    "cell", "orig.ident", "nCount_RNA", "nFeature_RNA", "n_sgrna_features",
    "sgrna", "n_sgrna_umis", "sg_target", "has_sgrna", "nCount_HTO",
    "nFeature_HTO", "nCount_ADT", "nFeature_ADT", "donor", "percent.mt",
    "percent.ribo", "HTO_maxID", "HTO_margin", "HTO_classification",
    "HTO_classification.global", "hash.ID", "nCount_TCRVJ", "nFeature_TCRVJ",
    "nCount_SCT", "nFeature_SCT", "S.Score", "G2M.Score", "Phase",
    "CC.Difference", "SCT.weight", "ADT.weight", "activation.score",
]
S14_SELECTED_COLUMNS = (
    "cell", "n_sgrna_features", "sgrna", "sg_target", "has_sgrna", "donor",
    "HTO_maxID", "HTO_classification.global", "activation.score",
)
S8_COLUMNS = [
    "sg_target", "median.activation.score", "mean.activation.score", "HTO_maxID",
    "statistic", "p.value", "method", "alternative", "stars", "padj",
]


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decode_strings(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        item.decode("utf-8") if isinstance(item, bytes) else str(item)
        for item in values
    )


def _read_h5_categorical(group: Any) -> np.ndarray:
    categories = _decode_strings(group["categories"][:])
    codes = np.asarray(group["codes"][:], dtype=np.int64)
    if np.any(codes < 0) or np.any(codes >= len(categories)):
        raise InputError("H5AD categorical codes contain missing or invalid values")
    return np.asarray([categories[index] for index in codes], dtype=object)


def _column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if letters is None:
        raise InputError(f"invalid XLSX cell reference: {reference!r}")
    value = 0
    for letter in letters.group(0):
        value = value * 26 + ord(letter) - ord("A") + 1
    return value - 1


def _xlsx_shared_strings(book: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in book.namelist():
        return []
    root = ET.fromstring(book.read("xl/sharedStrings.xml"))
    output = []
    for item in root.findall(f"{{{_MAIN_NS}}}si"):
        output.append("".join(node.text or "" for node in item.iter(f"{{{_MAIN_NS}}}t")))
    return output


def _xlsx_shared_strings_stream(book: zipfile.ZipFile) -> list[str]:
    """Read shared strings without retaining their expanded XML tree."""

    if "xl/sharedStrings.xml" not in book.namelist():
        return []
    output: list[str] = []
    with book.open("xl/sharedStrings.xml") as handle:
        for _, item in ET.iterparse(handle, events=("end",)):
            if item.tag == f"{{{_MAIN_NS}}}si":
                output.append(
                    "".join(node.text or "" for node in item.iter(f"{{{_MAIN_NS}}}t"))
                )
                item.clear()
    return output


def _xlsx_sheet_path(book: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(book.read("xl/workbook.xml"))
    relation_id = None
    for sheet in workbook.findall(f".//{{{_MAIN_NS}}}sheet"):
        if sheet.attrib.get("name") == sheet_name:
            relation_id = sheet.attrib.get(f"{{{_DOC_REL_NS}}}id")
            break
    if relation_id is None:
        raise InputError(f"XLSX is missing sheet {sheet_name!r}")
    relationships = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
    for relation in relationships.findall(f"{{{_PKG_REL_NS}}}Relationship"):
        if relation.attrib.get("Id") == relation_id:
            target = relation.attrib["Target"].lstrip("/")
            return target if target.startswith("xl/") else f"xl/{target}"
    raise InputError(f"XLSX has no relationship target for sheet {sheet_name!r}")


def _xlsx_value(cell: ET.Element, shared: list[str]) -> Any:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(f"{{{_MAIN_NS}}}t"))
    value_node = cell.find(f"{{{_MAIN_NS}}}v")
    if value_node is None or value_node.text is None:
        return None
    raw = value_node.text
    if cell_type == "s":
        try:
            return shared[int(raw)]
        except (IndexError, ValueError) as exc:
            raise InputError("XLSX contains an invalid shared-string index") from exc
    if cell_type in {"str", "e"}:
        return raw
    try:
        number = float(raw)
    except ValueError:
        return raw
    return int(number) if number.is_integer() else number


def parse_xlsx_table(payload: bytes, sheet_name: str) -> tuple[list[str], list[list[Any]]]:
    """Parse a simple tabular XLSX without adding an undeclared runtime dependency."""

    with zipfile.ZipFile(io.BytesIO(payload)) as book:
        shared = _xlsx_shared_strings(book)
        sheet_path = _xlsx_sheet_path(book, sheet_name)
        root = ET.fromstring(book.read(sheet_path))
    rows: list[list[Any]] = []
    for row in root.findall(f".//{{{_MAIN_NS}}}row"):
        values: list[Any] = []
        for cell in row.findall(f"{{{_MAIN_NS}}}c"):
            index = _column_index(cell.attrib.get("r", ""))
            if index >= len(values):
                values.extend([None] * (index + 1 - len(values)))
            values[index] = _xlsx_value(cell, shared)
        while values and values[-1] is None:
            values.pop()
        if values:
            rows.append(values)
    if not rows:
        raise InputError("XLSX sheet is empty")
    if any(value is None or not str(value).strip() for value in rows[0]):
        raise InputError("XLSX header contains an empty field")
    headers = [str(value).strip() for value in rows[0]]
    if len(set(headers)) != len(headers):
        raise InputError("XLSX header fields are not unique")
    return headers, rows[1:]


def iter_xlsx_selected_rows(
    payload: bytes,
    sheet_name: str,
    expected_headers: list[str],
    selected_fields: tuple[str, ...],
):
    """Stream selected XLSX columns after verifying the complete ordered header."""

    selected = set(selected_fields)
    if not selected <= set(expected_headers):
        raise ValueError("selected XLSX fields are absent from expected headers")
    with zipfile.ZipFile(io.BytesIO(payload)) as book:
        shared = _xlsx_shared_strings_stream(book)
        sheet_path = _xlsx_sheet_path(book, sheet_name)
        header_seen = False
        index_to_field = {index: field for index, field in enumerate(expected_headers)}
        with book.open(sheet_path) as handle:
            for _, row in ET.iterparse(handle, events=("end",)):
                if row.tag != f"{{{_MAIN_NS}}}row":
                    continue
                values: dict[int, Any] = {}
                for cell in row.findall(f"{{{_MAIN_NS}}}c"):
                    index = _column_index(cell.attrib.get("r", ""))
                    if not header_seen or index_to_field.get(index) in selected:
                        values[index] = _xlsx_value(cell, shared)
                if not header_seen:
                    width = max(values, default=-1) + 1
                    headers = [values.get(index) for index in range(width)]
                    if headers != expected_headers:
                        raise InputError("XLSX header/order differs from the frozen schema")
                    header_seen = True
                else:
                    yield {
                        field: values.get(expected_headers.index(field))
                        for field in selected_fields
                    }
                row.clear()
        if not header_seen:
            raise InputError("XLSX sheet is empty")


def verify_and_read_arce_members(
    config: dict[str, Any]
) -> tuple[dict[str, bytes], dict[str, Any]]:
    archive_spec = config["dataset"]["archive"]
    archive_path = ROOT / archive_spec["path"]
    if not archive_path.is_file():
        raise FileNotFoundError(archive_path)
    if archive_path.stat().st_size != archive_spec["bytes"]:
        raise InputError("Arce archive byte length differs from the frozen contract")
    archive_hash = sha256_file(archive_path)
    if archive_hash != archive_spec["sha256"]:
        raise InputError("Arce archive SHA-256 differs from the frozen contract")
    archive_md5 = md5_file(archive_path)
    if archive_md5 != archive_spec["md5"]:
        raise InputError("Arce archive MD5 differs from the frozen contract")
    member_keys = [
        key for key in (
            "screen_member", "activation_cells_member", "activation_summary_member"
        )
        if key in config["dataset"]
    ]
    payloads: dict[str, bytes] = {}
    identities: dict[str, Any] = {}
    with zipfile.ZipFile(archive_path) as archive:
        for key in member_keys:
            member_spec = config["dataset"][key]
            try:
                payload = archive.read(member_spec["path"])
            except KeyError as exc:
                raise InputError(f"Arce archive is missing frozen member {key}") from exc
            if len(payload) != member_spec["bytes"]:
                raise InputError(f"Arce {key} byte length differs from the frozen contract")
            member_hash = sha256_bytes(payload)
            if member_hash != member_spec["sha256"]:
                raise InputError(f"Arce {key} SHA-256 differs from the frozen contract")
            payloads[key] = payload
            identities[key] = {
                "path": member_spec["path"],
                "bytes": len(payload),
                "sha256": member_hash,
                "sheet": member_spec["sheet"],
            }
    return payloads, {
        "archive_path": archive_spec["path"],
        "archive_bytes": archive_spec["bytes"],
        "archive_md5": archive_md5,
        "archive_sha256": archive_hash,
        "members": identities,
    }


def _finite_number(value: Any, *, field: str, row_number: int) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(f"S1 row {row_number}: {field} is not numeric") from exc
    if not math.isfinite(result):
        raise InputError(f"S1 row {row_number}: {field} is not finite")
    return result


def expected_s1_columns(benchmark: dict[str, Any]) -> list[str]:
    columns = ["id"]
    suffixes = (
        "num.",
        "neg|score.",
        "neg|p-value.",
        "neg|fdr.",
        "neg|rank.",
        "neg|goodsgrna.",
        "neg|lfc.",
        "pos|score.",
        "pos|p-value.",
        "pos|fdr.",
        "pos|rank.",
        "pos|goodsgrna.",
        "pos|lfc.",
    )
    for context in benchmark["contexts"]:
        columns.extend(f"{suffix}{context}" for suffix in suffixes)
    return columns


def load_s1(payload: bytes, config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    benchmark = config["benchmark"]
    sheet_name = config["dataset"]["screen_member"]["sheet"]
    headers, raw_rows = parse_xlsx_table(payload, sheet_name)
    expected_headers = expected_s1_columns(benchmark)
    if headers != expected_headers:
        raise InputError("Arce S1 header/order differs from the frozen schema")
    if len(raw_rows) != config["expected"]["screen_rows"]:
        raise InputError("Arce S1 row count differs from the frozen contract")
    index = {field: position for position, field in enumerate(headers)}
    screen: dict[str, dict[str, Any]] = {}
    for row_number, raw in enumerate(raw_rows, start=2):
        if len(raw) < len(headers):
            raw = [*raw, *([None] * (len(headers) - len(raw)))]
        if len(raw) != len(headers) or any(value is None for value in raw):
            raise InputError(f"S1 row {row_number}: incomplete row")
        target = str(raw[index["id"]]).strip()
        if not target or target in screen:
            raise InputError(f"S1 row {row_number}: target is empty or duplicated")
        guide_counts: dict[str, int] = {}
        outcomes: dict[str, dict[str, float]] = {}
        for context in benchmark["contexts"]:
            count_value = _finite_number(
                raw[index[f"num.{context}"]], field=f"num.{context}", row_number=row_number
            )
            if not count_value.is_integer() or count_value < 0:
                raise InputError(f"S1 row {row_number}: invalid guide count")
            guide_counts[context] = int(count_value)
            observed = _finite_number(
                raw[index[f"{benchmark['observed_lfc_prefix']}{context}"]],
                field=f"{benchmark['observed_lfc_prefix']}{context}",
                row_number=row_number,
            )
            parallel = _finite_number(
                raw[index[f"{benchmark['parallel_lfc_prefix']}{context}"]],
                field=f"{benchmark['parallel_lfc_prefix']}{context}",
                row_number=row_number,
            )
            if observed != parallel:
                raise InputError(f"S1 row {row_number}: parallel LFC fields differ")
            outcomes[context] = {
                "lfc": observed,
                "positive_fdr": _finite_number(
                    raw[index[f"{benchmark['positive_fdr_prefix']}{context}"]],
                    field=f"{benchmark['positive_fdr_prefix']}{context}",
                    row_number=row_number,
                ),
                "negative_fdr": _finite_number(
                    raw[index[f"{benchmark['negative_fdr_prefix']}{context}"]],
                    field=f"{benchmark['negative_fdr_prefix']}{context}",
                    row_number=row_number,
                ),
            }
        required = benchmark["required_guides_per_context"]
        screen[target] = {
            "guide_counts": guide_counts,
            "four_guide_eligible": all(value == required for value in guide_counts.values()),
            "outcomes": outcomes,
        }
    eligible = sum(item["four_guide_eligible"] for item in screen.values())
    if eligible != config["expected"]["four_guide_eligible"]:
        raise InputError("Arce S1 four-guide eligibility differs from the frozen contract")
    return screen


def _truthy_xlsx(value: Any) -> bool:
    return value is True or value == 1 or str(value).strip().lower() == "true"


def load_s14(payload: bytes, config: dict[str, Any]) -> dict[str, Any]:
    """Stream and validate S14 while retaining only registered score strata."""

    benchmark = config["activation_benchmark"]
    member = config["dataset"]["activation_cells_member"]
    cells: set[str] = set()
    groups: dict[tuple[str, str, str, str], list[float]] = {}
    pooled: dict[tuple[str, str], list[float]] = {}
    target_guides: dict[str, set[str]] = {}
    guide_targets: dict[str, set[str]] = {}
    context_counts: dict[str, int] = {}
    donor_counts: dict[str, int] = {}
    rows = 0
    for row_number, row in enumerate(
        iter_xlsx_selected_rows(
            payload, member["sheet"], S14_COLUMNS, S14_SELECTED_COLUMNS
        ),
        start=2,
    ):
        rows += 1
        if any(row[field] is None for field in S14_SELECTED_COLUMNS):
            raise InputError(f"S14 row {row_number}: missing registered field")
        cell = str(row["cell"]).strip()
        guide = str(row["sgrna"]).strip()
        target = str(row["sg_target"]).strip()
        donor = str(row["donor"]).strip()
        context = str(row["HTO_maxID"]).strip()
        if not cell or cell in cells:
            raise InputError(f"S14 row {row_number}: cell is empty or duplicated")
        cells.add(cell)
        if row["n_sgrna_features"] != 1:
            raise InputError(f"S14 row {row_number}: expected one sgRNA feature")
        if not _truthy_xlsx(row["has_sgrna"]):
            raise InputError(f"S14 row {row_number}: has_sgrna is not true")
        if str(row["HTO_classification.global"]).strip() != "Singlet":
            raise InputError(f"S14 row {row_number}: cell is not an HTO singlet")
        if donor not in benchmark["donors"] or context not in benchmark["contexts"]:
            raise InputError(f"S14 row {row_number}: unknown donor or context")
        try:
            score = float(row["activation.score"])
        except (TypeError, ValueError) as exc:
            raise InputError(f"S14 row {row_number}: activation.score is not numeric") from exc
        if not math.isfinite(score):
            raise InputError(f"S14 row {row_number}: activation.score is not finite")
        groups.setdefault((target, guide, donor, context), []).append(score)
        pooled.setdefault((target, context), []).append(score)
        target_guides.setdefault(target, set()).add(guide)
        guide_targets.setdefault(guide, set()).add(target)
        context_counts[context] = context_counts.get(context, 0) + 1
        donor_counts[donor] = donor_counts.get(donor, 0) + 1

    expected = config["expected"]
    if rows != expected["activation_rows"] or len(cells) != expected["activation_unique_cells"]:
        raise InputError("Arce S14 row or unique-cell count differs from the frozen contract")
    if context_counts != expected["activation_context_counts"]:
        raise InputError("Arce S14 context counts differ from the frozen contract")
    if donor_counts != expected["activation_donor_counts"]:
        raise InputError("Arce S14 donor counts differ from the frozen contract")
    if len(target_guides) != expected["activation_targets"]:
        raise InputError("Arce S14 target count differs from the frozen contract")
    if len(guide_targets) != expected["activation_guides"]:
        raise InputError("Arce S14 guide count differs from the frozen contract")
    if any(len(targets) != 1 for targets in guide_targets.values()):
        raise InputError("Arce S14 guide maps to multiple targets")
    control = benchmark["control_target"]
    if target_guides.get(control) != set(benchmark["control_guides"]):
        raise InputError("Arce S14 Non-Targeting guide identities differ")
    regular_targets = sorted(set(target_guides) - {control})
    if len(regular_targets) != expected["activation_regular_targets"]:
        raise InputError("Arce S14 regular-target count differs")
    if any(
        len(target_guides[target]) != benchmark["regular_guides_per_target"]
        for target in regular_targets
    ):
        raise InputError("Arce S14 regular target does not have exactly two guides")
    if len(groups) != expected["activation_guide_groups"]:
        raise InputError("Arce S14 guide-stratum count differs from the frozen contract")
    target_group_count = len({(target, donor, context) for target, _, donor, context in groups})
    if target_group_count != expected["activation_target_groups"]:
        raise InputError("Arce S14 target-stratum count differs from the frozen contract")
    expected_groups = {
        (target, guide, donor, context)
        for target, guides in target_guides.items()
        for guide in guides
        for donor in benchmark["donors"]
        for context in benchmark["contexts"]
    }
    if set(groups) != expected_groups:
        raise InputError("Arce S14 donor/guide/context factorial is incomplete")
    sizes = [len(values) for values in groups.values()]
    if "activation_guide_group_cells_min" in expected:
        if (
            min(sizes) != expected["activation_guide_group_cells_min"]
            or float(np.median(sizes)) != expected["activation_guide_group_cells_median"]
            or max(sizes) != expected["activation_guide_group_cells_max"]
            or sum(size < 20 for size in sizes)
            != expected["activation_strata_below_20_cells"]
        ):
            raise InputError("Arce S14 guide-stratum size distribution differs")
    return {
        "groups": groups,
        "pooled": pooled,
        "target_guides": {key: tuple(sorted(value)) for key, value in target_guides.items()},
        "regular_targets": tuple(regular_targets),
        "data_quality": {
            "rows": rows,
            "unique_cells": len(cells),
            "targets": len(target_guides),
            "regular_targets": len(regular_targets),
            "guides": len(guide_targets),
            "guide_groups": len(groups),
            "target_groups": target_group_count,
            "context_counts": context_counts,
            "donor_counts": donor_counts,
            "guide_group_cells_min": min(sizes),
            "guide_group_cells_median": float(np.median(sizes)),
            "guide_group_cells_max": max(sizes),
            "missing_registered_fields": 0,
            "non_singlets": 0,
        },
    }


def load_s8_and_check_reproduction(
    payload: bytes, activation: dict[str, Any], config: dict[str, Any]
) -> dict[str, Any]:
    member = config["dataset"]["activation_summary_member"]
    headers, raw_rows = parse_xlsx_table(payload, member["sheet"])
    if headers != S8_COLUMNS:
        raise InputError("Arce S8 header/order differs from the frozen schema")
    if len(raw_rows) != config["expected"]["activation_summary_rows"]:
        raise InputError("Arce S8 row count differs from the frozen contract")
    index = {field: position for position, field in enumerate(headers)}
    observed: dict[tuple[str, str], tuple[float, float]] = {}
    control = config["activation_benchmark"]["control_target"]
    inferential_fields = ("statistic", "p.value", "method", "alternative", "padj")
    for row_number, raw in enumerate(raw_rows, start=2):
        raw = [*raw, *([None] * (len(headers) - len(raw)))]
        target = str(raw[index["sg_target"]]).strip()
        context = str(raw[index["HTO_maxID"]]).strip()
        key = (target, context)
        if key in observed:
            raise InputError(f"S8 row {row_number}: duplicate target/context")
        median = _finite_number(
            raw[index["median.activation.score"]],
            field="median.activation.score",
            row_number=row_number,
        )
        mean = _finite_number(
            raw[index["mean.activation.score"]],
            field="mean.activation.score",
            row_number=row_number,
        )
        missing_inference = [
            raw[index[field]] is None
            or str(raw[index[field]]).strip().upper() in {"", "NA"}
            for field in inferential_fields
        ]
        if target == control:
            if not all(missing_inference):
                raise InputError("S8 Non-Targeting row unexpectedly contains inference fields")
        elif any(missing_inference):
            raise InputError("S8 perturbation row is missing a published inference field")
        observed[key] = (median, mean)
    if set(observed) != set(activation["pooled"]):
        raise InputError("Arce S8 and S14 target/context axes differ")
    median_errors = []
    mean_errors = []
    for key, values in activation["pooled"].items():
        published_median, published_mean = observed[key]
        median_errors.append(abs(float(np.median(values)) - published_median))
        mean_errors.append(abs(float(np.mean(values)) - published_mean))
    tolerance = config["activation_benchmark"]["reproduction_tolerance"]
    if max(median_errors) > tolerance or max(mean_errors) > tolerance:
        raise InputError("Arce S8 aggregates do not reproduce from S14")
    return {
        "status": "PASS",
        "rows": len(observed),
        "maximum_absolute_median_error": max(median_errors),
        "maximum_absolute_mean_error": max(mean_errors),
        "tolerance": tolerance,
        "role": "archival aggregate reproduction only; S8 is derived from S14",
        "published_p_values_used_for_inference": False,
    }


def predictor_identity(predictors: dict[str, dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for target in sorted(predictors):
        item = predictors[target]
        digest.update(target.encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes([int(item["admitted"])]))
        raw_effect = item.get("raw_log_fc", item.get("effect"))
        digest.update(struct.pack("<d", float(raw_effect)))
    return digest.hexdigest()


def load_generator(config: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError("Arce external validation requires requirements-external.txt") from exc

    spec = config["generator"]
    path = ROOT / spec["path"]
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size != spec["bytes"]:
        raise InputError("Zhu H5AD byte length differs from the frozen contract")
    with h5py.File(path, "r") as handle:
        required_paths = (
            f"obs/{spec['condition_field']}",
            f"obs/{spec['perturbation_field']}",
            f"obs/{spec['admission_field']}",
            f"var/{spec['gene_field']}",
            f"layers/{spec['layer']}",
        )
        for required in required_paths:
            if required not in handle:
                raise InputError(f"Zhu H5AD is missing {required}")
        conditions = _read_h5_categorical(handle[f"obs/{spec['condition_field']}"])
        perturbations = _read_h5_categorical(handle[f"obs/{spec['perturbation_field']}"])
        admitted = np.asarray(handle[f"obs/{spec['admission_field']}"][:], dtype=bool)
        genes = _decode_strings(handle[f"var/{spec['gene_field']}"][:])
        if len(set(genes)) != len(genes):
            raise InputError("Zhu H5AD gene axis is not unique")
        try:
            target_column = genes.index(spec["target_gene"])
        except ValueError as exc:
            raise InputError(f"Zhu H5AD lacks target gene {spec['target_gene']}") from exc
        layer = handle[f"layers/{spec['layer']}"]
        expected_shape = (config["expected"]["effect_rows"], config["expected"]["effect_genes"])
        if layer.shape != expected_shape:
            raise InputError("Zhu H5AD effect layer shape differs from the frozen contract")
        if len(conditions) != layer.shape[0] or len(perturbations) != layer.shape[0]:
            raise InputError("Zhu H5AD observation fields do not match the effect layer")
        rows = np.flatnonzero(conditions == spec["condition"])
        effects = np.asarray(layer[rows, target_column], dtype=float)
    if len(rows) != config["expected"]["generator_condition_rows"]:
        raise InputError("Zhu generator condition row count differs from the frozen contract")
    if len(set(perturbations[rows])) != len(rows):
        raise InputError("Zhu generator perturbations are not unique in the selected condition")
    if not np.all(np.isfinite(effects)):
        raise InputError("Zhu IL2RA predictor contains non-finite values")
    raw_predictors = {
        str(perturbations[row]): {
            "effect": float(effect),
            "admitted": bool(admitted[row]),
        }
        for row, effect in zip(rows, effects, strict=True)
    }
    admitted_rows = sum(item["admitted"] for item in raw_predictors.values())
    if admitted_rows != config["expected"]["generator_admitted_rows"]:
        raise InputError("Zhu admitted generator count differs from the frozen contract")
    identity = predictor_identity(raw_predictors)
    if identity != spec["predictor_sha256"]:
        raise InputError("Zhu IL2RA predictor identity differs from the frozen contract")
    multiplier = float(config["benchmark"]["predictor_multiplier"])
    predictors = {
        target: {
            "raw_log_fc": item["effect"],
            "regulator_score": multiplier * item["effect"],
            "admitted": item["admitted"],
        }
        for target, item in raw_predictors.items()
    }
    return predictors, {
        "path": spec["path"],
        "bytes": spec["bytes"],
        "condition": spec["condition"],
        "condition_rows": len(rows),
        "admitted_rows": admitted_rows,
        "target_gene": spec["target_gene"],
        "layer": spec["layer"],
        "predictor_sha256": identity,
        "full_file_sha256_verified": False,
        "identity_note": "The exact used target/effect/admission projection is hash-bound; the 16.8 GB container is byte-bound.",
    }


def _exclusion(item: dict[str, Any], predictor: dict[str, Any] | None) -> str:
    if not item["four_guide_eligible"]:
        return "guide_count_not_four_all_contexts"
    if predictor is None:
        return "generator_missing"
    if not predictor["admitted"]:
        return "generator_not_admitted"
    return ""


def build_prediction_rows(
    screen: dict[str, dict[str, Any]],
    predictors: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    for target in sorted(screen):
        item = screen[target]
        predictor = predictors.get(target)
        exclusion = _exclusion(item, predictor)
        for context in config["benchmark"]["contexts"]:
            outcome = item["outcomes"][context]
            rows.append(
                {
                    "target": target,
                    "context": context,
                    "n_guides": item["guide_counts"][context],
                    "observed_lfc": outcome["lfc"],
                    "observed_positive_fdr": outcome["positive_fdr"],
                    "observed_negative_fdr": outcome["negative_fdr"],
                    "raw_zhu_il2ra_log_fc": None if predictor is None else predictor["raw_log_fc"],
                    "predicted_il2ra_regulator_score": (
                        None if predictor is None else predictor["regulator_score"]
                    ),
                    "four_guide_eligible": item["four_guide_eligible"],
                    "generator_available": predictor is not None,
                    "generator_admitted": False if predictor is None else predictor["admitted"],
                    "analysis_eligible": exclusion == "",
                    "exclusion_reason": exclusion,
                    "observed_rank_abs": None,
                    "predicted_rank_abs": None,
                }
            )
    for context in config["benchmark"]["contexts"]:
        selected = [row for row in rows if row["context"] == context and row["analysis_eligible"]]
        observed = np.asarray([row["observed_lfc"] for row in selected], dtype=float)
        predicted = np.asarray(
            [row["predicted_il2ra_regulator_score"] for row in selected], dtype=float
        )
        observed_ranks = rankdata(-np.abs(observed), method="average")
        predicted_ranks = rankdata(-np.abs(predicted), method="average")
        for row, observed_rank, predicted_rank in zip(
            selected, observed_ranks, predicted_ranks, strict=True
        ):
            row["observed_rank_abs"] = float(observed_rank)
            row["predicted_rank_abs"] = float(predicted_rank)
    return rows


def build_activation_guide_rows(
    activation: dict[str, Any],
    predictors: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    """Compute equally guide-weighted within-donor/context control contrasts."""

    benchmark = config["activation_benchmark"]
    groups = activation["groups"]
    control = benchmark["control_target"]
    medians = {key: float(np.median(values)) for key, values in groups.items()}
    rows: list[dict[str, Any]] = []
    for context in benchmark["contexts"]:
        for donor in benchmark["donors"]:
            control_medians = [
                medians[(control, guide, donor, context)]
                for guide in benchmark["control_guides"]
            ]
            baseline = float(np.median(control_medians))
            for target in activation["regular_targets"]:
                predictor = predictors.get(target)
                for guide in activation["target_guides"][target]:
                    values = groups[(target, guide, donor, context)]
                    guide_median = medians[(target, guide, donor, context)]
                    rows.append(
                        {
                            "target": target,
                            "guide": guide,
                            "donor": donor,
                            "context": context,
                            "n_cells": len(values),
                            "guide_median_supplied_activation_score": guide_median,
                            "ntc_baseline_median_of_guide_medians": baseline,
                            "supplied_activation_score_delta": guide_median - baseline,
                            "support_flag_lt_20_cells": len(values) < 20,
                            "generator_available": predictor is not None,
                            "generator_admitted": False if predictor is None else predictor["admitted"],
                        }
                    )
    return rows


def summarize_activation_robustness(
    rows: list[dict[str, Any]], config: dict[str, Any]
) -> dict[str, Any]:
    """Report guide/donor reliability without treating cells as replicates."""

    benchmark = config["activation_benchmark"]
    result: dict[str, Any] = {}
    for context in benchmark["contexts"]:
        selected = [row for row in rows if row["context"] == context]
        deltas = {
            (row["target"], row["donor"], row["guide"]): row[
                "supplied_activation_score_delta"
            ]
            for row in selected
        }
        target_donor: dict[tuple[str, str], float] = {}
        guide_pairs: dict[str, tuple[list[float], list[float]]] = {}
        guide_sign_agree = []
        all_four_sign_agree = []
        for donor in benchmark["donors"]:
            first_values: list[float] = []
            second_values: list[float] = []
            for target in sorted({row["target"] for row in selected}):
                guides = sorted(
                    row["guide"]
                    for row in selected
                    if row["target"] == target and row["donor"] == donor
                )
                values = [deltas[(target, donor, guide)] for guide in guides]
                target_donor[(target, donor)] = float(np.median(values))
                first_values.append(values[0])
                second_values.append(values[1])
                guide_sign_agree.append(np.sign(values[0]) == np.sign(values[1]))
            guide_pairs[donor] = (first_values, second_values)
        targets = sorted({row["target"] for row in selected})
        for target in targets:
            values = [
                row["supplied_activation_score_delta"]
                for row in selected
                if row["target"] == target
            ]
            signs = np.sign(values)
            all_four_sign_agree.append(bool(np.all(signs == signs[0]) and signs[0] != 0))
        donor_a = [target_donor[(target, benchmark["donors"][0])] for target in targets]
        donor_b = [target_donor[(target, benchmark["donors"][1])] for target in targets]
        guide_correlations = {}
        for donor, (first_values, second_values) in guide_pairs.items():
            guide_correlations[donor] = {
                "spearman": float(spearmanr(first_values, second_values).statistic),
                "kendall": float(kendalltau(first_values, second_values).statistic),
            }
        result[context] = {
            "targets": len(targets),
            "guide_pair_sign_agreement_fraction": float(np.mean(guide_sign_agree)),
            "guide_pair_sign_agreement_n": int(np.count_nonzero(guide_sign_agree)),
            "guide_pair_comparisons": len(guide_sign_agree),
            "all_two_guides_two_donors_same_nonzero_sign_fraction": float(
                np.mean(all_four_sign_agree)
            ),
            "all_two_guides_two_donors_same_nonzero_sign_n": int(
                np.count_nonzero(all_four_sign_agree)
            ),
            "guide_rank_concordance_by_donor": guide_correlations,
            "donor_rank_concordance": {
                "spearman": float(spearmanr(donor_a, donor_b).statistic),
                "kendall": float(kendalltau(donor_a, donor_b).statistic),
                "sign_agreement_fraction": float(
                    np.mean(np.sign(donor_a) == np.sign(donor_b))
                ),
            },
        }
    sparse = [row for row in rows if row["support_flag_lt_20_cells"]]
    return {
        "estimand": benchmark["estimand"],
        "inference": benchmark["inference"],
        "panel_provenance": benchmark["panel_provenance"],
        "score_claim_ceiling": benchmark["score_claim_ceiling"],
        "all_strata_retained": True,
        "outcome_based_filtering": False,
        "strata_below_20_cells": len(sparse),
        "minimum_cell_strata": [
            {
                key: row[key]
                for key in ("target", "guide", "donor", "context", "n_cells")
            }
            for row in sparse
        ],
        "contexts": result,
        "published_s8_p_values_used": False,
        "cell_level_inference_emitted": False,
    }


def _top_targets(targets: list[str], values: np.ndarray, k: int) -> set[str]:
    ordered = sorted(zip(targets, values, strict=True), key=lambda item: (-abs(item[1]), item[0]))
    return {target for target, _ in ordered[:k]}


def _point_metrics(
    targets: list[str], observed: np.ndarray, predicted: np.ndarray, top_k: list[int]
) -> dict[str, Any]:
    spearman = float(spearmanr(predicted, observed).statistic)
    kendall = float(kendalltau(predicted, observed).statistic)
    nonzero = (predicted != 0) & (observed != 0)
    directional = float(np.mean(np.sign(predicted[nonzero]) == np.sign(observed[nonzero])))
    result: dict[str, Any] = {
        "spearman": spearman,
        "kendall": kendall,
        "directional_agreement": directional,
        "directional_n": int(np.count_nonzero(nonzero)),
        "top_k": {},
    }
    for k in top_k:
        if k > len(targets):
            raise InputError(f"top-k={k} exceeds eligible target count")
        first = _top_targets(targets, observed, k)
        second = _top_targets(targets, predicted, k)
        overlap = len(first & second)
        result["top_k"][str(k)] = {
            "overlap": overlap,
            "jaccard": overlap / len(first | second),
        }
    return result


def _null_summary(
    values: list[float],
    observed: float,
    interval: list[float],
    *,
    center: float | None,
) -> dict[str, float]:
    null = np.asarray(values, dtype=float)
    lower, upper = np.quantile(null, interval)
    if center is None:
        center = float(np.mean(null))
    distance = abs(observed - center)
    p_value = (1 + int(np.count_nonzero(np.abs(null - center) >= distance))) / (len(null) + 1)
    return {
        "observed": observed,
        "null_center": center,
        "null_lower": float(lower),
        "null_upper": float(upper),
        "permutation_p": float(p_value),
    }


def _directional_null_center(observed: np.ndarray, predicted: np.ndarray) -> float:
    observed_positive = int(np.count_nonzero(observed > 0))
    observed_negative = int(np.count_nonzero(observed < 0))
    predicted_positive = int(np.count_nonzero(predicted > 0))
    predicted_negative = int(np.count_nonzero(predicted < 0))
    denominator = (observed_positive + observed_negative) * (
        predicted_positive + predicted_negative
    )
    if denominator == 0:
        raise InputError("directional agreement requires nonzero signs in both vectors")
    return (
        observed_positive * predicted_positive
        + observed_negative * predicted_negative
    ) / denominator


def evaluate_context(
    rows: list[dict[str, Any]], context: str, benchmark: dict[str, Any], *, seed: int
) -> dict[str, Any]:
    selected = [row for row in rows if row["context"] == context and row["analysis_eligible"]]
    targets = [row["target"] for row in selected]
    observed = np.asarray([row["observed_lfc"] for row in selected], dtype=float)
    predicted = np.asarray(
        [row["predicted_il2ra_regulator_score"] for row in selected], dtype=float
    )
    point = _point_metrics(targets, observed, predicted, benchmark["top_k"])
    null: dict[str, list[float]] = {
        "spearman": [],
        "kendall": [],
        "directional_agreement": [],
        **{f"top_{k}_overlap": [] for k in benchmark["top_k"]},
    }
    rng = np.random.default_rng(seed)
    for _ in range(benchmark["permutations"]):
        permuted = predicted[rng.permutation(len(predicted))]
        trial = _point_metrics(targets, observed, permuted, benchmark["top_k"])
        null["spearman"].append(trial["spearman"])
        null["kendall"].append(trial["kendall"])
        null["directional_agreement"].append(trial["directional_agreement"])
        for k in benchmark["top_k"]:
            null[f"top_{k}_overlap"].append(float(trial["top_k"][str(k)]["overlap"]))
    interval = benchmark["permutation_interval"]
    permutation = {
        "spearman": _null_summary(null["spearman"], point["spearman"], interval, center=0.0),
        "kendall": _null_summary(null["kendall"], point["kendall"], interval, center=0.0),
        "directional_agreement": _null_summary(
            null["directional_agreement"],
            point["directional_agreement"],
            interval,
            center=_directional_null_center(observed, predicted),
        ),
        "top_k": {},
    }
    for k in benchmark["top_k"]:
        values = np.asarray(null[f"top_{k}_overlap"], dtype=float)
        lower, upper = np.quantile(values, interval)
        observed_overlap = point["top_k"][str(k)]["overlap"]
        permutation["top_k"][str(k)] = {
            "observed_overlap": observed_overlap,
            "null_lower": float(lower),
            "null_upper": float(upper),
            "permutation_p": float(
                (1 + np.count_nonzero(values >= observed_overlap)) / (len(values) + 1)
            ),
        }
    return {
        "n_targets": len(targets),
        "permutation_seed": seed,
        "point": point,
        "target_label_permutation": permutation,
    }


def _check_attrition(
    screen: dict[str, dict[str, Any]], predictors: dict[str, dict[str, Any]], config: dict[str, Any]
) -> dict[str, int]:
    four = {target for target, item in screen.items() if item["four_guide_eligible"]}
    available = four & set(predictors)
    admitted = {target for target in available if predictors[target]["admitted"]}
    observed = {
        "screen_rows": len(screen),
        "four_guide_eligible": len(four),
        "generator_available_eligible": len(available),
        "analysis_eligible": len(admitted),
    }
    expected = config["expected"]
    for field, value in observed.items():
        if value != expected[field]:
            raise InputError(f"Arce benchmark attrition differs for {field}: {value} != {expected[field]}")
    return observed


def _format_csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".17g")
    return value


def render_predictions(rows: list[dict[str, Any]]) -> str:
    output = io.StringIO(newline="")
    fields = list(rows[0])
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: _format_csv_value(row[field]) for field in fields})
    return output.getvalue()


def run(
    config_path: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("schema_version") != "1.0.0":
        raise InputError("unsupported Arce external-validation config schema")
    payloads, arce_identity = verify_and_read_arce_members(config)
    screen = load_s1(payloads["screen_member"], config)
    activation = load_s14(payloads["activation_cells_member"], config)
    activation_reproduction = load_s8_and_check_reproduction(
        payloads["activation_summary_member"], activation, config
    )
    predictors, generator_identity = load_generator(config)
    attrition = _check_attrition(screen, predictors, config)
    rows = build_prediction_rows(screen, predictors, config)
    activation_rows = build_activation_guide_rows(activation, predictors, config)
    activation_robustness = summarize_activation_robustness(activation_rows, config)
    activation_available = sum(
        target in predictors for target in activation["regular_targets"]
    )
    activation_admitted = sum(
        target in predictors and predictors[target]["admitted"]
        for target in activation["regular_targets"]
    )
    if (
        activation_available != config["expected"]["activation_generator_available"]
        or activation_admitted != config["expected"]["activation_generator_admitted"]
    ):
        raise InputError("Arce S14/Zhu target attrition differs from the frozen contract")
    benchmark = config["benchmark"]
    context_metrics = {
        context: evaluate_context(
            rows,
            context,
            benchmark,
            seed=benchmark["permutation_seed"] + index,
        )
        for index, context in enumerate(benchmark["contexts"])
    }
    metadata = {
        "schema_version": "1.0.0",
        "generated_on": config["generated_on"],
        "benchmark": benchmark["name"],
        "claim_ceiling": config["claim_ceiling"],
        "config_sha256": sha256_file(config_path),
        "input_verification": {
            "arce": arce_identity,
            "zhu_generator": generator_identity,
        },
        "selection_contract": {
            "screen_eligibility": "exactly four S1 guides in every registered context",
            "generator_eligibility": "Rest row exists and ontarget_significant is true",
            "outcome_fields_used_for_selection": [],
            "rank_association": benchmark["rank_association"],
            "top_k_ranking": benchmark["top_k_ranking"],
            "predictor_orientation_multiplier": benchmark["predictor_multiplier"],
            "orientation_rationale": benchmark["orientation_rationale"],
        },
        "attrition": attrition,
        "contexts": context_metrics,
        "activation_data_quality": activation["data_quality"],
        "activation_summary_reproduction": activation_reproduction,
        "activation_attrition": {
            "regular_targets": len(activation["regular_targets"]),
            "zhu_rest_available": activation_available,
            "zhu_rest_admitted": activation_admitted,
            "role": "provenance annotation only; no IL2RA-to-activation.score association is computed",
        },
        "activation_robustness": activation_robustness,
        "limitations": [
            "S1 is an aggregate functional screen and does not expose donor-level effects.",
            "The Zhu predictor is donor-collapsed CRISPRi whereas Arce S1 is CRISPR-KO CD25/IL2RA screening.",
            "Permutation intervals quantify target-label exchangeability, not biological donor uncertainty.",
            "All context-by-metric permutation p-values are unadjusted exploratory diagnostics across correlated tests; they support no FWER- or FDR-controlled inference.",
            "S8 is a deterministic pooled-cell summary of S14, not independent replication; its cell-level tests are not used.",
            "The supplied activation.score has no frozen gene set, formula, normalization, or independence proof in the local tables.",
            "S14 is the authors' preselected 28-regulator panel; concordance is conditional on that panel and is neither genome-wide generality nor independent validation.",
            "Two donors permit descriptive concordance only, not donor-population inference or donor generalization.",
            "The Zhu IL2RA one-gene predictor and supplied global activation.score are endpoint-mismatched, so their correlation is intentionally not reported.",
        ],
        "permutation_inference": benchmark["multiplicity"],
    }
    return rows, activation_rows, metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--predictions", type=Path)
    parser.add_argument("--activation-guide-effects", type=Path)
    parser.add_argument("--metadata", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = args.config if args.config.is_absolute() else ROOT / args.config
    config = json.loads(config_path.read_text(encoding="utf-8"))
    prediction_path = args.predictions or ROOT / config["outputs"]["predictions"]
    activation_path = (
        args.activation_guide_effects
        or ROOT / config["outputs"]["activation_guide_effects"]
    )
    metadata_path = args.metadata or ROOT / config["outputs"]["metadata"]
    rows, activation_rows, metadata = run(config_path)
    prediction_text = render_predictions(rows)
    activation_text = render_predictions(activation_rows)
    metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    if args.write:
        prediction_path.parent.mkdir(parents=True, exist_ok=True)
        prediction_path.write_text(prediction_text, encoding="utf-8", newline="")
        activation_path.parent.mkdir(parents=True, exist_ok=True)
        activation_path.write_text(activation_text, encoding="utf-8", newline="")
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.write_text(metadata_text, encoding="utf-8")
    elif args.check:
        if prediction_path.read_text(encoding="utf-8") != prediction_text:
            raise AssertionError(f"prediction artifact differs: {prediction_path}")
        if activation_path.read_text(encoding="utf-8") != activation_text:
            raise AssertionError(f"activation artifact differs: {activation_path}")
        if metadata_path.read_text(encoding="utf-8") != metadata_text:
            raise AssertionError(f"metadata artifact differs: {metadata_path}")
    print(
        json.dumps(
            {
                "status": "PASS",
                "mode": "write" if args.write else "check" if args.check else "preview",
                "analysis_eligible": metadata["attrition"]["analysis_eligible"],
                "activation_guide_rows": len(activation_rows),
                "predictions": str(prediction_path),
                "activation_guide_effects": str(activation_path),
                "metadata": str(metadata_path),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
