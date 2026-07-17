import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pytest

from reachability import InputError
from scripts.run_arce_external_validation import (
    _directional_null_center,
    build_prediction_rows,
    evaluate_context,
    expected_s1_columns,
    load_generator,
    load_s1,
    parse_xlsx_table,
    predictor_identity,
    render_predictions,
    verify_and_read_screen_member,
)


CONTEXTS = ["Resting_Teff", "Stimulated_Teff", "Resting_Treg"]


def _column_name(index):
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(ord("A") + remainder) + result
    return result


def _xlsx_bytes(headers, rows, sheet_name="S1_all_screens_gene_summary"):
    def cell_xml(row_number, column, value):
        reference = f"{_column_name(column)}{row_number}"
        if isinstance(value, str):
            escaped = (
                value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            )
            return f'<c r="{reference}" t="inlineStr"><is><t>{escaped}</t></is></c>'
        return f'<c r="{reference}"><v>{value}</v></c>'

    xml_rows = []
    for row_number, values in enumerate([headers, *rows], start=1):
        cells = "".join(cell_xml(row_number, column, value) for column, value in enumerate(values))
        xml_rows.append(f'<row r="{row_number}">{cells}</row>')
    sheet_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>'
    )
    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<sheets><sheet name="{sheet_name}" sheetId="1" r:id="rId1"/></sheets></workbook>'
    )
    relationships = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/></Relationships>'
    )
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as book:
        book.writestr("xl/workbook.xml", workbook_xml)
        book.writestr("xl/_rels/workbook.xml.rels", relationships)
        book.writestr("xl/worksheets/sheet1.xml", sheet_xml)
    return output.getvalue()


def _benchmark_config(screen_rows):
    benchmark = {
        "contexts": CONTEXTS,
        "required_guides_per_context": 4,
        "observed_lfc_prefix": "pos|lfc.",
        "parallel_lfc_prefix": "neg|lfc.",
        "positive_fdr_prefix": "pos|fdr.",
        "negative_fdr_prefix": "neg|fdr.",
        "top_k": [1, 2],
        "permutations": 50,
        "permutation_interval": [0.025, 0.975],
    }
    return {
        "benchmark": benchmark,
        "dataset": {"screen_member": {"sheet": "S1_all_screens_gene_summary"}},
        "expected": {
            "screen_rows": screen_rows,
            "four_guide_eligible": screen_rows - 1,
        },
    }


def _screen_row(headers, target, guides, lfc):
    values = []
    for header in headers:
        if header == "id":
            values.append(target)
        elif header.startswith("num."):
            values.append(guides)
        elif "lfc." in header:
            values.append(lfc)
        elif "fdr." in header:
            values.append(0.01)
        elif "goodsgrna." in header:
            values.append(guides)
        elif "rank." in header:
            values.append(1)
        elif "p-value." in header:
            values.append(0.001)
        else:
            values.append(1.0)
    return values


def test_minimal_xlsx_parser_and_eligibility_ignore_outcomes():
    config = _benchmark_config(3)
    headers = expected_s1_columns(config["benchmark"])
    rows = [
        _screen_row(headers, "P1", 4, 1.0),
        _screen_row(headers, "P2", 3, -999.0),
        _screen_row(headers, "P3", 4, -2.0),
    ]
    payload = _xlsx_bytes(headers, rows)
    parsed_headers, parsed_rows = parse_xlsx_table(payload, "S1_all_screens_gene_summary")
    assert parsed_headers == headers
    assert len(parsed_rows) == 3
    first = load_s1(payload, config)
    rows[0] = _screen_row(headers, "P1", 4, -1e12)
    second = load_s1(_xlsx_bytes(headers, rows), config)
    assert {
        target: item["four_guide_eligible"] for target, item in first.items()
    } == {
        target: item["four_guide_eligible"] for target, item in second.items()
    }
    assert first["P2"]["four_guide_eligible"] is False


def test_s1_rejects_parallel_lfc_schema_disagreement():
    config = _benchmark_config(2)
    config["expected"]["four_guide_eligible"] = 2
    headers = expected_s1_columns(config["benchmark"])
    rows = [_screen_row(headers, "P1", 4, 1.0), _screen_row(headers, "P2", 4, 2.0)]
    rows[0][headers.index("neg|lfc.Resting_Teff")] = -1.0
    with pytest.raises(InputError, match="parallel LFC"):
        load_s1(_xlsx_bytes(headers, rows), config)


def test_archive_and_member_hashes_are_both_enforced(tmp_path, monkeypatch):
    member = _xlsx_bytes(["id"], [["P1"]])
    archive_buffer = io.BytesIO()
    member_path = "data_tables/S1_all_screens_gene_summary.xlsx"
    with zipfile.ZipFile(archive_buffer, "w") as archive:
        archive.writestr(member_path, member)
    archive_bytes = archive_buffer.getvalue()
    path = tmp_path / "arce.zip"
    path.write_bytes(archive_bytes)
    config = {
        "dataset": {
            "archive": {
                "path": "arce.zip",
                "bytes": len(archive_bytes),
                "md5": hashlib.md5(archive_bytes, usedforsecurity=False).hexdigest(),
                "sha256": hashlib.sha256(archive_bytes).hexdigest(),
            },
            "screen_member": {
                "path": member_path,
                "bytes": len(member),
                "sha256": hashlib.sha256(member).hexdigest(),
                "sheet": "S1_all_screens_gene_summary",
            },
        }
    }
    monkeypatch.setattr("scripts.run_arce_external_validation.ROOT", tmp_path)
    observed, identity = verify_and_read_screen_member(config)
    assert observed == member
    assert identity["archive_md5"] == config["dataset"]["archive"]["md5"]
    config["dataset"]["screen_member"]["sha256"] = "0" * 64
    with pytest.raises(InputError, match="member SHA-256"):
        verify_and_read_screen_member(config)


def _write_tiny_h5ad(path: Path):
    h5py = pytest.importorskip("h5py")
    string = h5py.string_dtype("utf-8")
    with h5py.File(path, "w") as handle:
        obs = handle.create_group("obs")
        condition = obs.create_group("culture_condition")
        condition.create_dataset("categories", data=np.asarray(["Rest", "Stim"], object), dtype=string)
        condition.create_dataset("codes", data=np.asarray([0, 0, 1], dtype=np.int8))
        perturbation = obs.create_group("target_contrast_gene_name")
        perturbation.create_dataset("categories", data=np.asarray(["P1", "P2", "P3"], object), dtype=string)
        perturbation.create_dataset("codes", data=np.asarray([0, 1, 2], dtype=np.int8))
        obs.create_dataset("ontarget_significant", data=np.asarray([True, False, True]))
        var = handle.create_group("var")
        var.create_dataset("gene_name", data=np.asarray(["IL2RA", "G2"], object), dtype=string)
        layers = handle.create_group("layers")
        layers.create_dataset("log_fc", data=np.asarray([[1.0, 0.0], [-2.0, 0.0], [3.0, 0.0]]))


def test_generator_projection_is_hash_bound_and_admission_is_explicit(tmp_path, monkeypatch):
    path = tmp_path / "tiny.h5ad"
    _write_tiny_h5ad(path)
    raw_predictors = {
        "P1": {"effect": 1.0, "admitted": True},
        "P2": {"effect": -2.0, "admitted": False},
    }
    config = {
        "generator": {
            "path": "tiny.h5ad",
            "bytes": path.stat().st_size,
            "condition": "Rest",
            "condition_field": "culture_condition",
            "perturbation_field": "target_contrast_gene_name",
            "admission_field": "ontarget_significant",
            "gene_field": "gene_name",
            "target_gene": "IL2RA",
            "layer": "log_fc",
            "predictor_sha256": predictor_identity(raw_predictors),
        },
        "benchmark": {"predictor_multiplier": 1.0},
        "expected": {
            "effect_rows": 3,
            "effect_genes": 2,
            "generator_condition_rows": 2,
            "generator_admitted_rows": 1,
        },
    }
    monkeypatch.setattr("scripts.run_arce_external_validation.ROOT", tmp_path)
    loaded, identity = load_generator(config)
    assert loaded == {
        "P1": {"raw_log_fc": 1.0, "regulator_score": 1.0, "admitted": True},
        "P2": {"raw_log_fc": -2.0, "regulator_score": -2.0, "admitted": False},
    }
    assert identity["predictor_sha256"] == predictor_identity(raw_predictors)
    config["generator"]["predictor_sha256"] = "0" * 64
    with pytest.raises(InputError, match="predictor identity"):
        load_generator(config)


def test_prediction_render_is_byte_deterministic():
    rows = [{"target": "P1", "score": 1.25, "eligible": True}]
    expected = "target,score,eligible\nP1,1.25,true\n"
    assert render_predictions(rows) == expected
    assert render_predictions(rows) == render_predictions(rows)


def test_directional_null_center_is_fixed_by_sign_margins():
    observed = np.asarray([1.0, 2.0, -1.0, -2.0])
    predicted = np.asarray([3.0, -1.0, -2.0, -4.0])
    # (2 observed-positive * 1 predicted-positive + 2 * 3 negatives) / 16
    assert _directional_null_center(observed, predicted) == pytest.approx(0.5)


def test_prediction_rows_preserve_all_exclusion_reasons():
    screen = {
        "P1": {"guide_counts": dict.fromkeys(CONTEXTS, 4), "four_guide_eligible": True, "outcomes": {c: {"lfc": 1.0, "positive_fdr": 0.1, "negative_fdr": 0.2} for c in CONTEXTS}},
        "P2": {"guide_counts": dict.fromkeys(CONTEXTS, 4), "four_guide_eligible": True, "outcomes": {c: {"lfc": -1.0, "positive_fdr": 0.1, "negative_fdr": 0.2} for c in CONTEXTS}},
        "P3": {"guide_counts": dict.fromkeys(CONTEXTS, 3), "four_guide_eligible": False, "outcomes": {c: {"lfc": 100.0, "positive_fdr": 0.1, "negative_fdr": 0.2} for c in CONTEXTS}},
    }
    predictors = {
        "P1": {"raw_log_fc": -0.5, "regulator_score": 0.5, "admitted": True},
        "P2": {"raw_log_fc": 0.5, "regulator_score": -0.5, "admitted": False},
        "P3": {"raw_log_fc": -100.0, "regulator_score": 100.0, "admitted": True},
    }
    rows = build_prediction_rows(screen, predictors, {"benchmark": {"contexts": CONTEXTS}})
    rest = {row["target"]: row for row in rows if row["context"] == CONTEXTS[0]}
    assert rest["P1"]["analysis_eligible"] is True
    assert rest["P2"]["exclusion_reason"] == "generator_not_admitted"
    assert rest["P3"]["exclusion_reason"] == "guide_count_not_four_all_contexts"


def test_target_label_permutation_metrics_are_deterministic():
    rows = []
    for index in range(8):
        rows.append(
            {
                "target": f"P{index}",
                "context": "C",
                "analysis_eligible": True,
                "observed_lfc": float(index + 1),
                "predicted_il2ra_regulator_score": float(index + 1),
            }
        )
    benchmark = {
        "top_k": [2, 4],
        "permutations": 100,
        "permutation_interval": [0.025, 0.975],
    }
    first = evaluate_context(rows, "C", benchmark, seed=19)
    second = evaluate_context(rows, "C", benchmark, seed=19)
    assert first == second
    assert first["point"]["spearman"] == pytest.approx(1.0)
    assert first["point"]["top_k"]["2"]["overlap"] == 2
