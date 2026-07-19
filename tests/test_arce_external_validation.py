import hashlib
import io
from pathlib import Path
import zipfile

import numpy as np
import pytest

from reachability import InputError
from scripts.run_arce_external_validation import (
    S14_COLUMNS,
    _directional_null_center,
    build_activation_guide_rows,
    build_prediction_rows,
    evaluate_context,
    expected_s1_columns,
    load_generator,
    load_s1,
    load_s14,
    load_s8_and_check_reproduction,
    parse_xlsx_table,
    predictor_identity,
    render_predictions,
    summarize_activation_robustness,
    verify_and_read_arce_members,
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
        if value is None:
            return ""
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
    observed, identity = verify_and_read_arce_members(config)
    assert observed["screen_member"] == member
    assert identity["archive_md5"] == config["dataset"]["archive"]["md5"]
    assert identity["members"]["screen_member"]["sha256"] == config["dataset"]["screen_member"]["sha256"]
    config["dataset"]["screen_member"]["sha256"] = "0" * 64
    with pytest.raises(InputError, match="member SHA-256"):
        verify_and_read_arce_members(config)


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


def _tiny_activation_config(rows):
    return {
        "dataset": {
            "activation_cells_member": {"sheet": "S14"},
            "activation_summary_member": {"sheet": "S8"},
        },
        "activation_benchmark": {
            "contexts": ["Resting-Teff"],
            "donors": ["A"],
            "control_target": "Non-Targeting",
            "control_guides": ["NTC1", "NTC2"],
            "regular_guides_per_target": 2,
            "reproduction_tolerance": 1e-10,
        },
        "expected": {
            "activation_rows": rows,
            "activation_unique_cells": rows,
            "activation_targets": 2,
            "activation_regular_targets": 1,
            "activation_guides": 4,
            "activation_guide_groups": 4,
            "activation_target_groups": 2,
            "activation_summary_rows": 2,
            "activation_context_counts": {"Resting-Teff": rows},
            "activation_donor_counts": {"A": rows},
        },
    }


def _s14_row(cell, guide, target, score):
    values = dict.fromkeys(S14_COLUMNS, 0)
    values.update(
        {
            "cell": cell,
            "n_sgrna_features": 1,
            "sgrna": guide,
            "sg_target": target,
            "has_sgrna": 1,
            "donor": "A",
            "HTO_maxID": "Resting-Teff",
            "HTO_classification.global": "Singlet",
            "activation.score": score,
        }
    )
    return [values[field] for field in S14_COLUMNS]


def test_s14_streaming_contract_and_equal_guide_weighted_control():
    rows = [
        _s14_row("c1", "NTC1", "Non-Targeting", 0.0),
        _s14_row("c2", "NTC2", "Non-Targeting", 100.0),
        _s14_row("c3", "NTC2", "Non-Targeting", 100.0),
        _s14_row("c4", "NTC2", "Non-Targeting", 100.0),
        _s14_row("c5", "g1", "P1", 60.0),
        _s14_row("c6", "g2", "P1", 70.0),
    ]
    config = _tiny_activation_config(len(rows))
    activation = load_s14(_xlsx_bytes(S14_COLUMNS, rows, "S14"), config)
    output = build_activation_guide_rows(
        activation,
        {"P1": {"admitted": True}},
        config,
    )
    assert activation["data_quality"]["guide_groups"] == 4
    assert {row["ntc_baseline_median_of_guide_medians"] for row in output} == {50.0}
    assert [row["supplied_activation_score_delta"] for row in output] == [10.0, 20.0]
    assert all(row["n_cells"] == 1 for row in output)


def test_s14_rejects_duplicate_cells_and_non_singlets():
    rows = [
        _s14_row("same", "NTC1", "Non-Targeting", 0.0),
        _s14_row("same", "NTC2", "Non-Targeting", 1.0),
        _s14_row("c3", "g1", "P1", 2.0),
        _s14_row("c4", "g2", "P1", 3.0),
    ]
    config = _tiny_activation_config(len(rows))
    with pytest.raises(InputError, match="duplicated"):
        load_s14(_xlsx_bytes(S14_COLUMNS, rows, "S14"), config)
    rows[1] = _s14_row("c2", "NTC2", "Non-Targeting", 1.0)
    rows[1][S14_COLUMNS.index("HTO_classification.global")] = "Doublet"
    with pytest.raises(InputError, match="not an HTO singlet"):
        load_s14(_xlsx_bytes(S14_COLUMNS, rows, "S14"), config)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("activation.score", None, "missing registered field"),
        ("activation.score", float("nan"), "not finite"),
        ("has_sgrna", 0, "has_sgrna is not true"),
        ("n_sgrna_features", 2, "expected one sgRNA feature"),
        ("donor", "C", "unknown donor or context"),
        ("HTO_maxID", "Unknown", "unknown donor or context"),
    ],
)
def test_s14_registered_qc_fields_fail_closed(field, value, message):
    rows = [
        _s14_row("c1", "NTC1", "Non-Targeting", 0.0),
        _s14_row("c2", "NTC2", "Non-Targeting", 1.0),
        _s14_row("c3", "g1", "P1", 2.0),
        _s14_row("c4", "g2", "P1", 3.0),
    ]
    rows[0][S14_COLUMNS.index(field)] = value
    config = _tiny_activation_config(len(rows))
    with pytest.raises(InputError, match=message):
        load_s14(_xlsx_bytes(S14_COLUMNS, rows, "S14"), config)


def test_s14_control_guide_identity_is_frozen():
    rows = [
        _s14_row("c1", "WRONG", "Non-Targeting", 0.0),
        _s14_row("c2", "NTC2", "Non-Targeting", 1.0),
        _s14_row("c3", "g1", "P1", 2.0),
        _s14_row("c4", "g2", "P1", 3.0),
    ]
    config = _tiny_activation_config(len(rows))
    with pytest.raises(InputError, match="guide identities differ"):
        load_s14(_xlsx_bytes(S14_COLUMNS, rows, "S14"), config)


def test_s8_is_only_an_exact_aggregate_reproduction_gate():
    rows = [
        _s14_row("c1", "NTC1", "Non-Targeting", 0.0),
        _s14_row("c2", "NTC2", "Non-Targeting", 2.0),
        _s14_row("c3", "g1", "P1", 4.0),
        _s14_row("c4", "g2", "P1", 6.0),
    ]
    config = _tiny_activation_config(len(rows))
    activation = load_s14(_xlsx_bytes(S14_COLUMNS, rows, "S14"), config)
    s8_headers = [
        "sg_target", "median.activation.score", "mean.activation.score", "HTO_maxID",
        "statistic", "p.value", "method", "alternative", "stars", "padj",
    ]
    s8_rows = [
        ["Non-Targeting", 1.0, 1.0, "Resting-Teff", None, None, None, None, None, None],
        ["P1", 5.0, 5.0, "Resting-Teff", 1.0, 0.1, "Wilcoxon", "two.sided", None, 2.8],
    ]
    report = load_s8_and_check_reproduction(
        _xlsx_bytes(s8_headers, s8_rows, "S8"), activation, config
    )
    assert report["status"] == "PASS"
    assert report["published_p_values_used_for_inference"] is False
    s8_rows[1][1] = 5.1
    with pytest.raises(InputError, match="do not reproduce"):
        load_s8_and_check_reproduction(
            _xlsx_bytes(s8_headers, s8_rows, "S8"), activation, config
        )


def test_robustness_reducer_uses_equal_weights_and_retains_low_support():
    config = {
        "activation_benchmark": {
            "contexts": ["C"],
            "donors": ["A", "B"],
            "estimand": "registered test estimand",
            "inference": "descriptive only",
            "panel_provenance": "preselected test panel",
            "score_claim_ceiling": "supplied score only",
        }
    }
    guide_values = {
        ("T1", "A"): (1.0, 1.0),
        ("T1", "B"): (2.0, 2.0),
        ("T2", "A"): (2.0, 2.0),
        ("T2", "B"): (9.0, -1.0),
        ("T3", "A"): (3.0, 3.0),
        ("T3", "B"): (6.0, 6.0),
    }
    rows = []
    for (target, donor), values in guide_values.items():
        for index, value in enumerate(values, start=1):
            low_support = target == "T2" and donor == "B" and index == 2
            rows.append(
                {
                    "target": target,
                    "guide": f"{target}_g{index}",
                    "donor": donor,
                    "context": "C",
                    "n_cells": 8 if low_support else 200,
                    "supplied_activation_score_delta": value,
                    "support_flag_lt_20_cells": low_support,
                }
            )
    report = summarize_activation_robustness(rows, config)
    context = report["contexts"]["C"]
    assert context["targets"] == 3
    assert context["donor_rank_concordance"]["spearman"] == pytest.approx(1.0)
    assert context["donor_rank_concordance"]["kendall"] == pytest.approx(1.0)
    assert context["guide_pair_sign_agreement_n"] == 5
    assert context["guide_pair_comparisons"] == 6
    assert context["all_two_guides_two_donors_same_nonzero_sign_n"] == 2
    assert context["all_two_guides_two_donors_same_nonzero_sign_fraction"] == pytest.approx(2 / 3)
    assert report["strata_below_20_cells"] == 1
    assert report["all_strata_retained"] is True

    # Cell counts never weight the fixed guide contrasts or remove a sparse stratum.
    duplicated_cells = [{**row, "n_cells": row["n_cells"] * 1000} for row in rows]
    duplicated = summarize_activation_robustness(duplicated_cells, config)
    assert duplicated["contexts"] == report["contexts"]
