import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from reachability import InputError
from scripts.run_goudy_combination_validation import (
    _canonical_json,
    _fold_gene_mask,
    fit_lodo,
    load_counts,
    parse_author_key,
    prepare_profiles,
    run,
    validate_author_key,
)


COMPONENTS = ("FAS", "RC3H1", "SUV39H1")
UNRELATED = ("PTPN2", "RASA2", "MED12")
GENES = (
    ("ENSG00000000001.1", "FAS"),
    ("ENSG00000000002.2", "RC3H1"),
    ("ENSG00000000003.3", "SUV39H1"),
    ("ENSG00000000004.1", "GENE1"),
    ("ENSG00000000005.1", "GENE2"),
    ("ENSG00000000006.1", "GENE3"),
    ("ENSG00000000007.1", "GENE4"),
)


def _identity(path: Path) -> dict[str, object]:
    payload = path.read_bytes()
    return {
        "path": str(path),
        "file": path.name,
        "url": f"https://example.invalid/{path.name}",
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "md5": hashlib.md5(payload).hexdigest(),  # noqa: S324 - fixture identity
    }


def _sample(
    serial: int,
    donor: str,
    role: str,
    genotype: str,
    *,
    conflict: bool = False,
    title_has_day: bool = True,
) -> dict[str, str]:
    title = f"CRISPRoff_{role}_Day7_{donor}" if title_has_day else f"CRISPRoff_{role}_{donor}"
    return {
        "gsm": f"GSM{1000 + serial}",
        "column": f"sample_{serial}",
        "title": title,
        "genotype": genotype,
        "time": "Day 18" if conflict else "Day 7",
    }


def _make_fixture(
    directory: Path,
    *,
    duplicate_gene: bool = False,
    bad_schema: bool = False,
    target_multiplier: float = 1.0,
    omit_soft_gsm: str | None = None,
    single_experiment_id: str = "FixtureSingleExperiment",
    multiplex_experiment_id: str = "FixtureMultiplexExperiment",
) -> tuple[dict, dict[str, dict[str, np.ndarray]]]:
    directory.mkdir(parents=True, exist_ok=True)
    donors: dict[str, dict] = {}
    serial = 1
    conflicts: list[str] = []
    for donor_index in range(1, 5):
        donor = f"D{donor_index}"
        donor_spec = {"component": {}, "unrelated": {}, "controls": {}}
        for role in COMPONENTS:
            spec = _sample(
                serial,
                donor,
                role,
                f"{role}_KD",
                conflict=donor == "D1",
            )
            donor_spec["component"][role] = spec
            if donor == "D1":
                conflicts.append(spec["gsm"])
            serial += 1
        for role in UNRELATED:
            spec = _sample(
                serial,
                donor,
                role,
                f"{role}_KD",
                conflict=donor == "D1",
            )
            donor_spec["unrelated"][role] = spec
            if donor == "D1":
                conflicts.append(spec["gsm"])
            serial += 1
        for role in ("AAVS1_Guide1", "AAVS1_Guide2"):
            spec = _sample(serial, donor, role, "WT", conflict=donor == "D1")
            donor_spec["controls"][role] = spec
            if donor == "D1":
                conflicts.append(spec["gsm"])
            serial += 1
        donor_spec["controls"]["NTC_multiplexing"] = _sample(
            serial, donor, "NTC_multiplexing", "WT", title_has_day=False
        )
        serial += 1
        donor_spec["triple"] = _sample(
            serial,
            donor,
            "FAS_SUV39H1_RC3H1",
            "FAS_SUV39H1_RC3H1_KD",
            title_has_day=False,
        )
        serial += 1
        donors[donor] = donor_spec

    all_specs = []
    for donor_spec in donors.values():
        for group in ("component", "unrelated", "controls"):
            all_specs.extend(donor_spec[group].values())
        all_specs.append(donor_spec["triple"])
    columns = [spec["column"] for spec in all_specs]
    values: dict[str, dict[str, np.ndarray]] = {}
    for donor_index, (donor, donor_spec) in enumerate(donors.items(), start=1):
        baseline = np.asarray(
            [100 + 8 * donor_index + 3 * gene_index for gene_index in range(len(GENES))],
            dtype=float,
        )
        values[donor] = {}
        for spec in donor_spec["controls"].values():
            values[donor][spec["column"]] = baseline.copy()
        for component_index, role in enumerate(COMPONENTS):
            profile = baseline.copy()
            profile[component_index] *= 0.20 * target_multiplier
            profile[3 + component_index] *= 1.5 + 0.1 * donor_index
            values[donor][donor_spec["component"][role]["column"]] = profile
        for unrelated_index, role in enumerate(UNRELATED):
            profile = baseline.copy()
            profile[3 + unrelated_index] *= 0.65
            values[donor][donor_spec["unrelated"][role]["column"]] = profile
        triple = baseline.copy()
        triple[:3] *= 0.12 * target_multiplier
        triple[3:6] *= 1.45 + 0.1 * donor_index
        triple[6] *= 2.2  # combination-specific residual keeps the cone non-trivial
        values[donor][donor_spec["triple"]["column"]] = triple

    by_column = {
        column: profile
        for donor_values in values.values()
        for column, profile in donor_values.items()
    }
    counts_path = directory / "tiny_counts.csv.gz"
    header = ["gene_id", *columns]
    if not bad_schema:
        header.append("gene_symbol")
    with gzip.open(counts_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for index, (gene_id, symbol) in enumerate(GENES):
            if duplicate_gene and index == 1:
                gene_id = GENES[0][0]
            row = [gene_id, *(by_column[column][index] for column in columns)]
            if not bad_schema:
                row.append(symbol)
            writer.writerow(row)

    soft_path = directory / "tiny_family.soft.gz"
    with gzip.open(soft_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("^SERIES = TEST1\n")
        for spec in all_specs:
            if spec["gsm"] == omit_soft_gsm:
                continue
            handle.write(f"^SAMPLE = {spec['gsm']}\n")
            handle.write(f"!Sample_title = {spec['title']}\n")
            handle.write(f"!Sample_geo_accession = {spec['gsm']}\n")
            handle.write("!Sample_characteristics_ch1 = cell type: Bulk T cells\n")
            handle.write(f"!Sample_characteristics_ch1 = genotype: {spec['genotype']}\n")
            handle.write("!Sample_characteristics_ch1 = treatment: CRISPRoff\n")
            handle.write(f"!Sample_characteristics_ch1 = time: {spec['time']}\n")
            handle.write(
                "!Sample_description = Library name: "
                f"{spec['column'].removeprefix('sample_')}\n"
            )
            handle.write("!Sample_description = tiny_counts.csv\n")

    author_key_path = directory / "rna_seq_meta_key.csv"
    author_rows: list[dict[str, str]] = []
    absent_declared: list[str] = []
    for donor_index, (donor, donor_spec) in enumerate(donors.items(), start=1):
        role_specs = [
            *(("component", role, donor_spec["component"][role]) for role in COMPONENTS),
            *(("unrelated", role, donor_spec["unrelated"][role]) for role in UNRELATED),
            *(
                ("controls", role, donor_spec["controls"][role])
                for role in ("AAVS1_Guide1", "AAVS1_Guide2")
            ),
        ]
        if donor_index <= 2:
            role_specs.extend(
                [
                    (
                        "controls",
                        "NTC_multiplexing",
                        donor_spec["controls"]["NTC_multiplexing"],
                    ),
                    ("triple", "triple", donor_spec["triple"]),
                ]
            )
        else:
            absent_declared.extend(
                [
                    donor_spec["controls"]["NTC_multiplexing"]["column"],
                    donor_spec["triple"]["column"],
                ]
            )
        for group, role, spec in role_specs:
            if group in {"component", "unrelated"}:
                target = role
            elif role.startswith("AAVS1_"):
                target = "AAVS1"
            elif role == "NTC_multiplexing":
                target = "NTC"
            else:
                target = "FAS, SUV39H1, RC3H1"
            author_rows.append(
                {
                    "sample": spec["column"],
                    "Experiment ID": (
                        multiplex_experiment_id
                        if group == "triple" or role == "NTC_multiplexing"
                        else single_experiment_id
                    ),
                    "Experimental Day Collected": "7",
                    "Donor": f"Donor {donor_index}",
                    "Cas": "CRISPRoff",
                    "Target": target,
                    "sgRNA notes": (
                        f"Frozen guide for {role}" if group == "component" else ""
                    ),
                }
            )
    author_only = ["sample_9001", "sample_9002", "sample_9003", "sample_9004"]
    for index, sample in enumerate(author_only, start=1):
        author_rows.append(
            {
                "sample": sample,
                "Experiment ID": multiplex_experiment_id,
                "Experimental Day Collected": "7",
                "Donor": f"Donor {1 if index <= 2 else 2}",
                "Cas": "CRISPRoff",
                "Target": "FAS, SUV39H1, RC3H1, MED12",
                "sgRNA notes": "",
            }
        )
    with author_key_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample",
                "Experiment ID",
                "Experimental Day Collected",
                "Donor",
                "Cas",
                "Target",
                "sgRNA notes",
            ],
        )
        writer.writeheader()
        writer.writerows(author_rows)

    config = {
        "schema_version": 2,
        "benchmark_id": "tiny_goudy",
        "source": {"accession": "TEST1", "citation": "synthetic"},
        "inputs": {
            "counts": _identity(counts_path),
            "soft": _identity(soft_path),
            "author_key": _identity(author_key_path),
        },
        "expected": {
            "soft_samples": len(all_specs) - int(omit_soft_gsm is not None),
            "count_samples": len(columns) - int(bad_schema),
            "gene_rows": len(GENES),
            "analysis_genes_before_mask": len(GENES),
            "analysis_genes_after_mask": len(GENES) - 3,
            "analysis_genes_before_mask_by_threshold": {
                key: len(GENES) for key in ("0.5", "1", "2", "5", "10")
            },
            "author_key_samples": len(author_rows),
            "author_key_count_overlap": len(author_rows) - len(author_only),
            "author_key_rows_absent_from_counts": author_only,
            "declared_rows_absent_from_author_key": absent_declared,
            "declared_analysis_samples": len(all_specs),
            "donors": 4,
            "controls_per_donor": 3,
        },
        "preprocessing": {
            "deposited_scale": "synthetic linear CPM",
            "transform": "log2(CPM + 1)",
            "gene_id_regex": r"^ENSG[0-9]{11}\.[0-9]+$",
            "control_mean_cpm_min": 1.0,
            "control_mean_cpm_sensitivity_thresholds": [0.5, 1, 2, 5, 10],
            "target_mask": [
                {"symbol": symbol, "gene_id": gene_id}
                for gene_id, symbol in GENES[:3]
            ],
        },
        "models": {
            "component_order": list(COMPONENTS),
            "unrelated_order": list(UNRELATED),
        },
        "metadata_conflict": {
            "status": "FLAGGED",
            "affected_gsms": conflicts,
            "computational_resolution": "Day7",
        },
        "claim_contract": {
            "modality": "CRISPRoff",
            "constituent_relationship": "target-matched constituent singles",
            "triple_guide_identity": "UNRESOLVED",
            "claim_ceiling": "descriptive synthetic cross-experiment alignment",
            "not_claimed": ["interaction", "prospective validation"],
        },
        "donors": donors,
    }
    return config, values


def _walk_keys(value):
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from _walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_keys(child)


def _prepare_fixture(config):
    counts = load_counts(Path(config["inputs"]["counts"]["path"]), config)
    author_key = parse_author_key(Path(config["inputs"]["author_key"]["path"]))
    author_contract = validate_author_key(author_key, counts, config)
    return prepare_profiles(counts, config, author_contract)


def test_tiny_gzip_fixture_is_deterministic_and_reports_required_contracts(tmp_path):
    config, _ = _make_fixture(tmp_path)
    first = run(config)
    second = run(config)
    assert first == second
    assert _canonical_json(first) == _canonical_json(second)
    assert first["execution_status"] == "PASS"
    assert first["geometric_model_status"] == "FAILS_DECLARED_GEOMETRIC_MODEL"
    assert first["biological_interpretation_status"] == (
        "INCONCLUSIVE_CROSS_EXPERIMENT_CONFOUNDING_LOW_RELIABILITY"
    )
    assert "passed" not in first and "status" not in first
    assert first["metadata_conflict"]["status"] == "FLAGGED"
    assert first["preprocessing"]["observed_transcriptome_genes_after_mask"] == 4
    assert "p_value" not in set(_walk_keys(first))
    on_target = first["on_target_controls"]["saturated_calibration_diagnostic"]
    assert on_target["warning_status"] == (
        "SATURATED_3X3_CALIBRATION_DIAGNOSTIC_NOT_VALIDATION"
    )
    assert set(on_target["per_donor"]) == set(config["donors"])
    for donor in config["donors"]:
        models = first["donor_results"][donor]["models"]
        assert {
            "component_cone",
            "equal_sum",
            "best_single_oracle",
            "unrelated_cone",
            "zero",
            "signed_least_squares",
            "on_target_calibrated_component_cone",
            "lodo_component_cone",
            "lodo_training_selected_best_single",
            "training_donor_mean_triple",
        } == set(models)
        geometry = models["component_cone"]["geometry"]
        assert geometry["strict_status"] in {"inside_tolerance", "outside_model_cone"}
        assert (geometry["certificate"] is None) == (
            geometry["strict_status"] == "inside_tolerance"
        )
        assert set(on_target["per_donor"][donor]["models"]) == {
            "component_cone",
            "equal_sum",
            "best_single_oracle",
        }
    assert set(first["filter_sensitivity"]["per_threshold"]) == {
        "0.5",
        "1",
        "2",
        "5",
        "10",
    }
    assert first["module_error"]["status"] == (
        "UNAVAILABLE_NO_PREREGISTERED_MODULE_SET"
    )
    assert first["additive_residual"]["per_donor"]


def test_author_crosswalk_reports_exact_omissions_and_unresolved_identity(tmp_path):
    config, _ = _make_fixture(
        tmp_path,
        single_experiment_id="AuthorSingles42",
        multiplex_experiment_id="AuthorMultiplex99",
    )
    report = run(config)
    author = report["sample_contract"]["author_key"]
    assert author["declared_analysis_role_count"] == 40
    assert author["author_key_confirmed_declared_role_count"] == 36
    assert author["author_key_rows_absent_from_count_matrix"] == [
        "sample_9001",
        "sample_9002",
        "sample_9003",
        "sample_9004",
    ]
    assert author["declared_rows_absent_from_author_key"] == sorted(
        config["expected"]["declared_rows_absent_from_author_key"]
    )
    assert author["experiment_role_summary"][
        "constituent_singles_unrelated_singles_and_aavs1"
    ]["experiment_ids"] == ["AuthorSingles42"]
    assert author["experiment_role_summary"]["multiplex_ntc_and_triple"][
        "experiment_ids"
    ] == ["AuthorMultiplex99"]
    assert author["experiment_role_summary"]["multiplex_ntc_and_triple"][
        "unresolved_donors"
    ] == ["D3", "D4"]
    assert author["triple_guide_identity"] == "UNRESOLVED"
    assert author["constituent_relationship"] == (
        "target-matched constituent singles"
    )
    assert author["same_guide_match_claimed"] is False
    assert report["preprocessing"]["control_pairing"]["D3"][
        "triple_role_experiment_id"
    ] == "UNRESOLVED"


def test_runner_has_no_hardcoded_experiment_ids():
    script = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "run_goudy_combination_validation.py"
    ).read_text(encoding="utf-8")
    assert "Co065" not in script
    assert "Co066" not in script


def test_hash_gate_precedes_gzip_parse(tmp_path):
    config, _ = _make_fixture(tmp_path)
    soft_path = Path(config["inputs"]["soft"]["path"])
    payload = bytearray(soft_path.read_bytes())
    payload[len(payload) // 2] ^= 0x01
    soft_path.write_bytes(payload)
    with pytest.raises(InputError, match="SHA-256"):
        run(config)


def test_author_key_hash_gate_precedes_csv_parse(tmp_path):
    config, _ = _make_fixture(tmp_path)
    key_path = Path(config["inputs"]["author_key"]["path"])
    key_path.write_bytes(b"not,the,registered,key\n")
    with pytest.raises(InputError, match="byte length"):
        run(config)


def test_rejects_count_schema_without_gene_symbol(tmp_path):
    config, _ = _make_fixture(tmp_path, bad_schema=True)
    with pytest.raises(InputError, match="schema"):
        run(config)


def test_rejects_missing_declared_soft_sample(tmp_path):
    probe, _ = _make_fixture(tmp_path / "probe")
    missing = probe["donors"]["D2"]["triple"]["gsm"]
    config, _ = _make_fixture(tmp_path / "missing", omit_soft_gsm=missing)
    with pytest.raises(InputError, match="missing declared sample"):
        run(config)


def test_rejects_duplicate_ensembl_version_rows(tmp_path):
    config, _ = _make_fixture(tmp_path, duplicate_gene=True)
    with pytest.raises(InputError, match="duplicate gene ID"):
        run(config)


def test_target_mask_separates_transcriptome_metrics_from_knockdown_controls(tmp_path):
    baseline_config, _ = _make_fixture(tmp_path / "baseline", target_multiplier=1.0)
    changed_config, _ = _make_fixture(tmp_path / "changed", target_multiplier=0.25)
    baseline = run(baseline_config)
    changed = run(changed_config)
    for donor in baseline["donor_results"]:
        baseline_models = baseline["donor_results"][donor]["models"]
        changed_models = changed["donor_results"][donor]["models"]
        for model in set(baseline_models) - {"on_target_calibrated_component_cone"}:
            assert baseline_models[model] == changed_models[model]
        assert (
            baseline_models["on_target_calibrated_component_cone"]
            != changed_models["on_target_calibrated_component_cone"]
        )
    assert baseline["on_target_controls"] != changed["on_target_controls"]


def test_effects_use_role_specific_controls_from_the_same_donor(tmp_path):
    config, raw_values = _make_fixture(tmp_path)
    prepared = _prepare_fixture(config)
    gene_index = 3  # GENE1 is the first coordinate after masking the three targets.
    sample_column = config["donors"]["D2"]["component"]["FAS"]["column"]
    single_controls = [
        raw_values["D2"][config["donors"]["D2"]["controls"][role]["column"]][
            gene_index
        ]
        for role in ("AAVS1_Guide1", "AAVS1_Guide2")
    ]
    expected = np.log2(raw_values["D2"][sample_column][gene_index] + 1.0) - np.mean(
        np.log2(np.asarray(single_controls) + 1.0)
    )
    assert prepared["profiles"]["D2"]["component"][0, 0] == pytest.approx(expected)
    triple_column = config["donors"]["D2"]["triple"]["column"]
    ntc_column = config["donors"]["D2"]["controls"]["NTC_multiplexing"]["column"]
    expected_triple = np.log2(
        raw_values["D2"][triple_column][gene_index] + 1.0
    ) - np.log2(raw_values["D2"][ntc_column][gene_index] + 1.0)
    assert prepared["profiles"]["D2"]["triple"][0] == pytest.approx(
        expected_triple
    )
    assert prepared["control_pairing"]["D2"] == {
        "constituent_relationship": "target-matched constituent singles",
        "single_role_experiment_ids": ["FixtureSingleExperiment"],
        "single_control_columns": [
            config["donors"]["D2"]["controls"]["AAVS1_Guide1"]["column"],
            config["donors"]["D2"]["controls"]["AAVS1_Guide2"]["column"],
        ],
        "triple_role_experiment_id": "FixtureMultiplexExperiment",
        "triple_control_experiment_id": "FixtureMultiplexExperiment",
        "triple_control_column": ntc_column,
        "cross_experiment_confounding": True,
    }


def test_lodo_coefficients_cannot_read_held_donor_target():
    components = {
        "D1": np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        "D2": np.asarray([[2.0, 0.0], [0.0, 1.0], [1.0, 2.0]]),
        "D3": np.asarray([[1.0, 1.0], [0.0, 2.0], [2.0, 0.0]]),
        "D4": np.asarray([[3.0, 0.0], [0.0, 3.0], [1.0, 1.0]]),
    }
    triples = {
        "D1": np.asarray([1.0, 0.5]),
        "D2": np.asarray([2.0, 1.0]),
        "D3": np.asarray([0.5, 2.0]),
        "D4": np.asarray([1.0, 1.0]),
    }
    first = fit_lodo(components, triples, "D4")
    changed = dict(triples)
    changed["D4"] = np.asarray([-1000.0, 700.0])
    second = fit_lodo(components, changed, "D4")
    np.testing.assert_array_equal(first["coefficients"], second["coefficients"])
    np.testing.assert_array_equal(first["prediction"], second["prediction"])
    assert first["training_selected_best_single"]["selected"] == second[
        "training_selected_best_single"
    ]["selected"]
    assert first["training_selected_best_single"]["nonnegative_scale"] == second[
        "training_selected_best_single"
    ]["nonnegative_scale"]
    np.testing.assert_array_equal(
        first["training_selected_best_single"]["prediction"],
        second["training_selected_best_single"]["prediction"],
    )
    assert first["training_donors"] == ("D1", "D2", "D3")


def test_lodo_filter_uses_training_donor_controls_only(tmp_path):
    config, _ = _make_fixture(tmp_path)
    prepared = _prepare_fixture(config)
    training = ("D1", "D2", "D3")
    baseline = _fold_gene_mask(prepared, training, threshold=1.0)
    changed = dict(prepared)
    changed["cpm"] = prepared["cpm"].copy()
    changed["cpm"][:, prepared["control_indices_by_donor"]["D4"]] = 0.0
    np.testing.assert_array_equal(
        baseline, _fold_gene_mask(changed, training, threshold=1.0)
    )
    report = run(config)
    lodo = report["donor_results"]["D4"]["models"]["lodo_component_cone"]
    comparator = report["donor_results"]["D4"]["models"][
        "lodo_training_selected_best_single"
    ]
    assert lodo["gene_filter"]["held_donor_controls_used"] is False
    assert lodo["design_status"] == "CONDITIONAL_TRANSDUCTIVE"
    assert comparator["selection_contract"].startswith("identity and scale selected")


def test_filter_reliability_metric_and_saturated_contracts(tmp_path):
    config, _ = _make_fixture(tmp_path)
    report = run(config)
    for threshold in report["filter_sensitivity"]["per_threshold"].values():
        assert threshold["retained_transcriptome_genes"] == len(GENES) - 3
        assert set(threshold["donor_ranges"]) == {"component_cone", "equal_sum"}
        assert "median component-cone cosine/nRMSE=" in threshold["conclusion"]
    pairwise = report["reliability_diagnostics"]["pairwise_donor_cosines"]
    assert pairwise["triple"]["pair_count"] == 6
    assert all(value["pair_count"] == 6 for value in pairwise["components"].values())
    assert set(report["reliability_diagnostics"]["aavs1_guide_replicate_noise"][
        "per_donor"
    ]) == set(config["donors"])
    assert "residual_fraction" not in set(_walk_keys(report))
    for donor in report["donor_results"].values():
        for model in donor["models"].values():
            assert set(model["metrics"]) == {
                "cosine",
                "normalized_rmse",
                "norm_ratio",
            }
    saturated = report["on_target_controls"]["saturated_calibration_diagnostic"]
    assert "not validation" in saturated["warning_status"].lower().replace("_", " ")


def test_metadata_conflict_set_is_fail_closed(tmp_path):
    config, _ = _make_fixture(tmp_path)
    config["metadata_conflict"]["affected_gsms"] = []
    with pytest.raises(InputError, match="metadata conflict set differs"):
        run(config)


def test_cli_write_then_exact_check(tmp_path):
    config, _ = _make_fixture(tmp_path / "inputs")
    config_path = tmp_path / "config.json"
    report_path = tmp_path / "report.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    script = Path(__file__).resolve().parents[1] / "scripts" / "run_goudy_combination_validation.py"
    subprocess.run(
        [sys.executable, str(script), "--config", str(config_path), "--write", str(report_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    completed = subprocess.run(
        [sys.executable, str(script), "--config", str(config_path), "--check", str(report_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "matches" in completed.stdout
