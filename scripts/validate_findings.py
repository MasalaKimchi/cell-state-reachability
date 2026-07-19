#!/usr/bin/env python3
"""Fail-closed validation of canonical findings and artifact lineage."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
from urllib.parse import unquote

import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
EVIDENCE = RESULTS / "evidence"
MANIFEST = RESULTS / "manifest.json"
MARKDOWN_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")

MANIFEST_PATHS = [
    ".gitignore",
    ".github/workflows/ci.yml",
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "pyproject.toml",
    "requirements.txt",
    "requirements-external.txt",
    "reachability.py",
    "validation.py",
    "effect_dictionary.py",
    "library_coverage.py",
    "demo_library_coverage.py",
    "reproduce.sh",
    "configs/validation_harness.json",
    "configs/source_reconstruction.json",
    "configs/arce_external_validation.json",
    "configs/zhu_arrayed_validation.json",
    "configs/donor_pair_transfer.json",
    "configs/guide_pair_transfer.json",
    "configs/goudy_combination_validation.json",
    "configs/schmidt_external_validation.json",
    "configs/library_coverage_crossdataset.json",
    "data/README.md",
    "data/fetch_de_stats.sh",
    "tests/test_reachability.py",
    "tests/test_validation.py",
    "tests/test_source_reconstruction.py",
    "tests/test_arce_external_validation.py",
    "tests/test_zhu_arrayed_validation.py",
    "tests/test_donor_pair_transfer.py",
    "tests/test_guide_pair_transfer.py",
    "tests/test_cache_reconstruction.py",
    "tests/test_effect_dictionary.py",
    "tests/test_goudy_combination_validation.py",
    "tests/test_schmidt_external_validation.py",
    "tests/test_library_coverage.py",
    "scripts/build_library_coverage_caches.py",
    "scripts/run_validation_harness.py",
    "scripts/run_goudy_combination_validation.py",
    "scripts/run_schmidt_external_validation.py",
    "scripts/run_library_coverage_crossdataset.py",
    "scripts/run_source_reconstruction.py",
    "scripts/run_arce_external_validation.py",
    "scripts/run_zhu_arrayed_validation.py",
    "scripts/run_donor_pair_transfer.py",
    "scripts/run_guide_pair_transfer.py",
    "scripts/validate_findings.py",
    "docs/FINDINGS.md",
    "docs/METHODS.md",
    "docs/SCIENTIFIC_VALIDATION_PLAN.md",
    "docs/VALIDATION_REPORT.md",
    "docs/figures/make_at_a_glance.py",
    "docs/figures/fig_at_a_glance.png",
    "docs/figures/fig_at_a_glance.pdf",
    "tutorial/README.md",
    "tutorial/tutorial.ipynb",
    "results/findings.json",
    "results/goudy_combination_validation.json",
    "results/schmidt_external_validation.json",
    "results/library_coverage_crossdataset.json",
    "results/validation_harness.json",
    "results/source_reconstruction.json",
    "results/donor_pair_transfer.json",
    "results/guide_pair_transfer.json",
    "results/README.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_paths() -> list[str]:
    evidence = [
        str(path.relative_to(ROOT))
        for path in sorted(EVIDENCE.iterdir())
        if path.is_file()
    ]
    return sorted(MANIFEST_PATHS + evidence)


def write_manifest() -> None:
    files = {}
    for relative in artifact_paths():
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(relative)
        files[relative] = {
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
            "executable": bool(path.stat().st_mode & 0o111),
        }
    findings = json.loads((RESULTS / "findings.json").read_text(encoding="utf-8"))
    payload = {
        "schema_version": "1.0.0",
        "generated_on": findings["generated_on"],
        "files": files,
    }
    MANIFEST.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {MANIFEST.relative_to(ROOT)} ({len(files)} files)")


def assert_close(actual: float, expected: float, tolerance: float = 1e-12) -> None:
    if (
        not math.isfinite(actual)
        or not math.isfinite(expected)
        or abs(actual - expected) > tolerance
    ):
        raise AssertionError(f"{actual} != {expected} within {tolerance}")


def _mean_metric(direction: dict, model: str, metric: str) -> float:
    return statistics.mean(
        split["metrics"][model][metric] for split in direction["splits"]
    )


def validate_activation_csv(ledger: dict, metadata: dict) -> None:
    path = EVIDENCE / "arce_activation_guide_effects.csv"
    expected_fields = [
        "target", "guide", "donor", "context", "n_cells",
        "guide_median_supplied_activation_score",
        "ntc_baseline_median_of_guide_medians",
        "supplied_activation_score_delta", "support_flag_lt_20_cells",
        "generator_available", "generator_admitted",
    ]
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != expected_fields:
            raise AssertionError("Arce activation CSV schema differs")
        rows = list(reader)
    if len(rows) != ledger["perturbation_guide_strata"]:
        raise AssertionError("Arce activation CSV row count differs")
    keys = {
        (row["target"], row["guide"], row["donor"], row["context"])
        for row in rows
    }
    if len(keys) != len(rows):
        raise AssertionError("Arce activation CSV has duplicate stratum keys")
    if len({row["target"] for row in rows}) != 28 or len({row["guide"] for row in rows}) != 56:
        raise AssertionError("Arce activation CSV target/guide axes differ")
    if {row["donor"] for row in rows} != {"A", "B"}:
        raise AssertionError("Arce activation CSV donor axis differs")
    contexts = set(ledger["donor_spearman"])
    if {row["context"] for row in rows} != contexts:
        raise AssertionError("Arce activation CSV context axis differs")
    if Counter(row["context"] for row in rows) != Counter({context: 112 for context in contexts}):
        raise AssertionError("Arce activation CSV context counts differ")
    if Counter(row["donor"] for row in rows) != Counter({"A": 224, "B": 224}):
        raise AssertionError("Arce activation CSV donor counts differ")

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    sparse = []
    for row in rows:
        n_cells = int(row["n_cells"])
        observed_delta = float(row["supplied_activation_score_delta"])
        expected_delta = (
            float(row["guide_median_supplied_activation_score"])
            - float(row["ntc_baseline_median_of_guide_medians"])
        )
        assert_close(observed_delta, expected_delta)
        flag = row["support_flag_lt_20_cells"] == "true"
        if flag != (n_cells < 20):
            raise AssertionError("Arce activation support flag differs from cell count")
        if flag:
            sparse.append((row["target"], row["guide"], row["donor"], row["context"], n_cells))
        grouped[(row["target"], row["donor"], row["context"])].append(row)
    if sparse != [("SOCS3", "SOCS3_2_CRISPRi", "A", "Resting-Treg", 8)]:
        raise AssertionError("Arce activation sparse-stratum identity differs")
    if len(grouped) != 224 or any(len(group) != 2 for group in grouped.values()):
        raise AssertionError("Arce activation CSV does not have two guides per target/donor/context")

    targets = sorted({row["target"] for row in rows})
    for context in sorted(contexts):
        donor_values = {}
        all_four_same = []
        for donor in ("A", "B"):
            donor_values[donor] = [
                statistics.median(
                    float(row["supplied_activation_score_delta"])
                    for row in grouped[(target, donor, context)]
                )
                for target in targets
            ]
        observed_rho = float(spearmanr(donor_values["A"], donor_values["B"]).statistic)
        assert_close(observed_rho, ledger["donor_spearman"][context])
        assert_close(
            observed_rho,
            metadata["contexts"][context]["donor_rank_concordance"]["spearman"],
        )
        for target in targets:
            values = [
                float(row["supplied_activation_score_delta"])
                for donor in ("A", "B")
                for row in grouped[(target, donor, context)]
            ]
            all_four_same.append(all(value > 0 for value in values) or all(value < 0 for value in values))
        fraction = statistics.mean(all_four_same)
        assert_close(fraction, ledger["all_two_guides_two_donors_same_sign_fraction"][context])
        assert_close(
            fraction,
            metadata["contexts"][context][
                "all_two_guides_two_donors_same_nonzero_sign_fraction"
            ],
        )


def validate_zhu_arrayed(ledger: dict) -> None:
    metadata_path = EVIDENCE / "zhu_arrayed_validation_meta.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    config_path = ROOT / "configs" / "zhu_arrayed_validation.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if metadata.get("status") != "PASS" or metadata["benchmark"] != ledger["benchmark"]:
        raise AssertionError("Zhu arrayed benchmark did not pass")
    if metadata["config_sha256"] != sha256(config_path):
        raise AssertionError("Zhu arrayed report/config identity differs")
    for name, identity in metadata["input_verification"].items():
        if (
            not identity["hash_verified"]
            or identity["sha256_actual"] != identity["sha256_expected"]
            or identity["sha256_expected"] != config["inputs"][name]["sha256"]
        ):
            raise AssertionError("Zhu arrayed input identity was not fully verified")
    design = metadata["design"]
    assert len(design["perturbations"]) == ledger["perturbations"]
    assert len(design["followup_donor_labels"]) == ledger["followup_donor_labels"]
    assert design["common_genes"] == ledger["common_genes"]
    if (
        not config["analysis"]["mask_all_panel_target_genes"]
        or not design["all_panel_target_genes_masked"]
        or ledger["scored_genes_per_profile"]
        != ledger["common_genes"] - ledger["perturbations"]
    ):
        raise AssertionError("Zhu arrayed panel-target masking contract differs")

    profile_path = EVIDENCE / "zhu_arrayed_profile_metrics.csv"
    with profile_path.open(newline="", encoding="utf-8") as handle:
        profile_rows = list(csv.DictReader(handle))
    if len(profile_rows) != 2 * ledger["perturbations"]:
        raise AssertionError("Zhu arrayed profile row count differs")
    if {(row["scale"], row["target"]) for row in profile_rows} != {
        (scale, target)
        for scale in ("raw", "panel_centered")
        for target in design["perturbations"]
    }:
        raise AssertionError("Zhu arrayed profile keys differ")
    if {int(row["genes_scored"]) for row in profile_rows} != {
        ledger["scored_genes_per_profile"]
    }:
        raise AssertionError("Zhu arrayed profile gene count differs")
    if any(int(row["retrieval_rank"]) != 1 for row in profile_rows):
        raise AssertionError("Zhu arrayed matching profile is not always rank one")

    profile = ledger["profile_replication"]
    raw = metadata["profile_summary"]["raw"]
    centered = metadata["profile_summary"]["panel_centered"]
    retrieval = metadata["conditional_permutation"]["retrieval"]
    assert_close(raw["cosine"], profile["raw_median_cosine"])
    assert_close(
        raw["median_cosine_gain_over_common_source"],
        profile["raw_median_cosine_gain_over_common_source"],
    )
    assert_close(raw["normalized_rmse"], profile["raw_median_normalized_rmse"])
    assert_close(centered["cosine"], profile["panel_centered_median_cosine"])
    assert_close(
        centered["normalized_rmse"], profile["panel_centered_median_normalized_rmse"]
    )
    assert_close(retrieval["raw"]["top1"], profile["raw_top1_retrieval"])
    assert_close(
        retrieval["panel_centered"]["top1"], profile["panel_centered_top1_retrieval"]
    )
    if metadata["conditional_permutation"]["permutations"] != math.factorial(
        ledger["perturbations"]
    ):
        raise AssertionError("Zhu arrayed permutation enumeration differs")
    for name, expected in ledger["rna_to_donor_median_flow_spearman"].items():
        assert_close(
            metadata["conditional_permutation"]["cytokine_rank_association"][name][
                "spearman"
            ],
            expected,
        )

    flow_path = EVIDENCE / "zhu_arrayed_flow_effects.csv"
    with flow_path.open(newline="", encoding="utf-8") as handle:
        flow_rows = list(csv.DictReader(handle))
    if len(flow_rows) != 78:
        raise AssertionError("Zhu arrayed donor/cytokine row count differs")
    keys = {
        (row["target"], row["donor"], row["cytokine"]) for row in flow_rows
    }
    if len(keys) != len(flow_rows):
        raise AssertionError("Zhu arrayed flow rows have duplicate keys")
    coverage = Counter((row["target"], row["cytokine"]) for row in flow_rows)
    if min(coverage.values()) != 3 or max(coverage.values()) != 6:
        raise AssertionError("Zhu arrayed donor coverage differs")
    for row in flow_rows:
        expected = math.log2(
            float(row["percent_positive"])
            / float(row["donor_ntc_mean_percent_positive"])
        )
        assert_close(float(row["log2_ratio_to_donor_ntc"]), expected)


def validate_schmidt(findings: dict) -> None:
    """Cross-bind the Schmidt functional-screen stress report and ledger."""

    ledger = findings.get("schmidt_external_validation")
    if not isinstance(ledger, dict):
        raise AssertionError("findings.schmidt_external_validation is missing")
    report_path = RESULTS / "schmidt_external_validation.json"
    config_path = ROOT / "configs" / "schmidt_external_validation.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))

    def require_keys(value: object, expected: set[str], label: str) -> dict:
        if not isinstance(value, dict) or set(value) != expected:
            raise AssertionError(f"{label} schema differs")
        return value

    def finite(value: object, label: str) -> float:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise AssertionError(f"{label} is not finite numeric data")
        return float(value)

    require_keys(
        ledger,
        {
            "benchmark", "report_version", "tier", "benchmark_role", "contexts",
            "assay", "modalities", "guide_libraries", "fixed_donor_labels",
            "orientation", "primary_contract", "whole_universe_donor_spearman",
            "whole_universe_modality_plus_library_spearman",
            "whole_universe_cross_context_spearman",
            "source_selected_top_k_transfer", "sensitivity", "claim_boundary",
            "interpretation", "source",
        },
        "Schmidt findings ledger",
    )

    require_keys(
        config,
        {
            "schema_version",
            "report_version",
            "benchmark_id",
            "source",
            "input",
            "archive_members",
            "author_script",
            "screens",
            "analysis",
            "claim_contract",
        },
        "Schmidt config",
    )
    require_keys(
        report,
        {
            "schema_version",
            "report_version",
            "benchmark_id",
            "execution_status",
            "config_sha256",
            "provenance",
            "analysis_contract",
            "eligibility",
            "whole_universe_by_minimum_guides",
            "source_selected_transfer",
            "claim_contract",
            "interpretation",
        },
        "Schmidt report",
    )
    benchmark = "schmidt_primary_tcell_crispra_crispri_transfer_v1"
    if (
        config["schema_version"] != 1
        or config["report_version"] != 1
        or report["schema_version"] != 1
        or report["report_version"] != 1
        or config["benchmark_id"] != benchmark
        or report["benchmark_id"] != benchmark
        or ledger["benchmark"] != benchmark
        or ledger["report_version"] != 1
    ):
        raise AssertionError("Schmidt benchmark/version identity differs")
    if (
        ledger["tier"] != "STRESS"
        or ledger["benchmark_role"]
        != "two-fixed-donor functional-screen rank-transfer stress test"
        or ledger["contexts"]
        != {"IL2": "primary human CD4 T cells", "IFNG": "primary human CD8 T cells"}
        or ledger["assay"] != "genome-wide cytokine-sort functional screen"
        or ledger["modalities"] != ["CRISPRa", "CRISPRi"]
        or ledger["guide_libraries"]
        != {"CRISPRa": "Calabrese", "CRISPRi": "Dolcetto"}
        or ledger["fixed_donor_labels"] != 2
        or ledger["orientation"] != {"CRISPRa": 1, "CRISPRi": -1}
        or ledger["source"] != "results/schmidt_external_validation.json"
    ):
        raise AssertionError("Schmidt findings context/claim grain differs")
    if report["execution_status"] != "PASS":
        raise AssertionError("Schmidt external validation did not pass")
    if report["config_sha256"] != sha256(config_path):
        raise AssertionError("Schmidt report/config byte identity differs")
    if sha256(report_path) != (
        "87745d31bcd08bf70d2a6c16287db9abe5fa836d04ed41fda632d7cc09da8810"
    ):
        raise AssertionError("Schmidt frozen report byte identity differs")

    provenance = require_keys(
        report["provenance"],
        {
            "source",
            "input",
            "archive_member_allowlist_count",
            "archive_member_allowlist_verified_exactly",
            "members",
            "donor_suffixes",
            "same_donors_across_modalities",
            "same_donors_basis",
        },
        "Schmidt provenance",
    )
    if provenance["source"] != config["source"]:
        raise AssertionError("Schmidt source metadata differs from config")
    identity = require_keys(
        provenance["input"],
        {
            "path",
            "url",
            "bytes",
            "sha256_expected",
            "sha256_actual",
            "md5_expected",
            "md5_actual",
            "hash_verified_before_zip_parse",
        },
        "Schmidt archive identity",
    )
    registered = config["input"]
    if (
        identity["path"] != registered["path"]
        or identity["url"] != registered["url"]
        or identity["bytes"] != registered["bytes"]
        or identity["sha256_expected"] != registered["sha256"]
        or identity["sha256_actual"] != registered["sha256"]
        or identity["md5_expected"] != registered["md5"]
        or identity["md5_actual"] != registered["md5"]
        or identity["hash_verified_before_zip_parse"] is not True
        or re.fullmatch(r"[0-9a-f]{64}", identity["sha256_actual"]) is None
        or re.fullmatch(r"[0-9a-f]{32}", identity["md5_actual"]) is None
    ):
        raise AssertionError("Schmidt archive was not verified before ZIP parsing")
    if (
        provenance["archive_member_allowlist_count"] != len(config["archive_members"])
        or provenance["archive_member_allowlist_verified_exactly"] is not True
        or len(config["archive_members"]) != len(set(config["archive_members"]))
    ):
        raise AssertionError("Schmidt ZIP member allow-list differs")
    if (
        provenance["donor_suffixes"] != {"r0": "Donor1", "r1": "Donor2"}
        or provenance["same_donors_across_modalities"] is not True
        or provenance["same_donors_basis"] != config["source"]["same_donors_basis"]
    ):
        raise AssertionError("Schmidt fixed-donor mapping differs")

    members = require_keys(
        provenance["members"],
        {"author_script", *config["screens"]},
        "Schmidt member identities",
    )
    author = require_keys(
        members["author_script"],
        {"member", "bytes", "sha256", "required_text_verified"},
        "Schmidt author script",
    )
    for field in ("member", "bytes", "sha256"):
        if author[field] != config["author_script"][field]:
            raise AssertionError(f"Schmidt author-script {field} differs")
    if author["required_text_verified"] is not True:
        raise AssertionError("Schmidt author transformation was not verified")
    schema = [
        "sgrna", "Gene", "control_count", "treatment_count", "control_mean",
        "treat_mean", "LFC", "control_var", "adj_var", "score", "p.low",
        "p.high", "p.twosided", "FDR", "high_in_treatment",
    ]
    for screen_id, spec in config["screens"].items():
        member = require_keys(
            members[screen_id],
            {"member", "rows", "bytes", "sha256", "schema_verified"},
            f"Schmidt {screen_id} member",
        )
        if any(
            member[field] != spec[field]
            for field in ("member", "rows", "bytes", "sha256")
        ) or member["schema_verified"] != schema:
            raise AssertionError(f"Schmidt {screen_id} member/schema identity differs")

    contract = require_keys(
        report["analysis_contract"],
        {
            "orientation",
            "aggregation",
            "eligibility",
            "excluded_gene",
            "minimum_guides_per_gene_per_donor",
            "primary_minimum_guides",
            "top_k",
            "primary_top_k",
            "sensitivity_status",
            "selection",
            "target_overlap",
        },
        "Schmidt analysis contract",
    )
    analysis = config["analysis"]
    for field in (
        "orientation",
        "aggregation",
        "eligibility",
        "excluded_gene",
        "minimum_guides_per_gene_per_donor",
        "primary_minimum_guides",
        "top_k",
        "primary_top_k",
        "sensitivity_status",
    ):
        if contract[field] != analysis[field]:
            raise AssertionError(f"Schmidt analysis contract {field} differs")
    if (
        contract["orientation"] != {"CRISPRa": 1, "CRISPRi": -1}
        or "source screen and training donor" not in contract["selection"]
        or "globally over the complete frozen universe" not in contract["target_overlap"]
    ):
        raise AssertionError("Schmidt orientation/training-only selection differs")
    eligibility = require_keys(
        report["eligibility"],
        {
            "common_gene_counts_by_minimum_guides",
            "outcome_fields_used",
            "identity_and_guide_coverage_only",
        },
        "Schmidt eligibility",
    )
    expected_counts = {"1": 18734, "3": 18568, "6": 18075}
    if (
        eligibility["common_gene_counts_by_minimum_guides"] != expected_counts
        or eligibility["outcome_fields_used"] != []
        or eligibility["identity_and_guide_coverage_only"] is not True
        or list(expected_counts.values()) != sorted(expected_counts.values(), reverse=True)
    ):
        raise AssertionError("Schmidt outcome-independent universe differs")
    if ledger["primary_contract"] != {
        "minimum_guides_per_gene_per_donor": 3,
        "top_k": 200,
        "common_complete_genes": 18568,
        "directions_per_transfer_class": 8,
    }:
        raise AssertionError("Schmidt primary ledger contract differs")

    whole = require_keys(
        report["whole_universe_by_minimum_guides"],
        {"1", "3", "6"},
        "Schmidt whole-universe sensitivity",
    )
    whole_metric_keys = {
        "signed_spearman", "signed_kendall", "sign_agreement",
        "absolute_effect_spearman", "cosine",
    }

    def validate_whole_metrics(value: object, label: str) -> dict:
        metrics = require_keys(value, whole_metric_keys, label)
        for name, observed in metrics.items():
            numeric = finite(observed, f"{label}.{name}")
            lower = 0.0 if name == "sign_agreement" else -1.0
            if not lower <= numeric <= 1.0:
                raise AssertionError(f"{label}.{name} is outside its range")
        return metrics

    whole_categories = {
        "donor_same_reagent",
        "modality_plus_library_same_context",
        "cross_context_cytokine_plus_cell_type",
    }
    for threshold, expected_n in expected_counts.items():
        threshold_rows = require_keys(
            whole[threshold], whole_categories,
            f"Schmidt threshold-{threshold} whole-universe metrics",
        )
        donor_rows = threshold_rows["donor_same_reagent"]
        if (
            not isinstance(donor_rows, list)
            or len(donor_rows) != 4
            or {row.get("screen") for row in donor_rows} != set(config["screens"])
        ):
            raise AssertionError("Schmidt donor whole-universe rows differ")
        for row in donor_rows:
            require_keys(
                row, {"screen", "n_genes", "donor_a", "donor_b", "metrics"},
                f"Schmidt threshold-{threshold} donor row",
            )
            if (
                row["n_genes"] != expected_n
                or {row["donor_a"], row["donor_b"]} != {"r0", "r1"}
            ):
                raise AssertionError("Schmidt donor whole-universe contract differs")
            metrics = validate_whole_metrics(
                row["metrics"], f"Schmidt threshold-{threshold} donor metrics"
            )
            if threshold == "3":
                assert_close(
                    metrics["signed_spearman"],
                    ledger["whole_universe_donor_spearman"][row["screen"]],
                )

        modality_rows = threshold_rows["modality_plus_library_same_context"]
        if (
            not isinstance(modality_rows, list)
            or len(modality_rows) != 2
            or {row.get("context") for row in modality_rows} != {"IFNG", "IL2"}
        ):
            raise AssertionError("Schmidt modality whole-universe rows differ")
        for row in modality_rows:
            require_keys(
                row, {"context", "screen_a", "screen_b", "n_genes", "metrics"},
                f"Schmidt threshold-{threshold} modality row",
            )
            spec_a = config["screens"][row["screen_a"]]
            spec_b = config["screens"][row["screen_b"]]
            if (
                row["n_genes"] != expected_n
                or spec_a["context"] != row["context"]
                or spec_b["context"] != row["context"]
                or {spec_a["modality"], spec_b["modality"]}
                != {"CRISPRa", "CRISPRi"}
            ):
                raise AssertionError("Schmidt modality whole-universe mapping differs")
            metrics = validate_whole_metrics(
                row["metrics"], f"Schmidt threshold-{threshold} modality metrics"
            )
            if threshold == "3":
                assert_close(
                    metrics["signed_spearman"],
                    ledger["whole_universe_modality_plus_library_spearman"][
                        row["context"]
                    ],
                )

        context_rows = threshold_rows["cross_context_cytokine_plus_cell_type"]
        if (
            not isinstance(context_rows, list)
            or len(context_rows) != 2
            or {row.get("modality") for row in context_rows}
            != {"CRISPRa", "CRISPRi"}
        ):
            raise AssertionError("Schmidt context whole-universe rows differ")
        for row in context_rows:
            require_keys(
                row, {"modality", "screen_a", "screen_b", "n_genes", "metrics"},
                f"Schmidt threshold-{threshold} context row",
            )
            spec_a = config["screens"][row["screen_a"]]
            spec_b = config["screens"][row["screen_b"]]
            if (
                row["n_genes"] != expected_n
                or spec_a["modality"] != row["modality"]
                or spec_b["modality"] != row["modality"]
                or {spec_a["context"], spec_b["context"]} != {"IFNG", "IL2"}
            ):
                raise AssertionError("Schmidt context whole-universe mapping differs")
            metrics = validate_whole_metrics(
                row["metrics"], f"Schmidt threshold-{threshold} context metrics"
            )
            if threshold == "3":
                assert_close(
                    metrics["signed_spearman"],
                    ledger["whole_universe_cross_context_spearman"][row["modality"]],
                )

    transfer = require_keys(
        report["source_selected_transfer"],
        {
            "rows",
            "primary_summary",
            "primary_directions_per_class",
            "correlated_descriptive_challenges",
        },
        "Schmidt source-selected transfer",
    )
    rows = transfer["rows"]
    thresholds = (1, 3, 6)
    top_values = (50, 100, 200, 500)
    classes = {
        "same_screen_held_donor",
        "donor_plus_modality_library_same_context",
        "donor_plus_cross_context_cytokine_cell_type_same_modality",
    }
    if (
        len(rows) != len(thresholds) * len(top_values) * 4 * 2 * len(classes)
        or transfer["primary_directions_per_class"] != 8
        or transfer["correlated_descriptive_challenges"] is not True
    ):
        raise AssertionError("Schmidt sensitivity grid size differs")
    expected_row_keys = {
        "minimum_guides", "top_k", "transfer_class", "source_screen",
        "target_screen", "training_donor", "held_donor", "universe_genes",
        "source_selected_gene_sha256", "held_target_global_top_k_sha256",
        "held_target_global_top_k_overlap_count",
        "held_target_global_top_k_overlap_fraction", "metrics",
    }
    expected_metric_keys = whole_metric_keys
    observed_keys = set()
    source_hashes: dict[tuple, set[str]] = defaultdict(set)
    expected_universe = {1: 18734, 3: 18568, 6: 18075}
    for row in rows:
        require_keys(row, expected_row_keys, "Schmidt transfer row")
        metrics = require_keys(row["metrics"], expected_metric_keys, "Schmidt metrics")
        key = (
            row["minimum_guides"], row["top_k"], row["transfer_class"],
            row["source_screen"], row["training_donor"],
        )
        if key in observed_keys:
            raise AssertionError("Schmidt transfer grid has duplicate rows")
        observed_keys.add(key)
        if (
            row["minimum_guides"] not in thresholds
            or row["top_k"] not in top_values
            or row["transfer_class"] not in classes
            or row["source_screen"] not in config["screens"]
            or row["target_screen"] not in config["screens"]
            or row["training_donor"] == row["held_donor"]
            or {row["training_donor"], row["held_donor"]} != {"r0", "r1"}
            or row["universe_genes"] != expected_universe[row["minimum_guides"]]
            or re.fullmatch(r"[0-9a-f]{64}", row["source_selected_gene_sha256"]) is None
            or re.fullmatch(r"[0-9a-f]{64}", row["held_target_global_top_k_sha256"]) is None
        ):
            raise AssertionError("Schmidt transfer row contract differs")
        source_spec = config["screens"][row["source_screen"]]
        target_spec = config["screens"][row["target_screen"]]
        if row["transfer_class"] == "same_screen_held_donor":
            mapping_valid = row["target_screen"] == row["source_screen"]
        elif row["transfer_class"] == "donor_plus_modality_library_same_context":
            mapping_valid = (
                target_spec["context"] == source_spec["context"]
                and target_spec["modality"] != source_spec["modality"]
            )
        else:
            mapping_valid = (
                target_spec["modality"] == source_spec["modality"]
                and target_spec["context"] != source_spec["context"]
            )
        if not mapping_valid:
            raise AssertionError("Schmidt transfer class/screen mapping differs")
        source_key = (
            row["minimum_guides"], row["top_k"], row["source_screen"],
            row["training_donor"],
        )
        source_hashes[source_key].add(row["source_selected_gene_sha256"])
        overlap_count = row["held_target_global_top_k_overlap_count"]
        if (
            isinstance(overlap_count, bool)
            or not isinstance(overlap_count, int)
            or not 0 <= overlap_count <= row["top_k"]
        ):
            raise AssertionError("Schmidt global top-k overlap count differs")
        assert_close(
            row["held_target_global_top_k_overlap_fraction"],
            overlap_count / row["top_k"],
        )
        for name, value in metrics.items():
            numeric = finite(value, f"Schmidt {name}")
            low = 0.0 if name == "sign_agreement" else -1.0
            if not low <= numeric <= 1.0:
                raise AssertionError(f"Schmidt {name} is outside its range")
    if len(observed_keys) != len(rows) or any(len(values) != 1 for values in source_hashes.values()):
        raise AssertionError("Schmidt source gene set changed across target evaluations")

    summary = require_keys(transfer["primary_summary"], classes, "Schmidt primary summary")
    ledger_summary = ledger["source_selected_top_k_transfer"]
    ledger_field_map = {
        "signed_spearman": "median_signed_spearman",
        "absolute_effect_spearman": "median_absolute_effect_spearman",
        "sign_agreement": "median_sign_agreement",
        "held_target_global_top_k_overlap_fraction": (
            "median_held_target_global_top_k_overlap_fraction"
        ),
    }
    summary_metrics = expected_metric_keys | {
        "held_target_global_top_k_overlap_fraction"
    }
    for transfer_class in classes:
        primary_rows = [
            row for row in rows
            if row["minimum_guides"] == 3
            and row["top_k"] == 200
            and row["transfer_class"] == transfer_class
        ]
        if len(primary_rows) != 8:
            raise AssertionError("Schmidt primary class lacks eight donor directions")
        class_summary = require_keys(
            summary[transfer_class], summary_metrics,
            f"Schmidt {transfer_class} primary summary",
        )
        for metric in summary_metrics:
            values = [
                row[metric] if metric in row else row["metrics"][metric]
                for row in primary_rows
            ]
            reported = require_keys(
                class_summary[metric], {"median", "minimum", "maximum"},
                f"Schmidt {transfer_class}.{metric} summary",
            )
            assert_close(reported["median"], statistics.median(values))
            assert_close(reported["minimum"], min(values))
            assert_close(reported["maximum"], max(values))
            if metric in ledger_field_map:
                assert_close(
                    reported["median"],
                    ledger_summary[transfer_class][ledger_field_map[metric]],
                )

    if report["claim_contract"] != config["claim_contract"]:
        raise AssertionError("Schmidt claim contract differs from config")
    claim = report["claim_contract"]
    if (
        claim["tier"] != "STRESS"
        or "donor-plus-modality-and-library" not in claim["cross_modality_label"]
        or "donor-plus-cross-context" not in claim["cross_context_label"]
        or "jointly change donor" not in claim["claim_ceiling"]
        or "No p-values, confidence intervals" not in claim["inference"]
        or "guide-held-out" not in claim["claim_ceiling"]
        or ledger["claim_boundary"] != claim["claim_ceiling"]
        or report["interpretation"]["status"]
        != "CONDITIONAL_TOP_EFFECT_CONCORDANCE_ONLY"
    ):
        raise AssertionError("Schmidt claim ceiling/status differs")
    if ledger["sensitivity"] != {
        "minimum_guides_per_gene_per_donor": [1, 3, 6],
        "top_k": [50, 100, 200, 500],
        "status": "EXPLORATORY_POST_HOC_DESCRIPTIVE",
    }:
        raise AssertionError("Schmidt exploratory sensitivity ledger differs")


def validate_goudy(findings: dict) -> None:
    """Cross-bind the hardened Goudy v2 stress report and findings ledger."""

    if "goudy_combination_validation" in findings:
        raise AssertionError("retired findings.goudy_combination_validation key remains")
    ledger = findings.get("goudy_cross_experiment_stress")
    if ledger is None:
        raise AssertionError("findings.goudy_cross_experiment_stress section is missing")

    report_path = RESULTS / "goudy_combination_validation.json"
    config_path = ROOT / "configs" / "goudy_combination_validation.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))

    def require_keys(value: object, expected: set[str], label: str) -> dict:
        if not isinstance(value, dict) or set(value) != expected:
            raise AssertionError(f"{label} schema differs")
        return value

    def finite_number(value: object, label: str) -> float:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise AssertionError(f"{label} is not finite numeric data")
        return float(value)

    def assert_range(actual: object, expected: list[float], label: str) -> None:
        if not isinstance(actual, list) or len(actual) != 2:
            raise AssertionError(f"{label} is not a two-value range")
        for observed, reference in zip(actual, expected, strict=True):
            assert_close(finite_number(observed, label), reference)

    def validate_summary(summary: object, values: list[float], label: str) -> None:
        item = require_keys(summary, {"mean", "median", "donor_range"}, label)
        assert_close(finite_number(item["mean"], label), statistics.fmean(values))
        assert_close(finite_number(item["median"], label), statistics.median(values))
        assert_range(item["donor_range"], [min(values), max(values)], label)

    require_keys(
        config,
        {
            "schema_version",
            "benchmark_id",
            "source",
            "inputs",
            "expected",
            "preprocessing",
            "models",
            "metadata_conflict",
            "claim_contract",
            "donors",
        },
        "Goudy config",
    )
    require_keys(
        report,
        {
            "report_version",
            "benchmark_id",
            "config_canonical_json_sha256",
            "execution_status",
            "geometric_model_status",
            "biological_interpretation_status",
            "source",
            "input_identity",
            "sample_contract",
            "metadata_conflict",
            "preprocessing",
            "models",
            "claim_contract",
            "donor_results",
            "donor_ranges",
            "on_target_controls",
            "filter_sensitivity",
            "reliability_diagnostics",
            "additive_residual",
            "module_error",
            "interpretation",
        },
        "Goudy report",
    )
    if config["schema_version"] != 2 or report["report_version"] != 2:
        raise AssertionError("unsupported Goudy config/report schema")
    benchmark = "goudy_gse306915_cross_experiment_stress"
    if config["benchmark_id"] != benchmark or report["benchmark_id"] != benchmark:
        raise AssertionError("Goudy v2 benchmark identity differs")
    canonical = json.dumps(
        config, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    if report["config_canonical_json_sha256"] != hashlib.sha256(canonical).hexdigest():
        raise AssertionError("Goudy report/config canonical identity differs")

    statuses = {
        "execution": "PASS",
        "geometric_model": "FAILS_DECLARED_GEOMETRIC_MODEL",
        "biological_interpretation": (
            "INCONCLUSIVE_CROSS_EXPERIMENT_CONFOUNDING_LOW_RELIABILITY"
        ),
    }
    if {
        "execution": report["execution_status"],
        "geometric_model": report["geometric_model_status"],
        "biological_interpretation": report["biological_interpretation_status"],
    } != statuses:
        raise AssertionError("Goudy execution/geometric/biological status differs")
    if report["source"] != config["source"]:
        raise AssertionError("Goudy source provenance differs from config")
    if report["models"] != config["models"]:
        raise AssertionError("Goudy model declaration differs from config")
    if report["claim_contract"] != config["claim_contract"]:
        raise AssertionError("Goudy claim contract differs from config")
    claim = report["claim_contract"]
    if (
        claim["tier"] != "STRESS"
        or claim["benchmark_role"] != "negative cross-experiment stress test"
        or claim["constituent_relationship"] != "target-matched constituent singles"
        or claim["triple_guide_identity"] != "UNRESOLVED"
        or "combination or additivity validation" not in claim["not_claimed"]
        or "same-guide matching between triple and constituent singles"
        not in claim["not_claimed"]
        or "inseparable from experiment, control type, and guide burden"
        not in claim["confounding"]
    ):
        raise AssertionError("Goudy negative-stress claim ceiling differs")

    inputs = require_keys(
        report["input_identity"], set(config["inputs"]), "Goudy input identity"
    )
    if set(inputs) != {"counts", "soft", "author_key"}:
        raise AssertionError("Goudy must bind exactly three registered inputs")
    for name, registered in config["inputs"].items():
        identity = require_keys(
            inputs[name],
            (set(registered) - {"sha256", "md5"})
            | {
                "sha256_expected",
                "sha256_actual",
                "md5_expected",
                "md5_actual",
                "hash_verified_before_parse",
            },
            f"Goudy {name} input identity",
        )
        for field in set(registered) - {"sha256", "md5"}:
            if identity[field] != registered[field]:
                raise AssertionError(f"Goudy {name}.{field} identity differs")
        if (
            identity["hash_verified_before_parse"] is not True
            or identity["sha256_expected"] != registered["sha256"]
            or identity["sha256_actual"] != registered["sha256"]
            or identity["md5_expected"] != registered["md5"]
            or identity["md5_actual"] != registered["md5"]
            or re.fullmatch(r"[0-9a-f]{64}", identity["sha256_actual"]) is None
            or re.fullmatch(r"[0-9a-f]{32}", identity["md5_actual"]) is None
        ):
            raise AssertionError(f"Goudy {name} input was not hash-verified")
    author_commit = config["inputs"]["author_key"]["commit"]
    if (
        re.fullmatch(r"[0-9a-f]{40}", author_commit) is None
        or f"/{author_commit}/" not in config["inputs"]["author_key"]["url"]
    ):
        raise AssertionError("Goudy author metadata is not commit-pinned")

    expected = config["expected"]
    donors = tuple(config["donors"])
    if donors != ("D1", "D2", "D3", "D4") or expected["donors"] != len(donors):
        raise AssertionError("Goudy fixed-donor contract differs")
    sample_contract = require_keys(
        report["sample_contract"], {"geo_soft", "author_key"}, "Goudy samples"
    )
    geo = require_keys(
        sample_contract["geo_soft"],
        {
            "series_accession",
            "soft_sample_count",
            "declared_roles",
            "observed_time_conflict_gsms",
        },
        "Goudy GEO sample contract",
    )
    author = require_keys(
        sample_contract["author_key"],
        {
            "author_key_row_count",
            "count_matrix_overlap_count",
            "cas_modality_row_counts",
            "author_key_rows_absent_from_count_matrix",
            "declared_analysis_role_count",
            "author_key_confirmed_declared_role_count",
            "declared_rows_absent_from_author_key",
            "unresolved_experiment_id_columns",
            "declared_roles",
            "experiment_role_summary",
            "constituent_relationship",
            "constituent_single_guide_notes",
            "triple_guide_identity",
            "same_guide_match_claimed",
        },
        "Goudy author-key contract",
    )
    if (
        geo["series_accession"] != config["source"]["accession"]
        or geo["soft_sample_count"] != expected["soft_samples"]
        or author["author_key_row_count"] != expected["author_key_samples"]
        or author["count_matrix_overlap_count"] != expected["author_key_count_overlap"]
        or author["declared_analysis_role_count"]
        != expected["declared_analysis_samples"]
        or author["author_key_confirmed_declared_role_count"] != 36
        or author["author_key_rows_absent_from_count_matrix"]
        != expected["author_key_rows_absent_from_counts"]
        or author["declared_rows_absent_from_author_key"]
        != expected["declared_rows_absent_from_author_key"]
        or author["unresolved_experiment_id_columns"]
        != expected["declared_rows_absent_from_author_key"]
        or author["cas_modality_row_counts"] != {"CRISPRoff": 40, "Cas9": 31}
        or author["same_guide_match_claimed"] is not False
        or author["triple_guide_identity"] != "UNRESOLVED"
    ):
        raise AssertionError("Goudy sample/author-key count or gap contract differs")

    declared_columns: list[str] = []
    confirmed_roles = 0
    for donor, donor_spec in config["donors"].items():
        expected_geo_roles: dict[str, dict[str, str]] = {}
        for group in ("component", "unrelated", "controls"):
            for role, spec in donor_spec[group].items():
                expected_geo_roles[f"{group}:{role}"] = {
                    "gsm": spec["gsm"],
                    "column": spec["column"],
                    "title": spec["title"],
                    "characteristic_time": spec["time"],
                }
                declared_columns.append(spec["column"])
        spec = donor_spec["triple"]
        expected_geo_roles["triple:triple"] = {
            "gsm": spec["gsm"],
            "column": spec["column"],
            "title": spec["title"],
            "characteristic_time": spec["time"],
        }
        declared_columns.append(spec["column"])
        if geo["declared_roles"][donor] != expected_geo_roles:
            raise AssertionError(f"Goudy {donor} GEO role crosswalk differs")
        author_roles = author["declared_roles"][donor]
        if set(author_roles) != set(expected_geo_roles):
            raise AssertionError(f"Goudy {donor} author role crosswalk differs")
        for role, role_row in author_roles.items():
            if role_row["column"] != expected_geo_roles[role]["column"]:
                raise AssertionError(f"Goudy {donor}.{role} column crosswalk differs")
            absent = role_row["column"] in expected["declared_rows_absent_from_author_key"]
            expected_status = "ABSENT_FROM_AUTHOR_KEY" if absent else "CONFIRMED"
            if role_row["author_key_status"] != expected_status:
                raise AssertionError(f"Goudy {donor}.{role} author-key status differs")
            confirmed_roles += not absent
            if absent and any(
                role_row[field] != "UNRESOLVED"
                for field in ("experiment_id", "experimental_day", "guide_identity")
            ):
                raise AssertionError(f"Goudy {donor}.{role} unresolved fields differ")
    if len(set(declared_columns)) != 40 or confirmed_roles != 36:
        raise AssertionError("Goudy declared sample crosswalk is not one-to-one")

    conflict = {**config["metadata_conflict"], "observed_affected_gsms": geo["observed_time_conflict_gsms"]}
    if (
        report["metadata_conflict"] != conflict
        or geo["observed_time_conflict_gsms"]
        != config["metadata_conflict"]["affected_gsms"]
        or config["metadata_conflict"]["status"] != "FLAGGED"
    ):
        raise AssertionError("Goudy Day18/Day7 conflict contract differs")

    preprocessing = report["preprocessing"]
    preprocessing_fields = set(config["preprocessing"]) | {
        "observed_gene_rows",
        "observed_count_samples",
        "observed_gene_universe_before_mask",
        "observed_gene_universe_before_mask_by_threshold",
        "observed_transcriptome_genes_after_mask",
        "universe_gene_ids_sha256",
        "masked_gene_ids_sha256",
        "control_columns",
        "control_pairing",
        "selected_sample_library_sums",
    }
    require_keys(preprocessing, preprocessing_fields, "Goudy preprocessing")
    for field, value in config["preprocessing"].items():
        if preprocessing[field] != value:
            raise AssertionError(f"Goudy preprocessing.{field} differs from config")
    if (
        preprocessing["observed_gene_rows"] != expected["gene_rows"]
        or preprocessing["observed_count_samples"] != expected["count_samples"]
        or preprocessing["observed_gene_universe_before_mask"]
        != expected["analysis_genes_before_mask"]
        or preprocessing["observed_transcriptome_genes_after_mask"]
        != expected["analysis_genes_after_mask"]
        or preprocessing["observed_gene_universe_before_mask_by_threshold"]
        != expected["analysis_genes_before_mask_by_threshold"]
        or len(preprocessing["control_columns"]) != 12
        or set(preprocessing["selected_sample_library_sums"]) != set(declared_columns)
        or any(
            finite_number(value, "Goudy selected sample library sum") <= 0
            for value in preprocessing["selected_sample_library_sums"].values()
        )
    ):
        raise AssertionError("Goudy preprocessing count/mask contract differs")
    for field in ("universe_gene_ids_sha256", "masked_gene_ids_sha256"):
        if re.fullmatch(r"[0-9a-f]{64}", preprocessing[field]) is None:
            raise AssertionError(f"Goudy preprocessing.{field} is not hash-bound")

    component_order = tuple(config["models"]["component_order"])
    model_names = tuple(
        field
        for field in config["models"]
        if field not in {"component_order", "unrelated_order"}
    )
    metric_names = {"cosine", "normalized_rmse", "norm_ratio"}
    donor_results = require_keys(
        report["donor_results"], set(donors), "Goudy donor results"
    )
    strict_inside = 0
    for donor in donors:
        models = require_keys(
            require_keys(donor_results[donor], {"models"}, f"Goudy {donor}")["models"],
            set(model_names),
            f"Goudy {donor} models",
        )
        for model_name in model_names:
            metrics = require_keys(
                models[model_name].get("metrics"),
                metric_names,
                f"Goudy {donor}.{model_name} metrics",
            )
            for metric, value in metrics.items():
                finite_number(value, f"Goudy {donor}.{model_name}.{metric}")
        geometry = models["component_cone"].get("geometry")
        if not isinstance(geometry, dict):
            raise AssertionError(f"Goudy {donor} component geometry is missing")
        status = geometry.get("strict_status")
        if status not in {"inside_tolerance", "outside_model_cone"}:
            raise AssertionError(f"Goudy {donor} component strict status differs")
        strict_inside += status == "inside_tolerance"
        training = [candidate for candidate in donors if candidate != donor]
        for model_name in (
            "lodo_component_cone",
            "lodo_training_selected_best_single",
            "training_donor_mean_triple",
        ):
            model = models[model_name]
            gene_filter = model.get("gene_filter")
            if (
                model.get("training_donors") != training
                or not isinstance(gene_filter, dict)
                or gene_filter.get("training_donors") != training
                or gene_filter.get("held_donor_controls_used") is not False
                or gene_filter.get("selection_scope") != "training-donor controls only"
            ):
                raise AssertionError(f"Goudy {donor}.{model_name} leaks held donor")
        frozen = models["on_target_calibrated_component_cone"]
        if (
            frozen.get("coefficient_order") != list(component_order)
            or "disjoint target-masked transcriptome" not in frozen.get("scoring_contract", "")
            or frozen.get("calibration_fit", {}).get("warning_status")
            != "SATURATED_3X3_CALIBRATION_DIAGNOSTIC_NOT_VALIDATION"
        ):
            raise AssertionError(f"Goudy {donor} on-target-frozen contract differs")
    if strict_inside != 0:
        raise AssertionError("Goudy primary component cone must be strict 0/4")

    donor_ranges = require_keys(
        report["donor_ranges"], set(model_names), "Goudy donor summaries"
    )
    for model_name in model_names:
        summaries = require_keys(
            donor_ranges[model_name], metric_names, f"Goudy {model_name} summaries"
        )
        for metric in metric_names:
            values = [
                finite_number(
                    donor_results[donor]["models"][model_name]["metrics"][metric],
                    f"Goudy {donor}.{model_name}.{metric}",
                )
                for donor in donors
            ]
            validate_summary(summaries[metric], values, f"Goudy {model_name}.{metric}")

    on_target = require_keys(
        report["on_target_controls"],
        {
            "per_donor",
            "summary",
            "saturated_calibration_diagnostic",
            "excluded_from_transcriptome_metrics",
        },
        "Goudy on-target controls",
    )
    if on_target["excluded_from_transcriptome_metrics"] is not True:
        raise AssertionError("Goudy on-target rows entered transcriptome metrics")
    target_rows = {row["symbol"]: row for row in preprocessing["target_mask"]}
    per_donor_targets = require_keys(
        on_target["per_donor"], set(donors), "Goudy per-donor on-target values"
    )
    all_negative_donors = 0
    for donor in donors:
        donor_targets = require_keys(
            per_donor_targets[donor], set(component_order), f"Goudy {donor} targets"
        )
        donor_negative = True
        for target in component_order:
            row = require_keys(
                donor_targets[target],
                {"gene_id", "single_effect", "triple_effect"},
                f"Goudy {donor}.{target} on-target row",
            )
            if row["gene_id"] != target_rows[target]["gene_id"]:
                raise AssertionError(f"Goudy {donor}.{target} target identity differs")
            single = finite_number(row["single_effect"], "Goudy target single effect")
            triple = finite_number(row["triple_effect"], "Goudy target triple effect")
            donor_negative &= single < 0 and triple < 0
        all_negative_donors += donor_negative
    target_summary = require_keys(
        on_target["summary"], set(component_order), "Goudy target summaries"
    )
    for target in component_order:
        summary = require_keys(
            target_summary[target],
            {
                "single_effect_mean",
                "single_effect_median",
                "single_effect_range",
                "single_negative_donors",
                "triple_effect_mean",
                "triple_effect_median",
                "triple_effect_range",
                "triple_negative_donors",
            },
            f"Goudy {target} target summary",
        )
        for effect in ("single", "triple"):
            values = [
                finite_number(
                    per_donor_targets[donor][target][f"{effect}_effect"],
                    f"Goudy {target} {effect} effect",
                )
                for donor in donors
            ]
            assert_close(summary[f"{effect}_effect_mean"], statistics.fmean(values))
            assert_close(summary[f"{effect}_effect_median"], statistics.median(values))
            assert_range(
                summary[f"{effect}_effect_range"],
                [min(values), max(values)],
                f"Goudy {target} {effect} range",
            )
            if summary[f"{effect}_negative_donors"] != sum(value < 0 for value in values):
                raise AssertionError(f"Goudy {target} {effect} sign count differs")

    warning = "SATURATED_3X3_CALIBRATION_DIAGNOSTIC_NOT_VALIDATION"
    saturated = require_keys(
        on_target["saturated_calibration_diagnostic"],
        {"scope", "warning_status", "interpretation", "per_donor", "donor_ranges"},
        "Goudy saturated calibration",
    )
    if (
        saturated["scope"] != "three on-target coordinates fitted by three atoms"
        or saturated["warning_status"] != warning
        or "not mechanistic, combination, or additivity validation"
        not in saturated["interpretation"]
    ):
        raise AssertionError("Goudy saturated 3x3 warning differs")
    saturated_per_donor = require_keys(
        saturated["per_donor"], set(donors), "Goudy saturated donor values"
    )
    saturated_models = {"component_cone", "equal_sum", "best_single_oracle"}
    for donor in donors:
        models = require_keys(
            require_keys(
                saturated_per_donor[donor], {"models"}, f"Goudy {donor} calibration"
            )["models"],
            saturated_models,
            f"Goudy {donor} calibration models",
        )
        for model_name in saturated_models:
            metrics = require_keys(
                models[model_name].get("metrics"),
                metric_names,
                f"Goudy {donor} calibration {model_name}",
            )
            for metric, value in metrics.items():
                finite_number(value, f"Goudy {donor} calibration {model_name}.{metric}")
    saturated_ranges = require_keys(
        saturated["donor_ranges"], saturated_models, "Goudy calibration summaries"
    )
    for model_name in saturated_models:
        summaries = require_keys(
            saturated_ranges[model_name],
            metric_names,
            f"Goudy calibration {model_name} summaries",
        )
        for metric in metric_names:
            values = [
                saturated_per_donor[donor]["models"][model_name]["metrics"][metric]
                for donor in donors
            ]
            validate_summary(
                summaries[metric], values, f"Goudy calibration {model_name}.{metric}"
            )

    filter_sensitivity = require_keys(
        report["filter_sensitivity"],
        {
            "selection_contract",
            "per_threshold",
            "component_cone_median_cosine_span",
            "filter_sensitivity_status",
            "conclusion",
        },
        "Goudy filter sensitivity",
    )
    if (
        filter_sensitivity["selection_contract"]
        != "Each mask uses only the 12 declared controls, never perturbation outcomes."
        or filter_sensitivity["filter_sensitivity_status"] != "FILTER_SENSITIVE"
    ):
        raise AssertionError("Goudy filter-selection contract differs")
    thresholds = config["preprocessing"]["control_mean_cpm_sensitivity_thresholds"]

    def threshold_key(value: float) -> str:
        return str(int(value)) if float(value).is_integer() else str(value)

    threshold_keys = [threshold_key(value) for value in thresholds]
    threshold_rows = require_keys(
        filter_sensitivity["per_threshold"],
        set(threshold_keys),
        "Goudy filter thresholds",
    )
    target_rows_removed: dict[str, int] = {}
    threshold_component_medians: dict[str, float] = {}
    for threshold, key in zip(thresholds, threshold_keys, strict=True):
        row = require_keys(
            threshold_rows[key],
            {
                "control_mean_cpm_min",
                "retained_genes_before_target_mask",
                "retained_transcriptome_genes",
                "gene_ids_sha256",
                "per_donor",
                "donor_ranges",
                "conclusion",
            },
            f"Goudy CPM {key} result",
        )
        before = expected["analysis_genes_before_mask_by_threshold"][key]
        if (
            row["control_mean_cpm_min"] != threshold
            or row["retained_genes_before_target_mask"] != before
            or re.fullmatch(r"[0-9a-f]{64}", row["gene_ids_sha256"]) is None
        ):
            raise AssertionError(f"Goudy CPM {key} mask identity differs")
        removed = before - row["retained_transcriptome_genes"]
        target_rows_removed[key] = removed
        if removed not in {2, 3}:
            raise AssertionError(f"Goudy CPM {key} target-mask count differs")
        threshold_donors = require_keys(
            row["per_donor"], set(donors), f"Goudy CPM {key} donors"
        )
        for donor in donors:
            threshold_models = require_keys(
                require_keys(
                    threshold_donors[donor],
                    {"models"},
                    f"Goudy CPM {key} {donor}",
                )["models"],
                {"component_cone", "equal_sum"},
                f"Goudy CPM {key} {donor} models",
            )
            for model_name in ("component_cone", "equal_sum"):
                metrics = require_keys(
                    threshold_models[model_name].get("metrics"),
                    metric_names,
                    f"Goudy CPM {key} {donor}.{model_name}",
                )
                for metric, value in metrics.items():
                    finite_number(value, f"Goudy CPM {key} {donor}.{model_name}.{metric}")
            if threshold_models["component_cone"].get("strict_status") != "outside_model_cone":
                raise AssertionError(f"Goudy CPM {key} component cone was not strict 0/4")
        summaries = require_keys(
            row["donor_ranges"],
            {"component_cone", "equal_sum"},
            f"Goudy CPM {key} summaries",
        )
        for model_name in ("component_cone", "equal_sum"):
            model_summaries = require_keys(
                summaries[model_name],
                metric_names,
                f"Goudy CPM {key} {model_name} summaries",
            )
            for metric in metric_names:
                values = [
                    threshold_donors[donor]["models"][model_name]["metrics"][metric]
                    for donor in donors
                ]
                validate_summary(
                    model_summaries[metric],
                    values,
                    f"Goudy CPM {key} {model_name}.{metric}",
                )
        threshold_component_medians[key] = summaries["component_cone"]["cosine"][
            "median"
        ]
    if target_rows_removed != {"0.5": 3, "1": 3, "2": 3, "5": 2, "10": 2}:
        raise AssertionError("Goudy threshold target-row removal contract differs")
    span = max(threshold_component_medians.values()) - min(
        threshold_component_medians.values()
    )
    assert_close(filter_sensitivity["component_cone_median_cosine_span"], span)

    reliability = require_keys(
        report["reliability_diagnostics"],
        {
            "pairwise_donor_cosines",
            "aavs1_guide_replicate_noise",
            "snr_status",
            "overall_status",
            "interpretation",
        },
        "Goudy reliability diagnostics",
    )
    if (
        reliability["snr_status"] != "SNR_LIMITED_LOW_DONOR_REPRODUCIBILITY"
        or reliability["overall_status"] != "SNR_LIMITED_FILTER_SENSITIVE"
    ):
        raise AssertionError("Goudy reliability status differs")
    pairwise = require_keys(
        reliability["pairwise_donor_cosines"],
        {"triple", "components"},
        "Goudy pairwise donor cosines",
    )
    components = require_keys(
        pairwise["components"], set(component_order), "Goudy component reliability"
    )
    pair_rows = {"triple": pairwise["triple"], **components}
    pair_medians: dict[str, float] = {}
    expected_pairs = {
        tuple(sorted((donor_a, donor_b)))
        for index, donor_a in enumerate(donors)
        for donor_b in donors[index + 1 :]
    }
    for name, raw in pair_rows.items():
        item = require_keys(
            raw, {"pair_count", "pairs", "summary"}, f"Goudy {name} reliability"
        )
        if item["pair_count"] != 6 or len(item["pairs"]) != 6:
            raise AssertionError(f"Goudy {name} pair count differs")
        observed_pairs: set[tuple[str, str]] = set()
        values: list[float] = []
        for pair in item["pairs"]:
            pair = require_keys(
                pair, {"donor_a", "donor_b", "cosine"}, f"Goudy {name} pair"
            )
            observed_pairs.add(tuple(sorted((pair["donor_a"], pair["donor_b"]))))
            values.append(finite_number(pair["cosine"], f"Goudy {name} pair cosine"))
        if observed_pairs != expected_pairs:
            raise AssertionError(f"Goudy {name} donor pairs differ")
        summary = require_keys(
            item["summary"], {"mean", "median", "range"}, f"Goudy {name} pair summary"
        )
        assert_close(summary["mean"], statistics.fmean(values))
        assert_close(summary["median"], statistics.median(values))
        assert_range(summary["range"], [min(values), max(values)], f"Goudy {name} pairs")
        pair_medians[name] = summary["median"]

    aavs1 = require_keys(
        reliability["aavs1_guide_replicate_noise"],
        {"definition", "per_donor", "donor_ranges"},
        "Goudy AAVS1 guide noise",
    )
    aavs1_donors = require_keys(
        aavs1["per_donor"], set(donors), "Goudy AAVS1 donor noise"
    )
    aavs1_metrics: dict[str, list[float]] = {
        "cosine": [],
        "normalized_rmse": [],
    }
    triple_noise_ratios: list[float] = []
    for donor in donors:
        row = require_keys(
            aavs1_donors[donor],
            {
                "metrics",
                "difference_l2_norm",
                "difference_norm_over_median_constituent_effect_norm",
                "difference_norm_over_triple_effect_norm",
                "prediction_role",
            },
            f"Goudy {donor} AAVS1 noise",
        )
        metrics = require_keys(row["metrics"], metric_names, f"Goudy {donor} AAVS1 metrics")
        for metric in aavs1_metrics:
            aavs1_metrics[metric].append(
                finite_number(metrics[metric], f"Goudy {donor} AAVS1 {metric}")
            )
        triple_noise_ratios.append(
            finite_number(
                row["difference_norm_over_triple_effect_norm"],
                f"Goudy {donor} AAVS1/triple ratio",
            )
        )
    aavs1_ranges = require_keys(
        aavs1["donor_ranges"],
        {"cosine", "normalized_rmse"},
        "Goudy AAVS1 summaries",
    )
    for metric, values in aavs1_metrics.items():
        summary = require_keys(
            aavs1_ranges[metric], {"median", "range"}, f"Goudy AAVS1 {metric}"
        )
        assert_close(summary["median"], statistics.median(values))
        assert_range(summary["range"], [min(values), max(values)], f"Goudy AAVS1 {metric}")

    residual = require_keys(
        report["additive_residual"],
        {"definition", "per_donor", "across_donor_mean_residual", "causal_interpretation"},
        "Goudy additive residual",
    )
    residual_fields = {
        "residual_l2_norm",
        "residual_norm_over_triple_norm",
        "mean_residual",
        "median_residual",
        "mean_absolute_residual",
        "median_absolute_residual",
        "positive_gene_count",
        "negative_gene_count",
        "top_absolute_residual_genes",
    }
    residual_donors = require_keys(
        residual["per_donor"], set(donors), "Goudy per-donor additive residual"
    )
    residual_ratios: list[float] = []
    for donor in donors:
        row = require_keys(
            residual_donors[donor], residual_fields, f"Goudy {donor} additive residual"
        )
        ratio = finite_number(
            row["residual_norm_over_triple_norm"], f"Goudy {donor} residual ratio"
        )
        residual_ratios.append(ratio)
        assert_close(
            ratio,
            donor_results[donor]["models"]["equal_sum"]["metrics"]["normalized_rmse"],
        )
    across_residual = require_keys(
        residual["across_donor_mean_residual"],
        residual_fields,
        "Goudy across-donor additive residual",
    )
    for field in (
        "residual_l2_norm",
        "residual_norm_over_triple_norm",
        "mean_residual",
        "median_residual",
        "mean_absolute_residual",
        "median_absolute_residual",
    ):
        finite_number(across_residual[field], f"Goudy mean residual {field}")
    module_error = require_keys(
        report["module_error"], {"status", "reason"}, "Goudy module error"
    )
    if module_error["status"] != "UNAVAILABLE_NO_PREREGISTERED_MODULE_SET":
        raise AssertionError("Goudy module-error status differs")

    ledger_fields = {
        "benchmark",
        "report_version",
        "tier",
        "benchmark_role",
        "cell_system",
        "assay",
        "modality",
        "fixed_donor_labels",
        "status",
        "primary_transcriptome",
        "on_target_controls",
        "filter_sensitivity",
        "reliability",
        "additive_residual",
        "provenance_limits",
        "claim_boundary",
        "interpretation",
        "source",
    }
    require_keys(ledger, ledger_fields, "Goudy findings ledger")
    require_keys(
        ledger["status"],
        {"execution", "geometric_model", "biological_interpretation"},
        "Goudy findings status",
    )
    primary = require_keys(
        ledger["primary_transcriptome"],
        {
            "analysis_genes_before_target_mask",
            "transcriptome_genes_after_target_mask",
            "component_cone_median_cosine",
            "component_cone_median_normalized_rmse",
            "component_cone_strict_inside_count",
            "equal_sum_median_cosine",
            "equal_sum_median_normalized_rmse",
            "best_single_in_sample_oracle_median_cosine",
            "best_single_in_sample_oracle_median_normalized_rmse",
            "lodo_component_cone_median_cosine",
            "lodo_component_cone_median_normalized_rmse",
            "lodo_training_selected_best_single_median_cosine",
            "lodo_training_selected_best_single_median_normalized_rmse",
            "on_target_frozen_disjoint_median_cosine",
            "on_target_frozen_disjoint_median_normalized_rmse",
        },
        "Goudy findings primary transcriptome",
    )
    expected_primary = {
        "analysis_genes_before_target_mask": preprocessing["observed_gene_universe_before_mask"],
        "transcriptome_genes_after_target_mask": preprocessing["observed_transcriptome_genes_after_mask"],
        "component_cone_median_cosine": donor_ranges["component_cone"]["cosine"]["median"],
        "component_cone_median_normalized_rmse": donor_ranges["component_cone"]["normalized_rmse"]["median"],
        "component_cone_strict_inside_count": strict_inside,
        "equal_sum_median_cosine": donor_ranges["equal_sum"]["cosine"]["median"],
        "equal_sum_median_normalized_rmse": donor_ranges["equal_sum"]["normalized_rmse"]["median"],
        "best_single_in_sample_oracle_median_cosine": donor_ranges["best_single_oracle"]["cosine"]["median"],
        "best_single_in_sample_oracle_median_normalized_rmse": donor_ranges["best_single_oracle"]["normalized_rmse"]["median"],
        "lodo_component_cone_median_cosine": donor_ranges["lodo_component_cone"]["cosine"]["median"],
        "lodo_component_cone_median_normalized_rmse": donor_ranges["lodo_component_cone"]["normalized_rmse"]["median"],
        "lodo_training_selected_best_single_median_cosine": donor_ranges["lodo_training_selected_best_single"]["cosine"]["median"],
        "lodo_training_selected_best_single_median_normalized_rmse": donor_ranges["lodo_training_selected_best_single"]["normalized_rmse"]["median"],
        "on_target_frozen_disjoint_median_cosine": donor_ranges["on_target_calibrated_component_cone"]["cosine"]["median"],
        "on_target_frozen_disjoint_median_normalized_rmse": donor_ranges["on_target_calibrated_component_cone"]["normalized_rmse"]["median"],
    }
    for field, expected_value in expected_primary.items():
        if isinstance(expected_value, int):
            if type(primary[field]) is not int or primary[field] != expected_value:
                raise AssertionError(f"Goudy findings primary.{field} differs")
        else:
            assert_close(finite_number(primary[field], f"Goudy findings primary.{field}"), expected_value)

    ledger_on_target = require_keys(
        ledger["on_target_controls"],
        {
            "all_three_targets_negative_in_single_and_triple_donors",
            "equal_sum_coordinate_cosine_median",
            "equal_sum_coordinate_cosine_range",
            "warning_status",
        },
        "Goudy findings on-target controls",
    )
    if (
        ledger_on_target["all_three_targets_negative_in_single_and_triple_donors"]
        != all_negative_donors
        or ledger_on_target["warning_status"] != warning
    ):
        raise AssertionError("Goudy findings on-target sign/warning differs")
    equal_sum_coordinates = saturated_ranges["equal_sum"]["cosine"]
    assert_close(
        ledger_on_target["equal_sum_coordinate_cosine_median"],
        equal_sum_coordinates["median"],
    )
    assert_range(
        ledger_on_target["equal_sum_coordinate_cosine_range"],
        equal_sum_coordinates["donor_range"],
        "Goudy findings on-target equal-sum range",
    )

    ledger_filter = require_keys(
        ledger["filter_sensitivity"],
        {
            "status",
            "control_mean_cpm_thresholds",
            "component_cone_median_cosine_by_threshold",
            "target_rows_removed_by_threshold",
            "conclusion",
        },
        "Goudy findings filter sensitivity",
    )
    if (
        ledger_filter["status"] != filter_sensitivity["filter_sensitivity_status"]
        or ledger_filter["control_mean_cpm_thresholds"] != thresholds
        or ledger_filter["target_rows_removed_by_threshold"] != target_rows_removed
        or ledger_filter["component_cone_median_cosine_by_threshold"]
        != threshold_component_medians
        or ledger_filter["conclusion"]
        != "No declared control-only expression threshold rescues the geometric model; SUV39H1 is already filtered at CPM 5 and 10, so only two target rows are removed there while scoring remains target-disjoint."
    ):
        raise AssertionError("Goudy findings filter-sensitivity contract differs")

    ledger_reliability = require_keys(
        ledger["reliability"],
        {
            "status",
            "triple_pairwise_donor_median_cosine",
            "component_pairwise_donor_median_cosine",
            "aavs1_guide_expression_median_cosine",
            "aavs1_guide_expression_median_normalized_rmse",
            "aavs1_difference_norm_over_triple_effect_norm_range",
        },
        "Goudy findings reliability",
    )
    if ledger_reliability["status"] != reliability["overall_status"]:
        raise AssertionError("Goudy findings reliability status differs")
    assert_close(
        ledger_reliability["triple_pairwise_donor_median_cosine"], pair_medians["triple"]
    )
    if ledger_reliability["component_pairwise_donor_median_cosine"] != {
        name: pair_medians[name] for name in component_order
    }:
        raise AssertionError("Goudy findings component reliability differs")
    assert_close(
        ledger_reliability["aavs1_guide_expression_median_cosine"],
        aavs1_ranges["cosine"]["median"],
    )
    assert_close(
        ledger_reliability["aavs1_guide_expression_median_normalized_rmse"],
        aavs1_ranges["normalized_rmse"]["median"],
    )
    assert_range(
        ledger_reliability["aavs1_difference_norm_over_triple_effect_norm_range"],
        [min(triple_noise_ratios), max(triple_noise_ratios)],
        "Goudy findings AAVS1/triple ratio",
    )

    ledger_residual = require_keys(
        ledger["additive_residual"],
        {
            "definition",
            "per_donor_relative_l2_median",
            "per_donor_relative_l2_range",
            "across_donor_mean_relative_l2",
            "across_donor_mean_absolute_residual",
            "module_error_status",
            "causal_interpretation",
        },
        "Goudy findings additive residual",
    )
    if (
        ledger_residual["definition"] != residual["definition"]
        or ledger_residual["module_error_status"] != module_error["status"]
        or ledger_residual["causal_interpretation"] != residual["causal_interpretation"]
    ):
        raise AssertionError("Goudy findings residual definition/status differs")
    assert_close(
        ledger_residual["per_donor_relative_l2_median"], statistics.median(residual_ratios)
    )
    assert_range(
        ledger_residual["per_donor_relative_l2_range"],
        [min(residual_ratios), max(residual_ratios)],
        "Goudy findings residual-ratio range",
    )
    assert_close(
        ledger_residual["across_donor_mean_relative_l2"],
        across_residual["residual_norm_over_triple_norm"],
    )
    assert_close(
        ledger_residual["across_donor_mean_absolute_residual"],
        across_residual["mean_absolute_residual"],
    )

    provenance = require_keys(
        ledger["provenance_limits"],
        {
            "author_key_confirmed_declared_roles",
            "declared_analysis_roles",
            "d3_d4_multiplex_control_and_triple_experiment_ids",
            "triple_guide_identity",
            "same_guide_match_claimed",
            "metadata_conflict",
        },
        "Goudy findings provenance limits",
    )
    if provenance != {
        "author_key_confirmed_declared_roles": author[
            "author_key_confirmed_declared_role_count"
        ],
        "declared_analysis_roles": author["declared_analysis_role_count"],
        "d3_d4_multiplex_control_and_triple_experiment_ids": "UNRESOLVED",
        "triple_guide_identity": author["triple_guide_identity"],
        "same_guide_match_claimed": author["same_guide_match_claimed"],
        "metadata_conflict": "Donor-1 GEO characteristics say Day18 while titles, paper methods, and pinned author metadata say Day7; the report uses Day7 and makes no duration claim.",
    }:
        raise AssertionError("Goudy findings provenance limits differ")

    expected_interpretation = (
        "Execution succeeds, but the declared geometric model fails. Weak transcriptome "
        "alignment, low donor reproducibility, filter sensitivity, and perfect "
        "single-versus-triple confounding by experiment, control type, and guide burden "
        "make the biological result inconclusive. This is not combination or additivity "
        "validation, same-guide matching, an interaction test, donor-population inference, "
        "prospective validation, or biological reachability evidence."
    )
    if (
        ledger["benchmark"] != benchmark
        or ledger["report_version"] != report["report_version"]
        or ledger["tier"] != claim["tier"]
        or ledger["benchmark_role"] != claim["benchmark_role"]
        or ledger["cell_system"] != claim["cell_system"]
        or ledger["assay"] != claim["assay"]
        or ledger["modality"] != claim["modality"]
        or ledger["fixed_donor_labels"] != len(donors)
        or ledger["status"] != statuses
        or ledger["claim_boundary"] != claim["claim_ceiling"]
        or ledger["interpretation"] != expected_interpretation
        or ledger["source"] != str(report_path.relative_to(ROOT))
    ):
        raise AssertionError("Goudy findings identity/claim/interpretation/source differs")


def validate_library_coverage(findings: dict) -> None:
    ledger = findings.get("library_coverage")
    if ledger is None:
        raise AssertionError("findings.library_coverage section is missing")

    report_path = RESULTS / "library_coverage_crossdataset.json"
    config_path = ROOT / "configs" / "library_coverage_crossdataset.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))

    def require_keys(value: object, expected: set[str], label: str) -> dict:
        if not isinstance(value, dict) or set(value) != expected:
            raise AssertionError(f"{label} schema differs")
        return value

    def finite_number(value: object, label: str) -> float:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise AssertionError(f"{label} is not finite numeric data")
        return float(value)

    def require_sha256(value: object, label: str) -> None:
        if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
            raise AssertionError(f"{label} is not a lowercase SHA-256")

    report_fields = {
        "schema_version",
        "generated_on",
        "status",
        "benchmark",
        "claim_ceiling",
        "provenance",
        "input_verification",
        "analysis_contract",
        "datasets",
        "config_sha256",
    }
    config_fields = {
        "schema_version",
        "generated_on",
        "benchmark",
        "claim_ceiling",
        "provenance",
        "cache_build",
        "inputs",
        "analysis",
        "datasets",
    }
    require_keys(report, report_fields, "library-coverage report")
    require_keys(config, config_fields, "library-coverage config")
    benchmark = "library_coverage_split_first_retrospective_v3"
    if (
        report["schema_version"] != "3.0.0"
        or config["schema_version"] != "3.0.0"
        or report["status"] != "PASS"
        or report["benchmark"] != benchmark
        or config["benchmark"] != benchmark
        or ledger.get("benchmark") != benchmark
    ):
        raise AssertionError("library-coverage v3 identity/status differs")
    if report["config_sha256"] != sha256(config_path):
        raise AssertionError("library-coverage report/config byte identity differs")
    for field, config_field in (
        ("generated_on", "generated_on"),
        ("claim_ceiling", "claim_ceiling"),
        ("provenance", "provenance"),
        ("analysis_contract", "analysis"),
    ):
        if report[field] != config[config_field]:
            raise AssertionError(f"library-coverage {field} differs from config")
    if (
        report["provenance"]["perturbation_modalities"] != ["CRISPRi", "CRISPRa"]
        or report["provenance"]["cell_systems"]
        != ["primary human CD4+ T cells", "K562"]
        or "does not require a minimum rank correlation or a top-1 match"
        not in report["analysis_contract"]["status_definition"]
        or "not sampling uncertainty" not in report["analysis_contract"]["inference"]
    ):
        raise AssertionError("library-coverage provenance/claim ceiling differs")

    # The frozen caches are portable outputs of an explicit, source-hash-gated builder.
    cache_build = require_keys(
        config["cache_build"],
        {"format", "writer_contract", "zhu", "norman", "replogle"},
        "library-coverage cache-build contract",
    )
    if (
        cache_build["format"] != "portable_effect_dictionary_v1"
        or cache_build["writer_contract"]
        != "pickle-free uncompressed NPZ with canonical little-endian C-contiguous arrays, fixed member order, ZIP timestamps, and permissions"
        or not (ROOT / "scripts" / "build_library_coverage_caches.py").is_file()
    ):
        raise AssertionError("library-coverage portable writer contract differs")
    zhu_build = require_keys(
        cache_build["zhu"],
        {"builder", "source_config", "layer", "gene_universe", "effect_rows"},
        "library-coverage Zhu cache builder",
    )
    if (
        zhu_build["builder"] != "zhu_de_layer"
        or zhu_build["source_config"] != "configs/source_reconstruction.json"
        or zhu_build["layer"] != "log_fc"
        or zhu_build["gene_universe"] != "shared_screen_genes"
        or "ontarget_significant=true" not in zhu_build["effect_rows"]
    ):
        raise AssertionError("library-coverage Zhu reconstruction contract differs")
    source_config_path = (ROOT / zhu_build["source_config"]).resolve()
    if not source_config_path.is_relative_to(ROOT.resolve()) or not source_config_path.is_file():
        raise AssertionError("library-coverage Zhu source config is unavailable or unsafe")
    source_config = json.loads(source_config_path.read_text(encoding="utf-8"))
    if (
        source_config["inputs"]["de_stats"]["sha256"]
        != config["provenance"]["zhu"]["upstream_object_sha256"]
    ):
        raise AssertionError("library-coverage Zhu upstream identity differs")

    raw_identities = {
        "norman": (
            "data/norman_perturbation.h5ad",
            2228849977,
            "a2a194c0eaa001d229e21a3f4f5b447db7c73f7a3a44c3d0464f317bda5f12a2",
        ),
        "replogle": (
            "data/replogle_k562_essential_perturbation.h5ad",
            2890786004,
            "04b7b3c28504ace115bb6ee192a0710f428aee3468cf0372b4fc0978cc05adb4",
        ),
    }
    cell_builder_fields = {
        "builder",
        "source",
        "condition_field",
        "control_field",
        "control_label",
        "gene_field",
        "gene_axis_note",
        "source_matrix_semantics",
    }
    source_fields = {
        "path",
        "bytes",
        "sha256",
        "url",
        "etag_observed",
        "last_modified_observed",
        "version_id_observed",
        "retrieved_on",
    }
    for name in ("norman", "replogle"):
        expected_fields = set(cell_builder_fields)
        if name == "replogle":
            expected_fields.add("feature_name_note")
        build = require_keys(
            cache_build[name], expected_fields, f"library-coverage {name} cache builder"
        )
        source = require_keys(
            build["source"], source_fields, f"library-coverage {name} raw source"
        )
        path, byte_count, digest = raw_identities[name]
        if (
            build["builder"] != "condition_minus_control_mean"
            or build["condition_field"] != "condition"
            or build["control_field"] != "control"
            or build["control_label"] != "ctrl"
            or build["gene_field"] != "gene_id"
            or source["path"] != path
            or source["bytes"] != byte_count
            or source["sha256"] != digest
            or "processed" not in build["source_matrix_semantics"]
            or "not raw counts" not in build["source_matrix_semantics"]
        ):
            raise AssertionError(f"library-coverage {name} raw reconstruction differs")
        require_sha256(source["sha256"], f"library-coverage {name} raw source")
        for field in source_fields - {"bytes"}:
            if not isinstance(source[field], str) or not source[field].strip():
                raise AssertionError(f"library-coverage {name} source.{field} is empty")
        relative_source = Path(source["path"])
        if relative_source.is_absolute() or ".." in relative_source.parts:
            raise AssertionError(f"library-coverage {name} source path is unsafe")

    portable_outputs = {
        "zhu": {
            "path": "slice/zhu_rest_effects_portable_v1.npz",
            "bytes": 246642804,
            "sha256": "75d7ce623ae8157c84a174b6dfc2a24baea34ce117e4ccb0ebceacfb174edba1",
            "arrays": {
                "E": {"shape": [6871, 8950], "dtype": "float32"},
                "perts": {"shape": [6871], "dtype": "<U11"},
                "genes": {"shape": [8950], "dtype": "<U10"},
            },
        },
        "norman": {
            "path": "slice/norman_effects_portable_v1.npz",
            "bytes": 6033788,
            "sha256": "636973133bb34b0bb7f0cc6c7e1ef1a799808c15cc1ea982b8097bb95a48d964",
            "arrays": {
                "E": {"shape": [283, 5045], "dtype": "float32"},
                "perts": {"shape": [283], "dtype": "<U15"},
                "genes": {"shape": [5045], "dtype": "<U15"},
                "ncells": {"shape": [283], "dtype": "int64"},
            },
        },
        "replogle": {
            "path": "slice/replogle_effects_portable_v1.npz",
            "bytes": 22110472,
            "sha256": "f48a92f7866ad1f055c9468c8d9ace289b5fc2e75b84ecf0aaf19d350aad25c3",
            "arrays": {
                "E": {"shape": [1087, 5000], "dtype": "float32"},
                "perts": {"shape": [1087], "dtype": "<U14"},
                "genes": {"shape": [5000], "dtype": "<U15"},
                "ncells": {"shape": [1087], "dtype": "int64"},
            },
        },
    }
    require_keys(config["inputs"], set(portable_outputs), "library-coverage inputs")
    require_keys(
        report["input_verification"],
        set(portable_outputs),
        "library-coverage input verification",
    )
    for name, expected_output in portable_outputs.items():
        registered = require_keys(
            config["inputs"][name],
            {"path", "bytes", "sha256", "effect_array", "label_array", "required_arrays"},
            f"library-coverage {name} portable cache",
        )
        if (
            registered["path"] != expected_output["path"]
            or registered["bytes"] != expected_output["bytes"]
            or registered["sha256"] != expected_output["sha256"]
            or registered["effect_array"] != "E"
            or registered["label_array"] != "perts"
            or registered["required_arrays"] != expected_output["arrays"]
        ):
            raise AssertionError(f"library-coverage {name} portable cache differs")
        identity = require_keys(
            report["input_verification"][name],
            {"path", "bytes", "sha256_expected", "sha256_actual", "hash_verified"},
            f"library-coverage {name} verified identity",
        )
        if (
            identity["hash_verified"] is not True
            or identity["path"] != registered["path"]
            or identity["bytes"] != registered["bytes"]
            or identity["sha256_expected"] != registered["sha256"]
            or identity["sha256_actual"] != registered["sha256"]
        ):
            raise AssertionError(f"library-coverage {name} input identity differs")
        require_sha256(identity["sha256_actual"], f"library-coverage {name} cache")

    analysis = report["analysis_contract"]
    threshold = finite_number(analysis["soft_cosine_threshold"], "soft threshold")
    sensitivity_contract = analysis["split_sensitivity"]
    seeds = sensitivity_contract["seeds"]
    thresholds = [float(value) for value in sensitivity_contract["soft_cosine_thresholds"]]
    if (
        len(seeds) != 12
        or len(set(seeds)) != 12
        or any(type(seed) is not int for seed in seeds)
        or thresholds != [0.3, 0.5, 0.7, 0.9]
    ):
        raise AssertionError("library-coverage sensitivity axes differ")

    coverage_fields = {
        "n_targets",
        "strict_inside_cone_count",
        "strict_inside_cone_fraction",
        "soft_cosine_threshold",
        "soft_coverage_count",
        "soft_coverage_fraction",
        "mean_cosine",
        "mean_residual_fraction",
    }
    alignment_fields = {
        "n_targets",
        "mean_cosine",
        "soft_cosine_threshold",
        "soft_coverage_count",
        "soft_coverage_fraction",
    }

    def validate_alignment(summary: dict, n_targets: int, label: str) -> None:
        if summary["n_targets"] != n_targets or summary["soft_cosine_threshold"] != threshold:
            raise AssertionError(f"{label} target/threshold contract differs")
        soft_count = summary["soft_coverage_count"]
        if type(soft_count) is not int or not 0 <= soft_count <= n_targets:
            raise AssertionError(f"{label} soft count differs")
        assert_close(
            finite_number(summary["soft_coverage_fraction"], f"{label} soft fraction"),
            soft_count / n_targets,
        )
        cosine = finite_number(summary["mean_cosine"], f"{label} mean cosine")
        if not -1.0 <= cosine <= 1.0:
            raise AssertionError(f"{label} mean cosine falls outside [-1, 1]")

    def validate_coverage(dataset: dict, label: str) -> int:
        coverage = require_keys(dataset["coverage"], coverage_fields, f"{label} coverage")
        n_targets = coverage["n_targets"]
        if type(n_targets) is not int or n_targets <= 0:
            raise AssertionError(f"{label} target count differs")
        validate_alignment(coverage, n_targets, f"{label} coverage")
        if (
            coverage["strict_inside_cone_count"] != 0
            or coverage["strict_inside_cone_fraction"] != 0.0
        ):
            raise AssertionError(f"{label} reference strict membership is not zero")
        finite_number(coverage["mean_residual_fraction"], f"{label} residual fraction")
        return n_targets

    acquisition_fields = {
        "candidate_status",
        "certificate_definition",
        "primary_realized_objective",
        "secondary_realized_objective",
        "candidate_pool_effect_atom_labels_unique",
        "candidate_pool_canonical_labels_unique",
        "top1_agreement_unit",
        "certificate_vs_realized_mean_cosine_gain_spearman",
        "certificate_top1",
        "certificate_top1_tie_count",
        "realized_mean_cosine_top1",
        "realized_mean_cosine_top1_tie_count",
        "certificate_matches_realized_mean_cosine_top1",
        "realized_soft_coverage_top1",
        "realized_soft_coverage_top1_tie_count",
    }
    full_top_fields = {
        "candidate_pool_position_zero_based",
        "source_effect_row_index_zero_based",
        "effect_atom_label",
        "canonical_perturbation_label",
        "certificate_score",
        "realized_mean_cosine_gain",
        "realized_soft_coverage_fraction_gain",
        "realized_strict_inside_cone_fraction_gain",
    }
    soft_top_fields = {
        "candidate_pool_position_zero_based",
        "source_effect_row_index_zero_based",
        "effect_atom_label",
        "canonical_perturbation_label",
        "realized_soft_coverage_fraction_gain",
        "realized_mean_cosine_gain",
    }

    def canonical_label(effect_label: str) -> str:
        parts = sorted(part for part in effect_label.split("+") if part != "ctrl")
        return "+".join(parts)

    def validate_acquisition(
        dataset: dict,
        candidate_count: int,
        source_rows: int,
        label: str,
    ) -> None:
        acquisition = require_keys(
            dataset["acquisition"], acquisition_fields, f"{label} acquisition"
        )
        if (
            acquisition["candidate_status"] != "supplied_already_measured_effect_atoms"
            or acquisition["candidate_pool_effect_atom_labels_unique"] is not True
            or acquisition["candidate_pool_canonical_labels_unique"] is not True
            or acquisition["top1_agreement_unit"]
            != "unique supplied effect atom; exactly one candidate exists per canonical label"
        ):
            raise AssertionError(f"{label} acquisition candidate contract differs")
        rho = finite_number(
            acquisition["certificate_vs_realized_mean_cosine_gain_spearman"],
            f"{label} acquisition Spearman",
        )
        if not -1.0 <= rho <= 1.0:
            raise AssertionError(f"{label} acquisition Spearman falls outside [-1, 1]")
        for name, fields in (
            ("certificate_top1", full_top_fields),
            ("realized_mean_cosine_top1", full_top_fields),
            ("realized_soft_coverage_top1", soft_top_fields),
        ):
            identity = require_keys(
                acquisition[name], fields, f"{label} acquisition.{name}"
            )
            position = identity["candidate_pool_position_zero_based"]
            source_index = identity["source_effect_row_index_zero_based"]
            if (
                type(position) is not int
                or not 0 <= position < candidate_count
                or type(source_index) is not int
                or not 0 <= source_index < source_rows
                or identity["canonical_perturbation_label"]
                != canonical_label(identity["effect_atom_label"])
            ):
                raise AssertionError(f"{label} acquisition.{name} identity differs")
            for field, value in identity.items():
                if field.endswith("_gain") and finite_number(
                    value, f"{label} acquisition.{name}.{field}"
                ) < -1e-8:
                    raise AssertionError(f"{label} acquisition realized gain is negative")
            if "certificate_score" in identity:
                finite_number(identity["certificate_score"], f"{label} certificate score")
        for tie_field in (
            "certificate_top1_tie_count",
            "realized_mean_cosine_top1_tie_count",
            "realized_soft_coverage_top1_tie_count",
        ):
            ties = acquisition[tie_field]
            if type(ties) is not int or not 1 <= ties <= candidate_count:
                raise AssertionError(f"{label} {tie_field} differs")
        observed_match = (
            acquisition["certificate_top1"]["candidate_pool_position_zero_based"]
            == acquisition["realized_mean_cosine_top1"][
                "candidate_pool_position_zero_based"
            ]
        )
        if acquisition["certificate_matches_realized_mean_cosine_top1"] is not observed_match:
            raise AssertionError(f"{label} acquisition top-1 agreement differs")

    datasets = require_keys(
        report["datasets"],
        {"zhu_crispri_tcell", "norman_k562_crispra", "replogle_k562_essential_crispri"},
        "library-coverage datasets",
    )
    standard_dataset_fields = {
        "effect_pool_rows",
        "design",
        "coverage",
        "comparators",
        "acquisition",
        "source_effect_rows",
        "preprocessing",
        "split_sensitivity",
    }
    standard_design_fields = {
        "rng",
        "seed",
        "library_atoms",
        "catalog_effect_rows",
        "supplied_candidate_effect_atoms",
        "unused_measured_effect_rows",
        "feature_count",
        "feature_selection",
        "feature_selection_rows",
        "feature_index_sha256",
        "library_source_row_index_sha256",
        "catalog_source_row_index_sha256",
        "candidate_source_row_index_sha256",
    }
    for name, input_name in (
        ("zhu_crispri_tcell", "zhu"),
        ("replogle_k562_essential_crispri", "replogle"),
    ):
        dataset = require_keys(datasets[name], standard_dataset_fields, name)
        design_config = config["datasets"][name]
        design = require_keys(dataset["design"], standard_design_fields, f"{name} design")
        source_rows = portable_outputs[input_name]["arrays"]["E"]["shape"][0]
        expected_unused = source_rows - (
            design_config["library_size"]
            + design_config["catalog_size"]
            + design_config["candidate_size"]
        )
        if (
            dataset["source_effect_rows"] != source_rows
            or dataset["effect_pool_rows"] != source_rows
            or design["rng"] != "numpy.default_rng(seed)"
            or design["seed"] != design_config["seed"]
            or design["library_atoms"] != design_config["library_size"]
            or design["catalog_effect_rows"] != design_config["catalog_size"]
            or design["supplied_candidate_effect_atoms"]
            != design_config["candidate_size"]
            or design["unused_measured_effect_rows"] != expected_unused
            or design["feature_count"] != design_config["n_features"]
            or design["feature_selection"] != design_config["feature_selection"]
            or design["feature_selection_rows"] != "current_library_only"
        ):
            raise AssertionError(f"{name} split-first reference design differs")
        for field in standard_design_fields:
            if field.endswith("sha256"):
                require_sha256(design[field], f"{name} design.{field}")
        expected_preprocessing = {"effect_definition": design_config["effect_definition"]}
        if name == "zhu_crispri_tcell":
            expected_preprocessing.update(
                {
                    "eligible_effect_pool": design_config["eligible_effect_pool"],
                    "outcome_based_effect_filtering": False,
                }
            )
        if dataset["preprocessing"] != expected_preprocessing:
            raise AssertionError(f"{name} preprocessing differs")
        n_targets = validate_coverage(dataset, name)
        if n_targets != design_config["catalog_size"]:
            raise AssertionError(f"{name} reference catalog size differs")
        comparators = require_keys(
            dataset["comparators"],
            {
                "best_single_library_atom_nonnegative_ray",
                "common_library_response_nonnegative_ray",
            },
            f"{name} comparators",
        )
        for comparator_name, comparator in comparators.items():
            require_keys(comparator, alignment_fields, f"{name} {comparator_name}")
            validate_alignment(comparator, n_targets, f"{name} {comparator_name}")
        validate_acquisition(
            dataset, design_config["candidate_size"], source_rows, name
        )

    norman_name = "norman_k562_crispra"
    norman = require_keys(
        datasets[norman_name],
        {
            "source_effect_rows",
            "single_effect_rows",
            "unique_single_genes",
            "double_effect_rows",
            "classification",
            "preprocessing",
            "coverage",
            "comparators",
            "acquisition_design",
            "acquisition",
            "split_sensitivity",
        },
        norman_name,
    )
    norman_config = config["datasets"][norman_name]
    norman_expected = norman_config["expected"]
    if (
        norman["source_effect_rows"] != 283
        or norman["single_effect_rows"] != norman_expected["single_effect_rows"]
        or norman["unique_single_genes"] != norman_expected["unique_single_genes"]
        or norman["double_effect_rows"] != norman_expected["double_effect_rows"]
        or norman["classification"]
        != "retrospective_single_to_double_additivity_and_coverage"
    ):
        raise AssertionError("Norman effect-row classification differs")
    norman_preprocessing = require_keys(
        norman["preprocessing"],
        {"effect_definition", "feature_count", "feature_selection", "feature_index_sha256"},
        "Norman preprocessing",
    )
    if (
        norman_preprocessing["effect_definition"] != norman_config["effect_definition"]
        or norman_preprocessing["feature_count"] != norman_config["n_features"]
        or norman_preprocessing["feature_selection"]
        != norman_config["feature_selection"]
    ):
        raise AssertionError("Norman full-single feature selection differs")
    require_sha256(norman_preprocessing["feature_index_sha256"], "Norman feature index")
    norman_targets = validate_coverage(norman, "Norman")
    if norman_targets != norman["double_effect_rows"]:
        raise AssertionError("Norman double catalog differs")
    norman_comparators = require_keys(
        norman["comparators"],
        {
            "paired_constituent_sum",
            "paired_constituent_two_atom_cone",
            "best_single_effect_nonnegative_ray",
            "common_single_response_nonnegative_ray",
        },
        "Norman comparators",
    )
    for name in ("best_single_effect_nonnegative_ray", "common_single_response_nonnegative_ray"):
        require_keys(norman_comparators[name], alignment_fields, f"Norman {name}")
        validate_alignment(norman_comparators[name], norman_targets, f"Norman {name}")
    pair_sum = require_keys(
        norman_comparators["paired_constituent_sum"],
        alignment_fields | {"constituent_selection", "opposite_position_fallback_count"},
        "Norman constituent sum",
    )
    pair_cone = require_keys(
        norman_comparators["paired_constituent_two_atom_cone"],
        alignment_fields
        | {
            "constituent_selection",
            "opposite_position_fallback_count",
            "strict_inside_cone_count",
            "strict_inside_cone_fraction",
        },
        "Norman paired cone",
    )
    for label, comparator in (("Norman constituent sum", pair_sum), ("Norman paired cone", pair_cone)):
        validate_alignment(comparator, norman_targets, label)
        if comparator["opposite_position_fallback_count"] != 1:
            raise AssertionError(f"{label} fallback count differs")
    if pair_cone["strict_inside_cone_count"] != 0 or pair_cone["strict_inside_cone_fraction"] != 0.0:
        raise AssertionError("Norman paired-cone strict membership is not zero")

    representative_fields = {
        "canonical_genes",
        "genes_with_both_cassette_positions",
        "opposite_position_fallbacks",
        "representative_rule",
    }
    expected_representatives = {
        "canonical_genes": 105,
        "genes_with_both_cassette_positions": 47,
        "opposite_position_fallbacks": 1,
        "representative_rule": (
            "choose measured GENE+ctrl; use ctrl+GENE only when GENE+ctrl is absent; "
            "choice is frozen before row-role assignment and does not use outcomes"
        ),
    }
    acquisition_design_fields = {
        "rng",
        "seed",
        "role_assignment_unit",
        "canonical_representative_selection",
        "library_candidate_canonical_gene_overlap_count",
        "current_library_canonical_genes",
        "supplied_candidate_canonical_genes",
        "unused_canonical_genes",
        "catalog_double_effect_rows",
        "feature_count",
        "feature_selection",
        "feature_selection_rows",
        "feature_index_sha256",
        "library_source_row_index_sha256",
        "candidate_source_row_index_sha256",
        "catalog_source_row_index_sha256",
        "library_canonical_labels_sha256",
        "candidate_canonical_labels_sha256",
    }
    acquisition_design = require_keys(
        norman["acquisition_design"], acquisition_design_fields, "Norman acquisition design"
    )
    representatives = require_keys(
        acquisition_design["canonical_representative_selection"],
        representative_fields,
        "Norman canonical representatives",
    )
    if representatives != expected_representatives:
        raise AssertionError("Norman fixed measured representatives differ")
    if (
        acquisition_design["rng"] != "numpy.default_rng(seed)"
        or acquisition_design["seed"] != norman_config["seed"]
        or acquisition_design["role_assignment_unit"]
        != "canonical_gene_with_one_fixed_measured_position_representative"
        or acquisition_design["library_candidate_canonical_gene_overlap_count"] != 0
        or acquisition_design["current_library_canonical_genes"] != 40
        or acquisition_design["current_library_canonical_genes"]
        != norman_config["acquisition_library_size"]
        or acquisition_design["supplied_candidate_canonical_genes"]
        != norman_config["acquisition_candidate_size"]
        or acquisition_design["unused_canonical_genes"] != 35
        or acquisition_design["catalog_double_effect_rows"] != norman["double_effect_rows"]
        or acquisition_design["feature_count"] != norman_config["n_features"]
        or acquisition_design["feature_selection"]
        != "stable top-variance coordinates from current canonical-gene acquisition library representatives only"
        or acquisition_design["feature_selection_rows"] != "current_library_only"
    ):
        raise AssertionError("Norman split-first acquisition roles differ")
    for field in acquisition_design_fields:
        if field.endswith("sha256"):
            require_sha256(acquisition_design[field], f"Norman acquisition {field}")
    validate_acquisition(
        norman,
        acquisition_design["supplied_candidate_canonical_genes"],
        norman["source_effect_rows"],
        "Norman",
    )

    distribution_fields = {
        "n_splits",
        "median",
        "q25",
        "q75",
        "p05",
        "p95",
        "minimum",
        "maximum",
        "positive_split_count",
    }
    metric_fields = (
        "cone_mean_cosine",
        "cone_strict_inside_fraction",
        "best_single_mean_cosine",
        "common_response_mean_cosine",
        "signed_span_mean_cosine",
        "cone_minus_best_single_mean_cosine",
        "cone_minus_common_response_mean_cosine",
        "signed_span_minus_cone_mean_cosine",
    )
    soft_keys = [f"cosine_at_least_{value:g}" for value in thresholds]
    aggregate_fields = set(metric_fields) | {
        f"cone_soft_coverage_{key}" for key in soft_keys
    }
    standard_protocol_fields = {
        "purpose",
        "rng",
        "seeds",
        "library_atoms",
        "catalog_effect_rows",
        "reserved_candidate_effect_atoms",
        "feature_count",
        "feature_selection_rows",
        "soft_cosine_thresholds",
        "acquisition_recomputed",
    }
    norman_protocol_fields = {
        "purpose",
        "rng",
        "seeds",
        "canonical_single_genes",
        "library_canonical_genes",
        "canonical_representative_selection",
        "eligible_double_effect_rows",
        "catalog_double_effect_rows",
        "catalog_seed",
        "catalog_role",
        "catalog_double_source_row_index_sha256",
        "feature_count",
        "feature_selection_rows",
        "soft_cosine_thresholds",
        "acquisition_recomputed",
    }
    base_split_fields = {
        "seed",
        "cone_mean_cosine",
        "cone_strict_inside_fraction",
        "cone_soft_coverage_by_threshold",
        "best_single_mean_cosine",
        "common_response_mean_cosine",
        "signed_span_mean_cosine",
        "cone_minus_best_single_mean_cosine",
        "cone_minus_common_response_mean_cosine",
        "signed_span_minus_cone_mean_cosine",
    }
    standard_hash_fields = {
        "library_source_row_index_sha256",
        "catalog_source_row_index_sha256",
        "candidate_source_row_index_sha256",
        "unused_source_row_index_sha256",
        "feature_index_sha256",
    }
    norman_hash_fields = {
        "library_representative_source_row_index_sha256",
        "unused_representative_source_row_index_sha256",
        "library_canonical_labels_sha256",
        "catalog_double_source_row_index_sha256",
        "feature_index_sha256",
    }

    def validate_sensitivity(dataset_name: str, *, norman_design: bool) -> dict:
        sensitivity = require_keys(
            datasets[dataset_name]["split_sensitivity"],
            {"protocol", "splits", "aggregate"},
            f"{dataset_name} split sensitivity",
        )
        expected_protocol_fields = (
            norman_protocol_fields if norman_design else standard_protocol_fields
        )
        protocol = require_keys(
            sensitivity["protocol"],
            expected_protocol_fields,
            f"{dataset_name} sensitivity protocol",
        )
        if (
            protocol["rng"] != "numpy.default_rng(seed)"
            or protocol["seeds"] != seeds
            or protocol["soft_cosine_thresholds"] != thresholds
            or protocol["acquisition_recomputed"] is not False
            or protocol["feature_count"] != 400
        ):
            raise AssertionError(f"{dataset_name} sensitivity protocol differs")
        if norman_design:
            if (
                protocol["canonical_single_genes"] != 105
                or protocol["library_canonical_genes"] != 40
                or protocol["canonical_representative_selection"]
                != expected_representatives
                or protocol["eligible_double_effect_rows"] != 131
                or protocol["catalog_double_effect_rows"] != 60
                or protocol["catalog_seed"]
                != config["datasets"][dataset_name]["sensitivity_catalog_seed"]
                or protocol["catalog_role"]
                != "fixed across every library-partition seed"
                or protocol["feature_selection_rows"]
                != "current 40-gene measured-representative library only"
            ):
                raise AssertionError("Norman fixed-catalog sensitivity design differs")
            require_sha256(
                protocol["catalog_double_source_row_index_sha256"],
                "Norman fixed sensitivity catalog",
            )
        else:
            design = config["datasets"][dataset_name]
            if (
                protocol["library_atoms"] != design["library_size"]
                or protocol["catalog_effect_rows"] != design["sensitivity_catalog_size"]
                or protocol["reserved_candidate_effect_atoms"]
                != design["sensitivity_candidate_reserve"]
                or protocol["feature_selection_rows"] != "current_library_only"
            ):
                raise AssertionError(f"{dataset_name} split-first sensitivity differs")

        rows = sensitivity["splits"]
        if not isinstance(rows, list) or len(rows) != len(seeds):
            raise AssertionError(f"{dataset_name} sensitivity split count differs")
        expected_hash_fields = norman_hash_fields if norman_design else standard_hash_fields
        for expected_seed, row in zip(seeds, rows, strict=True):
            require_keys(
                row,
                base_split_fields | expected_hash_fields,
                f"{dataset_name} sensitivity split {expected_seed}",
            )
            if row["seed"] != expected_seed:
                raise AssertionError(f"{dataset_name} sensitivity seed order differs")
            for field in expected_hash_fields:
                require_sha256(row[field], f"{dataset_name} split {expected_seed}.{field}")
            if norman_design and (
                row["catalog_double_source_row_index_sha256"]
                != protocol["catalog_double_source_row_index_sha256"]
            ):
                raise AssertionError("Norman sensitivity catalog is not fixed")
            values = {
                field: finite_number(
                    row[field], f"{dataset_name} split {expected_seed}.{field}"
                )
                for field in metric_fields
            }
            if values["cone_strict_inside_fraction"] != 0.0:
                raise AssertionError(
                    f"{dataset_name} split {expected_seed} has strict membership"
                )
            assert_close(
                values["cone_minus_best_single_mean_cosine"],
                values["cone_mean_cosine"] - values["best_single_mean_cosine"],
            )
            assert_close(
                values["cone_minus_common_response_mean_cosine"],
                values["cone_mean_cosine"] - values["common_response_mean_cosine"],
            )
            assert_close(
                values["signed_span_minus_cone_mean_cosine"],
                values["signed_span_mean_cosine"] - values["cone_mean_cosine"],
            )
            if (
                values["cone_minus_best_single_mean_cosine"] <= 0.0
                or values["cone_minus_common_response_mean_cosine"] <= 0.0
                or values["signed_span_minus_cone_mean_cosine"] <= 0.0
            ):
                raise AssertionError(
                    f"{dataset_name} split {expected_seed} violates nested capacity"
                )
            soft = require_keys(
                row["cone_soft_coverage_by_threshold"],
                set(soft_keys),
                f"{dataset_name} split {expected_seed} soft thresholds",
            )
            fractions = [
                finite_number(soft[key], f"{dataset_name} split {expected_seed}.{key}")
                for key in soft_keys
            ]
            if any(not 0.0 <= value <= 1.0 for value in fractions) or any(
                left < right for left, right in zip(fractions, fractions[1:])
            ):
                raise AssertionError(
                    f"{dataset_name} split {expected_seed} soft coverage is not monotone"
                )

        aggregate = require_keys(
            sensitivity["aggregate"], aggregate_fields, f"{dataset_name} aggregate"
        )
        for field in aggregate_fields:
            observed = require_keys(
                aggregate[field],
                distribution_fields,
                f"{dataset_name} aggregate.{field}",
            )
            if field.startswith("cone_soft_coverage_"):
                row_key = field.removeprefix("cone_soft_coverage_")
                values = [
                    float(row["cone_soft_coverage_by_threshold"][row_key]) for row in rows
                ]
            else:
                values = [float(row[field]) for row in rows]
            array = np.asarray(values, dtype=float)
            expected_distribution = {
                "n_splits": len(rows),
                "median": float(np.median(array)),
                "q25": float(np.quantile(array, 0.25)),
                "q75": float(np.quantile(array, 0.75)),
                "p05": float(np.quantile(array, 0.05)),
                "p95": float(np.quantile(array, 0.95)),
                "minimum": float(array.min()),
                "maximum": float(array.max()),
                "positive_split_count": int(np.count_nonzero(array > 0)),
            }
            if observed["n_splits"] != expected_distribution["n_splits"] or observed[
                "positive_split_count"
            ] != expected_distribution["positive_split_count"]:
                raise AssertionError(f"{dataset_name} aggregate.{field} counts differ")
            for summary_field in distribution_fields - {"n_splits", "positive_split_count"}:
                assert_close(
                    finite_number(
                        observed[summary_field],
                        f"{dataset_name} aggregate.{field}.{summary_field}",
                    ),
                    expected_distribution[summary_field],
                )
        strict = aggregate["cone_strict_inside_fraction"]
        if strict["positive_split_count"] != 0 or any(
            strict[field] != 0.0
            for field in distribution_fields - {"n_splits", "positive_split_count"}
        ):
            raise AssertionError(f"{dataset_name} strict aggregate is not identically zero")
        for field in (
            "cone_minus_best_single_mean_cosine",
            "cone_minus_common_response_mean_cosine",
            "signed_span_minus_cone_mean_cosine",
        ):
            if (
                aggregate[field]["positive_split_count"] != len(rows)
                or aggregate[field]["minimum"] <= 0.0
            ):
                raise AssertionError(f"{dataset_name} nested-capacity aggregate differs")
        return sensitivity

    sensitivities = {
        "zhu_crispri_tcell": validate_sensitivity(
            "zhu_crispri_tcell", norman_design=False
        ),
        "norman_k562_crispra": validate_sensitivity(
            "norman_k562_crispra", norman_design=True
        ),
        "replogle_k562_essential_crispri": validate_sensitivity(
            "replogle_k562_essential_crispri", norman_design=False
        ),
    }

    def sensitivity_projection(name: str, *, norman_design: bool) -> dict:
        sensitivity = sensitivities[name]
        protocol = sensitivity["protocol"]
        aggregate = sensitivity["aggregate"]
        projected = {
            "cone_mean_cosine_median": aggregate["cone_mean_cosine"]["median"],
            "cone_mean_cosine_range": [
                aggregate["cone_mean_cosine"]["minimum"],
                aggregate["cone_mean_cosine"]["maximum"],
            ],
            "soft_coverage_fraction_at_0_5_median": aggregate[
                "cone_soft_coverage_cosine_at_least_0.5"
            ]["median"],
            "soft_coverage_fraction_at_0_5_range": [
                aggregate["cone_soft_coverage_cosine_at_least_0.5"]["minimum"],
                aggregate["cone_soft_coverage_cosine_at_least_0.5"]["maximum"],
            ],
            "strict_inside_positive_splits": aggregate[
                "cone_strict_inside_fraction"
            ]["positive_split_count"],
            "cone_minus_best_single_mean_cosine_median": aggregate[
                "cone_minus_best_single_mean_cosine"
            ]["median"],
            "cone_minus_common_response_mean_cosine_median": aggregate[
                "cone_minus_common_response_mean_cosine"
            ]["median"],
            "signed_span_minus_cone_mean_cosine_median": aggregate[
                "signed_span_minus_cone_mean_cosine"
            ]["median"],
        }
        if norman_design:
            return {
                "fixed_catalog_double_effect_rows": protocol[
                    "catalog_double_effect_rows"
                ],
                "library_canonical_genes_per_split": protocol[
                    "library_canonical_genes"
                ],
                **projected,
            }
        return {"catalog_targets_per_split": protocol["catalog_effect_rows"], **projected}

    zhu = datasets["zhu_crispri_tcell"]
    replogle = datasets["replogle_k562_essential_crispri"]
    expected_dataset_ledger = {
        "zhu_crispri_tcell": {
            "source_effect_rows": zhu["source_effect_rows"],
            "effect_pool_rows": zhu["effect_pool_rows"],
            "catalog_targets": zhu["coverage"]["n_targets"],
            "coverage_mean_cosine": zhu["coverage"]["mean_cosine"],
            "strict_inside_cone_count": zhu["coverage"]["strict_inside_cone_count"],
            "soft_coverage_count": zhu["coverage"]["soft_coverage_count"],
            "best_single_mean_cosine": zhu["comparators"]
            ["best_single_library_atom_nonnegative_ray"]["mean_cosine"],
            "best_single_soft_coverage_count": zhu["comparators"]
            ["best_single_library_atom_nonnegative_ray"]["soft_coverage_count"],
            "certificate_vs_realized_mean_cosine_gain_spearman": zhu["acquisition"]
            ["certificate_vs_realized_mean_cosine_gain_spearman"],
            "certificate_top1": zhu["acquisition"]["certificate_top1"][
                "effect_atom_label"
            ],
            "realized_mean_cosine_top1": zhu["acquisition"]
            ["realized_mean_cosine_top1"]["effect_atom_label"],
            "certificate_matches_realized_mean_cosine_top1": zhu["acquisition"]
            ["certificate_matches_realized_mean_cosine_top1"],
            "sensitivity": sensitivity_projection(
                "zhu_crispri_tcell", norman_design=False
            ),
        },
        "norman_k562_crispra": {
            "source_effect_rows": norman["source_effect_rows"],
            "single_effect_rows": norman["single_effect_rows"],
            "canonical_single_genes": norman["unique_single_genes"],
            "double_effect_rows": norman["double_effect_rows"],
            "all_singles_cone_mean_cosine": norman["coverage"]["mean_cosine"],
            "strict_inside_cone_count": norman["coverage"]["strict_inside_cone_count"],
            "soft_coverage_count": norman["coverage"]["soft_coverage_count"],
            "paired_constituent_sum_mean_cosine": pair_sum["mean_cosine"],
            "paired_constituent_two_atom_cone_mean_cosine": pair_cone["mean_cosine"],
            "best_single_mean_cosine": norman_comparators[
                "best_single_effect_nonnegative_ray"
            ]["mean_cosine"],
            "best_single_soft_coverage_count": norman_comparators[
                "best_single_effect_nonnegative_ray"
            ]["soft_coverage_count"],
            "common_response_mean_cosine": norman_comparators[
                "common_single_response_nonnegative_ray"
            ]["mean_cosine"],
            "acquisition_library_canonical_genes": acquisition_design[
                "current_library_canonical_genes"
            ],
            "acquisition_library_candidate_gene_overlap_count": acquisition_design[
                "library_candidate_canonical_gene_overlap_count"
            ],
            "certificate_vs_realized_mean_cosine_gain_spearman": norman[
                "acquisition"
            ]["certificate_vs_realized_mean_cosine_gain_spearman"],
            "certificate_top1": norman["acquisition"]["certificate_top1"][
                "effect_atom_label"
            ],
            "realized_mean_cosine_top1": norman["acquisition"]
            ["realized_mean_cosine_top1"]["effect_atom_label"],
            "certificate_matches_realized_mean_cosine_top1": norman["acquisition"]
            ["certificate_matches_realized_mean_cosine_top1"],
            "sensitivity": sensitivity_projection(
                "norman_k562_crispra", norman_design=True
            ),
        },
        "replogle_k562_essential_crispri": {
            "source_effect_rows": replogle["source_effect_rows"],
            "catalog_targets": replogle["coverage"]["n_targets"],
            "coverage_mean_cosine": replogle["coverage"]["mean_cosine"],
            "strict_inside_cone_count": replogle["coverage"]
            ["strict_inside_cone_count"],
            "soft_coverage_count": replogle["coverage"]["soft_coverage_count"],
            "best_single_mean_cosine": replogle["comparators"]
            ["best_single_library_atom_nonnegative_ray"]["mean_cosine"],
            "best_single_soft_coverage_count": replogle["comparators"]
            ["best_single_library_atom_nonnegative_ray"]["soft_coverage_count"],
            "certificate_vs_realized_mean_cosine_gain_spearman": replogle[
                "acquisition"
            ]["certificate_vs_realized_mean_cosine_gain_spearman"],
            "certificate_top1": replogle["acquisition"]["certificate_top1"][
                "effect_atom_label"
            ],
            "realized_mean_cosine_top1": replogle["acquisition"]
            ["realized_mean_cosine_top1"]["effect_atom_label"],
            "certificate_matches_realized_mean_cosine_top1": replogle["acquisition"]
            ["certificate_matches_realized_mean_cosine_top1"],
            "sensitivity": sensitivity_projection(
                "replogle_k562_essential_crispri", norman_design=False
            ),
        },
    }

    def compare_projection(observed: object, expected: object, label: str) -> None:
        if isinstance(expected, dict):
            if not isinstance(observed, dict) or set(observed) != set(expected):
                raise AssertionError(f"{label} schema differs")
            for key, value in expected.items():
                compare_projection(observed[key], value, f"{label}.{key}")
        elif isinstance(expected, list):
            if not isinstance(observed, list) or len(observed) != len(expected):
                raise AssertionError(f"{label} list differs")
            for index, (actual, value) in enumerate(zip(observed, expected, strict=True)):
                compare_projection(actual, value, f"{label}[{index}]")
        elif isinstance(expected, bool):
            if observed is not expected:
                raise AssertionError(f"{label} differs")
        elif isinstance(expected, int):
            if type(observed) is not int or observed != expected:
                raise AssertionError(f"{label} differs")
        elif isinstance(expected, float):
            assert_close(finite_number(observed, label), expected)
        elif observed != expected:
            raise AssertionError(f"{label} differs")

    expected_ledger_fields = {
        "tier",
        "benchmark",
        "description",
        "reproducibility_boundary",
        "soft_cosine_threshold",
        "split_sensitivity",
        "datasets",
        "synthetic_contract",
        "interpretation",
        "source",
    }
    require_keys(ledger, expected_ledger_fields, "library-coverage findings ledger")
    expected_split_ledger = {
        "deterministic_splits": len(seeds),
        "purpose": (
            "algorithmic library/row-partition sensitivity over correlated measured "
            "effects; not biological sampling uncertainty"
        ),
        "acquisition_recomputed": False,
    }
    if (
        ledger["tier"] != "EXTENSION"
        or ledger["benchmark"] != benchmark
        or ledger["description"]
        != "Split-first catalog bookkeeping over reachability.project_cone: strict point-estimate cone membership, cosine-threshold alignment, simple rays, a signed-span capacity ceiling, and retrospective scoring of supplied, already-measured candidate effects."
        or ledger["reproducibility_boundary"]
        != "The three pickle-free cache artifacts are rebuilt byte-for-byte by a hash-gated deterministic builder. The registered Norman/Replogle URLs remain mutable and historical S3 versions are not anonymously retrievable, so a durable release archive is still required for manuscript-grade availability."
        or ledger["synthetic_contract"]
        != "Data-free tests distinguish strict from thresholded coverage, validate weighted gap normals, enforce nested signed-span/cone/ray capacity, require non-negative realized gains, and keep certificate order separate from realized order."
        or ledger["interpretation"]
        != "Strict membership is zero in every reference catalog and every one of the 12 sensitivity splits; the 0.5 cosine bar is descriptive alignment, not reachability. The cone beats best-single and common-response rays in every split, while the signed span is uniformly better, exposing the cost of non-negativity. Certificate scores correlate with realized mean-cosine gains (Spearman 0.861-0.921) but match the realized top supplied effect atom in 0/3 audits. Zhu and Replogle remain same-screen compression tests; Norman's full-single reference and 40-gene representative sensitivity designs are retrospective additivity/alignment diagnostics. No candidate is unmeasured, and no result is prospective design or a biological verdict."
        or ledger["source"] != str(report_path.relative_to(ROOT))
    ):
        raise AssertionError("library-coverage findings claim boundary differs")
    assert_close(
        finite_number(ledger["soft_cosine_threshold"], "findings soft threshold"),
        threshold,
    )
    compare_projection(
        ledger["split_sensitivity"], expected_split_ledger, "findings split sensitivity"
    )
    compare_projection(
        ledger["datasets"], expected_dataset_ledger, "findings library datasets"
    )

    zhu_coverage = zhu["coverage"]
    zhu_best = zhu["comparators"]["best_single_library_atom_nonnegative_ray"]
    norman_coverage = norman["coverage"]
    replogle_coverage = replogle["coverage"]
    replogle_best = replogle["comparators"]["best_single_library_atom_nonnegative_ray"]
    match_count = sum(
        dataset["acquisition"]["certificate_matches_realized_mean_cosine_top1"]
        for dataset in datasets.values()
    )
    zhu_sensitivity = expected_dataset_ledger["zhu_crispri_tcell"]["sensitivity"]
    norman_sensitivity = expected_dataset_ledger["norman_k562_crispra"]["sensitivity"]
    replogle_sensitivity = expected_dataset_ledger[
        "replogle_k562_essential_crispri"
    ]["sensitivity"]
    expected_summary = f"""<!-- BEGIN VALIDATED LIBRARY COVERAGE SUMMARY -->
| Audit | Cone mean cosine | Strict membership | Cosine ≥0.5 | Simple comparator | Certificate vs realized gain |
|---|---:|---:|---:|---|---|
| Zhu CD4 CRISPRi | {zhu_coverage['mean_cosine']:.3f} | {zhu_coverage['strict_inside_cone_count']}/{zhu_coverage['n_targets']} | {zhu_coverage['soft_coverage_count']}/{zhu_coverage['n_targets']} | Best atom: {zhu_best['mean_cosine']:.3f}; {zhu_best['soft_coverage_count']}/{zhu_best['n_targets']} | ρ={zhu['acquisition']['certificate_vs_realized_mean_cosine_gain_spearman']:.3f}; top-1 differs |
| Norman K562 CRISPRa | {norman_coverage['mean_cosine']:.3f} | {norman_coverage['strict_inside_cone_count']}/{norman_coverage['n_targets']} | {norman_coverage['soft_coverage_count']}/{norman_coverage['n_targets']} | Constituent sum: {pair_sum['mean_cosine']:.3f}; two-atom cone: {pair_cone['mean_cosine']:.3f} | ρ={norman['acquisition']['certificate_vs_realized_mean_cosine_gain_spearman']:.3f}; top-1 differs |
| Replogle K562 CRISPRi | {replogle_coverage['mean_cosine']:.3f} | {replogle_coverage['strict_inside_cone_count']}/{replogle_coverage['n_targets']} | {replogle_coverage['soft_coverage_count']}/{replogle_coverage['n_targets']} | Best atom: {replogle_best['mean_cosine']:.3f}; {replogle_best['soft_coverage_count']}/{replogle_best['n_targets']} | ρ={replogle['acquisition']['certificate_vs_realized_mean_cosine_gain_spearman']:.3f}; top-1 differs |

| 12-split sensitivity | Cone cosine median [range] | Cosine ≥0.5 median [range] | Strict-positive splits | Cone − best atom | Signed span − cone |
|---|---:|---:|---:|---:|---:|
| Zhu CD4 CRISPRi | {zhu_sensitivity['cone_mean_cosine_median']:.3f} [{zhu_sensitivity['cone_mean_cosine_range'][0]:.3f}, {zhu_sensitivity['cone_mean_cosine_range'][1]:.3f}] | {zhu_sensitivity['soft_coverage_fraction_at_0_5_median']:.3f} [{zhu_sensitivity['soft_coverage_fraction_at_0_5_range'][0]:.3f}, {zhu_sensitivity['soft_coverage_fraction_at_0_5_range'][1]:.3f}] | {zhu_sensitivity['strict_inside_positive_splits']}/{len(seeds)} | +{zhu_sensitivity['cone_minus_best_single_mean_cosine_median']:.3f} | +{zhu_sensitivity['signed_span_minus_cone_mean_cosine_median']:.3f} |
| Norman K562 CRISPRa | {norman_sensitivity['cone_mean_cosine_median']:.3f} [{norman_sensitivity['cone_mean_cosine_range'][0]:.3f}, {norman_sensitivity['cone_mean_cosine_range'][1]:.3f}] | {norman_sensitivity['soft_coverage_fraction_at_0_5_median']:.3f} [{norman_sensitivity['soft_coverage_fraction_at_0_5_range'][0]:.3f}, {norman_sensitivity['soft_coverage_fraction_at_0_5_range'][1]:.3f}] | {norman_sensitivity['strict_inside_positive_splits']}/{len(seeds)} | +{norman_sensitivity['cone_minus_best_single_mean_cosine_median']:.3f} | +{norman_sensitivity['signed_span_minus_cone_mean_cosine_median']:.3f} |
| Replogle K562 CRISPRi | {replogle_sensitivity['cone_mean_cosine_median']:.3f} [{replogle_sensitivity['cone_mean_cosine_range'][0]:.3f}, {replogle_sensitivity['cone_mean_cosine_range'][1]:.3f}] | {replogle_sensitivity['soft_coverage_fraction_at_0_5_median']:.3f} [{replogle_sensitivity['soft_coverage_fraction_at_0_5_range'][0]:.3f}, {replogle_sensitivity['soft_coverage_fraction_at_0_5_range'][1]:.3f}] | {replogle_sensitivity['strict_inside_positive_splits']}/{len(seeds)} | +{replogle_sensitivity['cone_minus_best_single_mean_cosine_median']:.3f} | +{replogle_sensitivity['signed_span_minus_cone_mean_cosine_median']:.3f} |

Strict membership is absent from every reference catalog and every sensitivity split, so a
soft directional bar cannot be renamed reachability. The cone beats best-single and
common-response rays in all {3 * len(seeds)} splits, while the signed span is uniformly better: useful
alignment exists, but non-negativity is a real capacity constraint. Rank correlations are
strong, yet certificate and realized top candidates disagree in {3 - match_count}/3 audits. Norman's
{norman_coverage['mean_cosine']:.3f} reference uses all {norman['single_effect_rows']} single rows; its {norman_sensitivity['cone_mean_cosine_median']:.3f} sensitivity median uses {norman_sensitivity['library_canonical_genes_per_split']} measured
canonical-gene representatives and is the appropriate partition-robustness result. Even in
the reference design, constituent-only baselines already clear the soft bar for {pair_sum['soft_coverage_count']}/{pair_sum['n_targets']}
doubles, so the full cone adds alignment rather than establishing a biological manifold.
<!-- END VALIDATED LIBRARY COVERAGE SUMMARY -->"""
    findings_text = (ROOT / "docs" / "FINDINGS.md").read_text(encoding="utf-8")
    if findings_text.count("<!-- BEGIN VALIDATED LIBRARY COVERAGE SUMMARY -->") != 1:
        raise AssertionError(
            "validated library-coverage findings block is missing or duplicated"
        )
    if expected_summary not in findings_text:
        raise AssertionError("library-coverage findings summary differs from v3 report")


def validate_donor_pair(ledger: dict) -> None:
    report = json.loads((RESULTS / "donor_pair_transfer.json").read_text(encoding="utf-8"))
    config_path = ROOT / "configs" / "donor_pair_transfer.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if report.get("status") != "PASS" or report["benchmark"] != ledger["benchmark"]:
        raise AssertionError("donor-pair benchmark did not pass")
    if report["config_sha256"] != sha256(config_path):
        raise AssertionError("donor-pair report/config identity differs")
    for name, identity in report["input_verification"].items():
        if (
            not identity["hash_verified"]
            or identity["sha256_actual"] != identity["sha256_expected"]
            or identity["sha256_expected"] != config["inputs"][name]["sha256"]
        ):
            raise AssertionError("donor-pair input identity was not fully verified")

    quality = report["data_quality"]
    if len(quality["modalities"]) != ledger["modalities"]:
        raise AssertionError("donor-pair modality count differs")
    if quality["rest_atoms_complete_all_modalities"] != ledger["rest_atoms_complete_all_modalities"]:
        raise AssertionError("donor-pair complete Rest atom count differs")
    if quality["shared_target_genes"] != ledger["shared_target_genes"]:
        raise AssertionError("donor-pair shared target universe differs")
    if "published two-donor DE-eligible" not in quality["eligibility_warning"]:
        raise AssertionError("donor-pair eligibility-selection warning is missing")

    challenges = report["challenges"]
    balanced = [row for row in challenges if not row["run_confounded"]]
    if len(challenges) != ledger["correlated_challenges"] or len(balanced) != ledger["run_balanced_challenges"]:
        raise AssertionError("donor-pair challenge counts differ")
    for row in challenges:
        for field in ("fit_gene_sha256", "score_gene_sha256"):
            if len(row[field]) != 64:
                raise AssertionError("donor-pair gene split identities are not hash-frozen")
    if "select atom/scalar only on training" not in report["protocol"]["baseline_application"]:
        raise AssertionError("donor-pair baselines are not declared training-only")
    if "no donor-population p-values" not in report["protocol"]["inference"]:
        raise AssertionError("donor-pair inference ceiling is missing")

    observed = report["summary"]["run_balanced_only"]
    expected = ledger["run_balanced"]
    if observed["challenges"] != ledger["run_balanced_challenges"]:
        raise AssertionError("donor-pair run-balanced summary count differs")
    for field in (
        "median_cone_cosine",
        "median_cosine_improvement_over_training_best_single",
        "fraction_cosine_improvement_positive_over_training_best_single",
        "median_cone_normalized_rmse",
        "median_training_best_single_normalized_rmse",
        "fraction_normalized_rmse_improvement_positive_over_training_best_single",
    ):
        assert_close(observed[field], expected[field])


def validate_guide_pair(ledger: dict) -> None:
    """Validate the frozen released guide-rank reciprocal-transfer stress."""

    report_path = RESULTS / "guide_pair_transfer.json"
    config_path = ROOT / "configs" / "guide_pair_transfer.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))

    benchmark = "zhu_released_guide_position_reciprocal_transfer_v1"
    identity_status = (
        "OFFICIAL_ALPHANUMERIC_GUIDE_RANKS_IDS_NOT_EMBEDDED_CROSSWALK_NOT_VERIFIED"
    )
    claim_ceiling = (
        "same-study reciprocal robustness across the released guide_1 and guide_2 "
        "alphanumeric guide-rank modalities under the authors' effectiveness-selected "
        "DE pipeline; not leakage-safe physical-guide generalization, named-sgRNA "
        "replication, donor-population inference, predictive utility, or state "
        "reachability"
    )
    modalities = ("guide_1", "guide_2")
    sources = ("hollbacker", "ota")
    seeds = (0, 1, 2)
    models = ("cone", "training_common_ray", "training_best_single", "zero")
    baselines = ("zero", "training_common_ray", "training_best_single")
    metrics = ("cosine", "normalized_rmse", "norm_ratio", "sign_agreement")

    def require_keys(value: object, expected: set[str], label: str) -> dict:
        if not isinstance(value, dict) or set(value) != expected:
            raise AssertionError(f"{label} schema differs")
        return value

    def finite(value: object, label: str) -> float:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
        ):
            raise AssertionError(f"{label} is not finite numeric data")
        return float(value)

    def hash_value(value: object, label: str) -> str:
        if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
            raise AssertionError(f"{label} is not a lowercase SHA-256")
        return value

    def canonical_hash(value: object) -> str:
        payload = json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    def scalar_vector_hash(value: float) -> str:
        array = np.asarray([value], dtype="<f8")
        return hashlib.sha256(array.tobytes(order="C")).hexdigest()

    def distribution(values: list[float]) -> dict[str, object]:
        if not values or not all(math.isfinite(value) for value in values):
            raise AssertionError("guide summary source values are not finite")
        return {
            "median": float(statistics.median(values)),
            "range": [float(min(values)), float(max(values))],
            "fraction_positive": float(
                statistics.mean(value > 0 for value in values)
            ),
        }

    def domain_summary(rows: list[dict], domain: str) -> dict[str, object]:
        return {
            "metrics": {
                model: {
                    metric: distribution(
                        [float(row[domain]["metrics"][model][metric]) for row in rows]
                    )
                    for metric in ("cosine", "normalized_rmse")
                }
                for model in models
            },
            "cone_improvement_over_baselines": {
                baseline: {
                    metric: distribution(
                        [
                            float(row[domain]["comparisons"][baseline][metric])
                            for row in rows
                        ]
                    )
                    for metric in (
                        "cosine_improvement",
                        "normalized_rmse_improvement",
                    )
                }
                for baseline in baselines
            },
        }

    def summary_scope(rows: list[dict]) -> dict[str, object]:
        return {
            "challenge_rows": len(rows),
            "unique_fits": len({row["fit_id"] for row in rows}),
            "within_training_guide": domain_summary(rows, "within_training_guide"),
            "reciprocal_held_guide": domain_summary(
                rows, "reciprocal_held_guide"
            ),
            "reciprocal_minus_within": {
                "metrics": {
                    model: {
                        metric: distribution(
                            [
                                float(
                                    row["reciprocal_minus_within"]["metrics"][
                                        model
                                    ][metric]
                                )
                                for row in rows
                            ]
                        )
                        for metric in metrics
                    }
                    for model in models
                },
                "comparisons": {
                    baseline: {
                        metric: distribution(
                            [
                                float(
                                    row["reciprocal_minus_within"]["comparisons"]
                                    [baseline][metric]
                                )
                                for row in rows
                            ]
                        )
                        for metric in (
                            "cosine_improvement",
                            "normalized_rmse_improvement",
                        )
                    }
                    for baseline in baselines
                },
            },
            "prediction_cosine": {
                model: distribution(
                    [float(row["prediction_cosine"][model]) for row in rows]
                )
                for model in (
                    "cone",
                    "training_common_ray",
                    "training_best_single",
                )
            },
        }

    require_keys(
        config,
        {
            "schema_version",
            "generated_on",
            "benchmark",
            "claim_ceiling",
            "source",
            "inputs",
            "author_provenance",
            "target",
            "h5mu",
            "analysis",
            "expected",
            "output",
        },
        "guide config",
    )
    require_keys(
        report,
        {
            "schema_version",
            "generated_on",
            "status",
            "benchmark",
            "claim_ceiling",
            "source",
            "config_sha256",
            "input_verification",
            "source_object_version",
            "author_provenance",
            "data_quality",
            "protocol",
            "solver_representation_sensitivity",
            "challenges",
            "summary",
            "limitations",
        },
        "guide report",
    )
    if (
        config["schema_version"] != "1.0.0"
        or report["schema_version"] != "1.0.0"
        or config["generated_on"] != "2026-07-19"
        or report["generated_on"] != config["generated_on"]
        or report["status"] != "PASS"
        or config["benchmark"] != benchmark
        or report["benchmark"] != benchmark
        or config["claim_ceiling"] != claim_ceiling
        or report["claim_ceiling"] != claim_ceiling
        or config["output"] != "results/guide_pair_transfer.json"
    ):
        raise AssertionError("guide benchmark version/status/claim identity differs")
    if report["config_sha256"] != sha256(config_path):
        raise AssertionError("guide report/config byte identity differs")
    if report["config_sha256"] != (
        "01d69a3090d5153d82c426e9b16cf728e9c222baa877d62a18bc0735b7ea85f4"
    ):
        raise AssertionError("guide frozen config identity differs")

    require_keys(
        ledger,
        {
            "benchmark",
            "report_schema_version",
            "report_sha256",
            "tier",
            "benchmark_role",
            "guide_identity_status",
            "official_modality_mapping",
            "h5mu_identity_field_status",
            "identity_crosswalk_status",
            "positional_modalities",
            "common_rest_atoms",
            "guide_1_only_rest_atoms",
            "shared_target_genes",
            "categorical_missing_rows",
            "missing_key_condition_suffix_counts",
            "fit_count",
            "representation_sensitivity_fits",
            "challenge_rows",
            "atom_order_sensitivity",
            "same_source_guide_position_transfer",
            "joint_guide_position_plus_target_source_transfer",
            "claim_boundary",
            "interpretation",
            "source",
        },
        "guide findings ledger",
    )
    if (
        ledger["benchmark"] != benchmark
        or ledger["report_schema_version"] != report["schema_version"]
        or hash_value(ledger["report_sha256"], "guide findings report hash")
        != sha256(report_path)
        or ledger["tier"] != "STRESS"
        or ledger["benchmark_role"]
        != "negative released first/second guide-rank reciprocal-transfer stress"
        or ledger["guide_identity_status"] != identity_status
        or ledger["claim_boundary"] != claim_ceiling
        or ledger["source"] != "results/guide_pair_transfer.json"
    ):
        raise AssertionError("guide findings identity/claim grain differs")

    expected_source = {
        "citation": (
            "Zhu R., Dann E. et al. (2025), Genome-scale perturb-seq in primary "
            "human CD4+ T cells maps context-specific regulators of T cell programs "
            "and human immune traits"
        ),
        "doi": "10.64898/2025.12.23.696273",
        "dataset_card": (
            "https://virtualcellmodels.cziscience.com/dataset/"
            "genome-scale-tcell-perturb-seq?access_dataset=true"
        ),
        "dataset_card_version": "v1.0",
        "dataset_card_accessed": "2026-07-19",
        "accessions": ["SRP643211", "GSE314342"],
        "data_license": (
            "Virtual Cells Platform metadata labels the dataset MIT plus its "
            "Acceptable Use Policy; the downloaded H5MU has no embedded license field"
        ),
        "license_boundary": (
            "Use is subject to the Virtual Cells Platform dataset metadata and "
            "Acceptable Use Policy; the author code repository is MIT-licensed."
        ),
        "redistribution": (
            "The 29.4 GB source object remains gitignored and is not redistributed."
        ),
    }
    if config["source"] != expected_source or report["source"] != expected_source:
        raise AssertionError("guide source/citation/license boundary differs")

    guide_input = {
        "path": "data/GWCD4i.DE_stats.by_guide.h5mu",
        "bytes": 29424424894,
        "sha256": "964eeafb3356a7322a1d5b1121802c6a1433456f3591e2d5797817df3bf9c2f6",
        "url": (
            "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
            "marson2025_data/GWCD4i.DE_stats.by_guide.h5mu"
        ),
        "s3_bucket": "genome-scale-tcell-perturb-seq",
        "s3_key": "marson2025_data/GWCD4i.DE_stats.by_guide.h5mu",
        "s3_version_id": "SQHf_ZhmdbCteM9f4HG1zs.2k25CL3gb",
        "etag_header": '"2e6705636ebaa276c7bc7c5a148ad096-3508"',
        "last_modified_header": "Thu, 28 May 2026 23:20:10 GMT",
        "last_modified_iso": "2026-05-28T23:20:10Z",
    }
    target_input = {
        "path": "data/Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv",
        "bytes": 6155771,
        "sha256": "c47d2df21414ca85e7aa255f4148904eec700fbcd9debc2f734ec97049698444",
        "url": (
            "https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/"
            "848d62fc2b7027f7218d6fc5f5b0c37255dc94af/metadata/suppl_tables/"
            "Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv"
        ),
    }
    if config["inputs"] != {
        "guide_h5mu": guide_input,
        "target_table": target_input,
    }:
        raise AssertionError("guide registered input contract differs")
    expected_input_verification = {
        name: {
            "path": spec["path"],
            "bytes": spec["bytes"],
            "sha256_expected": spec["sha256"],
            "sha256_actual": spec["sha256"],
            "hash_verified": True,
        }
        for name, spec in config["inputs"].items()
    }
    if report["input_verification"] != expected_input_verification:
        raise AssertionError("guide inputs were not exactly full-file-hash verified")
    object_fields = (
        "url",
        "s3_bucket",
        "s3_key",
        "s3_version_id",
        "etag_header",
        "last_modified_header",
        "last_modified_iso",
    )
    expected_object_version = {field: guide_input[field] for field in object_fields}
    if report["source_object_version"] != expected_object_version:
        raise AssertionError("guide source object/version identity differs")

    expected_author = {
        "repository": "https://github.com/emdann/GWT_perturbseq_analysis_2025",
        "commit": "848d62fc2b7027f7218d6fc5f5b0c37255dc94af",
        "initial_sample_selection_rule": (
            "keep_min_cells & keep_effective_guides & keep_total_counts"
        ),
        "selection_field": "keep_effective_guides",
        "minimum_passing_replicates_per_guide_condition": 3,
        "minimum_cells_per_guide_per_condition_sample": 5,
        "target_condition_membership_source": (
            "cond_targets from unreleased for_DE_by_guide.csv"
        ),
        "required_testable_guides_per_target_condition": 2,
        "final_sample_selection_rule": "initial keep_for_DE & keep_test_genes",
        "de_formula": "~ log10_n_cells + target",
        "donor_term_in_de_formula": False,
        "unreleased_intermediate": "for_DE_by_guide.csv",
        "official_modality_mapping": (
            "guide_1 is the lowest alphanumeric sgRNA ID and guide_2 the second "
            "within each perturbed-gene/culture-condition pair; targets with one "
            "passing guide appear only in guide_1"
        ),
        "h5mu_identity_field_status": (
            "no guide ID or sequence field is embedded in the guide-level H5MU"
        ),
        "public_identity_sources": [
            "GWCD4i.pseudobulk_merged.h5ad obs/guide_id",
            "sgrna_library_metadata.suppl_table.csv sgRNA",
        ],
        "identity_crosswalk_status": (
            "exact ranked-modality-to-sgRNA IDs not reconstructed or hash-verified "
            "in this benchmark"
        ),
        "guide_identity_status": identity_status,
        "interpretation": (
            "The official dataset card defines guide_1/guide_2 as the first/second "
            "alphanumeric sgRNA IDs within each target-condition pair, and public "
            "pseudobulk/library artifacts expose guide IDs. The H5MU does not embed "
            "those IDs, this benchmark does not cross-verify the exact mapping, and "
            "effectiveness-selected inclusion is not a structural-QC-only universe."
        ),
        "frozen_source_files": {
            "prep_DE_by_guide.ipynb": (
                "27cba79fdb96a6cecfa702759cd51a7b127d8faef9b08428d0c9bdd9c066c8ed"
            ),
            "run_guide_DE_chunk.py": (
                "fb068b7ef2d4413d55716cad3b21ce68c3f5a4819b712b6c8de63eb9258e58c7"
            ),
            "DE_config_by_guide.yaml": (
                "12b57d759d3a62c807aaa9dc646b15ce1aa5db939c315b3fc0602ba848801da6"
            ),
            "data_sharing_readme.md": (
                "1651f775007e3c8b544d0c4f069f52c25d2527bdde6d1f7feac9e88d37076ef8"
            ),
        },
    }
    if (
        config["author_provenance"] != expected_author
        or report["author_provenance"] != expected_author
    ):
        raise AssertionError("guide author-selection provenance differs")

    expected_counts = {
        "modalities": 2,
        "genes": 10282,
        "modality_rows": {"guide_1": 33488, "guide_2": 26078},
        "condition_counts": {
            "guide_1": {"Rest": 11075, "Stim8hr": 11231, "Stim48hr": 11068},
            "guide_2": {"Rest": 8323, "Stim8hr": 8998, "Stim48hr": 8719},
        },
        "categorical_missing_rows": {"guide_1": 114, "guide_2": 38},
        "missing_key_condition_suffix_counts": {
            "guide_1": {"Rest": 35, "Stim8hr": 18, "Stim48hr": 61},
            "guide_2": {"Rest": 0, "Stim8hr": 0, "Stim48hr": 38},
        },
        "rest_common_atoms": 8323,
        "guide_1_only_rest_atoms": 2752,
        "guide_1_rest_keys_sha256": (
            "06b97f432317c166e96d2abae7405f57d0106b3300d5e8cda5ec17c1e2169ada"
        ),
        "guide_2_rest_keys_sha256": (
            "a75887d0a38d0a3a800fbf9c28eae574e68a863dd549ef5b614a19621a808fed"
        ),
        "ordered_gene_ids_sha256": (
            "42a8a03a5cc9249c21da42fad070108398fba6977bb31d55c1de1be8c2b38faa"
        ),
        "ordered_gene_symbols_sha256": (
            "20b005b895d4ff83500c27b4db95f190cf16523664f902ffbf10aa44bb89212b"
        ),
        "shared_target_genes": 8950,
        "between_source_target_cosine": 0.7905978929805006,
        "fits": 12,
        "representation_sensitivity_fits": 1,
        "challenge_rows": 24,
    }
    if config["expected"] != expected_counts:
        raise AssertionError("guide expected data/design contract differs")
    for modality in modalities:
        counts = expected_counts["condition_counts"][modality]
        if (
            sum(counts.values())
            + expected_counts["categorical_missing_rows"][modality]
            != expected_counts["modality_rows"][modality]
        ):
            raise AssertionError("guide missing-category attrition arithmetic differs")
        if (
            sum(
                expected_counts["missing_key_condition_suffix_counts"][
                    modality
                ].values()
            )
            != expected_counts["categorical_missing_rows"][modality]
        ):
            raise AssertionError("guide missing-key suffix arithmetic differs")

    expected_quality_keys = {
        "h5mu_encoding",
        "modalities",
        "genes",
        "shared_source_safe_target_genes",
        "between_source_target_cosine",
        "guide_1_rest_atoms",
        "common_rest_atoms",
        "guide_1_only_rest_atoms",
        "guide_2_keys_exact_subset_of_guide_1",
        "selected_rest_target_mapping_exact",
        "metadata_hashes",
        "selected_matrix_sha256",
        "guide_identity_status",
        "eligibility_warning",
    }
    quality = require_keys(report["data_quality"], expected_quality_keys, "guide quality")
    expected_modalities = {
        modality: {
            "rows": expected_counts["modality_rows"][modality],
            "condition_counts": expected_counts["condition_counts"][modality],
            "categorical_missing_rows": expected_counts[
                "categorical_missing_rows"
            ][modality],
            "missing_key_condition_suffix_counts": expected_counts[
                "missing_key_condition_suffix_counts"
            ][modality],
        }
        for modality in modalities
    }
    expected_metadata_hashes = {
        "guide_1_rest_keys_sha256": expected_counts["guide_1_rest_keys_sha256"],
        "guide_2_common_rest_keys_sha256": expected_counts[
            "guide_2_rest_keys_sha256"
        ],
        "guide_1_only_rest_keys_sha256": (
            "6284d0e4cbe5a2d266b00589e03c5e12a3dee149632847d7ecb2dfb2d8ed84ab"
        ),
        "ordered_gene_ids_sha256": expected_counts["ordered_gene_ids_sha256"],
        "ordered_gene_symbols_sha256": expected_counts[
            "ordered_gene_symbols_sha256"
        ],
    }
    expected_matrix_hashes = {
        "guide_1": "e61544e2e260c50b916948496509806de7de51e83c2f83d6acc9e5e97042cf0e",
        "guide_2": "cb4c06a06b39e9498bfa680f90e2b39585aa208400c39b8972a051370e8b3844",
    }
    if (
        quality["h5mu_encoding"] != "MuData with two AnnData modalities"
        or quality["modalities"] != expected_modalities
        or quality["genes"] != expected_counts["genes"]
        or quality["shared_source_safe_target_genes"]
        != expected_counts["shared_target_genes"]
        or quality["between_source_target_cosine"]
        != expected_counts["between_source_target_cosine"]
        or quality["guide_1_rest_atoms"] != 11075
        or quality["common_rest_atoms"] != 8323
        or quality["guide_1_only_rest_atoms"] != 2752
        or quality["guide_2_keys_exact_subset_of_guide_1"] is not True
        or quality["selected_rest_target_mapping_exact"] is not True
        or quality["metadata_hashes"] != expected_metadata_hashes
        or quality["selected_matrix_sha256"] != expected_matrix_hashes
        or quality["guide_identity_status"] != identity_status
        or quality["eligibility_warning"]
        != (
            "The author pipeline uses keep_effective_guides, and the intermediate "
            "for_DE_by_guide.csv is absent. Presence in the released positional "
            "modalities can therefore select on perturbation effectiveness."
        )
    ):
        raise AssertionError("guide released-object quality/selection boundary differs")
    if quality["guide_1_rest_atoms"] != (
        quality["common_rest_atoms"] + quality["guide_1_only_rest_atoms"]
    ):
        raise AssertionError("guide Rest attrition arithmetic differs")

    expected_analysis = {
        "gene_order": "lexicographic_gene_symbol",
        "rng": "numpy.default_rng(seed)",
        "split_seeds": [0, 1, 2],
        "train_guide_slots": ["guide_1", "guide_2"],
        "train_target_sources": ["hollbacker", "ota"],
        "representation_sensitivity": {
            "train_guide_slot": "guide_1",
            "train_target_source": "hollbacker",
            "seed": 0,
            "alternative_atom_order": "reverse_lexicographic",
        },
        "score_target_sources": (
            "same source and the opposite source for every frozen fit"
        ),
        "baselines": ["zero", "training_common_ray", "training_best_single"],
        "inference": (
            "descriptive correlated challenges only; no p-values, confidence "
            "intervals, physical-guide generalization, or donor-population inference"
        ),
    }
    if config["analysis"] != expected_analysis:
        raise AssertionError("guide exhaustive analysis contract differs")
    protocol = report["protocol"]
    expected_protocol = {
        "condition": "Rest",
        "layer": "log_fc",
        "target_gene_universe": (
            "both target sources and guide-H5MU gene symbols; no held-out-source "
            "sign or magnitude filter"
        ),
        "gene_order": "lexicographic_gene_symbol",
        "split_rng": "numpy.default_rng(seed)",
        "split_seeds": [0, 1, 2],
        "fit_count": 12,
        "representation_sensitivity_fits": 1,
        "challenge_rows": 24,
        "opaque_key_contract": (
            "target_condition values are matched by exact string identity only; "
            "they are never parsed to reconstruct a guide identity. Cross-modality "
            "target ID/name equality is required only on the selected common Rest "
            "rows; no unsupported global cross-condition bijection is asserted"
        ),
        "coefficient_application": (
            "fit once on one positional modality, one target source, and fit genes; "
            "apply identical coefficients to both positional modalities on disjoint "
            "score genes"
        ),
        "baseline_application": (
            "select common-ray scale and best-single opaque key/scale only on the "
            "training modality, training source, and fit genes; freeze them for both "
            "score domains and both target-source scores"
        ),
        "whole_universe_reporting": (
            "all 8,323 released common Rest atoms are retained; no outcome-ranked "
            "top-k subset is selected"
        ),
        "prediction_cosine_zero_policy": (
            "two zero predictions map to 1 and exactly one zero prediction maps to 0"
        ),
        "inference": expected_analysis["inference"],
    }
    if protocol != expected_protocol:
        raise AssertionError("guide frozen-fit/inference protocol differs")

    expected_split_hashes = {
        0: (
            "a7e9b9765b708dd7c8fba2fe2bcfb44fc93b20b368defa5bf9c64f1df098f34b",
            "d281720b3185ca0fba2f07444473b90da33ba1fe93ea1144764f0f1782272a28",
        ),
        1: (
            "8c660049916a4c562662e5c4217cebb4c88c1bc19f5714b344df2c3d31544830",
            "ffcb4ba852727cedd068a9322a1b2a8ab88bc53e7e247ebfa606a01dc588c81b",
        ),
        2: (
            "a0673b2bddd34a50e19ef2f108aa372dfcf75e3d7d26803c2900e7a1c29ed8f5",
            "54bbb358363a31679aba03ab9a7bc00316eba5c0ab8d0688cbee32dcd1ee5590",
        ),
    }
    expected_design = []
    for train_slot in modalities:
        held_slot = "guide_2" if train_slot == "guide_1" else "guide_1"
        for train_source in sources:
            opposite = next(source for source in sources if source != train_source)
            for seed in seeds:
                for score_source in (train_source, opposite):
                    expected_design.append(
                        (
                            f"{train_slot}|{train_source}|seed={seed}",
                            train_slot,
                            held_slot,
                            train_source,
                            score_source,
                            (
                                "same_source_guide_position_transfer"
                                if score_source == train_source
                                else "joint_guide_position_plus_target_source_transfer"
                            ),
                            seed,
                        )
                    )

    challenges = report["challenges"]
    if not isinstance(challenges, list) or len(challenges) != 24:
        raise AssertionError("guide challenge row count differs")
    observed_design = [
        (
            row.get("fit_id"),
            row.get("train_guide_slot"),
            row.get("held_guide_slot"),
            row.get("train_target_source"),
            row.get("score_target_source"),
            row.get("source_transfer_class"),
            row.get("seed"),
        )
        for row in challenges
    ]
    if observed_design != expected_design:
        raise AssertionError("guide challenges are not the exact Cartesian design")

    row_keys = {
        "fit_id",
        "train_guide_slot",
        "held_guide_slot",
        "guide_identity_status",
        "train_target_source",
        "score_target_source",
        "source_transfer_class",
        "seed",
        "fit_genes",
        "score_genes",
        "fit_gene_sha256",
        "score_gene_sha256",
        "training",
        "within_training_guide",
        "reciprocal_held_guide",
        "reciprocal_minus_within",
        "prediction_cosine",
    }
    training_keys = {
        "cone_fit_cosine",
        "cone_kkt_violation",
        "cone_support",
        "coefficient_sha256",
        "baseline_sha256",
        "common_alpha",
        "best_single_opaque_key",
        "best_single_alpha",
        "atom_key_order_sha256",
        "training_model_sha256",
    }
    rows_by_fit: dict[str, list[dict]] = defaultdict(list)
    for index, row in enumerate(challenges):
        label = f"guide challenge {index}"
        require_keys(row, row_keys, label)
        if row["guide_identity_status"] != identity_status:
            raise AssertionError(f"{label} overstates guide identity")
        if row["fit_genes"] != 4475 or row["score_genes"] != 4475:
            raise AssertionError(f"{label} gene split sizes differ")
        expected_fit_hash, expected_score_hash = expected_split_hashes[row["seed"]]
        if (
            row["fit_gene_sha256"] != expected_fit_hash
            or row["score_gene_sha256"] != expected_score_hash
        ):
            raise AssertionError(f"{label} frozen gene identities differ")

        training = require_keys(row["training"], training_keys, f"{label} training")
        hashes = require_keys(
            training["baseline_sha256"], set(baselines), f"{label} baseline hashes"
        )
        coefficient_hash = hash_value(
            training["coefficient_sha256"], f"{label} coefficient hash"
        )
        for name, value in hashes.items():
            hash_value(value, f"{label} {name} hash")
        common_alpha = finite(training["common_alpha"], f"{label} common alpha")
        best_alpha = finite(training["best_single_alpha"], f"{label} best alpha")
        if common_alpha < 0 or best_alpha < 0:
            raise AssertionError(f"{label} training scalar is negative")
        opaque_key = training["best_single_opaque_key"]
        if not isinstance(opaque_key, str) or not opaque_key:
            raise AssertionError(f"{label} best-single opaque key is invalid")
        expected_baseline_hashes = {
            "training_common_ray": scalar_vector_hash(common_alpha),
            "training_best_single": canonical_hash(
                {"opaque_key": opaque_key, "alpha_hex": best_alpha.hex()}
            ),
            "zero": scalar_vector_hash(0.0),
        }
        if hashes != expected_baseline_hashes:
            raise AssertionError(f"{label} canonical baseline hashes differ")
        atom_hash = hash_value(
            training["atom_key_order_sha256"], f"{label} atom-order hash"
        )
        if atom_hash != expected_metadata_hashes["guide_2_common_rest_keys_sha256"]:
            raise AssertionError(f"{label} atom order differs from common Rest universe")
        expected_model_hash = canonical_hash(
            {
                "coefficient_sha256": coefficient_hash,
                "baseline_sha256": expected_baseline_hashes,
                "atom_key_order_sha256": atom_hash,
            }
        )
        if training["training_model_sha256"] != expected_model_hash:
            raise AssertionError(f"{label} canonical training-model hash differs")
        fit_cosine = finite(training["cone_fit_cosine"], f"{label} fit cosine")
        kkt = finite(training["cone_kkt_violation"], f"{label} KKT violation")
        support = training["cone_support"]
        if (
            not -1 <= fit_cosine <= 1
            or not 0 <= kkt <= 1e-8
            or type(support) is not int
            or not 0 < support <= 8323
        ):
            raise AssertionError(f"{label} training diagnostics are invalid")

        domain_values = {}
        for domain, expected_slot in (
            ("within_training_guide", row["train_guide_slot"]),
            ("reciprocal_held_guide", row["held_guide_slot"]),
        ):
            value = require_keys(
                row[domain], {"guide_slot", "metrics", "comparisons"}, f"{label} {domain}"
            )
            if value["guide_slot"] != expected_slot:
                raise AssertionError(f"{label} {domain} guide slot differs")
            observed_metrics = require_keys(
                value["metrics"], set(models), f"{label} {domain} metrics"
            )
            for model in models:
                model_metrics = require_keys(
                    observed_metrics[model], set(metrics), f"{label} {domain} {model}"
                )
                cosine = finite(
                    model_metrics["cosine"], f"{label} {domain} {model} cosine"
                )
                normalized_rmse = finite(
                    model_metrics["normalized_rmse"],
                    f"{label} {domain} {model} normalized RMSE",
                )
                norm_ratio = finite(
                    model_metrics["norm_ratio"],
                    f"{label} {domain} {model} norm ratio",
                )
                sign_agreement = finite(
                    model_metrics["sign_agreement"],
                    f"{label} {domain} {model} sign agreement",
                )
                if (
                    not -1 <= cosine <= 1
                    or normalized_rmse < 0
                    or norm_ratio < 0
                    or not 0 <= sign_agreement <= 1
                ):
                    raise AssertionError(f"{label} {domain} metric range differs")
            if observed_metrics["zero"] != {
                "cosine": 0.0,
                "normalized_rmse": 1.0,
                "norm_ratio": 0.0,
                "sign_agreement": 0.0,
            }:
                raise AssertionError(f"{label} zero baseline metrics differ")
            if common_alpha == 0 and observed_metrics["training_common_ray"] != (
                observed_metrics["zero"]
            ):
                raise AssertionError(f"{label} zero common-ray metrics differ")
            comparisons = require_keys(
                value["comparisons"], set(baselines), f"{label} {domain} comparisons"
            )
            for baseline in baselines:
                comparison = require_keys(
                    comparisons[baseline],
                    {"cosine_improvement", "normalized_rmse_improvement"},
                    f"{label} {domain} {baseline} comparison",
                )
                expected_cosine = (
                    observed_metrics["cone"]["cosine"]
                    - observed_metrics[baseline]["cosine"]
                )
                expected_rmse = (
                    observed_metrics[baseline]["normalized_rmse"]
                    - observed_metrics["cone"]["normalized_rmse"]
                )
                assert_close(
                    finite(
                        comparison["cosine_improvement"],
                        f"{label} {domain} cosine comparison",
                    ),
                    expected_cosine,
                )
                assert_close(
                    finite(
                        comparison["normalized_rmse_improvement"],
                        f"{label} {domain} RMSE comparison",
                    ),
                    expected_rmse,
                )
            domain_values[domain] = value

        deltas = require_keys(
            row["reciprocal_minus_within"],
            {"metrics", "comparisons"},
            f"{label} reciprocal-minus-within",
        )
        delta_metrics = require_keys(
            deltas["metrics"], set(models), f"{label} metric deltas"
        )
        for model in models:
            model_deltas = require_keys(
                delta_metrics[model], set(metrics), f"{label} {model} metric deltas"
            )
            for metric in metrics:
                expected_delta = (
                    domain_values["reciprocal_held_guide"]["metrics"][model][metric]
                    - domain_values["within_training_guide"]["metrics"][model][metric]
                )
                assert_close(
                    finite(model_deltas[metric], f"{label} {model} {metric} delta"),
                    expected_delta,
                )
        delta_comparisons = require_keys(
            deltas["comparisons"], set(baselines), f"{label} comparison deltas"
        )
        for baseline in baselines:
            comparison_deltas = require_keys(
                delta_comparisons[baseline],
                {"cosine_improvement", "normalized_rmse_improvement"},
                f"{label} {baseline} comparison deltas",
            )
            for metric in ("cosine_improvement", "normalized_rmse_improvement"):
                expected_delta = (
                    domain_values["reciprocal_held_guide"]["comparisons"][baseline][
                        metric
                    ]
                    - domain_values["within_training_guide"]["comparisons"][baseline][
                        metric
                    ]
                )
                assert_close(
                    finite(
                        comparison_deltas[metric],
                        f"{label} {baseline} {metric} delta",
                    ),
                    expected_delta,
                )
        prediction_cosines = require_keys(
            row["prediction_cosine"],
            {"cone", "training_common_ray", "training_best_single"},
            f"{label} prediction cosine",
        )
        for model, value in prediction_cosines.items():
            observed_cosine = finite(value, f"{label} {model} prediction cosine")
            if not -1 - 1e-12 <= observed_cosine <= 1 + 1e-12:
                raise AssertionError(f"{label} prediction cosine is out of range")
        if common_alpha == 0 and prediction_cosines["training_common_ray"] != 1.0:
            raise AssertionError(f"{label} prediction zero policy differs")
        rows_by_fit[row["fit_id"]].append(row)

    if len(rows_by_fit) != 12:
        raise AssertionError("guide unique fit count differs")
    model_hashes = set()
    coefficient_hashes = set()
    for fit_id, rows in rows_by_fit.items():
        if len(rows) != 2 or rows[0]["training"] != rows[1]["training"]:
            raise AssertionError(f"guide fit {fit_id} was not reused unchanged")
        for field in (
            "fit_genes",
            "score_genes",
            "fit_gene_sha256",
            "score_gene_sha256",
        ):
            if rows[0][field] != rows[1][field]:
                raise AssertionError(f"guide fit {fit_id} split identity was not reused")
        model_hashes.add(rows[0]["training"]["training_model_sha256"])
        coefficient_hashes.add(rows[0]["training"]["coefficient_sha256"])
    if len(model_hashes) != 12 or len(coefficient_hashes) != 12:
        raise AssertionError("guide frozen fit identities are not one-to-one")

    sensitivity = require_keys(
        report["solver_representation_sensitivity"],
        {
            "train_guide_slot",
            "held_guide_slot",
            "train_target_source",
            "seed",
            "fit_gene_sha256",
            "score_gene_sha256",
            "alternative_atom_order",
            "canonical_atom_order_sha256",
            "alternative_atom_order_sha256",
            "canonical_coefficient_sha256",
            "alternative_coefficient_mapped_to_canonical_order_sha256",
            "canonical_support",
            "alternative_support",
            "canonical_fit_cosine",
            "alternative_fit_cosine",
            "canonical_kkt_violation",
            "alternative_kkt_violation",
            "common_ray_alpha_absolute_difference",
            "best_single_key_match",
            "best_single_alpha_absolute_difference",
            "stability",
            "interpretation",
        },
        "guide solver representation sensitivity",
    )
    sensitivity_fit_id = "guide_1|hollbacker|seed=0"
    canonical_training = rows_by_fit[sensitivity_fit_id][0]["training"]
    if (
        sensitivity["train_guide_slot"] != "guide_1"
        or sensitivity["held_guide_slot"] != "guide_2"
        or sensitivity["train_target_source"] != "hollbacker"
        or sensitivity["seed"] != 0
        or sensitivity["fit_gene_sha256"] != expected_split_hashes[0][0]
        or sensitivity["score_gene_sha256"] != expected_split_hashes[0][1]
        or sensitivity["alternative_atom_order"] != "reverse_lexicographic"
        or sensitivity["canonical_atom_order_sha256"]
        != canonical_training["atom_key_order_sha256"]
        or sensitivity["canonical_coefficient_sha256"]
        != canonical_training["coefficient_sha256"]
        or sensitivity["canonical_support"] != canonical_training["cone_support"]
        or sensitivity["canonical_fit_cosine"]
        != canonical_training["cone_fit_cosine"]
        or sensitivity["canonical_kkt_violation"]
        != canonical_training["cone_kkt_violation"]
        or sensitivity["best_single_key_match"] is not True
        or sensitivity["interpretation"]
        != (
            "descriptive fixed-solver sensitivity only; one reversed atom order does "
            "not prove coefficient identifiability"
        )
    ):
        raise AssertionError("guide reverse-order sensitivity scope/binding differs")
    for field in (
        "alternative_atom_order_sha256",
        "alternative_coefficient_mapped_to_canonical_order_sha256",
    ):
        hash_value(sensitivity[field], f"guide sensitivity {field}")
    if sensitivity["alternative_atom_order_sha256"] == sensitivity[
        "canonical_atom_order_sha256"
    ]:
        raise AssertionError("guide reverse-order sensitivity did not change atom order")
    alternative_support = sensitivity["alternative_support"]
    alternative_fit = finite(
        sensitivity["alternative_fit_cosine"], "guide alternative fit cosine"
    )
    alternative_kkt = finite(
        sensitivity["alternative_kkt_violation"], "guide alternative KKT"
    )
    if (
        type(alternative_support) is not int
        or not 0 < alternative_support <= 8323
        or not -1 <= alternative_fit <= 1
        or not 0 <= alternative_kkt <= 1e-8
    ):
        raise AssertionError("guide alternative solver diagnostics are invalid")
    assert_close(alternative_fit, sensitivity["canonical_fit_cosine"])
    for field in (
        "common_ray_alpha_absolute_difference",
        "best_single_alpha_absolute_difference",
    ):
        value = finite(sensitivity[field], f"guide sensitivity {field}")
        if not 0 <= value <= 1e-12:
            raise AssertionError(f"guide sensitivity {field} is unstable")
    stability = require_keys(
        sensitivity["stability"],
        {
            "coefficients",
            "fit_prediction",
            "within_score_prediction",
            "held_score_prediction",
        },
        "guide solver stability",
    )
    for name, values in stability.items():
        require_keys(values, {"cosine", "relative_l2_difference"}, f"guide {name}")
        cosine = finite(values["cosine"], f"guide {name} cosine")
        difference = finite(
            values["relative_l2_difference"], f"guide {name} relative difference"
        )
        if not -1 - 1e-12 <= cosine <= 1 + 1e-12 or difference < 0:
            raise AssertionError(f"guide {name} solver sensitivity is invalid")

    summary = require_keys(
        report["summary"],
        {
            "same_source_guide_position_transfer",
            "joint_guide_position_plus_target_source_transfer",
        },
        "guide report summary",
    )
    for scope_name in summary:
        selected = [
            row for row in challenges if row["source_transfer_class"] == scope_name
        ]
        recomputed = summary_scope(selected)
        if summary[scope_name] != recomputed:
            raise AssertionError(f"guide {scope_name} summary is not reproducible")

    expected_ledger_direct = {
        "positional_modalities": 2,
        "common_rest_atoms": 8323,
        "guide_1_only_rest_atoms": 2752,
        "shared_target_genes": 8950,
        "categorical_missing_rows": {"guide_1": 114, "guide_2": 38},
        "fit_count": 12,
        "representation_sensitivity_fits": 1,
        "challenge_rows": 24,
    }
    for field, value in expected_ledger_direct.items():
        if ledger[field] != value:
            raise AssertionError(f"guide findings {field} differs")
    atom_ledger = require_keys(
        ledger["atom_order_sensitivity"],
        {
            "coefficient_cosine",
            "held_score_prediction_cosine",
            "maximum_relative_l2_difference",
        },
        "guide findings atom-order sensitivity",
    )
    assert_close(
        atom_ledger["coefficient_cosine"], stability["coefficients"]["cosine"]
    )
    assert_close(
        atom_ledger["held_score_prediction_cosine"],
        stability["held_score_prediction"]["cosine"],
    )
    assert_close(
        atom_ledger["maximum_relative_l2_difference"],
        max(value["relative_l2_difference"] for value in stability.values()),
    )

    scope_ledger_fields = {
        "challenge_rows",
        "unique_fits",
        "median_within_training_guide_cone_cosine",
        "median_reciprocal_held_guide_cone_cosine",
        "median_reciprocal_held_guide_training_best_single_cosine",
        "median_reciprocal_held_guide_cosine_improvement_over_training_best_single",
        "fraction_reciprocal_held_guide_cosine_improvement_positive_over_training_best_single",
        "median_reciprocal_held_guide_cone_normalized_rmse",
        "median_reciprocal_held_guide_training_best_single_normalized_rmse",
        "median_reciprocal_held_guide_normalized_rmse_improvement_over_training_best_single",
        "fraction_reciprocal_held_guide_normalized_rmse_improvement_positive_over_training_best_single",
        "median_reciprocal_minus_within_cone_cosine",
        "median_within_reciprocal_cone_prediction_cosine",
    }
    for scope_name in summary:
        observed_scope = summary[scope_name]
        ledger_scope = require_keys(
            ledger[scope_name], scope_ledger_fields, f"guide findings {scope_name}"
        )
        expected_compact = {
            "challenge_rows": observed_scope["challenge_rows"],
            "unique_fits": observed_scope["unique_fits"],
            "median_within_training_guide_cone_cosine": observed_scope[
                "within_training_guide"
            ]["metrics"]["cone"]["cosine"]["median"],
            "median_reciprocal_held_guide_cone_cosine": observed_scope[
                "reciprocal_held_guide"
            ]["metrics"]["cone"]["cosine"]["median"],
            "median_reciprocal_held_guide_training_best_single_cosine": observed_scope[
                "reciprocal_held_guide"
            ]["metrics"]["training_best_single"]["cosine"]["median"],
            "median_reciprocal_held_guide_cosine_improvement_over_training_best_single": observed_scope[
                "reciprocal_held_guide"
            ]["cone_improvement_over_baselines"]["training_best_single"][
                "cosine_improvement"
            ]["median"],
            "fraction_reciprocal_held_guide_cosine_improvement_positive_over_training_best_single": observed_scope[
                "reciprocal_held_guide"
            ]["cone_improvement_over_baselines"]["training_best_single"][
                "cosine_improvement"
            ]["fraction_positive"],
            "median_reciprocal_held_guide_cone_normalized_rmse": observed_scope[
                "reciprocal_held_guide"
            ]["metrics"]["cone"]["normalized_rmse"]["median"],
            "median_reciprocal_held_guide_training_best_single_normalized_rmse": observed_scope[
                "reciprocal_held_guide"
            ]["metrics"]["training_best_single"]["normalized_rmse"]["median"],
            "median_reciprocal_held_guide_normalized_rmse_improvement_over_training_best_single": observed_scope[
                "reciprocal_held_guide"
            ]["cone_improvement_over_baselines"]["training_best_single"][
                "normalized_rmse_improvement"
            ]["median"],
            "fraction_reciprocal_held_guide_normalized_rmse_improvement_positive_over_training_best_single": observed_scope[
                "reciprocal_held_guide"
            ]["cone_improvement_over_baselines"]["training_best_single"][
                "normalized_rmse_improvement"
            ]["fraction_positive"],
            "median_reciprocal_minus_within_cone_cosine": observed_scope[
                "reciprocal_minus_within"
            ]["metrics"]["cone"]["cosine"]["median"],
            "median_within_reciprocal_cone_prediction_cosine": observed_scope[
                "prediction_cosine"
            ]["cone"]["median"],
        }
        if ledger_scope != expected_compact:
            raise AssertionError(f"guide findings {scope_name} compact values differ")

    expected_limitations = [
        "guide_1 and guide_2 are positional modality labels; physical guide IDs are not present, so this is not physical guide-held-out generalization.",
        "guide_2 is a strict subset of guide_1 and the released universe is conditioned on the authors' keep_effective_guides field.",
        "The author pipeline also required keep_min_cells and keep_total_counts, retained two testable guide positions using at least three replicates and five cells per guide, and fit ~ log10_n_cells + target without a donor term.",
        "The object contains 114 guide_1 and 38 guide_2 rows with code -1 across every required categorical field; these are recorded as missing rather than imputed, and no selected common Rest row is missing its target mapping.",
        "The absent for_DE_by_guide.csv prevents reconstruction of a leakage-safe structural-QC-only universe from this object.",
        "Random gene splits contain correlated coordinates and are descriptive sensitivity challenges, not independent replicates.",
        "The same source study supplies both positional modalities; this is not independent external validation.",
        "Because the fit is underdetermined, transferred coefficients and their held-slot predictions are relative to the frozen SciPy NNLS solver and lexicographic opaque-key order; the fitted cone point need not identify a unique coefficient vector.",
        "No p-values or confidence intervals are emitted, and no donor, functional, state-conversion, or intervention claim is supported.",
    ]
    if report["limitations"] != expected_limitations:
        raise AssertionError("guide positional/effectiveness/solver limitations differ")
    expected_interpretation = (
        "The cone's within-position held-gene alignment does not survive reciprocal "
        "position transfer: same-source median cosine falls from 0.251 to -0.019 "
        "(median paired change -0.291), with positive gain over the training-selected "
        "best single in 3/12 rows and normalized-RMSE gain in 0/12. Joint "
        "target-source transfer is likewise negative. This falsifies robustness "
        "across the released positional summaries, but absent physical guide IDs and "
        "effectiveness-selected inclusion prevent a physical-guide generalization claim."
    )
    if ledger["interpretation"] != expected_interpretation:
        raise AssertionError("guide findings interpretation/claim ceiling differs")

    def inference_keys(value: object) -> set[str]:
        if isinstance(value, dict):
            return set(value) | set().union(
                *(inference_keys(item) for item in value.values()), set()
            )
        if isinstance(value, list):
            return set().union(*(inference_keys(item) for item in value), set())
        return set()

    forbidden = {
        "p_value",
        "p_values",
        "pvalue",
        "confidence_interval",
        "confidence_intervals",
    }
    if inference_keys({"challenges": challenges, "summary": summary}) & forbidden:
        raise AssertionError("guide report emitted undeclared inferential statistics")


def validate_values(findings: dict) -> None:
    if (
        findings.get("schema_version") != "3.2.0"
        or findings.get("generated_on") != "2026-07-19"
    ):
        raise AssertionError("unsupported findings schema")

    source = json.loads((RESULTS / "source_reconstruction.json").read_text())
    if source.get("status") != "PASS":
        raise AssertionError("source reconstruction did not pass")
    for status in source["data_quality"]["input_verification"].values():
        if not status["hash_verified"] or status["sha256_actual"] != status["sha256_expected"]:
            raise AssertionError("source input identity was not byte-hash verified")

    split_report = source["source_bound_splits"]
    if split_report["status"] != "PASS":
        raise AssertionError("source-bound split contract did not pass")
    ledger = findings["source_bound_alignment"]
    assert len(split_report["splits"]) == ledger["n_fixed_random_gene_splits"]
    observed_split_values = [row["held_out_cosine"] for row in split_report["splits"]]
    if observed_split_values != ledger["split_values"]:
        raise AssertionError("source-bound split values differ")
    assert_close(split_report["held_out_cosine_mean"], ledger["held_out_cosine_mean"])
    assert_close(split_report["held_out_cosine_sd"], ledger["held_out_cosine_sd"])
    assert_close(split_report["held_out_cosine_min"], ledger["held_out_cosine_min"])
    assert_close(split_report["held_out_cosine_max"], ledger["held_out_cosine_max"])
    if split_report["protocol"] != {
        "gene_order": ledger["gene_order"],
        "rng": ledger["split_rng"],
    }:
        raise AssertionError("source-bound split protocol differs")
    for split in split_report["splits"]:
        for field in ("fit_gene_sha256", "score_gene_sha256"):
            if len(split[field]) != 64:
                raise AssertionError("split gene identities are not hash-frozen")

    lineage = source["data_quality"]["target_lineage"]["zscore"]
    counts = lineage["counts"]
    ledger = findings["target_lineage"]
    assert counts["union"] == ledger["union_genes"]
    assert counts["shared"] == ledger["shared_source_genes"]
    assert counts["sign_concordant"] == ledger["sign_concordant_genes"]
    assert counts["shared_screen"] == ledger["shared_screen_genes"]
    assert counts["final"] == ledger["registered_genes"]
    assert_close(
        source["source_transfer"]["log_fc"]["between_source_screen_target_cosine"],
        ledger["between_source_log_fc_cosine_on_shared_screen"],
    )
    assert_close(
        source["source_transfer"]["zscore"]["between_source_screen_target_cosine"],
        ledger["between_source_zscore_cosine_on_shared_screen"],
    )

    transfer = findings["cross_source_directional_transfer"]
    directions = source["source_transfer"]["log_fc"]["directions"]
    for key in ("ota_to_hollbacker", "hollbacker_to_ota"):
        observed = directions[key]
        expected = transfer[key]
        assert len(observed["splits"]) == transfer["n_correlated_splits"]
        assert_close(_mean_metric(observed, "cone", "cosine"), expected["mean_cone_cosine"])
        assert_close(
            observed["cosine_improvement_mean"],
            expected["mean_cosine_improvement_over_test_selected_better_baseline"],
        )
        assert_close(
            _mean_metric(observed, "cone", "normalized_rmse"),
            expected["mean_cone_normalized_rmse"],
        )
        assert_close(
            _mean_metric(observed, "best_single_atom", "normalized_rmse"),
            expected["mean_best_single_normalized_rmse"],
        )

    arce = json.loads((EVIDENCE / "arce_external_validation_meta.json").read_text())
    arce_config_path = ROOT / "configs" / "arce_external_validation.json"
    arce_config = json.loads(arce_config_path.read_text())
    if arce["config_sha256"] != sha256(arce_config_path):
        raise AssertionError("Arce report/config identity differs")
    arce_input = arce["input_verification"]
    arce_members = arce_input["arce"]["members"]
    if (
        arce_input["arce"]["archive_sha256"]
        != arce_config["dataset"]["archive"]["sha256"]
        or arce_input["zhu_generator"]["predictor_sha256"]
        != arce_config["generator"]["predictor_sha256"]
    ):
        raise AssertionError("Arce input identities differ from the frozen contract")
    for key in (
        "screen_member", "activation_cells_member", "activation_summary_member"
    ):
        if arce_members[key]["sha256"] != arce_config["dataset"][key]["sha256"]:
            raise AssertionError(f"Arce {key} identity differs from the frozen contract")
    source_input = source["data_quality"]["input_verification"]["de_stats"]
    if (
        arce_input["zhu_generator"]["path"] != source_input["path"]
        or arce_input["zhu_generator"]["bytes"] != source_input["bytes"]
    ):
        raise AssertionError("Arce predictor is not cross-bound to the verified Zhu H5AD")
    if "unadjusted exploratory" not in arce["permutation_inference"]:
        raise AssertionError("Arce permutation multiplicity ceiling is missing")
    ledger = findings["arce_external_validation"]
    if arce["benchmark"] != ledger["benchmark"]:
        raise AssertionError("Arce benchmark identity differs")
    assert arce["attrition"]["analysis_eligible"] == ledger["analysis_targets"]
    for context, value in ledger["spearman"].items():
        assert_close(arce["contexts"][context]["point"]["spearman"], value)
    for context, value in ledger["top_25_overlap"].items():
        assert arce["contexts"][context]["point"]["top_k"]["25"]["overlap"] == value
    activation = ledger["activation_score_robustness"]
    data_quality = arce["activation_data_quality"]
    reproduction = arce["activation_summary_reproduction"]
    robustness = arce["activation_robustness"]
    assert data_quality["rows"] == activation["s14_cells"]
    assert data_quality["guide_groups"] == activation["guide_strata"]
    for field in ("min", "median", "max"):
        assert_close(
            data_quality[f"guide_group_cells_{field}"],
            arce_config["expected"][f"activation_guide_group_cells_{field}"],
        )
    assert reproduction["rows"] == activation["s8_aggregate_rows_reproduced"]
    if reproduction["published_p_values_used_for_inference"]:
        raise AssertionError("Arce S8 cell-level p-values entered current inference")
    assert_close(
        max(
            reproduction["maximum_absolute_mean_error"],
            reproduction["maximum_absolute_median_error"],
        ),
        activation["maximum_absolute_s8_error"],
    )
    if not robustness["all_strata_retained"] or robustness["cell_level_inference_emitted"]:
        raise AssertionError("Arce activation robustness violated its claim contract")
    assert (
        robustness["strata_below_20_cells"]
        == arce_config["expected"]["activation_strata_below_20_cells"]
    )
    for context, value in activation["donor_spearman"].items():
        assert_close(
            robustness["contexts"][context]["donor_rank_concordance"]["spearman"],
            value,
        )
    for context, value in activation[
        "all_two_guides_two_donors_same_sign_fraction"
    ].items():
        assert_close(
            robustness["contexts"][context][
                "all_two_guides_two_donors_same_nonzero_sign_fraction"
            ],
            value,
        )
    validate_activation_csv(activation, robustness)

    validate_zhu_arrayed(findings["zhu_arrayed_followup"])
    validate_donor_pair(findings["donor_pair_transfer"])
    validate_guide_pair(findings["guide_position_transfer"])
    validate_schmidt(findings)
    validate_goudy(findings)
    validate_library_coverage(findings)

    harness = json.loads((RESULTS / "validation_harness.json").read_text())
    ledger = findings["systemic_harness"]
    if harness["status"] != ledger["status"]:
        raise AssertionError("systemic harness did not pass")
    assert len(harness["scenarios"]) == ledger["scenarios"]
    max_t = harness["scenarios"]["maxT_null_calibration"]
    assert max_t["familywise_rejections"] == ledger["max_t_false_families"]
    assert max_t["trials"] == ledger["max_t_trials"]
    assert_close(max_t["familywise_error_upper_95"], ledger["max_t_exact_upper_95"])

    retired = split_report["retired_archived_comparison"]
    ledger = findings["retired_provenance"]
    assert_close(retired["mean"], ledger["archived_multisplit_mean"])
    assert_close(retired["sd"], ledger["archived_multisplit_sd"])
    assert_close(retired["fixed_split_cosine"], ledger["archived_fixed_split_cosine"])

    if findings.get("open_requirements") != [
        "leakage-safe pseudobulk donor- and guide-held-out evaluation",
        "module/pathway holdouts and calibrated structured nulls",
        "nested PCA, ridge, unconstrained and capacity-matched baselines",
        "same-experiment, guide-burden- and control-matched measured in-domain perturbation combinations",
        "paired CRISPRi and CRISPRa dictionaries",
        "independent whole-state protein, chromatin, function, fitness and durability",
        "prospective established-state validation",
    ]:
        raise AssertionError("open requirements changed without evidence")

    for section in findings.values():
        if not isinstance(section, dict):
            continue
        sources = section.get("sources", [section.get("source")])
        for relative in sources:
            if relative and not (ROOT / relative).is_file():
                raise FileNotFoundError(relative)


def validate_manifest() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    findings = json.loads((RESULTS / "findings.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.0.0":
        raise AssertionError("unsupported manifest schema")
    if manifest.get("generated_on") != findings.get("generated_on"):
        raise AssertionError("manifest/findings generation date differs")
    if sorted(manifest["files"]) != artifact_paths():
        raise AssertionError("manifest path set is stale")
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        if (
            path.stat().st_size != expected["bytes"]
            or sha256(path) != expected["sha256"]
            or bool(path.stat().st_mode & 0o111) != expected["executable"]
        ):
            raise AssertionError(f"artifact mismatch: {relative}")


def validate_local_markdown_links() -> None:
    """Fail when a maintained Markdown document points to a missing local path."""

    markdown_paths = [
        ROOT / relative
        for relative in artifact_paths()
        if relative.endswith(".md")
    ]
    failures: list[str] = []
    for document in markdown_paths:
        text = document.read_text(encoding="utf-8")
        for raw_target in MARKDOWN_LINK.findall(text):
            target = raw_target.strip().split(maxsplit=1)[0].strip("<>")
            if not target or target.startswith(("#", "http://", "https://", "mailto:")):
                continue
            local_part = unquote(target.split("#", 1)[0])
            if not local_part:
                continue
            resolved = (document.parent / local_part).resolve()
            if not resolved.exists():
                failures.append(
                    f"{document.relative_to(ROOT)} -> {target}"
                )
    if failures:
        raise AssertionError(
            "broken local Markdown links:\n  " + "\n  ".join(sorted(failures))
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    findings = json.loads((RESULTS / "findings.json").read_text(encoding="utf-8"))
    validate_values(findings)
    if args.write_manifest:
        write_manifest()
    validate_local_markdown_links()
    validate_manifest()
    print("findings and artifact lineage: OK")


if __name__ == "__main__":
    main()
