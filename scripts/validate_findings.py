#!/usr/bin/env python3
"""Fail-closed validation of canonical findings and artifact lineage."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
EVIDENCE = RESULTS / "evidence"
MANIFEST = RESULTS / "manifest.json"

MANIFEST_PATHS = [
    ".gitignore",
    "LICENSE",
    "README.md",
    "requirements.txt",
    "requirements-external.txt",
    "reachability.py",
    "validation.py",
    "reproduce.sh",
    "configs/validation_harness.json",
    "configs/source_reconstruction.json",
    "configs/arce_external_validation.json",
    "data/README.md",
    "data/fetch_de_stats.sh",
    "tests/test_reachability.py",
    "tests/test_validation.py",
    "tests/test_source_reconstruction.py",
    "tests/test_arce_external_validation.py",
    "scripts/run_validation_harness.py",
    "scripts/run_source_reconstruction.py",
    "scripts/run_arce_external_validation.py",
    "scripts/validate_findings.py",
    "docs/FINDINGS.md",
    "docs/METHODS.md",
    "docs/SCIENTIFIC_VALIDATION_PLAN.md",
    "docs/VALIDATION_REPORT.md",
    "docs/figures/make_at_a_glance.py",
    "docs/figures/fig_at_a_glance.png",
    "docs/figures/fig_at_a_glance.pdf",
    "results/findings.json",
    "results/validation_harness.json",
    "results/source_reconstruction.json",
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
    payload = {
        "schema_version": "1.0.0",
        "generated_on": "2026-07-17",
        "files": files,
    }
    MANIFEST.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {MANIFEST.relative_to(ROOT)} ({len(files)} files)")


def assert_close(actual: float, expected: float, tolerance: float = 1e-12) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{actual} != {expected} within {tolerance}")


def _mean_metric(direction: dict, model: str, metric: str) -> float:
    return statistics.mean(
        split["metrics"][model][metric] for split in direction["splits"]
    )


def validate_values(findings: dict) -> None:
    if findings.get("schema_version") != "2.0.0":
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
    if (
        arce_input["arce"]["archive_sha256"]
        != arce_config["dataset"]["archive"]["sha256"]
        or arce_input["arce"]["member_sha256"]
        != arce_config["dataset"]["screen_member"]["sha256"]
        or arce_input["zhu_generator"]["predictor_sha256"]
        != arce_config["generator"]["predictor_sha256"]
    ):
        raise AssertionError("Arce input identities differ from the frozen contract")
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

    for section in findings.values():
        if not isinstance(section, dict):
            continue
        sources = section.get("sources", [section.get("source")])
        for relative in sources:
            if relative and not (ROOT / relative).is_file():
                raise FileNotFoundError(relative)


def validate_manifest() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.0.0":
        raise AssertionError("unsupported manifest schema")
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    findings = json.loads((RESULTS / "findings.json").read_text(encoding="utf-8"))
    validate_values(findings)
    if args.write_manifest:
        write_manifest()
    validate_manifest()
    print("findings and artifact lineage: OK")


if __name__ == "__main__":
    main()
