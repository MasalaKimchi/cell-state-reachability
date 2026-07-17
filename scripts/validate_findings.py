#!/usr/bin/env python3
"""Validate the frozen updated findings and artifact lineage."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
EVIDENCE = RESULTS / "evidence"
MANIFEST = RESULTS / "manifest.json"

MANIFEST_PATHS = [
    ".gitignore",
    "LICENSE",
    "README.md",
    "ROADMAP.md",
    "requirements.txt",
    "reachability.py",
    "reproduce.sh",
    "data/README.md",
    "data/fetch_de_stats.sh",
    "tests/test_reachability.py",
    "scripts/validate_findings.py",
    "docs/FINDINGS.md",
    "docs/METHODS.md",
    "docs/VALIDATION_REPORT.md",
    "docs/figures/make_at_a_glance.py",
    "docs/figures/fig_at_a_glance.png",
    "docs/figures/fig_at_a_glance.pdf",
    "results/findings.json",
    "results/README.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_paths() -> list[str]:
    evidence = [str(path.relative_to(ROOT)) for path in sorted(EVIDENCE.iterdir()) if path.is_file()]
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
    MANIFEST.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {MANIFEST.relative_to(ROOT)} ({len(files)} files)")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def assert_close(actual: float, expected: float, tolerance: float = 1e-12) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"{actual} != {expected} within {tolerance}")


def validate_values(findings: dict) -> None:
    headline = findings["headline"]
    split = [float(row["heldout_cosine"]) for row in rows(EVIDENCE / "headline_heldout_split_stability.csv")]
    assert len(split) == headline["n_fixed_random_gene_splits"] == 12
    assert_close(statistics.mean(split), headline["held_out_cosine_mean"])
    assert_close(statistics.stdev(split), headline["held_out_cosine_sd"])
    assert_close(min(split), headline["held_out_cosine_min"])
    assert_close(max(split), headline["held_out_cosine_max"])
    if split != headline["split_values"]:
        raise AssertionError("split values differ from findings ledger")

    fixed_rows = rows(EVIDENCE / "historical_fixed_split_null.csv")
    if len(fixed_rows) != 1:
        raise AssertionError("historical fixed-split evidence must have one row")
    fixed = fixed_rows[0]
    assert_close(float(fixed["observed_heldout"]), headline["historical_fixed_split_cosine"])
    assert int(fixed["n_shuffles"]) == headline["diagnostic_target_shuffles"]
    below = int(fixed["n_shuffles"]) - int(fixed["exceedances"])
    assert below == headline["diagnostic_target_shuffles_below_observed"]
    assert_close(float(fixed["plus_one_p"]), headline["diagnostic_plus_one_p"])
    assert_close(
        (1 + int(fixed["exceedances"])) / (1 + int(fixed["n_shuffles"])),
        headline["diagnostic_plus_one_p"],
    )

    by_id = {entry["id"]: entry for entry in findings["updated_findings"]}
    baseline = {
        row["predictor"]: float(row["held_out_cosine"])
        for row in rows(EVIDENCE / "baseline_comparison.csv")
        if row["axis"] == "Th2->Th1_Rest(flagship)"
    }
    expected = by_id["baseline_position"]["values"]
    assert_close(baseline["cone_nnls"], expected["cone"])
    assert_close(baseline["B1a_predict_mean"], expected["predict_mean"])
    assert_close(baseline["B1b_pca_k1"], expected["pca_1"])
    assert_close(baseline["B1b_pca_k5"], expected["pca_5"])
    assert_close(baseline["B2_random_cone_mean"], expected["matched_random_cone_mean"])
    assert_close(baseline["B1b_pca_k20"], expected["pca_20"])
    assert_close(baseline["B3_unconstrained_ls"], expected["unconstrained_ls"])
    assert_close(baseline["B3_ridge_lam100"], expected["ridge_100"])

    metric = {row["weighting"]: row for row in rows(EVIDENCE / "metric_calibration.csv")}
    expected = by_id["metric_calibration"]["values"]
    assert_close(float(metric["uniform"]["dynamic_range_fraction"]), expected["uniform"])
    assert_close(float(metric["specificity"]["dynamic_range_fraction"]), expected["specificity"])

    context = {row["condition"]: row for row in rows(EVIDENCE / "context_condition_comparison.csv")}
    for condition, value in by_id["context"]["values"].items():
        assert_close(float(context[condition]["heldout_cosine"]), value)

    confound = json.loads((EVIDENCE / "confounder_robustness_summary.json").read_text())
    expected = by_id["specific_challenges"]["values"]
    cell_cycle_delta = (
        confound["cellcycle_ablation"]["ablated_heldout"]
        - confound["cellcycle_ablation"]["baseline_heldout"]
    )
    assert_close(cell_cycle_delta, expected["cell_cycle_delta"])
    assert_close(confound["magnitude_deconfound"]["norm_matched_null_z"], expected["norm_matched_null_z"])
    max_rescaling = max(row["abs_delta"] for row in confound["guide_efficiency"]["late_rescaling_invariance"])
    assert_close(max_rescaling, expected["max_rescaling_delta"])

    generator_summary = json.loads((EVIDENCE / "generator_significance_summary.json").read_text())
    generator_rows = {row["generator_set"]: row for row in rows(EVIDENCE / "generator_significance.csv")}
    expected = by_id["generator_filter"]["values"]
    assert generator_summary["n_generators_total"] == expected["n_generators"]
    assert_close(generator_summary["frac_generators_significant"], expected["fraction_source_significant"])
    assert_close(float(generator_rows["effect_ge_median"]["heldout_cosine"]), expected["top_half_held_out"])

    ranking = json.loads((EVIDENCE / "ranking_validation_summary.json").read_text())
    expected = by_id["ranking_scope"]["values"]
    assert ranking["n_panel"] == expected["n_panel"]
    assert_close(ranking["cone_directional_auroc"], expected["cone_directional_auroc"])
    assert_close(ranking["magnitude_directional_auroc"], expected["magnitude_directional_auroc"])

    combination = json.loads((EVIDENCE / "additivity_verdict_flip_summary.json").read_text())
    expected = by_id["combination_scope"]["values"]
    assert combination["n_doubles"] == expected["n"]
    assert combination["verdict_flip_count"] == expected["threshold_flips"]
    assert combination["modality_flip_count"] == expected["modality_flips"]
    assert_close(combination["median_cos_measured_vs_additive"], expected["median_measured_additive_cosine"])

    target_meta = json.loads((EVIDENCE / "reviewer2_ota_hollbacher_meta.json").read_text())
    expected = by_id["target_source_agreement"]["values"]
    assert_close(target_meta["between_study_cosine"], expected["between_source_cosine"])
    assert target_meta["genes_shared"] == expected["shared_genes"]
    assert_close(target_meta["pct_concordant"], expected["sign_concordance_percent"])

    coverage = {row["top_K"]: row for row in rows(EVIDENCE / "reviewer2_deg_survival.csv") if row["signature"].startswith("Th2/Th1")}
    source_scores = {row["target"]: row for row in rows(EVIDENCE / "reviewer2_ota_hollbacher_split.csv")}
    expected = by_id["target_observation_scope"]["values"]
    assert int(coverage["all(25672)"]["surviving"]) == expected["target_genes_in_screen"]
    assert int(coverage["50"]["surviving"]) == expected["top_50_surviving"]
    assert expected["target_genes_total"] == 25672
    assert_close(float(source_scores["Ota 2021 only"]["held_out_cosine"]), expected["ota_only_held_out"])
    assert_close(float(source_scores["Höllbacher 2021 only"]["held_out_cosine"]), expected["hollbacher_only_held_out"])
    assert_close(float(source_scores["merged (headline; sign-concordant)"]["held_out_cosine"]), expected["merged_held_out"])

    for entry in findings["updated_findings"]:
        sources = entry.get("sources", [entry.get("source")])
        for source in sources:
            if not source or not (ROOT / source).is_file():
                raise FileNotFoundError(source)


def validate_manifest() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.0.0":
        raise AssertionError("unsupported manifest schema")
    expected_paths = artifact_paths()
    if sorted(manifest["files"]) != expected_paths:
        raise AssertionError("manifest path set is stale")
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        executable = bool(path.stat().st_mode & 0o111)
        if (
            path.stat().st_size != expected["bytes"]
            or sha256(path) != expected["sha256"]
            or executable != expected["executable"]
        ):
            raise AssertionError(f"artifact mismatch: {relative}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    findings = json.loads((RESULTS / "findings.json").read_text(encoding="utf-8"))
    if findings.get("schema_version") != "1.0.0":
        raise AssertionError("unsupported findings schema")
    validate_values(findings)
    if args.write_manifest:
        write_manifest()
    validate_manifest()
    print("findings and artifact lineage: OK")


if __name__ == "__main__":
    main()
