import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
import pytest

from reachability import InputError
from scripts.run_schmidt_external_validation import (
    SCREEN_IDS,
    SGRNA_HEADER,
    _canonical_json,
    eligible_genes,
    parse_screen,
    run,
    select_top_k,
)


def _screen_payload(
    screen_index: int,
    *,
    bad_header: bool = False,
    nonfinite: bool = False,
    duplicate: bool = False,
    bad_suffix: bool = False,
) -> bytes:
    handle = io.StringIO(newline="")
    header = list(SGRNA_HEADER)
    if bad_header:
        header[-1] = "changed"
    writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    rows = []
    for gene_index in range(8):
        gene = f"GENE{gene_index:02d}"
        for guide_index in range(2):
            guide = f"{gene}_GUIDE{guide_index}"
            for donor_index, donor in enumerate(("r0", "r1")):
                if bad_suffix and gene_index == 0 and guide_index == 0 and donor_index == 0:
                    donor = "r2"
                # Screen, donor, and guide contributions keep every metric non-constant.
                lfc = (
                    (gene_index - 3.5) * (1.0 + 0.05 * screen_index)
                    + donor_index * (0.2 + gene_index * 0.03)
                    + guide_index * 0.07
                    + (0.15 if (gene_index + screen_index) % 3 == 0 else -0.05)
                )
                if nonfinite and gene_index == 0 and guide_index == 0 and donor_index == 0:
                    lfc_value = "nan"
                else:
                    lfc_value = f"{lfc:.8f}"
                row = {
                    "sgrna": f"{guide}_{donor}",
                    "Gene": gene,
                    "control_count": "10",
                    "treatment_count": "12",
                    "control_mean": "10",
                    "treat_mean": "12",
                    "LFC": lfc_value,
                    "control_var": "1",
                    "adj_var": "1",
                    "score": "2",
                    "p.low": "0.1",
                    "p.high": "0.9",
                    "p.twosided": "0.2",
                    "FDR": "0.3",
                    "high_in_treatment": "True",
                }
                if bad_header:
                    row["changed"] = row.pop("high_in_treatment")
                rows.append(row)
    if duplicate:
        rows.append(rows[0].copy())
    writer.writerows(rows)
    return handle.getvalue().encode()


def _identity(payload: bytes) -> dict[str, object]:
    return {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def _make_fixture(directory: Path) -> Path:
    archive = directory / "screens.zip"
    script = (
        "# r0 = Donor1, r1 = Donor2\n"
        "Donor1_LFC = median(Donor1_LFC)\n"
        "Donor2_LFC = median(Donor2_LFC)\n"
        "mutate(rev_lfc = -1 * LFC\n"
    ).encode()
    payloads = {
        screen_id: _screen_payload(index)
        for index, screen_id in enumerate(SCREEN_IDS)
    }
    script_member = "fixture/Generate-screen-analysis-table.R"
    members = [script_member]
    member_names = {}
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(script_member, script)
        for screen_id, payload in payloads.items():
            member = f"fixture/{screen_id}.txt"
            members.append(member)
            member_names[screen_id] = member
            zf.writestr(member, payload)
    archive_bytes = archive.read_bytes()

    screens = {}
    for screen_id in SCREEN_IDS:
        modality, context = screen_id.split("_")
        payload = payloads[screen_id]
        screens[screen_id] = {
            "modality": modality,
            "context": context,
            "cell_system": "fixture cells",
            "member": member_names[screen_id],
            "rows": 32,
            **_identity(payload),
        }
    config = {
        "schema_version": 1,
        "report_version": 1,
        "benchmark_id": "fixture_schmidt",
        "source": {
            "same_donors_basis": "fixture uses the same two labels",
        },
        "input": {
            "path": str(archive),
            "url": "https://example.invalid/screens.zip",
            "bytes": len(archive_bytes),
            "sha256": hashlib.sha256(archive_bytes).hexdigest(),
            "md5": hashlib.md5(archive_bytes).hexdigest(),  # noqa: S324 - fixture identity
        },
        "archive_members": members,
        "author_script": {
            "member": script_member,
            **_identity(script),
            "required_text": [
                "r0 = Donor1, r1 = Donor2",
                "Donor1_LFC = median(Donor1_LFC)",
                "Donor2_LFC = median(Donor2_LFC)",
                "mutate(rev_lfc = -1 * LFC",
            ],
        },
        "screens": screens,
        "analysis": {
            "donor_suffixes": {"r0": "Donor1", "r1": "Donor2"},
            "orientation": {"CRISPRa": 1, "CRISPRi": -1},
            "excluded_gene": "NO-TARGET",
            "eligibility": "identity/coverage only",
            "aggregation": "median guide LFC",
            "minimum_guides_per_gene_per_donor": [1, 2],
            "primary_minimum_guides": 2,
            "top_k": [3, 4],
            "primary_top_k": 4,
            "expected_common_genes": {"1": 8, "2": 8},
            "sensitivity_status": "EXPLORATORY_POST_HOC_DESCRIPTIVE",
        },
        "claim_contract": {
            "tier": "STRESS",
            "claim_ceiling": "fixture ceiling",
        },
    }
    config_path = directory / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path


def test_fixture_run_has_complete_transfer_grid_and_training_only_hashes(tmp_path: Path) -> None:
    report = run(_make_fixture(tmp_path))
    assert report["execution_status"] == "PASS"
    assert report["eligibility"]["outcome_fields_used"] == []
    rows = report["source_selected_transfer"]["rows"]
    assert len(rows) == 2 * 2 * 4 * 2 * 3
    assert report["source_selected_transfer"]["primary_directions_per_class"] == 8
    assert {row["transfer_class"] for row in rows} == {
        "same_screen_held_donor",
        "donor_plus_modality_library_same_context",
        "donor_plus_cross_context_cytokine_cell_type_same_modality",
    }
    grouped = {}
    for row in rows:
        key = (
            row["minimum_guides"],
            row["top_k"],
            row["source_screen"],
            row["training_donor"],
        )
        grouped.setdefault(key, set()).add(row["source_selected_gene_sha256"])
        assert row["training_donor"] != row["held_donor"]
        assert 0 <= row["held_target_global_top_k_overlap_fraction"] <= 1
    assert all(len(hashes) == 1 for hashes in grouped.values())


def test_fixture_is_deterministic_and_check_mode_is_byte_exact(tmp_path: Path) -> None:
    config_path = _make_fixture(tmp_path)
    assert _canonical_json(run(config_path)) == _canonical_json(run(config_path))
    output = tmp_path / "report.json"
    script = Path(__file__).parents[1] / "scripts" / "run_schmidt_external_validation.py"
    subprocess.run(
        [sys.executable, str(script), "--config", str(config_path), "--output", str(output), "--write"],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [sys.executable, str(script), "--config", str(config_path), "--output", str(output), "--check"],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(output.read_text())
    payload["execution_status"] = "DRIFT"
    output.write_text(json.dumps(payload), encoding="utf-8")
    failed = subprocess.run(
        [sys.executable, str(script), "--config", str(config_path), "--output", str(output), "--check"],
        capture_output=True,
        text=True,
    )
    assert failed.returncode != 0


def test_archive_hash_gate_precedes_zip_parse(tmp_path: Path) -> None:
    config_path = _make_fixture(tmp_path)
    config = json.loads(config_path.read_text())
    config["input"]["sha256"] = "0" * 64
    config_path.write_text(json.dumps(config))
    with pytest.raises(InputError, match="SHA-256"):
        run(config_path)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"bad_header": True}, "schema"),
        ({"nonfinite": True}, "non-finite"),
        ({"duplicate": True}, "duplicate"),
        ({"bad_suffix": True}, "donor suffix"),
    ],
)
def test_parser_fails_closed_on_corruption(kwargs: dict, message: str) -> None:
    payload = _screen_payload(0, **kwargs)
    expected_rows = 33 if kwargs.get("duplicate") else 32
    with pytest.raises(InputError, match=message):
        parse_screen(payload, "CRISPRa_IFNG", expected_rows, {"r0": "D1", "r1": "D2"})


def test_eligibility_does_not_read_effect_values() -> None:
    parsed = parse_screen(_screen_payload(0), "CRISPRa_IFNG", 32, {"r0": "D1", "r1": "D2"})
    screens = {screen_id: parsed for screen_id in SCREEN_IDS}
    before = eligible_genes(screens, 2, ("r0", "r1"), "NO-TARGET")
    for screen in screens.values():
        for donor_effects in screen["effects"].values():
            donor_effects["r0"] = 1e99
            donor_effects["r1"] = -1e99
    after = eligible_genes(screens, 2, ("r0", "r1"), "NO-TARGET")
    assert before == after


def test_top_k_uses_only_source_absolute_effect_and_lexical_ties() -> None:
    genes = ["B", "A", "C", "D"]
    source = np.asarray([2.0, -2.0, 1.0, 0.5])
    assert select_top_k(genes, source, 2) == ["A", "B"]
    target_a = np.asarray([100.0, 0.0, -100.0, 50.0])
    target_b = -target_a
    assert select_top_k(genes, source, 2) == select_top_k(genes, source, 2)
    assert not np.array_equal(target_a, target_b)  # target values never enter source selection
