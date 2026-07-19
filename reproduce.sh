#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

python - <<'PY'
from importlib.metadata import version
from pathlib import Path

expected = dict(
    line.strip().split("==", 1)
    for line in Path("requirements.txt").read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
)
actual = {package: version(package) for package in expected}
if actual != expected:
    raise SystemExit(f"environment mismatch: expected {expected}, found {actual}")
print("verified environment:", ", ".join(f"{key}={value}" for key, value in actual.items()))
PY
python -m pytest -q
python reachability.py
python demo_library_coverage.py
python scripts/run_validation_harness.py --check results/validation_harness.json
python scripts/validate_findings.py

echo "Reproduction complete: geometry, coverage layer, systemic harness, frozen findings, and lineage pass."
