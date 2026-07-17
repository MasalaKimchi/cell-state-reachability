#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

python -c "import numpy, scipy; print('numpy', numpy.__version__, '| scipy', scipy.__version__)"
python -m pytest -q
python reachability.py
python scripts/validate_findings.py

echo "Reproduction complete: strict geometry, frozen findings, and artifact lineage pass."
