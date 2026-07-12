#!/usr/bin/env bash
# ==========================================================================
# cell-state-reachability — one-command reproduction
# ==========================================================================
# Verifies the method end-to-end from a clean checkout:
#   1. builds a minimal environment (numpy + scipy + pytest),
#   2. runs the test suite (packaged 38-assert self-test + 10 property tests),
#   3. reruns the synthetic invariant battery with visible output.
#
# This verifies the software. It does not regenerate the real-data headline,
# because the 16.8 GB Tier-2 matrix is intentionally not committed.
#
# Usage:
#   ./reproduce.sh                 # uses `python` on PATH (needs numpy/scipy)
#   ./reproduce.sh --venv          # create an isolated .venv and install pins
#   ./reproduce.sh --conda         # create the conda env from environment.yml
#
# Exit code 0 == everything reproduced. Any failure aborts with a nonzero code.
# --------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# make reachability.py importable without installing it as a package
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"

MODE="${1:-plain}"
PY="python"

banner() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

case "$MODE" in
  --venv)
    banner "1/3  Building isolated venv (.venv) with validated pins"
    python -m venv .venv
    # shellcheck disable=SC1091
    source .venv/bin/activate
    pip install --quiet --upgrade pip
    pip install --quiet "numpy==2.4.6" "scipy==1.17.1" pytest
    PY="python"
    ;;
  --conda)
    banner "1/3  Building conda env from environment.yml"
    conda env create -f environment.yml -q || conda env update -f environment.yml -q
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate cell-state-reachability
    PY="python"
    ;;
  plain)
    banner "1/3  Using python on PATH (expecting numpy + scipy present)"
    $PY -c "import numpy, scipy; print('numpy', numpy.__version__, '| scipy', scipy.__version__)"
    ;;
  *)
    echo "unknown option: $MODE  (use --venv, --conda, or no argument)"; exit 2;;
esac

banner "2/3  Running the test suite (self-test + property tests)"
$PY -m pytest tests/test_reachability.py -q

banner "3/3  Re-running the synthetic invariant battery"
$PY - <<'PYEOF'
import reachability as rx
# The packaged self-test already re-derives every invariant on synthetic
# fixtures; here we surface its diagnostics so a human sees the method run.
print("Running reachability._selftest() end-to-end ...\n")
rx._selftest()
print("\nSoftware verification complete: method imports, all invariants hold,")
print("verdicts / KKT certification / signed decomposition all reproduce.")
PYEOF

banner "DONE — reproduction succeeded"
echo "All tests passed and the synthetic method invariants reproduced."
