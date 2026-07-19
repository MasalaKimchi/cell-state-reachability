#!/usr/bin/env bash
# One-command reproduction of the CombiCone manuscript's central result.
#
# Regenerates the headline emergence certificate for all measured Norman doubles
# from the frozen substrate, verifies it against the committed reference (geometry
# exact; verdict tiers stable), then runs the unit suite. See docs/ for the full
# evidence notes and scripts/reproduce_paper.py for the reproducibility contract.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

SUBSTRATE="${COMBICONE_SUBSTRATE:-combicone_substrate.npz}"
REFERENCE="${COMBICONE_CERT_REF:-docs/phase2/emergence_certificate.csv}"

if [[ ! -f "$SUBSTRATE" ]]; then
  echo "ERROR: substrate not found at $SUBSTRATE" >&2
  echo "  Set COMBICONE_SUBSTRATE=/path/to/combicone_substrate.npz" >&2
  exit 2
fi

echo "=== [1/3] Regenerate + verify the emergence certificate ==="
if [[ -f "$REFERENCE" ]]; then
  python scripts/reproduce_paper.py --substrate "$SUBSTRATE" \
    --check "$REFERENCE" --out results/emergence_certificate_reproduced.csv
else
  echo "(no committed reference at $REFERENCE — regenerating without --check)"
  python scripts/reproduce_paper.py --substrate "$SUBSTRATE" \
    --out results/emergence_certificate_reproduced.csv
fi

echo
echo "=== [2/3] Certificate trust dossier (negative controls + sensitivity) ==="
if [[ -f "$REFERENCE" ]]; then
  python scripts/certificate_dossier.py --substrate "$SUBSTRATE" \
    --certificate "$REFERENCE" --ledger docs/phase2/manuscript_numbers.json \
    --out results/certificate_dossier.json --n-neg 100 --n-boot 120
else
  echo "(skipped — needs the reference certificate CSV)"
fi

echo
echo "=== [3/3] Unit + property test suite ==="
python -m pytest -q

echo
echo "Reproduction complete. Regenerated certificate + dossier in results/;"
echo "the paper's numeric ledger is docs/phase2/manuscript_numbers.json."
