#!/usr/bin/env bash
# Build the static API documentation site with pdoc.
#
#   pip install -e '.[docs]'   # installs pdoc
#   bash build_docs.sh
#
# Output: site/api/*.html (module reference) alongside the hand-authored
# site/index.html landing page. Serve locally with:  python -m http.server -d site
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"

if ! python -c "import pdoc" 2>/dev/null; then
  echo "pdoc not installed. Run: pip install -e '.[docs]'" >&2
  exit 2
fi

rm -rf site/api
mkdir -p site/api site/figures
python -m pdoc -o site/api --docformat numpy --no-show-source \
  reachability combicone screen_ingest acquisition combicone_cli \
  effect_dictionary library_coverage neural_baseline

# keep the landing-page figure self-contained
cp -f docs/figures/fig_emergence_keystone.png site/figures/ 2>/dev/null || true

echo "API docs -> site/api/  (landing page: site/index.html)"
echo "Preview:  python -m http.server -d site 8000"
