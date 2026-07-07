#!/usr/bin/env bash
# OPTIONAL (Tier-2 only): fetch the gene-level effect matrix GWCD4i.DE_stats.h5ad.
#
# You do NOT need this to run the graded pipeline. Tier-1 (the no-auth CSVs already in
# data/) produces directional knockdown nominations on its own — see the ROADMAP
# "Tier-1 directional nominations" section and src/counterfactual.tier1_directional_nominations.
# This script is only for the full gene-space reachability solver.
#
# Reality check (from a real run): `vcp data search --exact` returns the 12 RAW
# per-donor/condition datasets (each millions of cells, hundreds of GB total). The
# derived `GWCD4i.DE_stats.h5ad` is NOT one of those hits, and current vcp-cli has no
# `--file` filter, so you cannot cherry-pick it by name. Options:
#
#   1) Check the dataset landing page for a direct derived-artifact download:
#        https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
#      and the bioRxiv "Data availability" section:
#        https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
#
#   2) If you locate a dataset ID that contains the derived DE artifact, download by ID
#      into ./data (this pulls the dataset's files; mind the size):
#        vcp login
#        vcp data download --id <DATASET_ID> -o ./data
#
#   3) The open, no-auth analysis repo (scripts + the supplementary tables we use):
#        https://github.com/emdann/GWT_perturbseq_analysis_2025
#
# After the file lands in ./data, verify:
#     python -m src.data_loader --check      # -> "DE_stats (h5ad, Tier 2) present: True"
set -euo pipefail
cd "$(dirname "$0")/.."

echo "Tier-2 is OPTIONAL. The graded pipeline runs on the no-auth Tier-1 CSVs."
echo "See the header of this script for how to locate GWCD4i.DE_stats.h5ad."
if command -v vcp >/dev/null 2>&1; then
  echo
  echo "Listing available datasets (note IDs; the derived DE artifact may be attached"
  echo "to the parent collection rather than a per-donor dataset):"
  vcp data search "Primary Human CD4+ T Cell Perturb-seq" --exact || true
fi
