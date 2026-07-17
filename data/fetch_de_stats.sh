#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FILE="${1:-GWCD4i.DE_stats.h5ad}"
PREFIX="s3://genome-scale-tcell-perturb-seq/marson2025_data"
TARGET_URL="https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/848d62fc2b7027f7218d6fc5f5b0c37255dc94af/metadata/suppl_tables/Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv"
OUT="$ROOT/data/$FILE"
PART="$OUT.part"

case "$FILE" in
  GWCD4i.DE_stats.h5ad) EXPECTED_BYTES=16786240107 ;;
  GWCD4i.DE_stats.by_donors.h5mu) EXPECTED_BYTES=16866278447 ;;
  GWCD4i.DE_stats.by_guide.h5mu) EXPECTED_BYTES=29424424894 ;;
  GWCD4i.pseudobulk_merged.h5ad) EXPECTED_BYTES=44566657140 ;;
  Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv) EXPECTED_BYTES=6155771 ;;
  *)
    echo "Unsupported or unresolved object: $FILE" >&2
    exit 2
    ;;
esac

if [[ -e "$OUT" ]]; then
  echo "Refusing to overwrite existing file: $OUT" >&2
  shasum -a 256 "$OUT"
  exit 2
fi
rm -f "$PART"

if [[ "$FILE" == Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv ]]; then
  command -v curl >/dev/null || { echo "curl is required" >&2; exit 2; }
  curl --fail --location --output "$PART" "$TARGET_URL"
else
  command -v aws >/dev/null || {
    echo "AWS CLI is required: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html" >&2
    exit 2
  }
  aws s3 cp --no-sign-request "$PREFIX/$FILE" "$PART"
fi

ACTUAL_BYTES="$(wc -c < "$PART" | tr -d ' ')"
if [[ "$ACTUAL_BYTES" != "$EXPECTED_BYTES" ]]; then
  echo "Byte-length mismatch: expected $EXPECTED_BYTES, found $ACTUAL_BYTES" >&2
  exit 1
fi
mv "$PART" "$OUT"
shasum -a 256 "$OUT"
