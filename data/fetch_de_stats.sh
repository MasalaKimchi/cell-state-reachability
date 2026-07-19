#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FILE="${1:-GWCD4i.DE_stats.h5ad}"
S3_PREFIX="https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
AUTHOR_PREFIX="https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/848d62fc2b7027f7218d6fc5f5b0c37255dc94af/metadata/suppl_tables"
VCP_PERTURB_PREFIX="https://cz-benchmarks-data.s3.amazonaws.com/datasets/v1/perturb/single_cell"
GEO_GSE306915_PREFIX="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306915"
GOUDY_AUTHOR_PREFIX="https://raw.githubusercontent.com/GilbertLabUCSF/T_Cell_CRISPRoff/53155c9207d8b4f70f0ae5d60e1f4c0513d41bd7/Cas9_KO_vs_CRISPRoff_KD"
SCHMIDT_ZENODO_PREFIX="https://zenodo.org/api/records/5784651/files"
OUT="$ROOT/data/$FILE"
PART="$OUT.part"

EXPECTED_SHA256=""
EXPECTED_MD5=""
EXPECTED_ETAG=""
case "$FILE" in
  GWCD4i.DE_stats.h5ad) EXPECTED_BYTES=16786240107; EXPECTED_SHA256=c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62; URL="$S3_PREFIX/$FILE" ;;
  GWCD4i.DE_stats.by_donors.h5mu) EXPECTED_BYTES=16866278447; EXPECTED_SHA256=2ee3cf90925600eb044619021da2bdd47d661f306a204586652256facf17af64; URL="$S3_PREFIX/$FILE" ;;
  GWCD4i.DE_stats.by_guide.h5mu) EXPECTED_BYTES=29424424894; EXPECTED_SHA256=964eeafb3356a7322a1d5b1121802c6a1433456f3591e2d5797817df3bf9c2f6; EXPECTED_ETAG='"2e6705636ebaa276c7bc7c5a148ad096-3508"'; URL="$S3_PREFIX/$FILE" ;;
  GWCD4i.pseudobulk_merged.h5ad) EXPECTED_BYTES=44566657140; URL="$S3_PREFIX/$FILE" ;;
  Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv) EXPECTED_BYTES=6155771; EXPECTED_SHA256=c47d2df21414ca85e7aa255f4148904eec700fbcd9debc2f734ec97049698444; URL="$AUTHOR_PREFIX/$FILE" ;;
  IL10IL21bulkRNAseq_DESeq2_results.csv) EXPECTED_BYTES=13952871; EXPECTED_SHA256=c20418a9285b10104dbae362b825971f86f97425800a92269e4433ce780e666d; URL="$AUTHOR_PREFIX/$FILE" ;;
  IL10_IL21_arrayed_validation.csv) EXPECTED_BYTES=2200; EXPECTED_SHA256=f60cdda392d6f29d10a539727ff7324b04d17e35c0512c889b733e00380b83dc; URL="$AUTHOR_PREFIX/$FILE" ;;
  norman_perturbation.h5ad) EXPECTED_BYTES=2228849977; EXPECTED_SHA256=a2a194c0eaa001d229e21a3f4f5b447db7c73f7a3a44c3d0464f317bda5f12a2; URL="$VCP_PERTURB_PREFIX/$FILE" ;;
  replogle_k562_essential_perturbation.h5ad) EXPECTED_BYTES=2890786004; EXPECTED_SHA256=04b7b3c28504ace115bb6ee192a0710f428aee3468cf0372b4fc0978cc05adb4; URL="$VCP_PERTURB_PREFIX/$FILE" ;;
  GSE306915_normalized_counts_CO065.csv.gz) EXPECTED_BYTES=9358021; EXPECTED_SHA256=02307f1019429530fae91d8da3d808a1c8e04241fe4657832205c94d01f43d42; URL="$GEO_GSE306915_PREFIX/suppl/$FILE" ;;
  GSE306915_family.soft.gz) EXPECTED_BYTES=10412; EXPECTED_SHA256=9059377ff91eee08ba71b52c787d4166baa0b2e29a9a3b02ba29566c63bbe5c4; URL="$GEO_GSE306915_PREFIX/soft/$FILE" ;;
  rna_seq_meta_key.csv) EXPECTED_BYTES=3671; EXPECTED_SHA256=ba27f8502a517dab3c25c7c8001e85e303659e42d84e96c0e48b20d51fbe3e2f; EXPECTED_MD5=ea121d9634ab466b77dbf6978762b93a; URL="$GOUDY_AUTHOR_PREFIX/$FILE" ;;
  Genome-wide-screens.zip) EXPECTED_BYTES=26152593; EXPECTED_SHA256=15571c41d76b2462d15f167f8920b0ec335f685b1582d18b0264f65f21b2fefd; EXPECTED_MD5=e0392eb7b2512720bb8cbf705ce9854f; URL="$SCHMIDT_ZENODO_PREFIX/$FILE/content" ;;
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
command -v curl >/dev/null || { echo "curl is required" >&2; exit 2; }
CURL_ARGS=(--fail --location --continue-at - --output "$PART")
if [[ -n "$EXPECTED_ETAG" ]]; then
  CURL_ARGS+=(--header "If-Match: $EXPECTED_ETAG")
fi
curl "${CURL_ARGS[@]}" "$URL"

ACTUAL_BYTES="$(wc -c < "$PART" | tr -d ' ')"
if [[ "$ACTUAL_BYTES" != "$EXPECTED_BYTES" ]]; then
  echo "Byte-length mismatch: expected $EXPECTED_BYTES, found $ACTUAL_BYTES" >&2
  exit 1
fi
ACTUAL_SHA256="$(shasum -a 256 "$PART" | awk '{print $1}')"
if [[ -n "$EXPECTED_SHA256" && "$ACTUAL_SHA256" != "$EXPECTED_SHA256" ]]; then
  echo "SHA-256 mismatch: expected $EXPECTED_SHA256, found $ACTUAL_SHA256" >&2
  exit 1
fi
if [[ -n "$EXPECTED_MD5" ]]; then
  if command -v md5 >/dev/null; then
    ACTUAL_MD5="$(md5 -q "$PART")"
  elif command -v md5sum >/dev/null; then
    ACTUAL_MD5="$(md5sum "$PART" | awk '{print $1}')"
  else
    echo "md5 or md5sum is required for this registered object" >&2
    exit 2
  fi
  if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "MD5 mismatch: expected $EXPECTED_MD5, found $ACTUAL_MD5" >&2
    exit 1
  fi
fi
mv "$PART" "$OUT"
echo "$ACTUAL_SHA256  $OUT"
