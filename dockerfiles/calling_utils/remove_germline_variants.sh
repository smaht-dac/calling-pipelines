#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
#     Remove germline variants from a somatic VCF using bcftools isec
#       - Produces variants present in INPUT_VCF but NOT in GERMLINE_VCF
#       - Output bgzipped VCF + tabix index
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -g germline.vcf.gz [-o prefix]

  -i   Input VCF to filter (bgzipped) with .tbi index (required)
  -g   Germline VCF to exclude (bgzipped) with .tbi index (required)
  -o   Output prefix (default: output) -> <prefix>.nogermline.vcf.gz
  -h   Help
EOF
  exit 1
}

# Defaults
INPUT_VCF=""
GERMLINE_VCF=""
OUTPUT_PRFX="output"
NTHREADS="$(nproc)"

# Parse args
while getopts "i:g:o:h" opt; do
  case $opt in
    i) INPUT_VCF="$OPTARG" ;;
    g) GERMLINE_VCF="$OPTARG" ;;
    o) OUTPUT_PRFX="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ -n "$INPUT_VCF" ]] || { echo "Error: input VCF (-i) is required"; usage; }
[[ -n "$GERMLINE_VCF" ]] || { echo "Error: germline VCF (-g) is required"; usage; }

# Tool checks
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }
command -v bgzip    >/dev/null 2>&1 || { echo "Error: bgzip not found in PATH"; exit 1; }
command -v tabix    >/dev/null 2>&1 || { echo "Error: tabix not found in PATH"; exit 1; }

# Input checks
[[ -f "$INPUT_VCF" ]]     || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "$GERMLINE_VCF" ]]  || { echo "Error: $GERMLINE_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]]    || { echo "Error: ${INPUT_VCF}.tbi index not found"; exit 1; }
[[ -f "${GERMLINE_VCF}.tbi" ]] || { echo "Error: ${GERMLINE_VCF}.tbi index not found"; exit 1; }

INT_GERMLINE_PASS_VCF="${OUTPUT_PRFX}_germline_snvs.vcf.gz"

# Cleanup
TMP_DIR=""
cleanup() { [[ -n "${TMP_DIR:-}" && -d "$TMP_DIR" ]] && rm -rf "$TMP_DIR"; }
trap cleanup EXIT

TMP_DIR="$(mktemp -d ./isec_tmp.XXXXXX)"

bcftools view -o "${INT_GERMLINE_PASS_VCF}" --write-index=tbi -f '.' -Ob \
	"$GERMLINE_VCF"

# Core logic (from script 1): isec -n~01, then take 0001.vcf
# Interpretation: output records unique to the 2nd file (INPUT_VCF) i.e., not in GERMLINE_VCF.
bcftools isec -n~01 -p "$TMP_DIR" \
  "${INT_GERMLINE_PASS_VCF}" \
  "$INPUT_VCF" \
  >/dev/null

OUT_VCF="${OUTPUT_PRFX}.nogermline.vcf"
OUT_GZ="${OUT_VCF}.gz"

if [[ ! -s "${TMP_DIR}/0001.vcf" ]]; then
  echo "Error: expected output not found or empty: ${TMP_DIR}/0001.vcf"
  exit 1
fi

mv "${TMP_DIR}/0001.vcf" "$OUT_VCF"

# Compress + index (preserving script 1 behavior)
bgzip -f "$OUT_VCF"
tabix -f -p vcf "$OUT_GZ"

echo "Done:"
echo "  VCF : $OUT_GZ"
echo "  TBI : ${OUT_GZ}.tbi"

