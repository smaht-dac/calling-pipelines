#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
#     Filter VCF by excluding multiple BED region lists
#       - Accept any number of -x <bed> files (plain or .gz)
#       - Exclude with bcftools (-T ^file) using all BED files
#       - Output bgzipped VCF + tabix index
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -x regions1.bed [-x regions2.bed ...] [-o prefix]

  -i   Input VCF (bgzipped) with .tbi index (required)
  -x   BED (or bgzipped BED) with regions to exclude (repeatable)
  -o   Output prefix (default: output) -> <prefix>.vcf.gz
  -h   Help
EOF
  exit 1
}

# Defaults
INPUT_VCF=""
OUTPUT_PRFX="output"
EXCLUDE_BEDS=()
NTHREADS="$(nproc)"

# Parse args
while getopts "i:x:o:h" opt; do
  case $opt in
    i) INPUT_VCF="$OPTARG" ;;
    x) EXCLUDE_BEDS+=("$OPTARG") ;;
    o) OUTPUT_PRFX="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ -n "$INPUT_VCF" ]] || { echo "Error: input VCF (-i) is required"; usage; }
(( ${#EXCLUDE_BEDS[@]} > 0 )) || { echo "Error: at least one BED (-x) is required"; usage; }

# Tool checks
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }

# Input checks
[[ -f "$INPUT_VCF" ]] || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]] || { echo "Error: ${INPUT_VCF}.tbi index not found"; exit 1; }
for bed in "${EXCLUDE_BEDS[@]}"; do
  [[ -f "$bed" ]] || { echo "Error: BED file not found: $bed"; exit 1; }
done

# Cleanup
TMP_FILES=()
cleanup() { rm -f "${TMP_FILES[@]}"; }
trap cleanup EXIT

# Apply each BED filter sequentially
CURRENT_IN="$INPUT_VCF"
for bed in "${EXCLUDE_BEDS[@]}"; do
  TMP_OUT="$(mktemp ./XXXXXX.vcf)"; TMP_FILES+=("$TMP_OUT")

  echo "Excluding regions in $bed ..."
  bcftools view --threads "$NTHREADS" -T "^$bed" -Ov -o "$TMP_OUT" "$CURRENT_IN"
  CURRENT_IN="$TMP_OUT"
done

# Compress-rename output and index
bcftools view --threads "$NTHREADS" -Oz -o "${OUTPUT_PRFX}.vcf.gz" "$CURRENT_IN"
bcftools index --threads "$NTHREADS" --tbi "${OUTPUT_PRFX}.vcf.gz"

# All done
echo "Done:"
echo "  VCF : ${OUTPUT_PRFX}.vcf.gz"
echo "  TBI : ${OUTPUT_PRFX}.vcf.gz.tbi"
