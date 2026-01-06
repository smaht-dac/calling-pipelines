#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
#     Preprocessing VCF files
#       - Filter PASS
#       - Normalize and atomize
#       - Remove duplicates
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -r reference.fasta [-o prefix]

  -i   Input VCF (bgzipped) with .tbi index (required)
  -r   Reference FASTA with .fai index (required)
  -o   Prefix for the output annotated VCF (default: output)
       The final output will be <prefix>.vcf.gz
  -h   Help
EOF
  exit 1
}

# Default
INPUT_VCF=""
REFERENCE_FASTA=""
OUTPUT_PRFX="output"
# Define threads
NTHREADS=$(nproc)

# Parse args
while getopts "i:r:o:h" opt; do
  case $opt in
    i) INPUT_VCF="$OPTARG" ;;
    r) REFERENCE_FASTA="$OPTARG" ;;
    o) OUTPUT_PRFX="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ -n "$INPUT_VCF" ]] || { echo "Error: input VCF (-i) is required"; usage; }
[[ -n "$REFERENCE_FASTA" ]] || { echo "Error: reference FASTA (-r) is required"; usage; }

# Tool checks
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }

# Input checks
[[ -f "$INPUT_VCF" ]] || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]] || { echo "Error: ${INPUT_VCF}.tbi index not found"; exit 1; }
[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

# Check header to set the filter correctly
TMP_HEAD="$(mktemp ./HEAD_XXXXXX)"
bcftools view --header-only "$INPUT_VCF" > "$TMP_HEAD"
trap 'rm -f "$TMP_HEAD"' EXIT

# Set filter based on header
if grep -q -m1 'FEX' "$TMP_HEAD"; then
    FILTER='FILTER=="PASS" || INFO/FEX=="PASS"'
else
    FILTER='FILTER=="PASS"'
fi

# Process VCF
echo "Processing VCF..."
bcftools view --threads "$NTHREADS" -i "$FILTER" -Ou "$INPUT_VCF" \
    | bcftools norm --threads "$NTHREADS" --check-ref x -m -any --atomize -f "$REFERENCE_FASTA" -Ou - \
    | bcftools norm --threads "$NTHREADS" -d exact -Ou - \
    | bcftools sort --threads "$NTHREADS" -T "tmp_bcftools.XXXXXX" -Oz -o "${OUTPUT_PRFX}.vcf.gz" \
    || { echo "Error: bcftools normalization failed"; exit 1; }

bcftools index --threads "$NTHREADS" --tbi "${OUTPUT_PRFX}.vcf.gz" || { echo "Error: bcftools index failed"; exit 1; }

# All done
echo "Done:"
echo "  VCF : ${OUTPUT_PRFX}.vcf.gz"
echo "  TBI : ${OUTPUT_PRFX}.vcf.gz.tbi"
