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
Usage: $0 -i input.vcf.gz -f additional.vcf.gz [-f additional.vcf.gz ...] -r reference.fasta [-o prefix]

  -i   Input VCF (bgzipped) with .tbi index (required)
  -f   Additional input VCF file to merge to main input(bgzipped with .tbi index, optional, can be specified multiple times)
  -r   Reference FASTA with .fai index (required)
  -o   Prefix for the output annotated VCF (default: output)
       The final output will be <prefix>.vcf.gz
  -h   Help
EOF
  exit 1
}

# Default
INPUT_VCF=""
ADDITIONAL_VCFS=()
REFERENCE_FASTA=""
OUTPUT_PRFX="output"
# Define threads
NTHREADS=$(nproc)

# Parse args
while getopts "i:f:r:o:h" opt; do
  case $opt in
    i) INPUT_VCF="$OPTARG" ;;
    f) ADDITIONAL_VCFS+=("$OPTARG") ;;
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

if (( ${#ADDITIONAL_VCFS[@]} > 0 )); then
    for vcf in "${ADDITIONAL_VCFS[@]}"; do
        [[ -f "$vcf" ]] || { echo "Error: additional VCF file not found: $vcf"; exit 1; }
        [[ -f "${vcf}.tbi" ]] || { echo "Error: index for additional VCF file not found: ${vcf}.tbi"; exit 1; }
    done
fi

# Cleanup function to remove temporary files on exit
TMP_FILES=()

cleanup() {
    rm -f "${TMP_FILES[@]}"
}
trap cleanup EXIT

# Concatenating additional VCF files if provided
if (( ${#ADDITIONAL_VCFS[@]} > 0 )); then
    echo "Merging additional VCF files..."
    ALL_VCFS=("$INPUT_VCF" "${ADDITIONAL_VCFS[@]}")
    MERGED_VCF=$(mktemp -u ./MERGED_XXXXXX.vcf.gz)
    TMP_FILES+=("$MERGED_VCF" "$MERGED_VCF.tbi")

    bcftools concat -a -Ou "${ALL_VCFS[@]}" \
    | bcftools sort -T "tmp_bcftools.XXXXXX" -Oz -o "$MERGED_VCF" \
    || { echo "Error: bcftools concat/sort step failed"; exit 1; }

    bcftools index --threads "$NTHREADS" --tbi "$MERGED_VCF" || { echo "Error: bcftools index failed for merged VCF"; exit 1; }
    INPUT_VCF="$MERGED_VCF"
fi

# Check header to set the filter correctly
TMP_HEAD="$(mktemp ./HEAD_XXXXXX)"
TMP_FILES+=("$TMP_HEAD")
bcftools view --header-only "$INPUT_VCF" > "$TMP_HEAD"

# Set filter based on header
if grep -q -m1 'ID=FEX' "$TMP_HEAD"; then
    FILTER='FILTER=="PASS" || INFO/FEX=="PASS"'
elif grep -q -m1 'longcallD' "$TMP_HEAD"; then
    # Keep PASS records that are NOT SVs (i.e., SVTYPE tag absent from INFO)
    FILTER='FILTER=="PASS" && INFO/SVTYPE="."'
else
    FILTER='FILTER=="PASS"'
fi

# Process VCF
echo "Processing VCF..."
bcftools view -v snps,indels --threads "$NTHREADS" -i "$FILTER" -Ou "$INPUT_VCF" \
    | bcftools norm --threads "$NTHREADS" --check-ref x -m -any --atomize -f "$REFERENCE_FASTA" -Ou - \
    | bcftools norm --threads "$NTHREADS" -d exact -Ou - \
    | bcftools sort -T "tmp_bcftools.XXXXXX" -Oz -o "${OUTPUT_PRFX}.vcf.gz" \
    || { echo "Error: bcftools normalization failed"; exit 1; }

bcftools index --threads "$NTHREADS" --tbi "${OUTPUT_PRFX}.vcf.gz" || { echo "Error: bcftools index failed"; exit 1; }

# All done
echo "Done:"
echo "  VCF : ${OUTPUT_PRFX}.vcf.gz"
echo "  TBI : ${OUTPUT_PRFX}.vcf.gz.tbi"
