#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
#     Preprocessing VCF files
#       - (NEW) If 2+ inputs: merge first
#       - Filter PASS (or PASS/FEX or longcallD special-case)
#       - Normalize and atomize
#       - Remove duplicates
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz [-i input2.vcf.gz ...] -r reference.fasta [-o prefix]

  -i   Input VCF (bgzipped) with .tbi index (required). Can be provided multiple times.
  -r   Reference FASTA with .fai index (required)
  -o   Prefix for the output annotated VCF (default: output)
       The final output will be <prefix>.vcf.gz
  -h   Help
EOF
  exit 1
}

# Default
REFERENCE_FASTA=""
OUTPUT_PRFX="output"
NTHREADS=$(nproc)

# Allow multiple -i
INPUT_VCFS=()

# Parse args (note: getopts supports repeated flags)
while getopts "i:r:o:h" opt; do
  case $opt in
    i) INPUT_VCFS+=("$OPTARG") ;;
    r) REFERENCE_FASTA="$OPTARG" ;;
    o) OUTPUT_PRFX="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ ${#INPUT_VCFS[@]} -ge 1 ]] || { echo "Error: at least one input VCF (-i) is required"; usage; }
[[ -n "$REFERENCE_FASTA" ]] || { echo "Error: reference FASTA (-r) is required"; usage; }

# Tool checks
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }

# Reference checks
[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai index not found"; exit 1; }

# Input checks
for vcf in "${INPUT_VCFS[@]}"; do
  [[ -f "$vcf" ]] || { echo "Error: $vcf not found"; exit 1; }
  [[ -f "${vcf}.tbi" ]] || { echo "Error: ${vcf}.tbi index not found"; exit 1; }
done

# Temp files
TMP_HEAD="$(mktemp ./HEAD_XXXXXX)"
TMP_MERGED="$(mktemp ./MERGED_XXXXXX.vcf.gz)"
TMP_MERGED_TBI="${TMP_MERGED}.tbi"

cleanup() {
  rm -f "$TMP_HEAD" "$TMP_MERGED" "$TMP_MERGED_TBI"
}
trap cleanup EXIT

# If 2+ inputs, merge them first; else use the single input as-is.
if [[ ${#INPUT_VCFS[@]} -ge 2 ]]; then
  echo "Merging ${#INPUT_VCFS[@]} VCFs..."
  # Merge samples/records across files. Assumes same reference contigs ordering.
  # -m all keeps all records; adjust if you want different multiallelic handling.
  bcftools merge --force-samples --threads "$NTHREADS" -m all -Oz -o "$TMP_MERGED" "${INPUT_VCFS[@]}" \
    || { echo "Error: bcftools merge failed"; exit 1; }
  bcftools index --threads "$NTHREADS" --tbi "$TMP_MERGED" \
    || { echo "Error: bcftools index of merged VCF failed"; exit 1; }
  INPUT_FOR_PIPELINE="$TMP_MERGED"
else
  INPUT_FOR_PIPELINE="${INPUT_VCFS[0]}"
fi

# Check header to set the filter correctly
bcftools view --header-only "$INPUT_FOR_PIPELINE" > "$TMP_HEAD"

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
bcftools view -v snps,indels --threads "$NTHREADS" -i "$FILTER" -Ou "$INPUT_FOR_PIPELINE" \
  | bcftools norm --threads "$NTHREADS" --check-ref x -m -any --atomize -f "$REFERENCE_FASTA" -Ou - \
  | bcftools norm --threads "$NTHREADS" -d exact -Ou - \
  | bcftools sort -T "tmp_bcftools.XXXXXX" -Oz -o "${OUTPUT_PRFX}.vcf.gz" \
  || { echo "Error: bcftools normalization failed"; exit 1; }

bcftools index --threads "$NTHREADS" --tbi "${OUTPUT_PRFX}.vcf.gz" \
  || { echo "Error: bcftools index failed"; exit 1; }

echo "Done:"
echo "  VCF : ${OUTPUT_PRFX}.vcf.gz"
echo "  TBI : ${OUTPUT_PRFX}.vcf.gz.tbi"
