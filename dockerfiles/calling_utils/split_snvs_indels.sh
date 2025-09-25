#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
#     Split VCF into SNVs and Indels
#       - Wraps bcftools view
#       - Creates compressed, indexed VCF files
#
# Requires bcftools in PATH
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz [-o prefix]

  -i   Input VCF (bgzipped) with .tbi index (required)
  -o   Output prefix (default: output)
  -h   Help
EOF
  exit 1
}

# Defaults
INPUT_VCF=""
OUTPUT_PRFX="output"
NTHREADS="$(nproc)"

# Parse args
while [[ $# -gt 0 ]]; do
  case $1 in
    -i) INPUT_VCF="$2"; shift 2 ;;
    -o) OUTPUT_PRFX="$2"; shift 2 ;;
    -h) usage ;;
    *)  echo "Unknown option: $1"; usage ;;
  esac
done

# Required checks
[[ -n "$INPUT_VCF" ]] || { echo "Error: input VCF (-i) is required"; usage; }

# Tool check
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }

# Input checks
[[ -f "$INPUT_VCF" ]] || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]] || { echo "Error: ${INPUT_VCF}.tbi index not found"; exit 1; }

# Output files
SNVS_OUT="${OUTPUT_PRFX}_snvs.vcf.gz"
INDELS_OUT="${OUTPUT_PRFX}_indels.vcf.gz"

# Run bcftools view for SNVs
echo "Extracting SNVs..."
bcftools view --threads "$NTHREADS" -v snps -Oz -o "$SNVS_OUT" "$INPUT_VCF" || { echo "Error: bcftools view for SNVs failed"; exit 1; }
bcftools index --threads "$NTHREADS" --tbi "$SNVS_OUT" || { echo "Error: bcftools index for SNVs failed"; exit 1; }

# Run bcftools view for indels
echo "Extracting indels..."
bcftools view --threads "$NTHREADS" -v indels -Oz -o "$INDELS_OUT" "$INPUT_VCF" || { echo "Error: bcftools view for Indels failed"; exit 1; }
bcftools index --threads "$NTHREADS" --tbi "$INDELS_OUT" || { echo "Error: bcftools index for Indels failed"; exit 1; }

# Done
echo "Done:"
echo "  SNVs   : $SNVS_OUT"
echo "  Indels : $INDELS_OUT"
