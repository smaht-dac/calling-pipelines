#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
#     Filter VCF by Panel-of-Errors (POE) FASTA
#       - Wraps filter_by_poe.py
#       - Compress and index final VCF
#
# Requires filter_by_poe.py that needs to be executable and in PATH
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -f error_panel.fa [-o prefix] [--extra-args ...]

  -i   Input VCF (bgzipped) with .tbi index (required)
  -f   Error panel FASTA with .fai (required)
  -o   Output prefix (default: output) -> <prefix>.vcf.gz
  -t   Number of threads for FASTA lookups (default: all available)
  -h   Help

Any additional arguments after options are passed directly to filter_by_poe.py
(e.g., --failed-out failed.vcf)
EOF
  exit 1
}

# Defaults
INPUT_VCF=""
PON_FASTA=""
OUTPUT_PRFX="output"
EXTRA_ARGS=()
NTHREADS="$(nproc)"

# Parse args
while [[ $# -gt 0 ]]; do
  case $1 in
    -i) INPUT_VCF="$2"; shift 2 ;;
    -f) PON_FASTA="$2"; shift 2 ;;
    -o) OUTPUT_PRFX="$2"; shift 2 ;;
    -t) NTHREADS="$2"; shift 2 ;;
    -h) usage ;;
    *)  EXTRA_ARGS+=("$1"); shift ;;
  esac
done

# Required checks
[[ -n "$INPUT_VCF" ]] || { echo "Error: input VCF (-i) is required"; usage; }
[[ -n "$PON_FASTA" ]] || { echo "Error: error panel FASTA (-f) is required"; usage; }

# Tool checks
command -v filter_by_poe.py >/dev/null 2>&1 || { echo "Error: filter_by_poe.py not found in PATH"; exit 1; }
command -v bcftools        >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }

# Input checks
[[ -f "$INPUT_VCF" ]] || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]] || { echo "Error: ${INPUT_VCF}.tbi index not found"; exit 1; }
[[ -f "$PON_FASTA" ]] || { echo "Error: $PON_FASTA not found"; exit 1; }
[[ -f "${PON_FASTA}.fai" ]] || { echo "Error: ${PON_FASTA}.fai index not found"; exit 1; }

# Run filter_by_poe.py
OUT_VCF="${OUTPUT_PRFX}.vcf"
trap 'rm -f "$OUT_VCF"' EXIT

echo "Running filter_by_poe.py..."
filter_by_poe.py --vcf "$INPUT_VCF" \
                 --fasta "$PON_FASTA" \
                 --out "$OUT_VCF" \
                 --threads "$NTHREADS" \
                 "${EXTRA_ARGS[@]}" \
                 || { echo "Error: filter_by_poe.py failed"; exit 1; }

# Compress and index
echo "Compressing and indexing..."
bcftools view --threads "$NTHREADS" -Oz -o "${OUT_VCF}.gz" "$OUT_VCF" || { echo "Error: bcftools view failed"; exit 1; }
bcftools index --threads "$NTHREADS" --tbi "${OUT_VCF}.gz" || { echo "Error: bcftools index failed"; exit 1; }

# All done
echo "Done:"
echo "  VCF : ${OUT_VCF}.gz"
echo "  TBI : ${OUT_VCF}.gz.tbi"
