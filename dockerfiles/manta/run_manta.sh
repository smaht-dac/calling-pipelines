#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

###############################################################################
# Run Manta in tumor-only mode (whole genome)
#
# - Accepts a single tumor CRAM file with index
# - Requires: Manta (configManta.py)
###############################################################################

usage() {
  cat <<EOF
Usage: $0 -t tumor.cram -r reference.fa [-o prefix] [-j threads] [-h]

  -t   Tumor CRAM with index (.crai) (required)
  -r   Reference FASTA with index (.fai) (required)
  -o   Output prefix (default: output) -> <prefix>.tumorSV.vcf.gz / <prefix>.candidateSV.vcf.gz
  -j   Threads/parallel jobs for runWorkflow.py (default: number of cores)
  -h   Help
EOF
  exit 1
}

# Defaults
TUMOR_CRAM=""
REFERENCE_FASTA=""
OUT_PREFIX="output"
JOBS="$(nproc)"

# Parse args
while getopts "t:r:o:j:h" opt; do
  case "$opt" in
    t) TUMOR_CRAM="$OPTARG" ;;
    r) REFERENCE_FASTA="$OPTARG" ;;
    o) OUT_PREFIX="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ -n "$TUMOR_CRAM" && -n "$REFERENCE_FASTA" ]] || usage

# Tool checks
command -v configManta.py >/dev/null 2>&1 || { echo "ERROR: configManta.py not found in PATH"; exit 1; }

# File checks
[[ -f "$TUMOR_CRAM" ]] || { echo "Error: $TUMOR_CRAM not found"; exit 1; }
[[ -f "${TUMOR_CRAM}.crai" ]] || { echo "Error: ${TUMOR_CRAM}.crai not found"; exit 1; }

[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

# Run directory
RUN_DIR="$(mktemp -d run.XXXXXX)"
# Cleanup
cleanup() { rm -rf "$RUN_DIR"; }
trap cleanup EXIT

MANTA_DIR="${RUN_DIR}/tmp.manta"

# **********************************************
# 1. Configure Manta tumor-only workflow
# **********************************************
echo "==> Configuring Manta tumor-only workflow"
configManta.py \
  --tumorBam "$TUMOR_CRAM" \
  --referenceFasta "$REFERENCE_FASTA" \
  --runDir "$MANTA_DIR" || { echo "Error: configManta.py failed"; exit 1; }

# **********************************************
# 2. Run Manta workflow
# **********************************************
echo "==> Running Manta workflow (-m local -j $JOBS)"
"${MANTA_DIR}/runWorkflow.py" -m local -j "$JOBS" || { echo "Error: runWorkflow.py failed"; exit 1; }

# **********************************************
# 3. Collect outputs
# **********************************************
echo "==> Collecting results"
TUMOR_SV_VCF="${MANTA_DIR}/results/variants/tumorSV.vcf.gz"
CANDIDATE_SV_VCF="${MANTA_DIR}/results/variants/candidateSV.vcf.gz"

TUMOR_SV_VCF_OUT="${OUT_PREFIX}.tumorSV.vcf.gz"
CANDIDATE_SV_VCF_OUT="${OUT_PREFIX}.candidateSV.vcf.gz"

cp "$TUMOR_SV_VCF" "$TUMOR_SV_VCF_OUT" || { echo "Error: tumorSV.vcf.gz not found"; exit 1; }
cp "${TUMOR_SV_VCF}.tbi" "${TUMOR_SV_VCF_OUT}.tbi" || { echo "Error: tumorSV.vcf.gz.tbi not found"; exit 1; }
cp "$CANDIDATE_SV_VCF" "$CANDIDATE_SV_VCF_OUT" || { echo "Error: candidateSV.vcf.gz not found"; exit 1; }
cp "${CANDIDATE_SV_VCF}.tbi" "${CANDIDATE_SV_VCF_OUT}.tbi" || { echo "Error: candidateSV.vcf.gz.tbi not found"; exit 1; }

echo "Done."
echo "  Results:"
echo "    - Tumor SV     : $TUMOR_SV_VCF_OUT"
echo "    - Candidate SV : $CANDIDATE_SV_VCF_OUT"
