#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

###############################################################################
# Run Manta in tumor-only mode (whole genome)
#
# - Accepts one or more tumor CRAM files with index (.crai)
# - When more than one CRAM is given, they are merged internally with
#   "samtools merge" before calling Manta; a single CRAM is passed directly
# - Requires: Manta (configManta.py), Samtools (samtools)
###############################################################################

usage() {
  cat <<EOF
Usage: $0 -t tumor.cram [-t tumor2.cram ...] -r reference.fa [-o prefix] [-j threads] [-h]

  -t   Tumor CRAM with index (.crai) (required, repeatable). If more than one is
       given, they are merged internally before calling Manta
  -r   Reference FASTA with index (.fai) (required)
  -o   Output prefix (default: output) -> <prefix>.tumorSV.vcf.gz / <prefix>.candidateSV.vcf.gz
  -j   Threads/parallel jobs for runWorkflow.py and samtools merge (default: number of cores)
  -h   Help
EOF
  exit 1
}

# Defaults
TUMOR_CRAMS=()
REFERENCE_FASTA=""
OUT_PREFIX="output"
JOBS="$(nproc)"

# Parse args
while getopts "t:r:o:j:h" opt; do
  case "$opt" in
    t) TUMOR_CRAMS+=("$OPTARG") ;;
    r) REFERENCE_FASTA="$OPTARG" ;;
    o) OUT_PREFIX="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ "${#TUMOR_CRAMS[@]}" -ge 1 && -n "$REFERENCE_FASTA" ]] || usage

# Tool checks
command -v configManta.py >/dev/null 2>&1 || { echo "ERROR: configManta.py not found in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH"; exit 1; }

# File checks
for cram in "${TUMOR_CRAMS[@]}"; do
  [[ -f "$cram" ]] || { echo "Error: $cram not found"; exit 1; }
  [[ -f "${cram}.crai" ]] || { echo "Error: ${cram}.crai not found"; exit 1; }
done

[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

# Run directory
RUN_DIR="$(mktemp -d run.XXXXXX)"
# Cleanup
cleanup() { rm -rf "$RUN_DIR"; }
trap cleanup EXIT

MANTA_DIR="${RUN_DIR}/tmp.manta"

# **********************************************
# 0. Merge input CRAMs (only when more than one)
# **********************************************
if [[ "${#TUMOR_CRAMS[@]}" -gt 1 ]]; then
  echo "==> Merging ${#TUMOR_CRAMS[@]} input CRAMs"
  TUMOR_INPUT="${RUN_DIR}/merged.cram"
  samtools merge -f --write-index -@ "$JOBS" \
    --reference "$REFERENCE_FASTA" \
    "$TUMOR_INPUT" "${TUMOR_CRAMS[@]}" || { echo "Error: samtools merge failed"; exit 1; }
else
  echo "==> Single input CRAM, skipping merge"
  TUMOR_INPUT="${TUMOR_CRAMS[0]}"
fi

# **********************************************
# 1. Configure Manta tumor-only workflow
# **********************************************
echo "==> Configuring Manta tumor-only workflow"
configManta.py \
  --tumorBam "$TUMOR_INPUT" \
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
