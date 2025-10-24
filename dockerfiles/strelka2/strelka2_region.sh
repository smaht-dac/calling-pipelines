#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

###############################################################################
# Run Strelka2 Somatic with a regions BED
#
# - Accepts normal and tumor BAM/CRAM files with indexes
# - Optional: run on a single BED line (-i) for array jobs
# - Requires: Strelka2 (configureStrelkaSomaticWorkflow.py)
# - Optional (only if -i is used): gunzip, bgzip, tabix
###############################################################################

usage() {
  cat <<EOF
Usage: $0 -n normal.bam -t tumor.bam -f reference.fa -b regions.bed.gz \\
          [-i index] [-o prefix] [-m local|sge] [-j threads] [-h]

  -n   Normal BAM/CRAM with index (.bai/.crai) (required)
  -t   Tumor  BAM/CRAM with index (.bai/.crai) (required)
  -f   Reference FASTA with index (.fai) (required)
  -b   Regions BED (bgzipped + tabix) to restrict calling (required)
  -i   1-based line index in the regions BED to run only that region (optional). Assumes no header lines.
  -o   Output prefix (default: output) -> <prefix>.snvs.vcf.gz / <prefix>.indels.vcf.gz
  -m   Run mode for runWorkflow.py (local|sge, default: local)
  -j   Threads/parallel jobs (default: number of cores)
  -h   Help
EOF
  exit 1
}

# Defaults
NORMAL_BAM=""
TUMOR_BAM=""
REFERENCE_FASTA=""
REGIONS_BED=""
LINE_INDEX=""
OUT_PREFIX="output"
RUN_MODE="local"
JOBS="$(nproc)"

# Parse args
while getopts "n:t:f:b:i:o:m:j:h" opt; do
  case "$opt" in
    n) NORMAL_BAM="$OPTARG" ;;
    t) TUMOR_BAM="$OPTARG" ;;
    f) REFERENCE_FASTA="$OPTARG" ;;
    b) REGIONS_BED="$OPTARG" ;;
    i) LINE_INDEX="$OPTARG" ;;
    o) OUT_PREFIX="$OPTARG" ;;
    m) RUN_MODE="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ -n "$NORMAL_BAM" && -n "$TUMOR_BAM" && -n "$REFERENCE_FASTA" && -n "$REGIONS_BED" ]] || usage

# Tool checks
command -v configureStrelkaSomaticWorkflow.py >/dev/null 2>&1 || { echo "ERROR: configureStrelkaSomaticWorkflow.py not found in PATH"; exit 1; }

# File checks
[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

[[ -f "$NORMAL_BAM" ]] || { echo "Error: $NORMAL_BAM not found"; exit 1; }
[[ -f "${NORMAL_BAM}.crai" || -f "${NORMAL_BAM}.bai" ]] || { echo "Error: index for $NORMAL_BAM not found (.crai or .bai)"; exit 1; }

[[ -f "$TUMOR_BAM" ]] || { echo "Error: $TUMOR_BAM not found"; exit 1; }
[[ -f "${TUMOR_BAM}.crai" || -f "${TUMOR_BAM}.bai" ]] || { echo "Error: index for $TUMOR_BAM not found (.crai or .bai)"; exit 1; }

[[ -f "$REGIONS_BED" ]] || { echo "Error: $REGIONS_BED not found"; exit 1; }
[[ -f "${REGIONS_BED}.tbi" ]] || { echo "Error: ${REGIONS_BED}.tbi not found"; exit 1; }

# Mode sanity
if [[ "$RUN_MODE" != "local" && "$RUN_MODE" != "sge" ]]; then
  echo "ERROR: -m must be local or sge"; exit 1
fi

# Run directory
RUN_DIR="$(mktemp -d run.XXXXXX)"
# Cleanup
cleanup() { rm -rf "$RUN_DIR"; }
trap cleanup EXIT

# If LINE_INDEX provided, ensure it's a positive integer and extract that region
if [[ -n "$LINE_INDEX" ]]; then
  # Checks - LINE_INDEX is a positive integer
  [[ "$LINE_INDEX" =~ ^[0-9]+$ && "$LINE_INDEX" -ge 1 ]] || { echo "ERROR: -i must be a positive integer (1-based)"; exit 1; }
  # Checks - bgzip, tabix, gunzip
  command -v bgzip >/dev/null 2>&1 || { echo "ERROR: bgzip not found in PATH"; exit 1; }
  command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found in PATH"; exit 1; }
  command -v gunzip >/dev/null 2>&1 || { echo "ERROR: gunzip not found in PATH"; exit 1; }

  # Extract region line to new BED
  REGIONS_BED_TMP="${RUN_DIR}/regions.bed"
  gunzip -c "$REGIONS_BED" > "$REGIONS_BED_TMP"

  REGION_LINE="$(sed -n "${LINE_INDEX}p" "$REGIONS_BED_TMP" | tr -d '\n\r')"
  if [[ -z "$REGION_LINE" ]]; then
    echo "ERROR: No region found at line $LINE_INDEX in $REGIONS_BED"; exit 1
  fi
  echo "==> Using region from line $LINE_INDEX: $REGION_LINE"

  LINE_BED_TMP="${RUN_DIR}/line.bed"
  echo "$REGION_LINE" > "$LINE_BED_TMP"

  bgzip -f "$LINE_BED_TMP"
  tabix -f -p bed "${LINE_BED_TMP}.gz"
  REGIONS_BED="${LINE_BED_TMP}.gz"
fi

# Strelka2
STRELKA_DIR="${RUN_DIR}/tmp.strelka"
echo "==> Configuring Strelka2 somatic workflow"
configureStrelkaSomaticWorkflow.py \
  --normalBam "$NORMAL_BAM" \
  --tumorBam "$TUMOR_BAM" \
  --referenceFasta "$REFERENCE_FASTA" \
  --callRegions "$REGIONS_BED" \
  --runDir "$STRELKA_DIR"

echo "==> Running Strelka2 workflow (-m $RUN_MODE, -j $JOBS)"
"${STRELKA_DIR}/runWorkflow.py" -m "$RUN_MODE" -j "$JOBS"

# Collect outputs
echo "==> Collecting results"
SNV_VCF="${STRELKA_DIR}/results/variants/somatic.snvs.vcf.gz"
INDEL_VCF="${STRELKA_DIR}/results/variants/somatic.indels.vcf.gz"

SNV_VCF_OUT="${OUT_PREFIX}.snvs.vcf.gz"
INDEL_VCF_OUT="${OUT_PREFIX}.indels.vcf.gz"

cp "$SNV_VCF" "$SNV_VCF_OUT"
cp "$SNV_VCF.tbi" "${SNV_VCF_OUT}.tbi"
cp "$INDEL_VCF" "$INDEL_VCF_OUT"
cp "$INDEL_VCF.tbi" "${INDEL_VCF_OUT}.tbi"

echo "Done."
echo "  Results:"
echo "    - SNVs  : $SNV_VCF_OUT"
echo "    - Indels: $INDEL_VCF_OUT"
