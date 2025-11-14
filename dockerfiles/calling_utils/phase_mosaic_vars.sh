#!/bin/bash
set -euo pipefail

#################################################################################
# PacBio Long-Read Phasing Pipeline (last of step2 filters)
#
# Requires:
#   - bcftools
#   - tabix
#   - bgzip
#   - python3 (including pysam)
#   - phasing_step1_get_closest_germline.sh
#   - phasing_step2_phase_mosaic.py
#
#################################################################################

usage() {
  cat <<EOF
Usage: $0 \\
    -i SAMPLE_ID \\
    -v GERMLINE_VCF \\
    -w INPUT_VCF \\
    -s SEX \\
    [--pb-cram CRAM ...] \\
    [-t THREADS]

Required arguments:
    -i          Sample ID
    -v          Germline VCF
    -w          Input mosaic VCF (SNV VCF)
    -s          Sex (M/F)
	--pb-cram   PacBio long-read CRAM with .crai file, also accept BAM with .bai or .csi index (repeatable)

    -t          Threads for phasing_step2_phase_mosaic.py

Example:
    $0 -i SMHT005 -v SMHT005_germline.vcf.gz -w SMHT005-3A_mosaic.vcf.gz --pb-cram a.bam --pb-cram b.bam -s M -t 8

EOF
  exit 1
}

#################################################################################
# Defaults
#################################################################################

SAMPLE_ID=""
GERMLINE_VCF=""
INPUT_VCF=""
SEX=""
THREADS="$(nproc)"

PB_CRAMS=()

#################################################################################
# Parse arguments
#################################################################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) SAMPLE_ID="$2"; shift 2;;
    -v) GERMLINE_VCF="$2"; shift 2;;
    -w) INPUT_VCF="$2"; shift 2;;
    -s) SEX="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;

	--pb-cram) PB_CRAMS+=("$2"); shift 2;;

    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

#################################################################################
# Required checks
#################################################################################

[[ -n "$SAMPLE_ID" ]]      || { echo "Error: -i SAMPLE_ID is required"; usage; }
[[ -n "$GERMLINE_VCF" ]]   || { echo "Error: -v GERMLINE_VCF is required"; usage; }
[[ -n "$INPUT_VCF" ]]      || { echo "Error: -w INPUT_VCF is required"; usage; }
[[ -n "$SEX" ]]            || { echo "Error: -s SEX is required"; usage; }

[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "Error: THREADS (-t) must be integer"; exit 1; }

[[ -f "$GERMLINE_VCF" ]] || { echo "Error: Germline VCF not found: $GERMLINE_VCF"; exit 1; }
[[ -f "$INPUT_VCF" ]]    || { echo "Error: Input VCF not found: $INPUT_VCF"; exit 1; }



# CRAM/BAM checks (index check too)
if (( ${#PB_CRAMS[@]} == 0 )); then
  echo "Error: no BAM/CRAM files provided. Please specify at least one with --pb-cram"
  exit 1
else
  for b in "${PB_CRAMS[@]}"; do
    [[ -f "$b" ]] || { echo "Error: file not found: $b"; exit 1; }

    # Determine file type based on extension
    if [[ "$b" == *.cram ]]; then
      if [[ ! -f "${b}.crai" ]]; then
        echo "Error: missing index for CRAM: ${b}. Expected ${b}.crai"
        exit 1
      fi
    elif [[ "$b" == *.bam ]]; then
      if [[ ! -f "${b}.bai" && ! -f "${b}.csi" ]]; then
        echo "Error: missing index for BAM: ${b}. Expected ${b}.bai or ${b}.csi"
        exit 1
      fi
    else
      echo "Error: unsupported file type (must be .bam or .cram): $b"
      exit 1
    fi
  done
fi

#################################################################################
# Tool checks
#################################################################################

command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found"; exit 1; }
command -v tabix     >/dev/null 2>&1 || { echo "Error: tabix not found"; exit 1; }
command -v bgzip     >/dev/null 2>&1 || { echo "Error: bgzip not found"; exit 1; }
command -v python3   >/dev/null 2>&1 || { echo "Error: python3 not found"; exit 1; }
command -v phasing_step2_phase_mosaic.py >/dev/null 2>&1 || { echo "Error: phasing_step2_phase_mosaic.py not found in PATH"; exit 1; }
command -v phasing_step1_get_closest_germline.sh >/dev/null 2>&1 || { echo "Error: phasing_step1_get_closest_germline.sh not found in PATH"; exit 1; }

#################################################################################
# Work directory
#################################################################################

WORKDIR="$(mktemp -d phasing_pipeline.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

echo "Workdir: $WORKDIR"

STEP4_TSV="${WORKDIR}/${SAMPLE_ID}.germline_map.tsv"
TAGS_TXT="${WORKDIR}/${SAMPLE_ID}_phasing_tags.tsv"
TAGS_GZ="${TAGS_TXT}.gz"
ANNOTATED_VCF="${WORKDIR}/annotated.vcf.gz"
FINAL_VCF="${SAMPLE_ID}.phased.vcf.gz"

#################################################################################
# STEP 4 – long-read phasing tag generation
#################################################################################

echo "------------------------------------------------------------"
echo "[STEP 4] Running phasing_step1_get_closest_germline.sh"
echo "------------------------------------------------------------"

phasing_step1_get_closest_germline.sh \
    -s "$SAMPLE_ID" \
    -v "$GERMLINE_VCF" \
    -w "$INPUT_VCF"

# generates "${WORKDIR}/${SAMPLE_ID}.germline_map.tsv" (STEP4_TSV)
[[ -s "$STEP4_TSV" ]] || { echo "Error: Step4 TSV is empty"; exit 1; }

#################################################################################
# STEP 5 – long-read haplotype classification
#################################################################################

echo "------------------------------------------------------------"
echo "[STEP 5] Running phasing_step2_phase_mosaic.py"
echo "------------------------------------------------------------"


# -b is a multi-arg, ie: -b bam1.bam  bam2.bam
phasing_step2_phase_mosaic.py \
    -w "$THREADS" \
    -t "$STEP4_TSV" \
    -b "$PB_CRAMS" \
    -s "$SEX" \
    -i "$SAMPLE_ID"

[[ -f "$TAGS_TXT" ]] || { echo "Error: phasing tags not produced"; exit 1; }

#################################################################################
# Compress and index tags
#################################################################################

echo "Compressing + indexing phasing tags..."

sort -k1,1V -k2,2n "$TAGS_TXT" | bgzip -c > "$TAGS_GZ"
tabix -s1 -b2 -e2 "$TAGS_GZ"

#################################################################################
# Annotate VCF
#################################################################################

echo "Annotating VCF with phasing information..."

bcftools annotate \
  -a "$TAGS_GZ" \
  -c CHROM,POS,PHASING \
  -H '##INFO=<ID=PHASING,Number=1,Type=String,Description="Phasing classification from long-read haplotyping">' \
  -Oz -o "$ANNOTATED_VCF" \
  "$INPUT_VCF"

tabix "$ANNOTATED_VCF"

#################################################################################
# Filter for passing variants
#################################################################################

echo "Filtering variants for TIER2 / interesting phasing categories..."

bcftools view \
  -i '(FILTER="TIER2") || (INFO/PHASING="MOSAIC_PHASED") || (INFO/PHASING="UNABLE_TO_PHASE")' \
  "$ANNOTATED_VCF" \
  -Oz -o "$FINAL_VCF"

tabix "$FINAL_VCF"

#################################################################################
# DONE
#################################################################################

echo "------------------------------------------------------------"
echo "Phasing pipeline complete"
echo "  Final annotated VCF:  $ANNOTATED_VCF"
echo "  Final filtered VCF:   $FINAL_VCF"
echo "------------------------------------------------------------"
