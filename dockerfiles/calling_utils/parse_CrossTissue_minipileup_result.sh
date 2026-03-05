#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
#   Parse SR-only minipileup VCF -> CrossTissue -> set filters -> final VCF
#
#   Steps (keeps Script 1 logic):
#     1) parse_minipileup_sr_only.py -> <prefix>.CrossTissue.vcf.gz (+tabix)
#     2) set_filter_vcf.py           -> <prefix>.final.vcf.gz      (+tabix)
#################################################################################

usage() {
  cat <<'EOF'
Usage: parse_CrossTissue_minipileup_result.sh -i orig.vcf.gz --mp mp.vcf.gz --tissue TISSUE [-o prefix]

  -i, --orig        Original/input VCF (bgzipped) (required)
  --mp              Minipileup SR-only VCF (bgzipped) (required)
  --tissue          Tissue label passed to parse_minipileup_sr_only.py (required)
  -o, --out-prefix  Output prefix (default: output)
  -h, --help        Help

Outputs:
  <prefix>.final.vcf.gz      (+.tbi)

Notes:
  - Assumes these are on PATH: tabix, python3,
      parse_minipileup_sr_only.py, set_filter_vcf.py
EOF
  exit 1
}

# Defaults
ORIG_VCF=""
MP_VCF=""
TISSUE=""
OUTPUT_PRFX="output"

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--orig)       ORIG_VCF="$2"; shift 2 ;;
    --mp)            MP_VCF="$2"; shift 2 ;;
    --tissue)        TISSUE="$2"; shift 2 ;;
    -o|--out-prefix) OUTPUT_PRFX="$2"; shift 2 ;;
    -h|--help)       usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Required checks
[[ -n "$ORIG_VCF" ]] || { echo "Error: -i/--orig is required"; usage; }
[[ -n "$MP_VCF" ]]   || { echo "Error: --mp is required"; usage; }
[[ -n "$TISSUE" ]]   || { echo "Error: --tissue is required"; usage; }

# Tool checks
command -v tabix >/dev/null 2>&1  || { echo "Error: tabix not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1|| { echo "Error: python3 not found in PATH"; exit 1; }
command -v parse_minipileup_sr_only.py >/dev/null 2>&1 || { echo "Error: parse_minipileup_sr_only.py not found in PATH"; exit 1; }
command -v set_filter_vcf.py >/dev/null 2>&1          || { echo "Error: set_filter_vcf.py not found in PATH"; exit 1; }

# Input checks
[[ -f "$ORIG_VCF" ]] || { echo "Error: original VCF not found: $ORIG_VCF"; exit 1; }
[[ -f "$MP_VCF" ]]   || { echo "Error: minipileup VCF not found: $MP_VCF"; exit 1; }

# (Optional but typical) index checks; warn rather than hard-fail if you prefer.
if [[ "$ORIG_VCF" == *.vcf.gz && ! -f "${ORIG_VCF}.tbi" ]]; then
  echo "Error: original VCF index not found: ${ORIG_VCF}.tbi"
  exit 1
fi
if [[ "$MP_VCF" == *.vcf.gz && ! -f "${MP_VCF}.tbi" ]]; then
  echo "Error: minipileup VCF index not found: ${MP_VCF}.tbi"
  exit 1
fi

# ------------------------------------------------------------------------------
# Step 1: Parse minipileup -> CrossTissue
# ------------------------------------------------------------------------------
CROSSTISSUE_VCF_GZ="${OUTPUT_PRFX}.CrossTissue.vcf.gz"

echo "Step 1/2: Parsing minipileup -> ${CROSSTISSUE_VCF_GZ}"
parse_minipileup_sr_only.py \
  --tissue "$TISSUE" \
  --orig_vcf "$ORIG_VCF" \
  --mp_vcf "$MP_VCF" \
  --out "$CROSSTISSUE_VCF_GZ" \
  || { echo "Error: parse_minipileup_sr_only.py failed"; exit 1; }

tabix -f -p vcf "$CROSSTISSUE_VCF_GZ" \
  || { echo "Error: tabix failed for $CROSSTISSUE_VCF_GZ"; exit 1; }

# ------------------------------------------------------------------------------
# Step 2: Set filters -> final VCF
# ------------------------------------------------------------------------------
FINAL_VCF_GZ="${OUTPUT_PRFX}.final.vcf.gz"

echo "Step 2/2: Setting filters -> ${FINAL_VCF_GZ}"
set_filter_vcf.py \
  -i "$CROSSTISSUE_VCF_GZ" \
  -o "$FINAL_VCF_GZ" \
  || { echo "Error: set_filter_vcf.py failed"; exit 1; }

tabix -f -p vcf "$FINAL_VCF_GZ" \
  || { echo "Error: tabix failed for $FINAL_VCF_GZ"; exit 1; }

echo "Done:"
echo "  final      : $FINAL_VCF_GZ"

