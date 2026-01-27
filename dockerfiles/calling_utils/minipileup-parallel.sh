#!/bin/bash
set -euo pipefail

#################################################################################
# Parallel minipileup on VCF intervals
# Requires:
#   - bcftools
#   - tabix
#   - minipileup
#   - samtools
#   - bgzip
#   - python3 with granite-suite package
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -r reference.fasta [-o prefix] \\
          [--sr-cram CRAM ...] [--pb-cram CRAM ...] [--ont-cram CRAM ...] \\
          [-t threads] [--args "additional minipileup args"] [--group INT]

  -i             Input VCF (bgzipped) with .tbi index (required)
  -r             Reference FASTA with .fai index (required)
  -o             Output prefix (default: output) -> produces <prefix>.vcf.gz
  -t             Parallel jobs/threads for xargs (default: number of CPUs)
  --args         Additional arguments to pass to minipileup (in quotes) (default: "-c -C -Q 20 -q 30 -s 0")

  --sr-cram      Short-read CRAM with .crai file, also accept BAM with .bai or .csi index (repeatable)
  --pb-cram      PacBio long-read CRAM with .crai file, also accept BAM with .bai or .csi index (repeatable)
  --pb-tissue    Tissue ID matching each --pb-cram (repeatable) (eg SMHT005-3AF)
  --ont-cram     ONT long-read CRAM with .crai file, also accept BAM with .bai or .csi index (repeatable)
  --ont-tissue   Tissue ID matching each --ont-cram (repeatable) (eg SMHT005-3AF)

  --group        Group intervals into batches of INT for processing (default: 100)
EOF
  exit 1
}

# Defaults
INPUT_VCF=""
REFERENCE_FASTA=""
OUTPUT_PRFX="output"
THREADS="$(nproc)"
MINPILEUP_ARGS="-c -C -Q 20 -q 30 -s 0"
GROUP=100

# mapping quality (-q)
# base quality (-Q)
# count alleles both strands (-C)
# vcf format (-c)
# drop alleles with depth<INT (-s)

SR_CRAMS=()
PB_CRAMS=()
PB_TISSUES=()
ONT_CRAMS=()
ONT_TISSUES=()


# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_VCF="$2"; shift 2;;
    -r) REFERENCE_FASTA="$2"; shift 2;;
    -o) OUTPUT_PRFX="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --args) MINPILEUP_ARGS="$2"; shift 2;;

    --sr-cram) SR_CRAMS+=("$2"); shift 2;;
    --pb-cram) PB_CRAMS+=("$2"); shift 2;;
    --pb-tissue) PB_TISSUES+=("$2"); shift 2;;
    --ont-cram) ONT_CRAMS+=("$2"); shift 2;;
    --ont-tissue) ONT_TISSUES+=("$2"); shift 2;;


    --group) GROUP="$2"; shift 2;;

    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

# Required checks
[[ -n "$INPUT_VCF" ]] || { echo "Error: input VCF (-i) is required"; usage; }
[[ -n "$REFERENCE_FASTA" ]] || { echo "Error: reference FASTA (-r) is required"; usage; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "Error: -t must be an integer"; exit 1; }
[[ "$GROUP" =~ ^[0-9]+$ ]] || { echo "Error: --group must be an integer"; exit 1; }
(( GROUP > 0 )) || { echo "Error: --group must be > 0"; exit 1; } 

# Input checks
[[ -f "$INPUT_VCF" ]] || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]] || { echo "Error: ${INPUT_VCF}.tbi index not found"; exit 1; }
[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

# CRAM/BAM checks (index check too)
ALL_CRAMS=("${SR_CRAMS[@]}" "${PB_CRAMS[@]}" "${ONT_CRAMS[@]}")
if (( ${#ALL_CRAMS[@]} == 0 )); then
  echo "Error: no BAM/CRAM files provided. Please specify at least one with --sr-cram/--pb-cram/--ont-cram"
  exit 1
else
  for b in "${ALL_CRAMS[@]}"; do
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

(( ${#PB_CRAMS[@]} > 0 )) || { echo "Error: at least one --pb-cram required"; exit 1; }
(( ${#PB_CRAMS[@]} == ${#PB_TISSUES[@]} )) || {
  echo "Error: number of --pb-cram entries must equal number of --pb-tissue entries"; exit 1;
}

if (( ${#ONT_CRAMS[@]} > 0 )); then
  (( ${#ONT_CRAMS[@]} == ${#ONT_TISSUES[@]} )) || {
    echo "Error: number of --ont-cram entries must equal number of --ont-tissue entries"; exit 1;
  }
fi

# Tool checks
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }
command -v tabix >/dev/null 2>&1     || { echo "Error: tabix not found in PATH"; exit 1; }
command -v minipileup >/dev/null 2>&1|| { echo "Error: minipileup not found in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1  || { echo "Error: samtools not found in PATH"; exit 1; }
command -v bgzip   >/dev/null 2>&1 || { echo "Error: bgzip not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "Error: python3 not found in PATH"; exit 1; }
python3 - <<'PY' || { echo "Error: granite Python package not importable"; exit 1; }
try:
    from granite.lib import vcf_parser  # noqa
except Exception as e:
    import sys
    sys.stderr.write(f"granite import failed: {e}\n")
    sys.exit(1)
PY

# Workdirs
WORKDIR="$(mktemp -d minipileup_shards.XXXXXX)"
cleanup() {
  rm -rf "$WORKDIR"
}
trap cleanup EXIT

# File paths
INTERVALS="${WORKDIR}/intervals.txt"
HEADER_VCF="${WORKDIR}/header.vcf"

# Extract intervals from VCF
echo "Extracting intervals to file..."
bcftools query -f '%CHROM:%POS-%END\n' "$INPUT_VCF" > "$INTERVALS"

# Group intervals into batches (to reduce overhead)
GROUPED="${WORKDIR}/intervals_grouped.txt"
awk -v n="$GROUP" '
{
  if (NR % n == 1) line = $0;
  else             line = line "\t" $0;
  if (NR % n == 0) { print line; line = "" }
}
END { if (line != "") print line }
' "$INTERVALS" > "$GROUPED"


# Create run function
run_region() {
  local REGION_LINE="$1" # one line containing ~100 "chr:start-end" tokens
  local HEADER="$2" # 1 = write header, 0 = only body
  shift 2
  local CRAMS=( "$@" )

  # Split the line into an array of regions
  # (word-splitting is intentional here because REGION_LINE is TAB-delimited)
  local OLDIFS=$IFS
  IFS=$'\t' read -r -a REG_ARR <<< "$REGION_LINE"
  IFS=$OLDIFS

  # Use the first region to derive a safe ID for temp files
  local key="${REG_ARR[0]}"
  local safe="${key//:/_}"; safe="${safe//-/__}"

  # If HEADER==1, grab the header once from the first region
  if [[ "$HEADER" -eq 1 ]]; then
    # Extract BAMs for that first region (just for header generation)
    local HBAMS=()
    for cram in "${CRAMS[@]}"; do
      local bam="${WORKDIR}/$(basename "${cram%.*}")_${safe}.bam"
      samtools view --reference "$REFERENCE_FASTA" --write-index -b -o "$bam" "$cram" chr1:10001-10001
      HBAMS+=("$bam")
    done
    minipileup -f "$REFERENCE_FASTA" \
      $MINPILEUP_ARGS \
      -r "$key" \
      "${HBAMS[@]}" | grep '^#' > "$HEADER_VCF"
    # Cleanup header-only BAMs
    for b in "${HBAMS[@]}"; do rm -f "$b" "${b}.csi"; done

    return 0 # avoid processing the body for the header run
  fi

  # For the group body: loop regions; for each region build BAMs, run minipileup, append
  local group_vcf="${WORKDIR}/${safe}.group.vcf"
  : > "$group_vcf"

  for region in "${REG_ARR[@]}"; do
    local rsafe="${region//:/_}"; rsafe="${rsafe//-/__}"

    # Extract region to BAM files
    local BAMS=()
    for cram in "${CRAMS[@]}"; do
      local bam="${WORKDIR}/$(basename "${cram%.*}")_${rsafe}.bam"
      samtools view --reference "$REFERENCE_FASTA" --write-index -b -o "$bam" "$cram" "$region"
      BAMS+=("$bam")
    done

    # Run minipileup for this region (avoid piping so segfault can't break the shell pipeline)
    tmp_out="${WORKDIR}/minipileup_${rsafe}.$$.$RANDOM.out"

    minipileup -f "$REFERENCE_FASTA" \
      $MINPILEUP_ARGS \
      -r "$region" \
      "${BAMS[@]}" > "$tmp_out"
    mp_status=$?

    if [[ $mp_status -eq 139 ]]; then
      echo "Seg fault at $region" >&2
      rm -f "$tmp_out"
      for b in "${BAMS[@]}"; do rm -f "$b" "${b}.csi"; done
      continue
    elif [[ $mp_status -ne 0 ]]; then
      # keep prior behavior: don't kill the overall job on other minipileup failures
      rm -f "$tmp_out"
      true
    else
      grep -v '^#' "$tmp_out" >> "$group_vcf" || true
      rm -f "$tmp_out"
    fi

    ## Run minipileup for this region
    #minipileup -f "$REFERENCE_FASTA" \
    #  $MINPILEUP_ARGS \
    #  -r "$region" \
    #  "${BAMS[@]}" | grep -v '^#' >> "$group_vcf" || true

    # Cleanup per-region BAMs
    for b in "${BAMS[@]}"; do rm -f "$b" "${b}.csi"; done
  done
}

# Make the function available to xargs bash -c
export -f run_region
export WORKDIR REFERENCE_FASTA MINPILEUP_ARGS HEADER_VCF

export MERGED_VCF="${WORKDIR}/merged.vcf"
export SORTED_VCF="${WORKDIR}/sorted.vcf"

# Safety check
if [[ ! -s "$GROUPED" ]]; then
    echo "Error: grouped intervals file is empty"
    run_region chr1:10001-10001 1 "${ALL_CRAMS[@]}"

    cat "$HEADER_VCF" > "$MERGED_VCF"
else
    # Create temporary VCF file with header only
    echo "Creating temporary VCF with header..."
    read -r first_line < "$GROUPED"
    run_region "$first_line" 1 "${ALL_CRAMS[@]}"
    
    [[ -s "$HEADER_VCF" ]] || { echo "Error: header not generated"; exit 1; }
    
    # Parallel execution
    echo "Running minipileup in parallel on intervals..."
    grep -Ev '^[[:space:]]*($|#)' "$GROUPED" | \
    xargs -P "$THREADS" -I{} bash -c 'run_region "$1" "$2" "${@:3}"' _ \
      "{}" 0 "${ALL_CRAMS[@]}" || { echo "Error: parallel minipileup execution failed"; exit 1; }
    
    # Collect shard outputs (sorted, safe if none)
    shopt -s nullglob
    mapfile -t SHARDS < <(printf "%s\n" "$WORKDIR"/*.group.vcf)
    shopt -u nullglob
    (( ${#SHARDS[@]} > 0 )) || { echo "Error: no grouped outputs found in $WORKDIR"; exit 1; }
    
    echo "Merging ${#SHARDS[@]} annotated groups..."
    
    # Write header once
    cat "$HEADER_VCF" > "$MERGED_VCF"
    # Append all shard bodies
    cat "${SHARDS[@]}" >> "$MERGED_VCF"

fi

echo "-- BCFTools ---------------------------"
bcftools sort -T "tmp_bcftools.XXXXXX" -O v -o "$SORTED_VCF" "$MERGED_VCF" \
  || { echo "Error: bcftools sort failed"; exit 1; }
echo "-----------------------------------------"

[[ -s "$SORTED_VCF" ]] || { echo "Error: sorted VCF is empty"; exit 1; }

# Write files to rename sample and add file type information (+ PB tissue metadata)
: > "${WORKDIR}/map.txt"

for cram in "${SR_CRAMS[@]}"; do
  printf "%s\tSR\t.\n" "$(basename "${cram%.*}")" >> "${WORKDIR}/map.txt"
done

for i in "${!PB_CRAMS[@]}"; do
  cram="${PB_CRAMS[$i]}"
  tissue="${PB_TISSUES[$i]}"
  printf "%s\tPB\t%s\n" "$(basename "${cram%.*}")" "$tissue" >> "${WORKDIR}/map.txt"
done


for i in "${!ONT_CRAMS[@]}"; do
  cram="${ONT_CRAMS[$i]}"
  tissue="${ONT_TISSUES[$i]}"
  printf "%s\tONT\t%s\n" "$(basename "${cram%.*}")" "$tissue" >> "${WORKDIR}/map.txt"
done



# Granite to rename the sample and add file type information
py_script="
from granite.lib import vcf_parser
import os, re, sys

# Load mapping: sample_basename -> (SR/PB/ONT, pb_tissue or '')
sample_map = {}
with open('${WORKDIR}/map.txt', 'r') as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        sample, ftype, tissue = parts[0], parts[1], parts[2]
        if tissue == '.':
            sample_map[sample] = ftype
        else:
            sample_map[sample] = '-'.join([ftype, tissue])

# Load the VCF
vcf_obj = vcf_parser.Vcf('${SORTED_VCF}')

# Current header columns (e.g. '#CHROM ... FORMAT  path/sample1_chr1_123__456  path/sample2_chr7_...')
cols = vcf_obj.header.columns.rstrip('\n').split('\t')
fixed = cols[:9]

for col in cols[9:]:
    base = os.path.basename(col)

    # Try to remove the region suffix we introduced ('_<something derived from region>').
    # 1) If suffix starts with '_chr...' (common case)
    m = re.match(r'^(.*)_chr.*$', base)
    if m:
        stem = m.group(1)
    else:
        # 2) Otherwise, try to strip the last _<token> that looks like a region blob
        #    (contains digits/underscores), keeping it conservative.
        m2 = re.match(r'^(.*)_[0-9A-Za-z_]+$', base)
        stem = m2.group(1) if m2 else base

    ftype = sample_map[stem]
    fixed.append(f'{stem}-{ftype}')

vcf_obj.header.columns = '\t'.join(fixed) + '\n'

with open('${OUTPUT_PRFX}.vcf', 'w') as fo:
    vcf_obj.write_header(fo)
    for vnt_obj in vcf_obj.parse_variants():
        vcf_obj.write_variant(fo, vnt_obj)
"

echo "Renaming samples and adding file type information..."
python3 -c "$py_script" || { echo "Error: Renaming script failed"; exit 1; }

# Compress and index final VCF
echo "Compressing and indexing final VCF..."
bgzip -f "${OUTPUT_PRFX}.vcf" || { echo "Error: bgzip failed"; exit 1; }
tabix -f -p vcf "${OUTPUT_PRFX}.vcf.gz" || { echo "Error: tabix indexing failed"; exit 1; }

echo "Final output: ${OUTPUT_PRFX}.vcf.gz"
