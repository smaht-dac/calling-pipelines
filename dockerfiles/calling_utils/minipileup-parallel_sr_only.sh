#!/bin/bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -r reference.fasta [-o prefix] \\
          --sr-cram FILE --sr-tissue ID [repeatable] \\
          [-t threads] [--args "minipileup args"] [--group INT]

Required:
  -i             Input VCF (bgzipped, indexed)
  -r             Reference FASTA (.fai indexed)
  --sr-cram      Short-read CRAM/BAM file (repeatable)
  --sr-tissue    Tissue ID matching each --sr-cram (repeatable)

Optional:
  -o             Output prefix (default: output)
  -t             Threads for parallel xargs (default: nproc)
  --args         Extra minipileup arguments (default: "-c -C -Q 30 -q 30 -s 0")
  --group        Number of VCF intervals per batch (default: 100)
EOF
  exit 1
}

INPUT_VCF=""
REFERENCE_FASTA=""
OUTPUT_PRFX="output"
THREADS="$(nproc)"
MINPILEUP_ARGS="-c -C -Q 30 -q 30 -s 0"
GROUP=100

SR_CRAMS=()
SR_TISSUES=()

# ---------------- ARG PARSING --------------------

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_VCF="$2"; shift 2;;
    -r) REFERENCE_FASTA="$2"; shift 2;;
    -o) OUTPUT_PRFX="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --args) MINPILEUP_ARGS="$2"; shift 2;;

    --sr-cram) SR_CRAMS+=("$2"); shift 2;;
    --sr-tissue) SR_TISSUES+=("$2"); shift 2;;

    --group) GROUP="$2"; shift 2;;

    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

# ---------------- VALIDATION --------------------

[[ -n "$INPUT_VCF" ]] || { echo "Error: -i required"; exit 1; }
[[ -n "$REFERENCE_FASTA" ]] || { echo "Error: -r required"; exit 1; }

(( ${#SR_CRAMS[@]} > 0 )) || { echo "Error: at least one --sr-cram required"; exit 1; }
(( ${#SR_CRAMS[@]} == ${#SR_TISSUES[@]} )) || {
  echo "Error: number of --sr-cram entries must equal number of --sr-tissue entries"; exit 1;
}

# Index checks
[[ -f "$INPUT_VCF" ]] || exit 1
[[ -f "${INPUT_VCF}.tbi" ]] || exit 1
[[ -f "$REFERENCE_FASTA" ]] || exit 1
[[ -f "${REFERENCE_FASTA}.fai" ]] || exit 1

for b in "${SR_CRAMS[@]}"; do
  [[ -f "$b" ]] || { echo "Missing CRAM/BAM: $b"; exit 1; }
  if [[ "$b" == *.cram && ! -f "${b}.crai" ]]; then
      echo "Missing .crai for $b"; exit 1;
  fi
  if [[ "$b" == *.bam && ! -f "${b}.bai" && ! -f "${b}.csi" ]]; then
      echo "Missing BAM index for $b"; exit 1;
  fi
done

command -v bcftools >/dev/null || exit 1
command -v minipileup >/dev/null || exit 1
command -v samtools >/dev/null || exit 1
command -v python3 >/dev/null || exit 1
command -v bgzip >/dev/null || exit 1
command -v tabix >/dev/null || exit 1

python3 - <<'PY' || { echo "Missing granite suite"; exit 1; }
from granite.lib import vcf_parser
PY

# ---------------- WORK DIRECTORIES --------------------

WORKDIR="$(mktemp -d minipileup_shards.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

INTERVALS="$WORKDIR/intervals.txt"
HEADER_VCF="$WORKDIR/header.vcf"

echo "Extracting intervals..."
bcftools query -f '%CHROM:%POS-%END\n' "$INPUT_VCF" > "$INTERVALS"

GROUPED="$WORKDIR/intervals_grouped.txt"
awk -v n="$GROUP" '
{
  if (NR % n == 1) line = $0;
  else             line = line "\t" $0;
  if (NR % n == 0) { print line; line = "" }
}
END { if (line != "") print line }
' "$INTERVALS" > "$GROUPED"

[[ -s "$GROUPED" ]] || { echo "Error: no intervals"; exit 1; }

# ---------------- REGION HANDLER --------------------

run_region() {
  local REGION_LINE="$1"
  local HEADER="$2"
  shift 2
  local CRAMS=( "$@" )

  IFS=$'\t' read -r -a REG_ARR <<< "$REGION_LINE"

  local first="${REG_ARR[0]}"
  local safe="${first//:/_}"; safe="${safe//-/__}"

  if [[ "$HEADER" -eq 1 ]]; then
    local HBAMS=()
    for cram in "${CRAMS[@]}"; do
      local bam="$WORKDIR/$(basename "${cram%.*}")_${safe}.bam"
      samtools view --reference "$REFERENCE_FASTA" --write-index -b -o "$bam" "$cram" chr1:10001-10001
      HBAMS+=("$bam")
    done

    minipileup -f "$REFERENCE_FASTA" \
      $MINPILEUP_ARGS \
      -r "$first" \
      "${HBAMS[@]}" | grep '^#' > "$HEADER_VCF"

	for b in "${HBAMS[@]}"; do rm -f "$b" "${b}.csi"; done

    return 0
  fi

  local outvcf="$WORKDIR/${safe}.group.vcf"
  : > "$outvcf"

  for region in "${REG_ARR[@]}"; do
    local rsafe="${region//:/_}"; rsafe="${rsafe//-/__}"
    local BAMS=()

    for cram in "${CRAMS[@]}"; do
      local bam="$WORKDIR/$(basename "${cram%.*}")_${rsafe}.bam"
      samtools view --reference "$REFERENCE_FASTA" --write-index -b -o "$bam" "$cram" "$region"
      BAMS+=("$bam")
    done

    minipileup -f "$REFERENCE_FASTA" \
       $MINPILEUP_ARGS \
       -r "$region" \
       "${BAMS[@]}" | grep -v '^#' >> "$outvcf" || true

	for b in "${BAMS[@]}"; do rm -f "$b" "${b}.csi"; done
  done
}

export -f run_region
export WORKDIR REFERENCE_FASTA MINPILEUP_ARGS HEADER_VCF

# ---------------- RUN HEADER --------------------

read -r firstline < "$GROUPED"
run_region "$firstline" 1 "${SR_CRAMS[@]}"

# ---------------- RUN IN PARALLEL --------------------

echo "Running minipileup..."
grep -Ev '^\s*$' "$GROUPED" | \
xargs -P "$THREADS" -I{} bash -c 'run_region "$1" 0 "${@:2}"' _ \
  "{}" "${SR_CRAMS[@]}"

# ---------------- MERGING --------------------

SHARDS=( "$WORKDIR"/*.group.vcf )
(( ${#SHARDS[@]} > 0 )) || { echo "No shard outputs"; exit 1; }

MERGED="$WORKDIR/merged.vcf"
SORTED="$WORKDIR/sorted.vcf"

cat "$HEADER_VCF" > "$MERGED"
cat "${SHARDS[@]}" >> "$MERGED"

bcftools sort -T "$WORKDIR/tmp" -O v -o "$SORTED" "$MERGED"

# ---------------- SAMPLE RENAMING BY TISSUE --------------------

echo "Building map for sample→tissue"
MAPFILE="$WORKDIR/map.txt"
: > "$MAPFILE"

for i in "${!SR_CRAMS[@]}"; do
    stem="$(basename "${SR_CRAMS[$i]%.*}")"
    tissue="${SR_TISSUES[$i]}"
    echo -e "${stem}\t${tissue}" >> "$MAPFILE"
done

py_script="
from granite.lib import vcf_parser
import os, re

# stem -> tissue
tmap = {}
with open('${MAPFILE}', 'r') as f:
    for line in f:
        s, t = line.strip().split('\t')
        tmap[s] = t

vcf = vcf_parser.Vcf('${SORTED}')
cols = vcf.header.columns.rstrip('\n').split('\t')
fixed = cols[:9]

for c in cols[9:]:
    base = os.path.basename(c)
    # Remove region suffix
    m = re.match(r'^(.*)_chr.*$', base)
    stem = m.group(1) if m else base.split('_')[0]

    tissue = tmap.get(stem, 'UNK')
    fixed.append(f'{stem}-{tissue}')

vcf.header.columns = '\t'.join(fixed) + '\n'

with open('${OUTPUT_PRFX}.vcf', 'w') as fo:
    vcf.write_header(fo)
    for v in vcf.parse_variants():
        vcf.write_variant(fo, v)
"

python3 -c "$py_script"

bgzip -f "${OUTPUT_PRFX}.vcf"
tabix -f -p vcf "${OUTPUT_PRFX}.vcf.gz"

echo "Done. Output: ${OUTPUT_PRFX}.vcf.gz"

