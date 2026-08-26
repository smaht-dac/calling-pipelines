#!/bin/bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -r reference.fasta [-o prefix] \\
          --sr-cram CRAM --sr-tissue ID [repeatable]

Required:
  -i             Input VCF (bgzipped, indexed)
  -r             Reference FASTA (.fai indexed)
  --sr-cram      Short-read CRAM/BAM file (repeatable)
  --sr-tissue    Tissue ID matching each --sr-cram (repeatable)

Optional:
  -o             Output prefix (default: output)
  --args         Additional arguments to pass to minipileup2 (in quotes)
                 (default: "-D -c -C -Q 30 -q 30 -s 0")
EOF
  exit 1
}

INPUT_VCF=""
REFERENCE_FASTA=""
OUTPUT_PRFX="output"
MINPILEUP_ARGS="-D -c -C -Q 30 -q 30 -s 0"

SR_CRAMS=()
SR_TISSUES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_VCF="$2"; shift 2;;
    -r) REFERENCE_FASTA="$2"; shift 2;;
    -o) OUTPUT_PRFX="$2"; shift 2;;
    --sr-cram)   SR_CRAMS+=("$2"); shift 2;;
    --sr-tissue) SR_TISSUES+=("$2"); shift 2;;
    --args)      MINPILEUP_ARGS="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

[[ -n "$INPUT_VCF" ]]        || { echo "Error: -i required"; exit 1; }
[[ -n "$REFERENCE_FASTA" ]]  || { echo "Error: -r required"; exit 1; }

(( ${#SR_CRAMS[@]} > 0 )) || { echo "Error: at least one --sr-cram required"; exit 1; }
(( ${#SR_CRAMS[@]} == ${#SR_TISSUES[@]} )) || {
  echo "Error: number of --sr-cram entries must equal number of --sr-tissue entries"; exit 1;
}

[[ -f "$INPUT_VCF" ]]             || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]]       || { echo "Error: ${INPUT_VCF}.tbi not found"; exit 1; }
[[ -f "$REFERENCE_FASTA" ]]       || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

for b in "${SR_CRAMS[@]}"; do
  [[ -f "$b" ]] || { echo "Error: missing CRAM/BAM: $b"; exit 1; }
  if [[ "$b" == *.cram ]]; then
    [[ -f "${b}.crai" ]] || { echo "Error: missing .crai for $b"; exit 1; }
  elif [[ "$b" == *.bam ]]; then
    [[ -f "${b}.bai" || -f "${b}.csi" ]] || { echo "Error: missing BAM index for $b"; exit 1; }
  else
    echo "Error: unsupported file type (must be .bam or .cram): $b"; exit 1
  fi
done

command -v minipileup2 >/dev/null 2>&1 || { echo "Error: minipileup2 not in PATH"; exit 1; }
command -v bcftools    >/dev/null 2>&1 || { echo "Error: bcftools not in PATH"; exit 1; }
command -v bgzip       >/dev/null 2>&1 || { echo "Error: bgzip not in PATH"; exit 1; }
command -v tabix       >/dev/null 2>&1 || { echo "Error: tabix not in PATH"; exit 1; }
command -v python3     >/dev/null 2>&1 || { echo "Error: python3 not in PATH"; exit 1; }
python3 - <<'PY' || { echo "Error: granite not importable"; exit 1; }
from granite.lib import vcf_parser
PY

WORKDIR="$(mktemp -d minipileup2_sr_work.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

echo "Running minipileup2 (SR only)..."
minipileup2 -t 2 -f "$REFERENCE_FASTA" -x "$INPUT_VCF" \
  $MINPILEUP_ARGS \
  "${SR_CRAMS[@]}" > "$WORKDIR/merged.vcf" \
  || { echo "Error: minipileup2 failed"; exit 1; }

[[ -s "$WORKDIR/merged.vcf" ]] || { echo "Error: minipileup2 produced empty output"; exit 1; }

bcftools sort -T "$WORKDIR/tmp" -O v -o "$WORKDIR/sorted.vcf" "$WORKDIR/merged.vcf" \
  || { echo "Error: bcftools sort failed"; exit 1; }

# Build sample map: stem -> tissue
MAPFILE="$WORKDIR/map.txt"
: > "$MAPFILE"
for i in "${!SR_CRAMS[@]}"; do
  printf "%s\t%s\n" "$(basename "${SR_CRAMS[$i]%.*}")" "${SR_TISSUES[$i]}" >> "$MAPFILE"
done

py_script="
from granite.lib import vcf_parser
import os, re

tmap = {}
with open('$MAPFILE') as f:
    for line in f:
        stem, tissue = line.strip().split('\t')
        tmap[stem] = tissue

vcf = vcf_parser.Vcf('$WORKDIR/sorted.vcf')
cols = vcf.header.columns.rstrip('\n').split('\t')
fixed = cols[:9]

for c in cols[9:]:
    stem = re.sub(r'\.(cram|bam)$', '', os.path.basename(c))
    tissue = tmap.get(stem, 'UNK')
    fixed.append(f'{stem}-{tissue}')

vcf.header.columns = '\t'.join(fixed) + '\n'

with open('${OUTPUT_PRFX}.vcf', 'w') as fo:
    vcf.write_header(fo)
    for v in vcf.parse_variants():
        vcf.write_variant(fo, v)
"

echo "Renaming samples..."
python3 -c "$py_script" || { echo "Error: sample renaming failed"; exit 1; }

bgzip -f "${OUTPUT_PRFX}.vcf"            || { echo "Error: bgzip failed"; exit 1; }
tabix -f -p vcf "${OUTPUT_PRFX}.vcf.gz"  || { echo "Error: tabix failed"; exit 1; }

echo "Done: ${OUTPUT_PRFX}.vcf.gz"
