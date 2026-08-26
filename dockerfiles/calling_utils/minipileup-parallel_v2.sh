#!/bin/bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -r reference.fasta [-o prefix] \\
          [--sr-cram CRAM ... --sr-core CORE ...] \\
          [--lr-cram CRAM ... --lr-core CORE ... --lr-tissue TISSUE --lr-type PB|ONT ...] \\
          [--args "additional minipileup2 args"]

  -i             Input VCF (bgzipped) with .tbi index (required)
  -r             Reference FASTA with .fai index (required)
  -o             Output prefix (default: output)

  --sr-cram      Short-read CRAM with .crai index (repeatable)
  --sr-core      Sequencing core ID matching each --sr-cram (repeatable)
  --lr-cram      Long-read CRAM (repeatable)
  --lr-core      Sequencing core ID matching each --lr-cram (repeatable)
  --lr-tissue    Tissue ID matching each --lr-cram (repeatable)
  --lr-type      Sequencing type matching each --lr-cram (repeatable, PB or ONT)
  --args         Additional arguments to pass to minipileup2 (in quotes)
                 (default: "-D -c -C -Q 20 -q 30 -s 0")

Sample columns are named:

  with cores:     {core}__{cram_stem}-SR
                  {core}__{cram_stem}-PB-{tissue}
                  {core}__{cram_stem}-ONT-{tissue}
  without cores:  {cram_stem}-SR, {cram_stem}-PB-{tissue}, ...

Core IDs are matched to their CRAMs BY POSITION, so --sr-core must be given in the
same relative order as --sr-cram (likewise --lr-core with --lr-cram). Supply cores
for every long-read CRAM, including CRAMs from other tissues -- the donor-pooled
long-read set is passed through whole and tier_filter_variants_SR_PB_ONT.py selects
the current tissue's cores itself. Omit both core options entirely to get the
original unprefixed column names.
EOF
  exit 1
}

INPUT_VCF=""
REFERENCE_FASTA=""
OUTPUT_PRFX="output"
MINPILEUP_ARGS="-D -c -C -Q 20 -q 30 -s 0"

SR_CRAMS=()
SR_CORES=()
LR_CRAMS=()
LR_CORES=()
LR_TISSUES=()
LR_TYPES=()
PB_CRAMS=()
PB_CORES=()
PB_TISSUES=()
ONT_CRAMS=()
ONT_CORES=()
ONT_TISSUES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_VCF="$2"; shift 2;;
    -r) REFERENCE_FASTA="$2"; shift 2;;
    -o) OUTPUT_PRFX="$2"; shift 2;;
    --sr-cram)   SR_CRAMS+=("$2"); shift 2;;
    --sr-core)   SR_CORES+=("$2"); shift 2;;
    --lr-cram)   LR_CRAMS+=("$2"); shift 2;;
    --lr-core)   LR_CORES+=("$2"); shift 2;;
    --lr-tissue) LR_TISSUES+=("$2"); shift 2;;
    --lr-type)   LR_TYPES+=("$2"); shift 2;;
    --args)      MINPILEUP_ARGS="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

[[ -n "$INPUT_VCF" ]]        || { echo "Error: -i required"; usage; }
[[ -n "$REFERENCE_FASTA" ]]  || { echo "Error: -r required"; usage; }

[[ -f "$INPUT_VCF" ]]             || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]]       || { echo "Error: ${INPUT_VCF}.tbi not found"; exit 1; }
[[ -f "$REFERENCE_FASTA" ]]       || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }

(( ${#LR_CRAMS[@]} > 0 )) || { echo "Error: at least one --lr-cram required"; exit 1; }
(( ${#LR_CRAMS[@]} == ${#LR_TISSUES[@]} )) || { echo "Error: --lr-cram/--lr-tissue count mismatch"; exit 1; }
(( ${#LR_CRAMS[@]} == ${#LR_TYPES[@]} ))   || { echo "Error: --lr-cram/--lr-type count mismatch"; exit 1; }

# Core IDs are optional. Supplying any core switches on prefixed sample columns and
# then EVERY CRAM needs one, matched by position -- a short array would silently
# shift cores onto the wrong CRAMs.
USE_CORES=0
if (( ${#SR_CORES[@]} > 0 || ${#LR_CORES[@]} > 0 )); then
  USE_CORES=1
  (( ${#SR_CRAMS[@]} == ${#SR_CORES[@]} )) || {
    echo "Error: --sr-cram was given ${#SR_CRAMS[@]} time(s) but --sr-core was given ${#SR_CORES[@]} time(s)."
    echo "       Each --sr-cram needs exactly one --sr-core, matched by position."
    exit 1
  }
  (( ${#LR_CRAMS[@]} == ${#LR_CORES[@]} )) || {
    echo "Error: --lr-cram was given ${#LR_CRAMS[@]} time(s) but --lr-core was given ${#LR_CORES[@]} time(s)."
    echo "       Each --lr-cram needs exactly one --lr-core, matched by position."
    exit 1
  }
fi

for i in "${!LR_CRAMS[@]}"; do
  t_up="$(printf "%s" "${LR_TYPES[$i]}" | tr '[:lower:]' '[:upper:]')"
  lr_core="."
  (( USE_CORES )) && lr_core="${LR_CORES[$i]}"
  case "$t_up" in
    PB)  PB_CRAMS+=("${LR_CRAMS[$i]}");  PB_TISSUES+=("${LR_TISSUES[$i]}");  PB_CORES+=("$lr_core");;
    ONT) ONT_CRAMS+=("${LR_CRAMS[$i]}"); ONT_TISSUES+=("${LR_TISSUES[$i]}"); ONT_CORES+=("$lr_core");;
    *) echo "Error: invalid --lr-type '${LR_TYPES[$i]}' for '${LR_CRAMS[$i]}'. Expected PB or ONT."; exit 1;;
  esac
done

# Misalignment check: most CRAM basenames are {tissue}-{core}-{...}, so a core that is
# not field 3 of its own basename suggests the positional arrays have drifted apart.
#
# Only a heuristic, so it only warns, and only when the convention is demonstrably in
# use for this run -- some donors key cores by accession instead (e.g. SMHT005-3A's
# cores are SMAANBIDMBAC/SMAANTB3AIM2 while every basename carries a literal 'XX'), and
# there a per-CRAM warning would fire on every column while telling us nothing.
if (( USE_CORES )); then
  core_field3() { printf "%s" "$(basename "${1%.*}")" | cut -d- -f3; }

  matches=0; total=0
  for i in "${!SR_CRAMS[@]}"; do
    [[ "${SR_CORES[$i]}" == "." ]] && continue
    total=$((total + 1))
    [[ "$(core_field3 "${SR_CRAMS[$i]}")" == "${SR_CORES[$i]}" ]] && matches=$((matches + 1))
  done
  for i in "${!LR_CRAMS[@]}"; do
    [[ "${LR_CORES[$i]}" == "." ]] && continue
    total=$((total + 1))
    [[ "$(core_field3 "${LR_CRAMS[$i]}")" == "${LR_CORES[$i]}" ]] && matches=$((matches + 1))
  done

  if (( matches > 0 && matches < total )); then
    echo "Warning: $((total - matches))/$total core IDs do not match field 3 of their own CRAM basename," >&2
    echo "         while the rest do. Check that each core array is in the same order as its CRAM array:" >&2
    for i in "${!SR_CRAMS[@]}"; do
      [[ "${SR_CORES[$i]}" == "." ]] && continue
      f3="$(core_field3 "${SR_CRAMS[$i]}")"
      [[ "$f3" != "${SR_CORES[$i]}" ]] && \
        echo "           --sr-core '${SR_CORES[$i]}' vs basename field 3 '$f3' ($(basename "${SR_CRAMS[$i]}"))" >&2
    done
    for i in "${!LR_CRAMS[@]}"; do
      [[ "${LR_CORES[$i]}" == "." ]] && continue
      f3="$(core_field3 "${LR_CRAMS[$i]}")"
      [[ "$f3" != "${LR_CORES[$i]}" ]] && \
        echo "           --lr-core '${LR_CORES[$i]}' vs basename field 3 '$f3' ($(basename "${LR_CRAMS[$i]}"))" >&2
    done
  fi
fi

(( ${#PB_CRAMS[@]} > 0 )) || { echo "Error: at least one PB long-read (--lr-type PB) required"; exit 1; }

ALL_CRAMS=("${SR_CRAMS[@]}" "${PB_CRAMS[@]}" "${ONT_CRAMS[@]}")
(( ${#ALL_CRAMS[@]} > 0 )) || { echo "Error: no CRAM/BAM files provided"; exit 1; }

for b in "${ALL_CRAMS[@]}"; do
  [[ -f "$b" ]] || { echo "Error: file not found: $b"; exit 1; }
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

WORKDIR="$(mktemp -d minipileup2_work.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

echo "Running minipileup2..."
minipileup2 -t 2 -f "$REFERENCE_FASTA" -x "$INPUT_VCF" \
  $MINPILEUP_ARGS \
  "${ALL_CRAMS[@]}" > "$WORKDIR/merged.vcf" \
  || { echo "Error: minipileup2 failed"; exit 1; }

[[ -s "$WORKDIR/merged.vcf" ]] || { echo "Error: minipileup2 produced empty output"; exit 1; }

bcftools sort -T "tmp_bcftools.XXXXXX" -O v -o "$WORKDIR/sorted.vcf" "$WORKDIR/merged.vcf" \
  || { echo "Error: bcftools sort failed"; exit 1; }

# Build sample map: stem -> type, tissue, core ('.' when absent)
: > "$WORKDIR/map.txt"
for i in "${!SR_CRAMS[@]}"; do
  sr_core="."
  (( USE_CORES )) && sr_core="${SR_CORES[$i]}"
  printf "%s\tSR\t.\t%s\n" "$(basename "${SR_CRAMS[$i]%.*}")" "$sr_core" >> "$WORKDIR/map.txt"
done
for i in "${!PB_CRAMS[@]}"; do
  printf "%s\tPB\t%s\t%s\n" "$(basename "${PB_CRAMS[$i]%.*}")" "${PB_TISSUES[$i]}" "${PB_CORES[$i]}" >> "$WORKDIR/map.txt"
done
for i in "${!ONT_CRAMS[@]}"; do
  printf "%s\tONT\t%s\t%s\n" "$(basename "${ONT_CRAMS[$i]%.*}")" "${ONT_TISSUES[$i]}" "${ONT_CORES[$i]}" >> "$WORKDIR/map.txt"
done

py_script="
from granite.lib import vcf_parser
import os, re, sys

use_cores = bool(int('$USE_CORES'))

sample_map = {}
with open('$WORKDIR/map.txt') as f:
    for line in f:
        stem, ftype, tissue, core = line.rstrip('\n').split('\t')
        sample_map[stem] = (ftype if tissue == '.' else f'{ftype}-{tissue}', core)

vcf_obj = vcf_parser.Vcf('$WORKDIR/sorted.vcf')
cols = vcf_obj.header.columns.rstrip('\n').split('\t')
fixed = cols[:9]

for col in cols[9:]:
    stem = re.sub(r'\.(cram|bam)\$', '', os.path.basename(col))
    if stem not in sample_map:
        sys.exit(
            f'Error: minipileup2 emitted sample column {col!r} (stem {stem!r}) '
            'that is not in the generated sample map. Known stems: '
            + ', '.join(sorted(sample_map)) + '. This means the CRAM list passed to '
            'minipileup2 and the --sr-cram/--lr-cram lists disagree.'
        )
    ftype, core = sample_map[stem]
    # Prefix every column when cores are in use, so a stem containing '__' can never
    # be mistaken for a core boundary. Core IDs contain no '_', so the first '__' in
    # a column name is always the core/stem separator.
    fixed.append(f'{core}__{stem}-{ftype}' if use_cores else f'{stem}-{ftype}')

vcf_obj.header.columns = '\t'.join(fixed) + '\n'

with open('${OUTPUT_PRFX}.vcf', 'w') as fo:
    vcf_obj.write_header(fo)
    for vnt in vcf_obj.parse_variants():
        vcf_obj.write_variant(fo, vnt)
"

echo "Renaming samples..."
python3 -c "$py_script" || { echo "Error: sample renaming failed"; exit 1; }

bgzip -f "${OUTPUT_PRFX}.vcf"   || { echo "Error: bgzip failed"; exit 1; }
tabix -f -p vcf "${OUTPUT_PRFX}.vcf.gz" || { echo "Error: tabix failed"; exit 1; }

echo "Done: ${OUTPUT_PRFX}.vcf.gz"
