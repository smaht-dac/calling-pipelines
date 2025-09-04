#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#################################################################################
# This script is used to run the VEP (Variant Effect Predictor) in parallel.
# Requires:
#   - bcftools
#   - tabix
#   - vep
#   - bgzip
#################################################################################

usage() {
  cat <<EOF
Usage: $0 -i input.vcf.gz -r reference.fasta -l regions.txt -v vep.tar.gz
          [-g gnomad.vcf.gz] [-o prefix] [-a assembly] [-s species]

  -i   Input VCF (bgzipped) with .tbi index (required)
  -r   Reference FASTA with .fai index (required)
  -l   Regions list file (one region per line, e.g. chr1:1-1000000) (required)
  -v   VEP database archive (tar.gz) (required)
  -g   gnomAD VCF file (bgzipped) with .tbi index (optional)
  -o   Prefix for the output annotated VCF (default: output)
       The final output will be <prefix>.vep.vcf.gz
  -a   Assembly (default: GRCh38)
  -s   Species (default: homo_sapiens)
  -h   Help
EOF
  exit 1
}

# Defaults
INPUT_VCF=""
REFERENCE_FASTA=""
REGIONFILE=""
VEP_DB=""
GNOMAD_VCF=""
OUTPUT_PRFX="output"
ASSEMBLY="GRCh38"
NTHREADS="$(nproc)"
SPECIES="homo_sapiens"

# Parse args
while getopts "i:r:l:v:g:o:a:s:h" opt; do
  case $opt in
    i) INPUT_VCF="$OPTARG" ;;
    r) REFERENCE_FASTA="$OPTARG" ;;
    l) REGIONFILE="$OPTARG" ;;
    v) VEP_DB="$OPTARG" ;;
    g) GNOMAD_VCF="$OPTARG" ;;
    o) OUTPUT_PRFX="$OPTARG" ;;
    a) ASSEMBLY="$OPTARG" ;;
    s) SPECIES="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Required checks
[[ -n "$INPUT_VCF" ]] || { echo "Error: input VCF (-i) is required"; usage; }
[[ -n "$REFERENCE_FASTA" ]] || { echo "Error: reference FASTA (-r) is required"; usage; }
[[ -n "$REGIONFILE" ]] || { echo "Error: regions list file (-l) is required"; usage; }
[[ -n "$VEP_DB" ]] || { echo "Error: VEP database archive (-v) is required"; usage; }

# Tool checks
command -v vep >/dev/null 2>&1 || { echo "Error: vep not found in PATH"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "Error: tabix not found in PATH"; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "Error: bgzip not found in PATH"; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found in PATH"; exit 1; }

# Input checks
[[ -f "$INPUT_VCF" ]] || { echo "Error: $INPUT_VCF not found"; exit 1; }
[[ -f "${INPUT_VCF}.tbi" ]] || { echo "Error: ${INPUT_VCF}.tbi index not found"; exit 1; }
[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }
[[ -f "$VEP_DB" ]] || { echo "Error: VEP database archive $VEP_DB not found"; exit 1; }
[[ -s "$REGIONFILE" ]] || { echo "Error: regions file $REGIONFILE is empty"; exit 1; }

# Optional input checks
if [[ -n "$GNOMAD_VCF" ]]; then
  [[ -f "$GNOMAD_VCF" ]] || { echo "Error: $GNOMAD_VCF not found"; exit 1; }
  [[ -f "${GNOMAD_VCF}.tbi" ]] || { echo "Error: ${GNOMAD_VCF}.tbi index not found"; exit 1; }
fi

# Unpack VEP
EXTRACT_DIR="$(mktemp -d vep_cache.XXXXXX)"
echo "Unpacking VEP database in $EXTRACT_DIR ..."
tar -xzf "$VEP_DB" -C "$EXTRACT_DIR"

# Build VEP flags/options
FLAGS=(
    --hgvs --symbol --uniprot --pubmed --var_synonyms
    --biotype --variant_class --canonical --regulatory
    --protein --domains --gene_phenotype
)

OPTIONS=(
    --fasta "$REFERENCE_FASTA"
    --assembly "$ASSEMBLY"
    --use_given_ref
    --offline
    --cache
    --dir_cache "$EXTRACT_DIR"
    --force_overwrite
    --vcf
    --compress_output bgzip
    --no_stats
    --species "$SPECIES"
)

# Add gnomAD options if provided
if [[ -n "$GNOMAD_VCF" ]]; then
    OPTIONS+=( --custom file="$GNOMAD_VCF",short_name=gnomADc,format=vcf,type=exact,fields=AC_joint%AN_joint%AF_joint%AC_exomes%AN_exomes%AF_exomes%AC_genomes%AN_genomes%AF_genomes%AC_grpmax_joint%AF_grpmax_joint%AN_grpmax_joint%AC_grpmax_exomes%AF_grpmax_exomes%AN_grpmax_exomes%AC_grpmax_genomes%AF_grpmax_genomes%AN_grpmax_genomes )
fi

# Temp workspace
WORKDIR="$(mktemp -d vep_shards.XXXXXX)"
cleanup() {
  rm -rf "$WORKDIR"
  rm -rf "$EXTRACT_DIR"
}
trap cleanup EXIT

# Create a per-shard runner used by xargs
# Expects {} to expand to a region string like chr1:1-1000000
run_shard() {
  local region="$1"
  local input_vcf="$2"
  local workdir="$3"
  shift 3
  local vep_args=( "$@" ) 

  # Sanitize region string
  local safe="${region//:/_}"; safe="${safe//-/__}"
  local shard_vcf="$workdir/${safe}.vcf"
  local shard_out="$workdir/${safe}.vep.vcf.gz"

  # Extract region
  tabix -h "$input_vcf" "$region" > "$shard_vcf" || { echo "Error: tabix failed for $region"; exit 1; }

  # Run VEP
  if grep -qv '^#' "$shard_vcf"; then
    vep -i "$shard_vcf" -o "$shard_out" "${vep_args[@]}" \
      || { echo "Error: VEP failed for $region"; exit 1; }
    tabix -f -p vcf "$shard_out" || { echo "Error: tabix failed for $region"; exit 1; }
  fi
}

# Make the function available to xargs’ bash -c
export -f run_shard

# Parallel execution
echo "Sharding and running VEP in $WORKDIR directory..."
grep -Ev '^[[:space:]]*($|#)' "$REGIONFILE" | \
xargs -P "$NTHREADS" -I{} bash -c 'run_shard "$@"' _ \
  {} "$INPUT_VCF" "$WORKDIR" "${OPTIONS[@]}" "${FLAGS[@]}" \
  || { echo "Error: one or more shards failed"; exit 1; }

# Collect shard outputs (sorted, safe if none)
shopt -s nullglob
mapfile -t SHARDS < <(printf "%s\n" "$WORKDIR"/*.vep.vcf.gz | sort -V)
shopt -u nullglob
(( ${#SHARDS[@]} > 0 )) || { echo "Error: no shard outputs found in $WORKDIR"; exit 1; }

# Merging annotated shards using bcftools
echo "Merging ${#SHARDS[@]} annotated shards..."
MERGED="$WORKDIR/merged.vep.vcf"
FINAL_GZ="${OUTPUT_PRFX}.vep.vcf.gz"

echo "-- BCFTools ---------------------------"
bcftools concat -a -D -O v -o "$MERGED" "${SHARDS[@]}" \
  || { echo "Error: bcftools concat failed"; exit 1; }

bcftools sort -T "tmp_bcftools.XXXXXX" -O z -o "$FINAL_GZ" "$MERGED" \
  || { echo "Error: bcftools sort failed"; exit 1; }

tabix -f -p vcf "$FINAL_GZ" \
  || { echo "Error: tabix failed"; exit 1; }
echo "-----------------------------------------"

echo "Done: $FINAL_GZ"
