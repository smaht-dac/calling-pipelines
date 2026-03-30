#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -i input.cram -r reference.fa [-o prefix] [-t threads] [-h]

  -i   Long read normal BAM/CRAM with index (.bai/.crai) (required)
  -r   Reference FASTA with index (.fai) (required)
  -o   Output file prefix (default: output) -> <prefix>.plup.gz
  -t   Threads (default: number of cores)
  -h   Help
EOF
  exit 1
}

# Defaults
LONG_READ_INPUT=""
REFERENCE_FASTA=""
OUTPUT_PREFIX="output"
THREADS="$(nproc)"

## Command line arguments
while getopts ":i:r:l:t:" opt; do
  case $opt in
    i) LONG_READ_INPUT="$OPTARG" ;; # Input cram file
    r) REFERENCE_FASTA="$OPTARG" ;; # Reference file
    o) OUTPUT_PREFIX="$OPTARG" ;; # Output file prefix
    t) THREADS="$OPTARG" ;; # Threads
    h) usage ;;
    *) usage ;;
  esac
done

# Check required arguments
[[ -n "$LONG_READ_INPUT" && -n "$REFERENCE_FASTA" ]] || usage


# File checks
[[ -f "$LONG_READ_INPUT" ]] || { echo "Error: $LONG_READ_INPUT not found"; exit 1; }
[[ -f "${LONG_READ_INPUT}.crai" || -f "${LONG_READ_INPUT}.bai" ]] || { echo "Error: index for $LONG_READ_INPUT not found (.crai or .bai)"; exit 1; }

[[ -f "$REFERENCE_FASTA" ]] || { echo "Error: $REFERENCE_FASTA not found"; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "Error: ${REFERENCE_FASTA}.fai not found"; exit 1; }


# Run kanpig plup
/kanpig/target/release/kanpig plup --bam $LONG_READ_INPUT -r $reference --threads $thread | bedtools sort -header | bgzip > ${outputFile}.plup.gz