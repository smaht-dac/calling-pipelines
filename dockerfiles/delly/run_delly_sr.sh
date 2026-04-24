#!/usr/bin/env bash
set -euo pipefail

## Command line arguments
while getopts ":n:r:x:s:" opt; do
  case $opt in
    # Output file prefix
    n) output_file_prefix="$OPTARG" ;;
    # Reference file
    r) reference_fasta="$OPTARG" ;;
    # TSV regions to exclude
    x) exclude_tsv="$OPTARG" ;;
    # Input cram file
    s) short_read_input="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

[[ -z "$output_file_prefix" || -z "$reference_fasta" || -z "$short_read_input" || -z "$exclude_tsv" ]] && {
  echo "Error: missing required arguments"
  exit 1
}

# **********************************************
# 1. Run Delly short-read command line
# **********************************************

delly call -g "$reference_fasta" -x "$exclude_tsv" -o "${output_file_prefix}.bcf" "$short_read_input" || { echo "Error: delly call failed"; exit 1; }

# **********************************************
# 2. Convert Delly short-read output to compressed vcf and index
# **********************************************

bcftools view -Oz -o "${output_file_prefix}.vcf.gz" "${output_file_prefix}.bcf" || { echo "Error: bcftools view failed"; exit 1; }
tabix -p vcf "${output_file_prefix}.vcf.gz" || { echo "Error: tabix failed"; exit 1; }
