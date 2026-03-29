#!/usr/bin/env bash
set -euo pipefail

## Command line arguments
while getopts ":n:r:l:" opt; do
  case $opt in
    # Sample name
    n) output_file_prefix="$OPTARG" ;;
    # Reference file
    r) reference_fasta="$OPTARG" ;;
    # Input cram file
    l) long_read_input="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

[[ -z "$output_file_prefix" || -z "$reference_fasta" || -z "$long_read_input" ]] && {
  echo "Error: missing required arguments"
  exit 1
}