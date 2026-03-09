#!/usr/bin/env bash

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

# **********************************************
# 1. Create Delly long-read command line
# **********************************************

command="delly lr -g $reference_fasta"
command+=" -o $sample_name-Delly-LR.bcf"
command+=" $long_read_input"

# **********************************************
# 2. Run Delly long-read command line
# **********************************************

eval $command || exit 1

# **********************************************
# 3. Convert Delly long-read output to compressed vcf and index
# **********************************************

bcftools view -Oz -o $output_file_prefix-Delly-LR.vcf.gz $output_file_prefix-Delly-LR.bcf
tabix -p vcf $output_file_prefix-Delly-LR.vcf.gz
