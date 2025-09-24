#!/bin/bash
set -euo pipefail

input_file=$1
output_prefix=$2

bgzip -c "$input_file" > "${output_prefix}.vcf.gz" || exit 1
tabix -p vcf "${output_prefix}.vcf.gz" || exit 1
