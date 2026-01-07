#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------
#  Map filtered germline VCF to CrossTech sites
#   Input: filtered germline VCF (AF > 0.001) and somatic VCF
#   Output: final mapping table linking each CrossTech variant
#            to the *closest* germline variant (within ±5 kb)
#
#  Note: this script assumes that the input somatic VCF only contains SNVs
#        Indels should be removed prior to running this step
# ---------------------------------------------------------------

usage() {
  echo "Usage: phasing_step1_get_closest_germline.sh -s sample -v germline.vcf.gz -w somatic.vcf.gz"
}

while getopts ":s:v:w:" opt; do
  case $opt in
    s) SAMPLE=$OPTARG ;;
    v) GERMLINE_VCF=$OPTARG ;;
    w) SOMATIC_VCF=$OPTARG ;;
    *) usage; exit 1 ;;
  esac
done

[[ -z ${SAMPLE-} || -z ${GERMLINE_VCF-} || -z ${SOMATIC_VCF-} ]] && { usage; exit 1; }

# ---- File paths ----
CROSSTECH_VCF="CROSSTECH.norm.vcf.gz"
SOMATIC_BED="${SAMPLE}.somatic_sites.bed"
SOMATIC_5KB="${SAMPLE}.somatic_5kb.bed"
GERM_SUBVCF="${SAMPLE}.germline.subset.vcf.gz"
GERM_BED="${SAMPLE}.germline.subset.bed"
CLOSEST="${SAMPLE}.closest.tsv"
MAP_TSV="${SAMPLE}.germline_map.tsv"
FAILED_TSV="${SAMPLE}.no_close_germline.tsv"

# ---- 1. Extract somatic CrossTech variants ----
echo "Extracting CrossTech somatic variants..."
bcftools view -i 'INFO/CrossTech=1' -Oz -o "$CROSSTECH_VCF" "$SOMATIC_VCF"
tabix -f "$CROSSTECH_VCF"

# ---- 2. Convert CrossTech VCF to BED ----
echo "Converting CrossTech VCF to BED..."
bcftools query -f'%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' "$CROSSTECH_VCF" \
  | sort -k1,1 -k2,2n -S1G > "$SOMATIC_BED"

# ---- 3. Make ±5 kb windows around CrossTech variants ----
echo "Creating ±5 kb windows around CrossTech sites..."
awk 'BEGIN{OFS="\t"}{s=$2-5000; if(s<0)s=0; e=$3+5000; print $1,s,e}' "$SOMATIC_BED" \
  | bedtools merge -i - > "$SOMATIC_5KB"

# ---- 4. Subset germline VCF to these regions ----
echo "Subsetting germline VCF to ±5 kb around CrossTech sites..."
bcftools view \
  -v snps -f '.' \
  -R "$SOMATIC_5KB" \
  -i 'GT~"0[|/]1" || GT~"1[|/]0"' \
  -Oz -o "$GERM_SUBVCF" "$GERMLINE_VCF"

tabix -f "$GERM_SUBVCF"

# ---- 5. Convert subsetted germline VCF to BED ----
echo "Converting subsetted germline VCF to BED..."
bcftools query -f'%CHROM\t%POS0\t%POS\t%REF\t%ALT\n' "$GERM_SUBVCF" \
  | sort -k1,1 -k2,2n -S1G > "$GERM_BED"

# Step 6a. Exclude identical variants from the germline subset
echo "Excluding identical variants from germline subset..."
awk 'NR==FNR{a[$1":"$3":"$4":"$5]; next}
     !($1":"$3":"$4":"$5 in a)' "$SOMATIC_BED" "$GERM_BED" > "${GERM_BED%.bed}.filtered.bed"

# Step 6b. Find closest (now guaranteed non-identical)
echo "Finding closest germline variants to each CrossTech somatic variant..."
bedtools closest -a "$SOMATIC_BED" -b "${GERM_BED%.bed}.filtered.bed" -d > "$CLOSEST"

# ---- 7. Split pass / fail by distance ----
echo "Splitting results into within-5kb and failed variants..."
awk '$11 != -1 && $11 <= 5000' "$CLOSEST" > tmp.pass.tsv
awk '$11 == -1 || $11 > 5000' "$CLOSEST" > "$FAILED_TSV"

# ---- 8. Format outputs ----
echo "Formatting output mapping table..."
awk '
BEGIN {OFS="\t"}
{
  var_chrom=$1; var_pos=$3; var_ref=$4; var_alt=$5;
  germ_pos=$8; germ_ref=$9; germ_alt=$10;
  if (germ_ref == "" || germ_ref == ".") {
    germ_pos="NA"; germ_ref="NA"; germ_alt="NA";
  }
  print var_chrom, var_pos, var_ref, var_alt, germ_pos, germ_ref, germ_alt;
}' tmp.pass.tsv | sort -k1,1 -k2,2n -S1G > tmp.pass.clean.tsv

awk '
BEGIN {OFS="\t"}
{
  var_chrom=$1; var_pos=$3; var_ref=$4; var_alt=$5;
  print var_chrom, var_pos, var_ref, var_alt, "NA", "NA", "NA";
}' "$FAILED_TSV" | sort -k1,1 -k2,2n -S1G > tmp.fail.clean.tsv

# ---- 9. Combine and finalize ----
{
  echo -e "chrom\tVar_pos\tVar_ref\tVar_alt\tGerm_pos\tGerm_ref\tGerm_alt"
  cat tmp.pass.clean.tsv tmp.fail.clean.tsv
} > "$MAP_TSV"

rm -f tmp.pass.tsv tmp.pass.clean.tsv tmp.fail.clean.tsv

echo "[Complete]"
echo "  Mapped variants  : $MAP_TSV"
echo "  No close germline: $FAILED_TSV"
echo "  Germline subset  : $GERM_SUBVCF (±5 kb windows)"
