#!/usr/bin/env python3
import argparse
import pysam


import pysam

def filter_clustered(input_vcf, output_vcf, window=50):
    vcf_in = pysam.VariantFile(input_vcf)
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    group = []            # records in the current proximity group
    group_clustered = False
    prev = None

    for rec in vcf_in.fetch():
        if prev is None or rec.chrom != prev.chrom or rec.pos - prev.pos > window:
            # group ended; flush previous group
            if group and not group_clustered:
                vcf_out.write(group[0])   # singleton group only
            # start new group
            group = [rec]
            group_clustered = False
        else:
            # still in the same group (adjacent within window)
            group.append(rec)
            group_clustered = True

        prev = rec

    # flush last group
    if group and not group_clustered:
        vcf_out.write(group[0])

    vcf_in.close()
    vcf_out.close()

    if output_vcf.endswith(".gz"): pysam.tabix_index(output_vcf,preset="vcf", force=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove clustered variants (SNVs + indels) within N bp")
    parser.add_argument("input_vcf", help="Input VCF file (bgzipped or not)")
    parser.add_argument("output_vcf", help="Output VCF file (will be bgzipped if ends with .gz)")
    parser.add_argument("--window", type=int, default=50, help="Window size in bp [default=50]")
    args = parser.parse_args()

    filter_clustered(args.input_vcf, args.output_vcf, args.window)

