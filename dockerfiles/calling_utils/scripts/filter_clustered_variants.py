#!/usr/bin/env python3
import argparse
import pysam

def filter_clustered(input_vcf, output_vcf, window=100):
    vcf_in = pysam.VariantFile(input_vcf)
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    last_rec = None
    last_clustered = False

    for rec in vcf_in.fetch():
        if last_rec and rec.chrom == last_rec.chrom and rec.pos - last_rec.pos <= window:
            # Cluster detected: drop last and current
            last_clustered = True
            continue
        else:
            # If the previous one was not clustered, write it out
            if last_rec and not last_clustered:
                vcf_out.write(last_rec)
            last_rec = rec
            last_clustered = False

    # Flush last record if safe
    if last_rec and not last_clustered:
        vcf_out.write(last_rec)

    vcf_in.close()
    vcf_out.close()

    if output_vcf.endswith(".gz"):
        pysam.tabix_index(output_vcf, preset="vcf", force=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove clustered variants (SNVs + indels) within N bp")
    parser.add_argument("input_vcf", help="Input VCF file (bgzipped or not)")
    parser.add_argument("output_vcf", help="Output VCF file (will be bgzipped if ends with .gz)")
    parser.add_argument("--window", type=int, default=100, help="Window size in bp [default=100]")
    args = parser.parse_args()

    filter_clustered(args.input_vcf, args.output_vcf, args.window)

