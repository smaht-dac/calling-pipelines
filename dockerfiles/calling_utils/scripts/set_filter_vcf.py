#!/usr/bin/env python3
import argparse
import pysam

###############################################################################
# Clone a pysam VariantRecord into a new record tied to output header
###############################################################################
def clone_record(rec, out_header):
    """
    Create a new record in the output VCF that duplicates the input record.
    This ensures FILTER + INFO fields added to the output header are valid.
    """

    key_order = [
        "CrossTech",
        "CrossCaller",
        "CrossTissue",
        "CALLERS",
        "SR_VAF",
        "POOLED_PB_VAF",
        "POOLED_ONT_VAF",
        "TISSUE_SR_VAFS",
        "SR_ADF",
        "SR_ADR",
        "PB_ADF",
        "PB_ADR",
        "ONT_ADF",
        "ONT_ADR",
        "SB_SRC",
        "SB_PVAL",
        "GERMLINE_PVAL",
        "GERMLINE_PVAL_SR",
        "GERMLINE_PVAL_PB",
        "GERMLINE_PVAL_ONT",
        "PB_PHASING"
    ]

    info = dict(rec.info)
    ordered_info = {
        k: info[k]
        for k in key_order
        if k in info
    }


    new_rec = out_header.new_record(
        contig=rec.contig,
        start=rec.start,
        stop=rec.stop,
        id=rec.id,
        qual=rec.qual,
        alleles=rec.alleles,
        info=ordered_info,
        filter=None               # we overwrite FILTER
    )

    # copy FORMAT/sample fields
    for sample in rec.samples:
        new_rec.samples[sample].update(rec.samples[sample].items())

    return new_rec

def fix_header(header):
    """
    Create final output header
    """

    final_headers= [
            '##INFO=<ID=CrossTech,Number=0,Type=Flag,Description="Alt supported in both short read and PacBio data at or above their combined thresholds">',
            '##INFO=<ID=CrossCaller,Number=0,Type=Flag,Description="Alt found in more than one variant caller">',
            '##INFO=<ID=CrossTissue,Number=0,Type=Flag,Description="Variant has VAF > 0 in another short read tissue">',
            '##INFO=<ID=CALLERS,Number=.,Type=String,Description="List of variant callers that reported this variant">',

            '##INFO=<ID=SR_VAF,Number=1,Type=Float,Description="VAF for short read in current tissue">',
            '##INFO=<ID=POOLED_PB_VAF,Number=1,Type=Float,Description="VAF for PacBio in current donor pooled tissues">',
            '##INFO=<ID=POOLED_ONT_VAF,Number=1,Type=Float,Description="VAF for ONT in current donor pooled tissues">',
            '##INFO=<ID=TISSUE_SR_VAFS,Number=.,Type=String,Description="VAFs for all tissues with short read nonzero VAF, based on pileup (BQ≥30)">',

            '##INFO=<ID=SR_ADF,Number=2,Type=Integer,Description="tissue short-read forward depths (REF,ALT)">',
            '##INFO=<ID=SR_ADR,Number=2,Type=Integer,Description="tissue short-read reverse depths (REF,ALT)">',
            '##INFO=<ID=PB_ADF,Number=2,Type=Integer,Description="donor pooled Long-read forward depths (REF,ALT)">',
            '##INFO=<ID=PB_ADR,Number=2,Type=Integer,Description="donor pooled Long-read reverse depths (REF,ALT)">',
            '##INFO=<ID=ONT_ADF,Number=2,Type=Integer,Description="donor pooled ONT forward depths (REF,ALT)">',
            '##INFO=<ID=ONT_ADR,Number=2,Type=Integer,Description="donor pooled ONT reverse depths (REF,ALT)">',

            '##INFO=<ID=SB_SRC,Number=1,Type=String,Description="Counts source used for Fisher strand test: PB, ONT, or SR">',
            '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher p-value for strand balance on chosen sample">',
            '##INFO=<ID=GERMLINE_PVAL,Number=1,Type=Float,Description="Minimum binomial p-value for germline deviation across all platforms tested">',
            '##INFO=<ID=GERMLINE_PVAL_SR,Number=1,Type=Float,Description="Binomial p-value for germline deviation in short-read data">',
            '##INFO=<ID=GERMLINE_PVAL_PB,Number=1,Type=Float,Description="Binomial p-value for germline deviation in PacBio data">',
            '##INFO=<ID=GERMLINE_PVAL_ONT,Number=1,Type=Float,Description="Binomial p-value for germline deviation in Oxford Nanopore (long-read) data">',

            '##INFO=<ID=PB_PHASING,Number=1,Type=String,Description="Phasing classification from pooled PacBio haplotyping">',

            '##FILTER=<ID=HighConf,Description="High confidence variant">',
            '##FILTER=<ID=LowConf,Description="Low confidence variant">',
            '##FILTER=<ID=.,Description="Variants passing all filters but with no CrossTech, CrossCaller, or CrossTissue evidence, lowest confidence variants">'
    ]

    header_list = str(header).split('\n')

    new_header = pysam.VariantHeader()
    for header_line in header_list:
        if "fileformat" in header_line or "contig" in header_line:
            new_header.add_line(header_line)
        if "SAMPLE" in header_line:
            sample_line = header_line

    for header_line in final_headers:
        new_header.add_line(header_line)

    new_header.add_line(sample_line)

    return new_header


###############################################################################
# Main
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Assign HighConf / LowConf / . FILTERs to variants"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input VCF (.vcf or .vcf.gz)"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output VCF (.vcf or .vcf.gz). Default: <input>.confidence.vcf.gz"
    )
    args = parser.parse_args()

    in_path = args.input

    # Determine default output path
    if args.output:
        out_path = args.output
    else:
        if in_path.endswith(".vcf.gz"):
            out_path = in_path.replace(".vcf.gz", ".confidence.vcf.gz")
        elif in_path.endswith(".vcf"):
            out_path = in_path.replace(".vcf", ".confidence.vcf")
        else:
            out_path = in_path + ".confidence.vcf.gz"

    ###########################################################################
    # Open input VCF and modify header to add filter definitions
    ###########################################################################
    vcf_in = pysam.VariantFile(in_path)
    header = vcf_in.header.copy()

    new_header = fix_header(header)

    ###########################################################################
    # Create output VCF with updated header
    ###########################################################################
    vcf_out = pysam.VariantFile(out_path, "w", header=new_header)

    ###########################################################################
    # Process variants
    ###########################################################################
    for rec in vcf_in:

        # detection flags
        crossTech   = "CrossTech"   in rec.info
        crossCaller = "CrossCaller" in rec.info
        crossTissue = "CrossTissue" in rec.info

        # clone record so filters can be added freely
        new = clone_record(rec, vcf_out.header)

        # Decision tree:
        # If CrossTech OR (CrossCaller AND CrossTissue) → HighConf
        # Else if CrossCaller OR CrossTissue → LowConf
        # Else → .

        new.filter.clear()

        if crossTech or (crossCaller and crossTissue):
            new.filter.add("HighConf")
        elif crossCaller or crossTissue:
            new.filter.add("LowConf")
        else:
            pass

        vcf_out.write(new)

    vcf_out.close()
    print(f"[set_filter_vcf] Wrote output: {out_path}")


###############################################################################
# Entrypoint
###############################################################################
if __name__ == "__main__":
    main()

