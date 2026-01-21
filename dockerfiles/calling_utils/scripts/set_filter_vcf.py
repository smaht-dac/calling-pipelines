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

    new_rec = out_header.new_record(
        contig=rec.contig,
        start=rec.start,
        stop=rec.stop,
        id=rec.id,
        qual=rec.qual,
        alleles=rec.alleles,
        info=dict(rec.info),      # copy INFO
        filter=None               # we overwrite FILTER
    )

    # copy FORMAT/sample fields
    for sample in rec.samples:
        new_rec.samples[sample].update(rec.samples[sample].items())

    return new_rec


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

    if "HighConf" not in header.filters:
        header.filters.add("HighConf", None, None, "High confidence variant")
    if "LowConf" not in header.filters:
        header.filters.add("LowConf", None, None, "Low confidence variant")

    ###########################################################################
    # Create output VCF with updated header
    ###########################################################################
    vcf_out = pysam.VariantFile(out_path, "w", header=header)

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

