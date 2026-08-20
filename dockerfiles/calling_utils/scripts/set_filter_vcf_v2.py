#!/usr/bin/env python3
import argparse
import pysam


FORMAT_DEFS = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 0/1=called+passed all filters, 0/0=did not pass filtering">',
    '##FORMAT=<ID=TECH,Number=1,Type=String,Description="Sequencing technology of this core: SR=short-read, PB=PacBio/long-read">',
    '##FORMAT=<ID=SR_ADF,Number=2,Type=Integer,Description="Per-core short-read forward depths (REF,ALT)">',
    '##FORMAT=<ID=SR_ADR,Number=2,Type=Integer,Description="Per-core short-read reverse depths (REF,ALT)">',
    '##FORMAT=<ID=PB_ADF,Number=2,Type=Integer,Description="Per-core tissue-matched PacBio forward depths (REF,ALT)">',
    '##FORMAT=<ID=PB_ADR,Number=2,Type=Integer,Description="Per-core tissue-matched PacBio reverse depths (REF,ALT)">',
    '##FORMAT=<ID=SR_VAF,Number=1,Type=Float,Description="Per-core short-read variant allele fraction">',
    '##FORMAT=<ID=PB_VAF,Number=1,Type=Float,Description="Per-core tissue-matched PacBio variant allele fraction">',
    '##FORMAT=<ID=ALT_SUPPORT_PASS,Number=1,Type=Integer,Description="1 if this core\'s ALT reads meet the standalone read cutoff for its coverage, 0 otherwise">',
    '##FORMAT=<ID=CALLERS,Number=.,Type=String,Description="Callers that reported this variant for this core">',
]

INFO_DEFS = [
    '##INFO=<ID=EvidenceScore,Number=1,Type=Integer,Description="Count of cross-evidence labels present (CrossTech, CrossCaller, CrossTissue, CrossCore)">',
    '##INFO=<ID=CrossCaller,Number=0,Type=Flag,Description="Alt found in 2+ unique callers across all cores">',
    '##INFO=<ID=CrossCore,Number=0,Type=Flag,Description="Variant has GT=0/1 in more than one core">',
    '##INFO=<ID=CrossTech,Number=0,Type=Flag,Description="Alt supported in both short read and tissue-matched PacBio at or above combined thresholds">',
    '##INFO=<ID=CrossTissue,Number=0,Type=Flag,Description="Variant has VAF > 0 in another tissue">',
    '##INFO=<ID=POOLED_PB_VAF,Number=1,Type=Float,Description="PacBio VAF pooled across all of a donor\'s tissues">',
    '##INFO=<ID=POOLED_ONT_VAF,Number=1,Type=Float,Description="ONT VAF pooled across all of a donor\'s tissues">',
    '##INFO=<ID=POOLED_PB_ADF,Number=2,Type=Integer,Description="PacBio forward depths pooled across all of a donor\'s tissues (REF,ALT)">',
    '##INFO=<ID=POOLED_PB_ADR,Number=2,Type=Integer,Description="PacBio reverse depths pooled across all of a donor\'s tissues (REF,ALT)">',
    '##INFO=<ID=POOLED_ONT_ADF,Number=2,Type=Integer,Description="ONT forward depths pooled across all of a donor\'s tissues (REF,ALT)">',
    '##INFO=<ID=POOLED_ONT_ADR,Number=2,Type=Integer,Description="ONT reverse depths pooled across all of a donor\'s tissues (REF,ALT)">',
    '##INFO=<ID=GERMLINE_PVAL,Number=1,Type=Float,Description="Min binomial p-value for germline deviation across pooled platforms (tissue-level SR; donor-level PB and ONT)">',
    '##INFO=<ID=GERMLINE_PVAL_SR,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled SR data (tissue-level)">',
    '##INFO=<ID=GERMLINE_PVAL_PB,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled PacBio data across all of a donor\'s tissues">',
    '##INFO=<ID=GERMLINE_PVAL_ONT,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled ONT data across all of a donor\'s tissues">',
    '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher exact test p-value for strand balance on pooled counts">',
    '##INFO=<ID=SB_SRC,Number=1,Type=String,Description="Platform used for pooled Fisher strand test: SR or PB">',
    '##INFO=<ID=SR_TISSUE_PRESENCE,Number=.,Type=String,Description="Tissues from this donor with short-read support for this variant based on pileups (BQ>30)">',
]

# '##INFO=<ID=REGION,Number=1,Type=String,Description="SMaHT region classification: easy, diff, or ext">',
# removed because not in final vcf

FILTER_DEFS = [
    '##FILTER=<ID=PASS,Description="At least one cross-evidence label present (CrossTech, CrossCaller, CrossTissue, or CrossCore)">',
    '##FILTER=<ID=LowEvidence,Description="No cross-evidence: lacks CrossTech, CrossCaller, CrossTissue, and CrossCore">',
]


def build_ordered_header(header):
    """
    Build the canonical final-VCF header from scratch.
    Contigs and sample names are taken from the input header;
    all FORMAT/INFO/FILTER definitions are the source of truth here.
    """
    new_header = pysam.VariantHeader()

    for line in str(header).split('\n'):
        if line.startswith('##contig='):
            new_header.add_line(line)

    for d in FORMAT_DEFS + INFO_DEFS + FILTER_DEFS:
        new_header.add_line(d)

    for sample in header.samples:
        new_header.add_sample(sample)

    return new_header


def _set_info(rec, key, val):
    try:
        rec.info[key] = val
    except Exception:
        # pysam reads Number=1 String fields containing commas as a tuple;
        # rejoin to reconstruct the original string.
        if isinstance(val, tuple):
            try:
                rec.info[key] = ','.join(str(v) for v in val)
            except Exception:
                pass


def write_record(rec, evidence_score, out_header, vcf_out):
    """
    Write a variant to vcf_out with INFO fields ordered:
    EvidenceScore, Cross* tags (alphabetically), then all others.
    """
    new = out_header.new_record(
        contig=rec.contig,
        start=rec.start,
        stop=rec.stop,
        id=rec.id,
        qual=rec.qual,
        alleles=rec.alleles,
        filter=None,
    )

    # Collect existing INFO, dropping TIER and any stale EvidenceScore
    info = {k: v for k, v in rec.info.items() if k not in ('TIER', 'EvidenceScore')}

    # Set INFO in desired order: EvidenceScore → Cross* → others
    _set_info(new, 'EvidenceScore', evidence_score)
    for key in sorted(k for k in list(info) if k.startswith('Cross')):
        _set_info(new, key, info.pop(key))
    for key, val in info.items():
        _set_info(new, key, val)

    # Copy FORMAT/sample fields
    for sample in rec.samples:
        for key, val in rec.samples[sample].items():
            try:
                new.samples[sample][key] = val
            except Exception:
                pass

    new.filter.clear()
    if evidence_score == 0:
        new.filter.add("LowEvidence")
    else:
        new.filter.add("PASS")

    vcf_out.write(new)


###############################################################################
# Main
###############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Assign PASS / LowEvidence FILTERs and EvidenceScore to variants"
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

    if args.output:
        out_path = args.output
    else:
        if in_path.endswith(".vcf.gz"):
            out_path = in_path.replace(".vcf.gz", ".confidence.vcf.gz")
        elif in_path.endswith(".vcf"):
            out_path = in_path.replace(".vcf", ".confidence.vcf")
        else:
            out_path = in_path + ".confidence.vcf.gz"

    vcf_in = pysam.VariantFile(in_path)
    new_header = build_ordered_header(vcf_in.header)
    vcf_out = pysam.VariantFile(out_path, "w", header=new_header)

    for rec in vcf_in:
        crossTech   = "CrossTech"   in rec.info
        crossTissue = "CrossTissue" in rec.info
        crossCaller = "CrossCaller" in rec.info
        crossCore   = "CrossCore"   in rec.info
        evidence_score = sum([crossTech, crossTissue, crossCaller, crossCore])

        write_record(rec, evidence_score, vcf_out.header, vcf_out)

    vcf_out.close()
    print(f"[set_filter_vcf] Wrote output: {out_path}")


###############################################################################
# Entrypoint
###############################################################################
if __name__ == "__main__":
    main()
