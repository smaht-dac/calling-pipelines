#!/usr/bin/env python3

################################################################################
### Libraries
################################################################################
import argparse, subprocess, re
from granite.lib import vcf_parser


################################################################################
### Globals
################################################################################
contig_re = re.compile(
    r'##contig=<ID=(?P<ID>[^,>]+)(?:,length=(?P<length>\d+))?'
)

natural_num_re = re.compile(r'(\d+)')

CANONICAL_CONTIGS = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}

VALID_CALLERS = [
    "TNhaplotyper2",
    "Strelka2",
    "longcallD",
    "RUFUS"
]

DEFINITIONS_TO_ADD = [
    '##FILTER=<ID=PASS,Description="Passed filters in at least one caller">',
    '##INFO=<ID=CALLERS,Number=.,Type=String,Description="List of variant callers that reported this variant">',
    ## RUFUS
    '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
    '##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">',
    '##ALT=<ID=INS:ME:MOB,Description="Insertion of ALU or L1 element">'
]


################################################################################
### Objects
################################################################################
#*******************************************************************************
# VcfHandler
#   This is a base class for loading VCF files using the granite VCF parser
#*******************************************************************************
class VcfHandler:

    def __init__(self, vcf_path):
        self.input_vcf = vcf_path
        self.vcf_reader = vcf_parser.Vcf(self.input_vcf)
        self.caller_name = "BaseVcfHandler"
        self.HEADER = self.vcf_reader.header
        self.RECORDS = {}

        # Load VCF
        for record in self.vcf_reader.parse_variants():
            self.RECORDS[record.repr()] = record

    def _create_empty_record(self, record_to_import):
        # Empty record with 8 fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
        record_placeholder = "\t".join([".", "1", ".", ".", ".", ".", ".", "."])
        empty_record = vcf_parser.Vcf.Variant(record_placeholder, [])

        # Complete CHROM, POS, REF, ALT from the record to import
        empty_record.CHROM = record_to_import.CHROM
        empty_record.POS = record_to_import.POS
        empty_record.REF = record_to_import.REF
        empty_record.ALT = record_to_import.ALT

        return empty_record

    def add_caller_specific_info(self, shared_record, source_record):
        """
        Hook to let subclasses add INFO tags or other annotations
        into the shared record, based on the per-caller source record.
        """
        pass

    def update_dict(self, target_dict):
        """
        Update shared records in a target dictionary
        with records from this VCF handler.
        """
        for record_repr, record in self.RECORDS.items():
            if record_repr not in target_dict:
                empty_record = self._create_empty_record(record)
                target_dict.setdefault(record_repr, {
                    'record': empty_record,
                    'callers': set()
                })

            # Add caller information
            target_dict[record_repr]['callers'].add(self.caller_name)

            # Add caller-specific INFO, if any
            shared_record = target_dict[record_repr]['record']
            self.add_caller_specific_info(shared_record, record)

#*******************************************************************************
# TNhaplotyper2
#*******************************************************************************
class TNhaplotyper2Vcf(VcfHandler):

    def __init__(self, vcf_path):
        super().__init__(vcf_path)
        self.caller_name = "TNhaplotyper2"

#*******************************************************************************
# Strelka2
#*******************************************************************************
class Strelka2Vcf(VcfHandler):

    def __init__(self, vcf_path):
        super().__init__(vcf_path)
        self.caller_name = "Strelka2"

#*******************************************************************************
# longcallD
#*******************************************************************************
class longcallDVcf(VcfHandler):

    def __init__(self, vcf_path):
        super().__init__(vcf_path)
        self.caller_name = "longcallD"

#*******************************************************************************
# RUFUS
#*******************************************************************************
class RUFUSVcf(VcfHandler):

    def __init__(self, vcf_path):
        super().__init__(vcf_path)
        self.caller_name = "RUFUS"


################################################################################
### Functions
################################################################################
def validate_input_vcf(value):
    """
    Validate and return a string of the form CALLER:VCF.
    """
    if ":" not in value:
        raise argparse.ArgumentTypeError(
            f"Invalid format '{value}'. Expected CALLER:VCF"
        )

    # Enforce allowed CALLER names
    caller, _ = value.split(":", 1)
    if caller not in VALID_CALLERS:
        allowed = ", ".join(VALID_CALLERS)
        raise argparse.ArgumentTypeError(
            f"Caller '{caller}' does not match allowed callers: {allowed}"
        )

    return value

def validate_output_vcf(path):
    """
    Enforce output VCF file extension.
    """
    if path.endswith(".vcf") or path.endswith(".vcf.gz"):
        return path
    raise argparse.ArgumentTypeError(
        "Output must end with .vcf or .vcf.gz"
    )

def natural_key(s):
    """
    Split a string into text and integer chunks so that
    'KI270706v1' < 'KI270707v1' in numeric order.
    Returns a tuple suitable for sorting.
    """
    parts = natural_num_re.split(s)
    key = []
    for p in parts:
        if p.isdigit():
            key.append(int(p))
        else:
            key.append(p)
    return tuple(key)

def contig_sort_key(contig_id):
    """
    Sort key for contigs.

    Canonical (chr1-chr22, chrX, chrY, chrM/MT *without* suffix) first,
    ordered by chromosome.
    Non-canonical (random, alt, decoy, etc.) after, ordered by natural
    sort of the ID (ignoring 'chr' prefix).
    """
    if contig_id.startswith("chr"):
        cid_nochr = contig_id[3:]
    else:
        cid_nochr = contig_id

    # base: part before the first '_', e.g. '1' in '1_KI270706v1_random'
    base = cid_nochr.split("_", 1)[0]
    has_suffix = "_" in cid_nochr

    # canonical = plain chr1..22, chrX, chrY, chrM/MT (no suffix)
    canonical = (not has_suffix) and base in CANONICAL_CONTIGS

    if canonical:
        if base == "X":
            chrom_order = 23
        elif base == "Y":
            chrom_order = 24
        elif base in ("M", "MT"):
            chrom_order = 25
        else:
            chrom_order = int(base)
        # group 0 = canonical
        return (0, chrom_order)

    # Non-canonical: group 1, then natural sort of the ID without 'chr'
    return (1, natural_key(cid_nochr))

def get_chrom_pos(record_repr):
    """
    Sort key for variants.
    """
    chrom, right_repr = record_repr.split(":")

    # Extract POS (all leading digits)
    pos_str = ""
    for c in right_repr:
        if c.isdigit():
            pos_str += c
        else:
            break
    if not pos_str:
        raise ValueError(f"Could not parse POS from record representation: {record_repr}")

    return contig_sort_key(chrom) + (int(pos_str),)

def parse_contig(line):
    """
    Parse a VCF ##contig= line and return (ID, length).
    length is an int if present, otherwise None.
    """
    match = contig_re.match(line)
    if match:
        contig_id = match.group("ID")
        length_str = match.group("length")
        length = int(length_str) if length_str is not None else None
        return contig_id, length
    else:
        raise ValueError(f"Invalid contig line: {line}")

def bgzip_and_tabix(path):
    """
    Compress a VCF with bgzip and create a tabix index.
    Produces path + '.gz' and path + '.gz.tbi'.
    """
    # Compress (overwrites existing .gz if present)
    subprocess.run(["bgzip", "-f", path], check=True)

    # Index with tabix (VCF preset)
    gz_path = path + ".gz"
    subprocess.run(["tabix", "-f", "-p", "vcf", gz_path], check=True)

################################################################################
### Main
################################################################################
def main(args):

    # Validate input VCF and create handlers
    caller_handlers = []
    for caller_vcf in args.input_vcf:
        caller, vcf_path = caller_vcf.split(":", 1)

        if caller == "TNhaplotyper2":
            handler = TNhaplotyper2Vcf(vcf_path)
        elif caller == "Strelka2":
            handler = Strelka2Vcf(vcf_path)
        elif caller == "longcallD":
            handler = longcallDVcf(vcf_path)
        elif caller == "RUFUS":
            handler = RUFUSVcf(vcf_path)

        caller_handlers.append(handler)

    # Merge VCF records
    merged_records = {}
    for handler in caller_handlers:
        handler.update_dict(merged_records)

    # Decide where to write the uncompressed file
    if args.output_vcf.endswith(".vcf.gz"):
        write_path = args.output_vcf[:-3]  # strip .gz -> .vcf
    else:
        write_path = args.output_vcf

    # Collect contigs and check lengths are consistent across callers
    contigs = {}          # contig_id -> length
    contigs_no_length = set()
    saw_chr_prefixed = False
    saw_unprefixed_canonical = False
    saw_MT = False
    saw_M = False

    for handler in caller_handlers:
        for header_line in handler.HEADER.definitions.splitlines():
            if header_line.startswith("##contig="):
                contig_id, length = parse_contig(header_line)

                # Check canonical contig naming style consistency (chr1 vs 1)
                if contig_id.startswith("chr"):
                    base = contig_id[3:]
                    if base in CANONICAL_CONTIGS:
                        saw_chr_prefixed = True
                else:
                    base = contig_id
                    if base in CANONICAL_CONTIGS:
                        saw_unprefixed_canonical = True

                # Check MT vs M naming consistency
                if base == "MT": saw_MT = True
                elif base == "M": saw_M = True

                if saw_chr_prefixed and saw_unprefixed_canonical:
                    raise ValueError(
                        "Inconsistent canonical contig naming across callers: "
                        "mix of 'chr1'-style and '1'-style contigs"
                    )
                if saw_MT and saw_M:
                    raise ValueError(
                        "Inconsistent mitochondrial contig naming across callers: "
                        "mix of 'chrM' and 'chrMT' (or 'M' and 'MT')"
                    )

                # Check length consistency
                if length is None:
                    contigs_no_length.add(contig_id)
                else:
                    # If length, check for consistency
                    if contig_id in contigs and contigs[contig_id] != length:
                        raise ValueError(
                            f"Conflicting lengths for contig {contig_id}: "
                            f"{contigs[contig_id]} vs {length}"
                        )
                    contigs.setdefault(contig_id, length)

    # If a contig appears without length in some callers, make sure at least one
    # caller provides its length
    for contig_id in contigs_no_length:
        if contig_id not in contigs:
            raise ValueError(
                f"Contig {contig_id} missing length information from other callers"
            )

    # Sort contigs and generate contig lines
    contig_lines = ""
    for contig_id in sorted(contigs.keys(), key=contig_sort_key):
        length = contigs[contig_id]
        contig_lines += f"##contig=<ID={contig_id},length={length}>\n"

    # Generate new definitions (FILTER/INFO/ALT + SAMPLE)
    new_definitions = "##fileformat=VCFv4.2\n" + contig_lines
    for def_line in DEFINITIONS_TO_ADD:
        new_definitions += def_line + "\n"

    new_definitions += f'##SAMPLE=<ID={args.sample_name}>\n'

    # Update columns line
    columns_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

    # Sort merged records by CHROM and POS and write to output VCF
    with open(write_path, "w") as out_vcf:
        # Write cleaned header
        out_vcf.write(new_definitions)
        out_vcf.write(columns_line)
        # Write variant records
        for record_repr in sorted(merged_records.keys(), key=get_chrom_pos):
            entry = merged_records[record_repr]
            # Update FILTER field
            if entry['record'].FILTER == ".":
                entry['record'].FILTER = "PASS"

            # Update CALLERS INFO field
            callers_info = ','.join(sorted(entry['callers']))
            entry['record'].add_tag_info(f"CALLERS={callers_info}")

            # Write record
            out_vcf.write(entry['record'].to_string())

    # If output is .vcf.gz, bgzip and tabix
    if args.output_vcf.endswith(".vcf.gz"):
        bgzip_and_tabix(write_path)


################################################################################
### Entry Point
################################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Merge VCF files from multiple variant callers into a single VCF.",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "-i", "--input_vcf",
        action="append",
        required=True,
        metavar="CALLER:VCF",
        type=validate_input_vcf,
        help=(
        "Caller and its VCF file, separated by a colon.\n"
        "Format: CALLER:VCF\n"
        "Example:\n"
        "    -i TNhaplotyper2:file1.vcf -i Strelka2:file2.vcf.gz\n"
        "Caller must be one of: TNhaplotyper2, Strelka2, longcallD, RUFUS."
        )
    )
    parser.add_argument(
        "-o", "--output_vcf",
        default="output.vcf",
        type=validate_output_vcf,
        help=(
        "Output VCF file path (must be .vcf or .vcf.gz).\n"
        "If .vcf.gz is requested, bgzip and tabix must be available in PATH."
        )
    )
    parser.add_argument(
        '-s', '--sample_name',
        type=str,
        required=True,
        help='Sample name to use in the output VCF.')

    args = parser.parse_args()
    main(args)
