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
    '##INFO=<ID=CORE_CALLS,Number=1,Type=String,Description="Per-core caller presence: core1:caller1,caller2|core2:caller1. Used by tier script to assign GT per sample column.">',
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

    def update_dict(self, target_dict, core):
        """
        Update shared records in a target dictionary
        with records from this VCF handler.
        core: sequencing core ID (e.g. '001C1' or the merged core
              '001B2-001C5'), supplied by -c/--vcf_core_id. Per-core caller
              presence is tracked for the CORE_CALLS INFO field.
        """
        for record_repr, record in self.RECORDS.items():
            if record_repr not in target_dict:
                empty_record = self._create_empty_record(record)
                target_dict.setdefault(record_repr, {
                    'record': empty_record,
                    'callers': set(),
                    'core_callers': {}   # core -> set(callers)
                })

            # Add caller information
            target_dict[record_repr]['callers'].add(self.caller_name)

            # Track per-core caller presence
            target_dict[record_repr]['core_callers'].setdefault(core, set()).add(self.caller_name)

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

def validate_core_id(value):
    """
    Validate a --vcf_core_id value of the form CALLER:CORE_ID.

    First, split on : to separate CALLER and CORE_ID. Then validate CALLER and CORE_ID separately.

    Two classes of character are rejected in CORE_ID validation:

    CORE_CALLS delimiters — the field uses '|' to separate cores, ':' to
    separate a core from its caller list, and ',' to separate callers within
    that list. A core ID containing any of these would make CORE_CALLS
    ambiguous to parse.

    VCF INFO structure — ';' separates INFO fields, '=' separates an INFO key
    from its value, and whitespace terminates the INFO column. A core ID
    containing any of these would not merely confuse a parser, it would
    produce a structurally invalid VCF record.

    Note that '-' is deliberately allowed: downstream, a hyphen marks a merged
    core (e.g. '001B2-001C5') and is split on to recover the constituents.
    """
    if ":" not in value:
        raise argparse.ArgumentTypeError(
            f"Invalid format '{value}'. Expected CALLER:CORE_ID"
        )

    caller, core_id = value.split(":", 1)
    
    # Enforce allowed CALLER names
    if caller not in VALID_CALLERS:
        allowed = ", ".join(VALID_CALLERS)
        raise argparse.ArgumentTypeError(
            f"Caller '{caller}' does not match allowed callers: {allowed}"
        )

    # Enforce non-empty core ID and disallow rejected characters
    if not core_id:
        raise argparse.ArgumentTypeError("Core ID must not be empty")
    for bad_char in (":", "|", ","):
        if bad_char in core_id:
            raise argparse.ArgumentTypeError(
                f"Invalid core ID '{core_id}': must not contain '{bad_char}' "
                "(reserved as a CORE_CALLS delimiter)"
            )
    for bad_char in (";", "="):
        if bad_char in core_id:
            raise argparse.ArgumentTypeError(
                f"Invalid core ID '{core_id}': must not contain '{bad_char}' "
                "(reserved by the VCF INFO column)"
            )
    if any(c.isspace() for c in core_id):
        raise argparse.ArgumentTypeError(
            f"Invalid core ID '{core_id}': must not contain whitespace "
            "(whitespace terminates the VCF INFO column)"
        )
    return value

def count_by_caller(caller_prefixed_values):
    """
    Count the number of entries per caller in a list of CALLER:VALUE strings (eg. VCF paths or core IDs).
    """
    counts = {}
    for value in caller_prefixed_values:
        caller, _rest = value.split(":", 1)
        counts[caller] = counts.get(caller, 0) + 1
    return counts

def normalize_core_id(core):
    """
    Normalize a (possibly merged) core ID by sorting its hyphen-separated
    constituents, e.g. '002B3-002B1' -> '002B1-002B3'.

    This ensures a merged core is represented identically regardless of the
    order its constituents were passed in on the command line, so that:
      - runs using '002B3-002B1' and '002B1-002B3' merge into the SAME
        core_callers dict entry (their caller sets get unioned) rather than
        two separate entries for what is semantically the same merged core
      - the CORE_CALLS INFO field always reports merged cores in a single,
        consistent constituent order

    Single (non-hyphenated) core IDs are returned unchanged.
    """
    if "-" not in core:
        return core
    return "-".join(sorted(core.split("-")))

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

    # Validate input VCF, normalize core IDs, and create handlers.
    # Cores are matched to -i by position: the Nth -c applies to
    # the Nth -i. Splitting on the first colon only keeps VCF paths that
    # contain colons (e.g. 's3://bucket/key.vcf.gz') intact.
    vcfs_by_caller = {}    # caller -> list of vcf_path, in order of appearance
    for caller_vcf in args.input_vcf:
        caller, vcf_path = caller_vcf.split(":", 1)
        vcfs_by_caller.setdefault(caller, []).append(vcf_path)

    cores_by_caller = {}   # caller -> list of normalized core, in order of appearance
    for caller_core in args.vcf_core_id:
        caller, core = caller_core.split(":", 1)
        cores_by_caller.setdefault(caller, []).append(normalize_core_id(core))

    caller_handlers = []   # list of (core, handler)
    for caller, vcf_paths in vcfs_by_caller.items():
        cores = cores_by_caller.get(caller, [])
        for vcf_path, core in zip(vcf_paths, cores):
            if caller == "TNhaplotyper2":
                handler = TNhaplotyper2Vcf(vcf_path)
            elif caller == "Strelka2":
                handler = Strelka2Vcf(vcf_path)
            elif caller == "longcallD":
                handler = longcallDVcf(vcf_path)
            elif caller == "RUFUS":
                handler = RUFUSVcf(vcf_path)

            caller_handlers.append((core, handler))

    # Merge VCF records
    merged_records = {}
    for core, handler in caller_handlers:
        handler.update_dict(merged_records, core=core)

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

    for _, handler in caller_handlers:
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

            # Write CORE_CALLS INFO field
            core_callers = entry['core_callers']
            core_calls_str = '|'.join(
                f"{core}:{','.join(sorted(callers))}"
                for core, callers in sorted(
                    core_callers.items(),
                    key=lambda x: ('-' in x[0], x[0])
                )
            )
            entry['record'].add_tag_info(f"CORE_CALLS={core_calls_str}")

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
        "-c", "--vcf_core_id",
        action="append",
        required=True,
        metavar="CALLER:CORE",
        type=validate_core_id,
        help=(
        "Sequencing core ID for the -i VCF at the SAME POSITION on the\n"
        "command line. Required, and repeated once per -i: the 1st -c applies\n"
        "to the 1st -i, the 2nd -c to the 2nd -i, and so on. The number of -c\n"
        "values must equal the number of -i values.\n"
        "\n"
        "Cores are used to build the CORE_CALLS INFO field, which records\n"
        "which callers reported each variant in each core, in the form\n"
        "    core1:caller1,caller2|core2:caller1\n"
        "The tier script reads CORE_CALLS to assign a GT per sample column.\n"
        "A core ID may therefore not contain ':', '|' or ','.\n"
        "\n"
        "The same core may be repeated across several -i entries (one per\n"
        "caller run on that core); its caller sets are unioned.\n"
        "\n"
        "Merged cores are written through unchanged aside from having their\n"
        "constituents sorted, so a hyphenated ID such as '001C5-001B2' is\n"
        "normalized to '001B2-001C5' (or a three-way merge is sorted the same\n"
        "way) and stays a single core here, distinct from its constituents.\n" 
        "Downstream, the tier script treats any core containing '-' as a merged\n"
        "core and splits it on '-' to find those constituents, so a hyphen is\n"
        "meaningful and must not be used in an ordinary core name.\n"
        "\n"
        "Example (two callers on core 001B2, one on merged core 001B2-001C5):\n"
        "    -i RUFUS:a.vcf.gz         -c RUFUS:001B2 \\\n"
        "    -i Strelka2:b.vcf.gz      -c Strelka2:001B2 \\\n"
        "    -i TNhaplotyper2:c.vcf.gz -c TNhaplotyper2:001B2-001C5\n"
        "\n"
        "Ordering note: -i and -c are collected into two independent lists in\n"
        "command-line order, so they need not be interleaved. They only need\n"
        "to appear in the same relative order as each other."
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

    # -c is matched to -i per caller (both carry a CALLER: prefix), so tally
    # -i and -c counts per caller and report any caller whose counts differ.
    vcf_counts = count_by_caller(args.input_vcf)
    core_counts = count_by_caller(args.vcf_core_id)

    mismatches = []
    for caller in sorted(set(vcf_counts) | set(core_counts)):
        n_vcf = vcf_counts.get(caller, 0)
        n_core = core_counts.get(caller, 0)
        if n_vcf != n_core:
            mismatches.append(
                f"  {caller}: {n_vcf} -i/--input_vcf entr{'y' if n_vcf == 1 else 'ies'} "
                f"vs. {n_core} -c/--vcf_core_id entr{'y' if n_core == 1 else 'ies'}"
            )

    if mismatches:
        parser.error(
            "-i/--input_vcf and -c/--vcf_core_id counts do not match for the "
            "following caller(s) (each -i needs exactly one -c for that same "
            "caller):\n" + "\n".join(mismatches)
        )

    main(args)
