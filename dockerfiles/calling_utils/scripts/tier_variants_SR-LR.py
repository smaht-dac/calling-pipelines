#!/usr/bin/env python3

################################################################################
### Libraries
################################################################################
import argparse, sys
import pysam
import math
from datetime import datetime

################################################################################
### Objects
################################################################################
#*******************************************************************************
# Statistical functions and filtering logic
#*******************************************************************************
def get_read_cutoffs(SR_cov: int, PB_cov: int,
                     error_rate: float = 0.001,
                     target_p: float = 1e-3):
    """
    Compute read count cutoffs for SR (Short-reads), PB (PacBio),
    and combined coverage using Poisson approximation to the binomial distribution.

    Parameters
    ----------
    SR_cov : int
        Short-read coverage
    PB_cov : int
        Long-read coverage
    error_rate : float, optional
        Sequencing error rate (default = 0.001)
    target_p : float, optional
        Target probability for random error (default = 1e-3)

    Returns
    -------
        {
            'SR': int,
            'combined_SR': int,
            'combined_PB': int,
        }
    """

    def poisson_tail(lmbda: float, k: int):
        """Compute P(X >= k) for Poisson(λ)."""
        # To avoid underflow for large λ, accumulate tail from k to λ+10√λ
        # but cap at reasonable bound
        max_k = int(lmbda + 10 * math.sqrt(lmbda)) + 50
        tail = 0.0
        term = math.exp(-lmbda)  # P(X=0)
        for i in range(1, max_k + 1):
            term *= lmbda / i
            if i >= k:
                tail += term
        return tail

    def find_cutoff(cov: int, error_rate: float, target_p: float):
        """Find minimal count where Poisson tail < target_p."""
        lmbda = cov * error_rate
        for r in range(1, cov + 1):
            if poisson_tail(lmbda, r) < target_p:
                return r
        return cov

    # independent cutoffs
    SR_cutoff = find_cutoff(SR_cov, error_rate, target_p)

    # combined cutoff using total coverage
    total_cov = SR_cov + PB_cov
    combined_total = find_cutoff(total_cov, error_rate, target_p)

    # distribute proportionally to coverage, ensure ≥1 of each
    combined_SR = max(1, round(SR_cov / total_cov * combined_total))
    combined_PB = max(1, combined_total - combined_SR)

    return dict(
        SR=SR_cutoff,
        combined_SR=combined_SR,
        combined_PB=combined_PB,
    )

#*******************************************************************************
# Parsing variants and storing minipileup counts
#*******************************************************************************
class SampleCounts:
    """Store counts for a sample.
    """
    def __init__(self, name: str):
        self.name = name
        self.REF_ADF = 0
        self.REF_ADR = 0
        self.ALT_ADF = 0
        self.ALT_ADR = 0

class OriginalVCF:
    """Store SNV ONLY from original VCF for lookup.
    """
    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.snvs = dict() # (chrom, pos, ref, alt) -> pysam.VariantRecord

        self.load_records()

    def load_records(self):
        """Load original VCF and populate snvs dict.
        Assumes a file normalized with bcftools norm -m -any.
        The function will only store mono allelic SNV (1bp REF and 1bp ALT).
        """
        with pysam.VariantFile(self.vcf_path) as vf:
            # for filter in self.filter_definitions:
            #     vf.header.add_line(filter)
            it = vf.fetch() if vf.index is not None else vf
            for record in it:
                if record.alts is None or len(record.alts) != 1:
                    continue
                REF = record.ref
                ALT = record.alts[0]
                if REF and ALT and len(REF) == 1 and len(ALT) == 1:
                    self.snvs[(record.chrom, record.pos, REF, ALT)] = record

class MinipileupVCF:
    """Store counts from minipileup for lookup.
    ONLY stores counts for records in OriginalVCF.snvs dict.
    Counts are stored per sample as a SampleCounts object.
    Also aggregate counts by group: SR (Short-reads), PB (PacBio), ONT (Oxford Nanopore).
    """
    def __init__(self, vcf_path: str, original_vcf: OriginalVCF):
        self.vcf_path = vcf_path
        self.original_vcf = original_vcf
        self.counts = dict()  # (chrom, pos, ref, alt) -> {sample: SampleCounts, ...}
        self.aggregate_counts = dict() # (chrom, pos, ref, alt) -> {PB: SampleCounts, SR: SampleCounts, ONT: SampleCounts}

        self.load_records()
        self.aggregate_by_group()

    def select_alt_index(self, record: pysam.VariantRecord, ALT: str):
        """Given a pysam.VariantRecord and an ALT allele,
        return the index of the ALT allele in the record.alts tuple.
        Try to match exactly first, then by first base only.
        If not found, return None.
        """
        for i, alt in enumerate(record.alts):
            if alt == ALT: return i
        for i, alt in enumerate(record.alts):
            if alt[0] == ALT: return i
        return None

    def add_counts(self, record: pysam.VariantRecord, sample: str, ALT_index: int):
        """Add counts from a pysam.VariantRecord sample call to a SampleCounts object
        ALT_index: index of the ALT allele in record.alts
        If there are multiple alternative alleles, only the one at ALT_index is counted as ALT,
        all others are counted as reference.
        """
        call = record.samples[sample]
        sc = SampleCounts(sample)
        ADF = call.get("ADF")
        ADR = call.get("ADR")
        # Reference counts
        sc.REF_ADF += ADF[0]
        sc.REF_ADR += ADR[0]
        # Alternate counts
        for i, _ in enumerate(record.alts):
            if i == ALT_index:
                # This is the ALT allele we count as alternate
                sc.ALT_ADF += ADF[i + 1]
                sc.ALT_ADR += ADR[i + 1]
            else:
                # This is an alternate allele we count as reference
                sc.REF_ADF += ADF[i + 1]
                sc.REF_ADR += ADR[i + 1]
        return sc

    def load_records(self):
        """Load minipileup VCF and populate counts.
        """
        with pysam.VariantFile(self.vcf_path) as vf:
            it = vf.fetch() if vf.index is not None else vf
            for record in it:
                # Not sure if this can happen, but just in case - skip records with no ALT
                if record.alts is None: continue
                # For each ALT allele, check if in original_vcf.snvs
                for ALT in {alt[0] for alt in record.alts}:
                    # This is because minipileup adds extra bases sometimes to REF/ALT
                    REF = record.ref[0]
                    key = (record.chrom, record.pos, REF, ALT)
                    if key in self.original_vcf.snvs: # We store this record
                        # Select correct ALT index
                        i = self.select_alt_index(record, ALT)
                        if i is None:
                            sys.exit(f"ERROR: Could not find ALT {ALT} in record at {record.chrom}:{record.pos} with ALTs {record.alts}")
                        # Get counts for all samples
                        for sample in record.samples:
                            sc = self.add_counts(record, sample, i)
                            # Store counts
                            self.counts.setdefault(key, dict())[sample] = sc

    def aggregate_by_group(self):
        """Aggregate counts by group:
        SR (Short-reads), PB (PacBio), ONT (Oxford Nanopore).
        """
        for key, counts_ in self.counts.items():
            agg = dict(
                    SR=SampleCounts("SR"),
                    PB=SampleCounts("PB"),
                    ONT=SampleCounts("ONT")
            )
            for sample, sc in counts_.items():
                if sample.endswith("-SR"): group = "SR"
                elif sample.endswith("-PB"): group = "PB"
                elif sample.endswith("-ONT"): group = "ONT"
                else:
                    sys.exit(f"ERROR: Sample {sample} does not end with -SR, -PB, or -ONT. Cannot determine group")
                agg[group].REF_ADF += sc.REF_ADF
                agg[group].REF_ADR += sc.REF_ADR
                agg[group].ALT_ADF += sc.ALT_ADF
                agg[group].ALT_ADR += sc.ALT_ADR
            self.aggregate_counts[key] = agg

#*******************************************************************************
# Tiering and filtering variants
#*******************************************************************************
class TieredVCF:
    """Tier variants based on counts from MinipileupVCF and cutoffs.
    """
    def __init__(self, original_vcf: OriginalVCF, minipileup_vcf: MinipileupVCF):
        self.original_vcf = original_vcf
        self.minipileup_vcf = minipileup_vcf
        self.snvs = original_vcf.snvs  # (chrom, pos, ref, alt) -> pysam.VariantRecord
        self.tiers = dict()  # (chrom, pos, ref, alt) -> tier
        # Only store TIER1 and TIER2 variants, others are not in dict
        self.filter_definitions = [
            '##FILTER=<ID=TIER1,Description="Alt present in long-read group at or above threshold">',
            '##FILTER=<ID=TIER2,Description="Alt present in short-read group at or above threshold; absent in long-read">'
        ]

        self.tier_variants()

    def tier_variants(self):
        """Tier variants in self.snvs based on counts from minipileup
        for  SR (Short-reads), PB (PacBio) and read cutoffs.
        """
        for key, record in self.snvs.items():
            # Retrieve aggregated counts
            agg = self.minipileup_vcf.aggregate_counts.get(key)
            if agg is None:
                continue
            # Short-read counts
            SR_REF_ADF, SR_REF_ADR = agg["SR"].REF_ADF, agg["SR"].REF_ADR
            SR_ALT_ADF, SR_ALT_ADR = agg["SR"].ALT_ADF, agg["SR"].ALT_ADR
            # PacBio counts
            PB_REF_ADF, PB_REF_ADR = agg["PB"].REF_ADF, agg["PB"].REF_ADR
            PB_ALT_ADF, PB_ALT_ADR = agg["PB"].ALT_ADF, agg["PB"].ALT_ADR
            # Totals
            SR_ALT_TOTAL = SR_ALT_ADF + SR_ALT_ADR
            PB_ALT_TOTAL = PB_ALT_ADF + PB_ALT_ADR
            SR_TOTAL = SR_REF_ADF + SR_REF_ADR + SR_ALT_TOTAL
            PB_TOTAL = PB_REF_ADF + PB_REF_ADR + PB_ALT_TOTAL
            # Determine tier
            if SR_TOTAL != 0:
                thresholds = get_read_cutoffs(SR_TOTAL, PB_TOTAL)
                # Tier classification
                if SR_ALT_TOTAL >= thresholds["combined_SR"] and PB_ALT_TOTAL >= thresholds["combined_PB"]:
                    self.tiers[key] = "TIER1"
                elif SR_ALT_TOTAL >= thresholds["SR"]:
                    self.tiers[key] = "TIER2"

    # def chrom_order(self, chrom):
    #     chrom = chrom.replace("chr", "")
    #     if chrom == "X": return 23
    #     elif chrom == "Y": return 24
    #     elif chrom in ("M", "MT"): return 25
    #     else: return int(chrom) if chrom.isdigit() else 26

    # def write_tiered_vcf(self, out_vcf_path: str):
    #     """Write tiered VCF to out_vcf_path.
    #     """
    #     with pysam.VariantFile(self.original_vcf.vcf_path) as vf_in:
    #         header = vf_in.header.copy()
    #         for definition in self.filter_definitions:
    #             header.add_line(definition)
    #         with pysam.VariantFile(out_vcf_path, "w", header=header) as vf_out:
    #             for key in sorted(self.snvs, key=lambda k: (self.chrom_order(k[0]), k[1])):
    #                 record = self.snvs[key]
    #                 # Make record compatible with the writer header
    #                 record.translate(header)
    #                 # Reset to PASS, then add tier if present
    #                 record.filter.clear()
    #                 tier = self.tiers.get(key)
    #                 if tier:
    #                     record.filter.add(tier)
    #                 vf_out.write(record)
