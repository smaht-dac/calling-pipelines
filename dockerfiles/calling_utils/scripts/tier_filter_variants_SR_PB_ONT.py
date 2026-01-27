#!/usr/bin/env python3

################################################################################
### Libraries
################################################################################
import argparse, sys
import pysam
import math
import re
from datetime import datetime
from scipy.stats import fisher_exact

# Scipy compatibility for binomtest/binom_test
try:
    from scipy.stats import binomtest as _binomtest
    def binom_pvalue(k, n, p, alternative="two-sided"):
        import re
        return _binomtest(k, n, p, alternative=alternative).pvalue
except Exception:
    from scipy.stats import binom_test as _binomtest
    def binom_pvalue(k, n, p, alternative="two-sided"):
        return _binomtest(k, n, p, alternative=alternative)

################################################################################
### Objects
################################################################################
#*******************************************************************************
# Statistical functions and filtering logic
#*******************************************************************************
def get_read_cutoffs(SR_cov: int, PB_cov: int,
                     error_rate: float = 0.001,
                     target_p: float = 1e-2):
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
        Target probability for random error (default = 1e-2)

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

    # independent cutoffs (hard code strict cutoff!)
    SR_cutoff = find_cutoff(SR_cov, error_rate, 1e-5)

    # combined cutoff using total coverage
    total_cov = SR_cov + PB_cov
    combined_total = find_cutoff(total_cov, error_rate, target_p)

    # distribute proportionally to coverage, ensure ≥1 of each
    combined_SR = max(1, round(SR_cov / total_cov * float(combined_total)))
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
    def __init__(self, vcf_path: str, original_vcf: OriginalVCF, current_tissue: str):
        self.vcf_path = vcf_path
        self.original_vcf = original_vcf
        self.current_tissue = current_tissue

        self.counts = dict()  # (chrom, pos, ref, alt) -> {sample: SampleCounts, ...}
        self.aggregate_counts = dict() # (chrom, pos, ref, alt) -> {PB: SampleCounts, SR: SampleCounts, ONT: SampleCounts}\        self.tissue_pb_counts = dict()  # (chrom,pos,ref,alt) -> SampleCounts("PB_TISSUE")
        self.tissue_pb_counts = dict()  # (chrom,pos,ref,alt) -> SampleCounts("PB_TISSUE")
        self.tissue_ont_counts = dict()  # (chrom,pos,ref,alt) -> SampleCounts("ONT_TISSUE")

        self.load_records()
        self.aggregate_by_group()

    def select_alt_index(self, record: pysam.VariantRecord, ALT: str):
        """Given a pysam.VariantRecord and an ALT allele,
        return the index of the ALT allele in the record.alts tuple.
        Try to match exactly first, then by first base only.
        If not found, return None.

        Check length of ref to make sure alt match is accurate. Example:
            orig vcf:
                G>A
            minipileup vcf:
                GNN -> AN, ANN, A
            Should select second alt allele
        """
        l_ref = len(record.ref)
        for i, alt in enumerate(record.alts):
            if alt == ALT: return i
        for i, alt in enumerate(record.alts):
            l_alt = len(alt)
            if l_ref == l_alt:
                alt_no_N = alt.strip('N') 
                # takes care of cases with TNN > ACN,ANN
                if alt_no_N == ALT: return i
        return None

    def add_counts(self, record: pysam.VariantRecord, sample: str, ALT_index: int):
        """Add counts from a pysam.VariantRecord sample call to a SampleCounts object
        ALT_index: index of the ALT allele in record.alts
        If there are multiple alternative alleles, only the one at ALT_index is counted as ALT,
        all others are counted as reference.
        """
        call = record.samples[sample]
        sc = SampleCounts(sample)
        ADF = call.get("ADF"); ADR = call.get("ADR")
        if ADF is None or ADR is None:
            return sc
        if len(ADF) < len(record.alts) + 1 or len(ADR) < len(record.alts) + 1:
            return sc
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
                            # not found in minipileup alts, remove variant
                            print(f"ERROR: Could not find ALT {ALT} in record at {record.chrom}:{record.pos}")
                            continue
                        # Get counts for all samples
                        for sample in record.samples:
                            sc = self.add_counts(record, sample, i)
                            # Store counts
                            self.counts.setdefault(key, dict())[sample] = sc

    def aggregate_by_group(self):
        """
        Aggregate counts by group:
        SR (Short-reads), PB (PacBio), ONT (Oxford Nanopore).

        Also compute tissue-specific PB counts (PB_TISSUE) using only samples
        that end with: -PB-<current_tissue> when current_tissue is provided.
        """
        pb_pat = re.compile(r'-PB(?:-|$)')  # matches ...-PB or ...-PB-...
        ont_pat = re.compile(r'-ONT(?:-|$)')  # matches ...-ONT or ...-ONT-...

        for key, counts_ in self.counts.items():
            agg = dict(
                    SR=SampleCounts("SR"),
                    PB=SampleCounts("PB"),
                    ONT=SampleCounts("ONT")
            )

            tissue_pb = SampleCounts("PB_TISSUE")
            tissue_ont = SampleCounts("ONT_TISSUE")

            for sample, sc in counts_.items():
                if sample.endswith("-SR"): group = "SR"
                elif pb_pat.search(sample): group = "PB"
                elif ont_pat.search(sample): group = "ONT"
                else:
                    sys.exit(f"ERROR: Sample {sample} does not end with -SR, -PB, or -ONT. Cannot determine group")
                agg[group].REF_ADF += sc.REF_ADF
                agg[group].REF_ADR += sc.REF_ADR
                agg[group].ALT_ADF += sc.ALT_ADF
                agg[group].ALT_ADR += sc.ALT_ADR

                # Additionally aggregate tissue-matched PB only
                if (
                    group == "PB"
                    and self.current_tissue
                    and sample.endswith(f"-PB-{self.current_tissue}")
                ):
                    tissue_pb.REF_ADF += sc.REF_ADF
                    tissue_pb.REF_ADR += sc.REF_ADR
                    tissue_pb.ALT_ADF += sc.ALT_ADF
                    tissue_pb.ALT_ADR += sc.ALT_ADR

                if (
                    group == "ONT"
                    and self.current_tissue
                    and sample.endswith(f"-ONT-{self.current_tissue}")
                ):
                    tissue_ont.REF_ADF += sc.REF_ADF
                    tissue_ont.REF_ADR += sc.REF_ADR
                    tissue_ont.ALT_ADF += sc.ALT_ADF
                    tissue_ont.ALT_ADR += sc.ALT_ADR

            self.aggregate_counts[key] = agg
            self.tissue_pb_counts[key] = tissue_pb
            self.tissue_ont_counts[key] = tissue_ont


#*******************************************************************************
# Tiering and filtering variants
#*******************************************************************************
class FisherTestResult:
    """Store Fisher's exact test result.
    """
    def __init__(self, group: str, p_value: float):
        self.group = group
        self.p_value = p_value

    def is_pass(self, alpha: float):
        """Check if p-value passes threshold at alpha.
        """
        if self.p_value >= alpha: return True
        return False

class BinomialTestResult:
    """Store binomial test result.
    """
    def __init__(self):
        self.p_value_SR = None
        self.p_value_PB = None
        self.p_value_ONT = None

    def is_pass(self, group: str, alpha: float):
        """Check if p-value passes threshold at alpha for read type group.
        Return None if test was not run (p_value is None).
        """
        if group == "SR":
            if self.p_value_SR is None: return None
            if self.p_value_SR <= alpha: return True
        if group == "PB":
            if self.p_value_PB is None: return None
            if self.p_value_PB <= alpha: return True
        if group == "ONT":
            if self.p_value_ONT is None: return None
            if self.p_value_ONT <= alpha: return True
        return False

class TieredVCF:
    """Tier variants based on counts from MinipileupVCF and cutoffs.
    """
    def __init__(self, original_vcf: OriginalVCF, minipileup_vcf: MinipileupVCF,
                 strand_alpha: float, germline_alpha: float, germline_alpha_SR: float, min_alt_PB: int, min_alt_binom: int):
        self.original_vcf = original_vcf
        self.minipileup_vcf = minipileup_vcf
        self.strand_alpha = strand_alpha
        self.germline_alpha = germline_alpha
        self.germline_alpha_SR = germline_alpha_SR
        self.min_alt_PB = min_alt_PB
        self.min_alt_binom = min_alt_binom
        self.snvs = original_vcf.snvs  # (chrom, pos, ref, alt) -> pysam.VariantRecord
        self.tiers = dict()  # (chrom, pos, ref, alt) -> tier
        self.tests = dict()  # (chrom, pos, ref, alt) -> {fisher: FisherTestResult, binomial: BinomialTestResult}
        # Only store TIER1 and TIER2 variants, others are not in dict
        self.definitions = [
            # CrossTech classification
            '##INFO=<ID=CrossTech,Number=0,Type=Flag,Description="Alt supported in both short read and PacBio data at or above their combined thresholds">',
            '##INFO=<ID=CrossCaller,Number=0,Type=Flag,Description="Alt found in more than one variant caller">',
            # CALLERS
            '##INFO=<ID=CALLERS,Number=.,Type=String,Description="List of variant callers that reported this variant">',
            # Fisher strand bias
            '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher exact test p-value for strand balance on the selected platform">',
            '##INFO=<ID=SB_SRC,Number=1,Type=String,Description="Platform used for Fisher strand test: SR (short read) or PB (PacBio long-read)">',
            # Binomial germline deviation
            '##INFO=<ID=GERMLINE_PVAL,Number=1,Type=Float,Description="Minimum binomial p-value for germline deviation across all platforms tested">',
            '##INFO=<ID=GERMLINE_PVAL_SR,Number=1,Type=Float,Description="Binomial p-value for germline deviation in short read data">',
            '##INFO=<ID=GERMLINE_PVAL_PB,Number=1,Type=Float,Description="Binomial p-value for germline deviation in PacBio (long-read) data">',
            '##INFO=<ID=GERMLINE_PVAL_ONT,Number=1,Type=Float,Description="Binomial p-value for germline deviation in Oxford Nanopore (long-read) data">',
            # Raw strand-specific counts
            '##INFO=<ID=SR_ADF,Number=2,Type=Integer,Description="Short-read forward depths (REF, ALT)">',
            '##INFO=<ID=SR_ADR,Number=2,Type=Integer,Description="Short-read reverse depths (REF, ALT)">',
            '##INFO=<ID=PB_ADF,Number=2,Type=Integer,Description="PacBio (long-read) forward depths (REF, ALT)">',
            '##INFO=<ID=PB_ADR,Number=2,Type=Integer,Description="PacBio (long-read) reverse depths (REF, ALT)">',
            '##INFO=<ID=ONT_ADF,Number=2,Type=Integer,Description="Oxford Nanopore (long-read) forward depths (REF, ALT)">',
            '##INFO=<ID=ONT_ADR,Number=2,Type=Integer,Description="Oxford Nanopore (long-read) reverse depths (REF, ALT)">',
            '##INFO=<ID=TISSUE_PB_VAF,Number=1,Type=Float,Description="PacBio VAF computed using only tissue-matched PB sample(s)">',
            '##INFO=<ID=TISSUE_ONT_VAF,Number=1,Type=Float,Description="ONT VAF computed using only tissue-matched ONT sample(s)">'
        ]

        self.filter_variants()

    def tier_variant(self, key: tuple):
        """Tier variants in self.snvs based on counts from minipileup
        for SR (Short-reads), PB (PacBio) and read cutoffs.
        """
        agg = self.minipileup_vcf.aggregate_counts[key]
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

    def fisher_strand_bias(self, key: tuple):
        """Compute Fisher's exact test p-value for strand bias.
        Use counts from minipileup for SR (Short-reads) or PB (PacBio) depending on total ALT reads.
        If enough PB ALT reads, use PB counts; otherwise use SR counts.
        The threshold for "enough" PB ALT reads is min_alt_PB.
        """
        agg = self.minipileup_vcf.aggregate_counts[key]
        # Short-read counts
        SR_REF_ADF, SR_REF_ADR = agg["SR"].REF_ADF, agg["SR"].REF_ADR
        SR_ALT_ADF, SR_ALT_ADR = agg["SR"].ALT_ADF, agg["SR"].ALT_ADR
        # PacBio counts
        PB_REF_ADF, PB_REF_ADR = agg["PB"].REF_ADF, agg["PB"].REF_ADR
        PB_ALT_ADF, PB_ALT_ADR = agg["PB"].ALT_ADF, agg["PB"].ALT_ADR
        # Totals
        PB_ALT_TOTAL = PB_ALT_ADF + PB_ALT_ADR
        # Determine which counts to use
        if PB_ALT_TOTAL >= self.min_alt_PB:
            group = "PB"
            REF_ADF, REF_ADR = PB_REF_ADF, PB_REF_ADR
            ALT_ADF, ALT_ADR = PB_ALT_ADF, PB_ALT_ADR
        else:
            group = "SR"
            REF_ADF, REF_ADR = SR_REF_ADF, SR_REF_ADR
            ALT_ADF, ALT_ADR = SR_ALT_ADF, SR_ALT_ADR
        # Fisher's exact test
        table = [[REF_ADF, REF_ADR],
                 [ALT_ADF, ALT_ADR]]
        try:
            _, p_value = fisher_exact(table, alternative="two-sided")
        except Exception as e:
            p_value = 0.0  # On error, force an extreme low p so the variant FAILS Fisher (conservative).
        # Store result
        self.tests.setdefault(key, dict())
        self.tests[key]["fisher"] = FisherTestResult(group, p_value)

    def binomial_germline_deviation(self, key: tuple):
        """Compute binomial test p-value for germline deviation.
        Use counts from minipileup for SR (Short-reads), PB (PacBio), and ONT (Oxford Nanopore).
        Store p-values for each group.
        """
        agg = self.minipileup_vcf.aggregate_counts[key]
        # For each group, compute p-value
        binom_result = BinomialTestResult()
        for group in ["SR", "PB", "ONT"]:
            REF_ADF = agg[group].REF_ADF
            REF_ADR = agg[group].REF_ADR
            ALT_ADF = agg[group].ALT_ADF
            ALT_ADR = agg[group].ALT_ADR
            ALT_TOTAL = ALT_ADF + ALT_ADR
            REF_TOTAL = REF_ADF + REF_ADR
            TOTAL = REF_TOTAL + ALT_TOTAL
            if ALT_TOTAL >= self.min_alt_binom and TOTAL > 0:
                # Binomial test against 0.5 (heterozygous expectation)
                p_value = binom_pvalue(ALT_TOTAL, TOTAL, 0.5, alternative='less')
            else:
                p_value = None  # We do not have enough ALT reads to run the test
            # Store p-value
            if group == "SR":
                binom_result.p_value_SR = p_value
            elif group == "PB":
                binom_result.p_value_PB = p_value
            elif group == "ONT":
                binom_result.p_value_ONT = p_value
        # Store result
        self.tests.setdefault(key, dict())
        self.tests[key]["binomial"] = binom_result

    def filter_variants(self):
        """
        Filter variants in self.snvs based on counts from minipileup
        for SR (Short-reads), PB (PacBio) and read cutoffs.
        Filter is based on Fisher's exact test for strand bias and binomial test for germline deviation.
        Add tier information to FILTER.
        """
        for key, _ in self.snvs.items():
            if key not in self.minipileup_vcf.aggregate_counts:
                continue  # No counts available, skip
            self.tier_variant(key)
            self.fisher_strand_bias(key)
            self.binomial_germline_deviation(key)

    def chrom_order(self, chrom):
        chrom = chrom.replace("chr", "")
        if chrom == "X": return 23
        elif chrom == "Y": return 24
        elif chrom in ("M", "MT"): return 25
        else: return int(chrom) if chrom.isdigit() else 26

    def write_tiered_vcf(self, out_vcf_path: str, keep_info: bool=False):
        """Write tiered VCF to out_vcf_path.
        """
        no_pileup_counts, fail_filters = 0, 0
        written, t1, t2 = 0, 0, 0
        with pysam.VariantFile(self.original_vcf.vcf_path) as vf_in:

            header = vf_in.header.copy()

            for definition in self.definitions:
                header.add_line(definition)

            with pysam.VariantFile(out_vcf_path, "w", header=header) as vf_out:
                for key in sorted(self.snvs, key=lambda k: (self.chrom_order(k[0]), k[1])):
                    if key not in self.minipileup_vcf.aggregate_counts:
                        no_pileup_counts += 1
                        continue  # No counts available, skip
                    record = self.snvs[key]

                    # Extract CALLERS (if present) then remove all INFO fields from original vcf
                    callers_value = record.info.get("CALLERS")

                    if keep_info == False:
                        # Clear all INFO fields
                        record.info.clear()
  
                    # Make record compatible with the writer header
                    record.translate(header)
                    # Reset to PASS, then add tier if present
                    record.filter.clear()
                    tier = self.tiers.get(key)
                    if tier not in {"TIER1", "TIER2"}:
                        continue
                    # Consider fisher and binomial results
                    fisher_result = self.tests[key]["fisher"]
                    binom_result = self.tests[key]["binomial"]
                    fisher_pass = fisher_result.is_pass(self.strand_alpha)
                    # Binomial gating per tier
                    if tier == "TIER2":
                        # require SR only
                        binom_pass = binom_result.is_pass("SR", self.germline_alpha_SR)
                        binom_pass_long = True  # ignore PB/ONT for TIER2
                    else:  # TIER1
                        # ignore SR here to match spec
                        binom_pass_pb = binom_result.is_pass("PB", self.germline_alpha)
                        binom_pass_ont = binom_result.is_pass("ONT", self.germline_alpha)
                        # Only tests that *ran* matter; is_pass returns None when not run, which should not fail
                        binom_pass, binom_pass_long = True, True
                        if binom_pass_pb is False or binom_pass_ont is False:
                            binom_pass_long = False

                    if fisher_pass is False or binom_pass is False or binom_pass_long is False:
                      fail_filters += 1
                      continue  # Variant fails filters, do not write


                    # Add flag for CrossTech (old Tier1 classification) and CrossCaller
                    if tier == 'TIER1':
                        record.info['CrossTech'] = True
                    if len(callers_value) > 1:
                        record.info['CrossCaller'] = True

                    if keep_info == False:
                        record.info["CALLERS"] = callers_value
                    # Add raw counts to INFO fields
                    agg = self.minipileup_vcf.aggregate_counts[key]
                    record.info["SR_ADF"] = [agg["SR"].REF_ADF, agg["SR"].ALT_ADF]
                    record.info["SR_ADR"] = [agg["SR"].REF_ADR, agg["SR"].ALT_ADR]
                    record.info["PB_ADF"] = [agg["PB"].REF_ADF, agg["PB"].ALT_ADF]
                    record.info["PB_ADR"] = [agg["PB"].REF_ADR, agg["PB"].ALT_ADR]
                    record.info["ONT_ADF"] = [agg["ONT"].REF_ADF, agg["ONT"].ALT_ADF]
                    record.info["ONT_ADR"] = [agg["ONT"].REF_ADR, agg["ONT"].ALT_ADR]

                    # Add tissue-matched PB VAF if available
                    tpb = self.minipileup_vcf.tissue_pb_counts.get(key)
                    if tpb is not None:
                        tpb_ref = tpb.REF_ADF + tpb.REF_ADR
                        tpb_alt = tpb.ALT_ADF + tpb.ALT_ADR
                        tpb_total = tpb_ref + tpb_alt
                        if tpb_total > 0:
                            record.info["TISSUE_PB_VAF"] = float(tpb_alt) / float(tpb_total)

                    # Add tissue-matched ONT VAF if available
                    tont = self.minipileup_vcf.tissue_ont_counts.get(key)
                    if tont is not None:
                        tont_ref = tont.REF_ADF + tont.REF_ADR
                        tont_alt = tont.ALT_ADF + tont.ALT_ADR
                        tont_total = tont_ref + tont_alt
                        if tont_total > 0:
                            record.info["TISSUE_ONT_VAF"] = float(tont_alt) / float(tont_total)


                    # Add Fisher test results to INFO fields
                    record.info["SB_PVAL"] = fisher_result.p_value
                    record.info["SB_SRC"] = fisher_result.group
                    # Add Binomial test results to INFO fields
                    glm_pvals = []
                    if binom_result.p_value_SR is not None:
                        if tier == "TIER2":
                            record.info["GERMLINE_PVAL_SR"] = binom_result.p_value_SR
                            glm_pvals.append(binom_result.p_value_SR)
                    if binom_result.p_value_PB is not None:
                        if tier == "TIER1":
                            record.info["GERMLINE_PVAL_PB"] = binom_result.p_value_PB
                            glm_pvals.append(binom_result.p_value_PB)
                    if binom_result.p_value_ONT is not None:
                        if tier == "TIER1":
                            record.info["GERMLINE_PVAL_ONT"] = binom_result.p_value_ONT
                            glm_pvals.append(binom_result.p_value_ONT)
                    if glm_pvals:
                        record.info["GERMLINE_PVAL"] = min(glm_pvals)
                    # Write record
                    vf_out.write(record)
                    written += 1
                    t1 += (tier == "TIER1")
                    t2 += (tier == "TIER2")
        # Report
        print(f"INFO: Wrote tiered VCF to {out_vcf_path}.")
        print(f"REPORT: wrote {written} records (TIER1={t1}, TIER2={t2}).")
        print(f"REPORT: {no_pileup_counts} variants with no pileup counts.")
        print(f"REPORT: {fail_filters} variants failing filters.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Tier variants and filter by strand balance (Fisher) and germline deviation (binomial). Only SNVs are considered for tiering.")

    parser.add_argument("-i", "--input_vcf", required=True, help="Input VCF with variants to tier. Assumes normalized with bcftools norm -m -any. Compressed (.vcf.gz) or uncompressed (.vcf) VCF")
    parser.add_argument("-m", "--minipileup_vcf", required=True, help="Minipileup VCF with ADF/ADR counts. Compressed (.vcf.gz) or uncompressed (.vcf) VCF")
    parser.add_argument("-o", "--output_vcf", required=True, help="Output VCF with tiered and filtered variants. Compressed (.vcf.gz) or uncompressed (.vcf) VCF")

    parser.add_argument("--current_tissue",default=None,help="Tissue ID for this run (e.g. SMHT005-3AF). Used to compute TISSUE_PB_VAF from samples named *-PB-<current_tissue>.")

    parser.add_argument("--strand_alpha", type=float, default=0.01,
                    help="Keep if Fisher p >= this (default: 0.01)")
    parser.add_argument("--germline_alpha", type=float, default=0.01,
                    help="Keep if Binomial p <= this (default: 0.01)")
    parser.add_argument("--germline_alpha_SR", type=float, default=0.01,
                    help="Keep if Binomial p <= this. SR (Short-reads) only (default: 0.01)")
    parser.add_argument("--min_alt_PB", type=int, default=2,
                    help="Min ALT reads for PB (PacBio) to use PB counts for the strand test (default: 2)")
    parser.add_argument("--min_alt_binom", type=int, default=1,
                    help="Min ALT reads required to run the binomial test (default: 1)")
    parser.add_argument("--keep_info", action='store_true',
                    help="Keep existing INFO annotations")

    args = parser.parse_args()

    ovcf  = OriginalVCF(args.input_vcf)
    mpvcf = MinipileupVCF(args.minipileup_vcf, ovcf, current_tissue=args.current_tissue)
    tvcf  = TieredVCF(
        ovcf, mpvcf,
        strand_alpha=args.strand_alpha,
        germline_alpha=args.germline_alpha,
        germline_alpha_SR=args.germline_alpha_SR,
        min_alt_PB=args.min_alt_PB,
        min_alt_binom=args.min_alt_binom,
    )
    tvcf.write_tiered_vcf(args.output_vcf, args.keep_info)

    if args.output_vcf.endswith(".vcf.gz"):
        pysam.tabix_index(args.output_vcf, preset="vcf", force=True)
