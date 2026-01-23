#!/usr/bin/env python3
import argparse
import pysam
from collections import defaultdict
import math

############################################################
# Helper functions (from parse_minipileup-style logic)
############################################################
def get_short_read_cutoff(short_cov: int,
                          error_rate: float = 0.001,
                          target_p: float = 1e-5):
    """
    Compute a read count cutoff for short-read coverage using a
    Poisson approximation to the binomial distribution.

    | Short-read depth | Typical cutoff |
    | ---------------- | -------------- |
    | 30×              | 2–3 reads      |
    | 100×             | 3–4 reads      |
    | 300×             | 4–6 reads      |
    | 1000×            | 7–9 reads      |


    Parameters
    ----------
    short_cov : int
        Short-read coverage
    error_rate : float, optional
        Sequencing error rate (default = 0.001)
    target_p : float, optional
        Target probability for random error (default = 1e-2)

    Returns
    -------
    int
        Minimum supporting read count such that
        P(error >= cutoff) < target_p
    """

    def poisson_tail(lmbda, k):
        """Compute P(X >= k) for Poisson(λ)."""
        max_k = int(lmbda + 10 * math.sqrt(lmbda)) + 50
        tail = 0.0
        term = math.exp(-lmbda)  # P(X = 0)
        for i in range(1, max_k + 1):
            term *= lmbda / i
            if i >= k:
                tail += term
        return tail

    lam = short_cov * error_rate

    for r in range(1, short_cov + 1):
        if poisson_tail(lam, r) < target_p:
            return r

    return short_cov




def first_base(allele: str):
    if allele is None:
        return None
    return allele[0] if len(allele) > 0 else None


def choose_alt_index(alleles, target_alt_base: str):
    """
    Select the ALT index matching the original SNV ALT base, applying special
    rules when REF/ALT in minipileup are padded (multi-base).
    """
    if target_alt_base is None or len(alleles) <= 1:
        return None

    ref_len = len(alleles[0])

    # 1. First: padded-case rule -> same-length ALT with same first base
    padded_matches = [
        i for i in range(1, len(alleles))
        if first_base(alleles[i]) == target_alt_base
        and len(alleles[i]) == ref_len
    ]
    if padded_matches:
        return padded_matches[0]  # exact padded match wins

    # 2. Otherwise fall back to single-base / prefix matching
    candidates = [
        i for i in range(1, len(alleles))
        if first_base(alleles[i]) == target_alt_base
    ]
    if not candidates:
        return None

    # Prefer single-base alt if present
    singles = [i for i in candidates if len(alleles[i]) == 1]
    return singles[0] if singles else candidates[0]



def get_numberR(sample, key):
    """
    Safely get Number=R array (ADF or ADR) from a pysam Call.
    Returns list of ints or None.
    """
    val = sample.get(key, None)
    if val is None:
        return None
    try:
        return list(val)
    except Exception:
        return None

def load_alt_base_map(original_vcf_path):
    """
    From the original VCF, build a map:
        (chrom, pos) -> alt_base (first base of ALT)
    Restrict to SNVs (single-base ref/alt) to mirror parse_minipileup.
    """
    alt_map = {}
    with pysam.VariantFile(original_vcf_path) as vf:
        it = vf.fetch() if vf.index is not None else vf
        for rec in it:
            if rec.alts is None or len(rec.alts) == 0:
                continue
            ref = rec.ref
            alt = rec.alts[0]
            if ref and alt and len(ref) == 1 and len(alt) == 1:
                alt_map[(rec.chrom, rec.pos)] = first_base(alt)
    return alt_map

############################################################
# Parse arguments
############################################################
parser = argparse.ArgumentParser(description="Annotate VCF using Minipileup results")
parser.add_argument("--tissue", required=True, help="Current tissue (e.g. 3A or SMHT004-3A)")
parser.add_argument("--orig_vcf", required=True, help="Original input VCF (bgzipped)")
parser.add_argument("--mp_vcf", required=True, help="Minipileup VCF (bgzipped)")
parser.add_argument("--out", required=True, help="Output annotated VCF file")
args = parser.parse_args()

current_tissue = args.tissue
# Normalize current tissue to the short code (e.g. '3A')
current_tissue_short = current_tissue.split("-")[-1]

############################################################
# 0. Build map from (chrom,pos) -> original ALT base
############################################################
alt_base_map = load_alt_base_map(args.orig_vcf)

############################################################
# 1. Load mp_vcf and compute sample → tissue VAFs
############################################################
mp = pysam.VariantFile(args.mp_vcf)

# Sample name -> tissue label (e.g. SMHT004-3A -> 3A)
sample2tissue = {s: s.split("-")[-1] for s in mp.header.samples}

# Store VAFs for each sample, then aggregate per tissue
# key = (chrom,pos,tissue) → list of (alt_count, depth)
tissue_vafs = defaultdict(list)

for rec in mp:
    chrom = rec.chrom
    pos = rec.pos
    key_pos = (chrom, pos)

    # Only process sites that exist in the original VCF (and are SNVs)
    target_alt_base = alt_base_map.get(key_pos)
    if target_alt_base is None:
        continue

    alleles = [rec.ref] + list(rec.alts or [])
    alt_idx = choose_alt_index(alleles, target_alt_base)
    if alt_idx is None:
        # No ALT index in this record that matches the original ALT
        continue

    for sample in rec.samples:
        call = rec.samples[sample]

        adf = get_numberR(call, "ADF")
        adr = get_numberR(call, "ADR")
        if adf is None or adr is None:
            continue

        # ALT count = chosen ALT index
        alt_count = 0
        if alt_idx < len(adf):
            alt_count += adf[alt_idx] or 0
        if alt_idx < len(adr):
            alt_count += adr[alt_idx] or 0

        # REF count = sum of all alleles except the chosen ALT
        ref_count = 0
        for i, v in enumerate(adf):
            if i != alt_idx and v is not None:
                ref_count += v
        for i, v in enumerate(adr):
            if i != alt_idx and v is not None:
                ref_count += v

        depth = ref_count + alt_count
        if depth == 0:
            continue

        vaf = alt_count / depth
        tissue = sample2tissue[sample].strip("[").strip("]").strip(" ").strip(",")

        tissue_vafs[(chrom, pos, tissue)].append((alt_count, depth))

############################################################
# 2. Aggregate per tissue and REMOVE 0-VAF tissues
############################################################
aggregated_vaf = defaultdict(dict)   # (chrom,pos) → {tissue: vaf}

for (chrom, pos, tissue), vafs in tissue_vafs.items():
    total_alt = 0
    total_depth = 0
    for (alt_count, depth) in vafs:
        total_alt += alt_count
        total_depth += depth

    if total_depth == 0:
        continue

    if total_alt >= get_short_read_cutoff(total_depth):
        mean_vaf = total_alt / total_depth
        aggregated_vaf[(chrom, pos)][tissue] = mean_vaf

############################################################
# 3. Build TISSUE_SR_VAFS and CrossTissue annotations
############################################################
summary = {}           # (chrom,pos) → "3A,0.05|3I,0.2"
crosstissue_flag = {}  # variants shared >1 tissue

for key_pos, tissue_dict in aggregated_vaf.items():
    chrom, pos = key_pos

    # skip if no tissues remain after filtering
    if len(tissue_dict) == 0:
        continue

    # Build tissue strings
    parts = [f"{t}:{vaf:.6f}" for t, vaf in tissue_dict.items()]
    summary[key_pos] = "|".join(parts)

    # CrossTissue if more than one tissue has VAF
    if len(tissue_dict) > 1:
        crosstissue_flag[key_pos] = True

############################################################
# 4. Annotate the original VCF
############################################################
orig = pysam.VariantFile(args.orig_vcf)

orig.header.info.add("SR_VAF", number=1, type="Float",
                     description="VAF for short read in current tissue")
orig.header.info.add("POOLED_PB_VAF", number=1, type="Float",
                     description="VAF for PacBio in current donor pooled tissues")
orig.header.info.add("POOLED_ONT_VAF", number=1, type="Float",
                     description="VAF for ONT in current donor pooled tissues")
orig.header.info.add("TISSUE_SR_VAFS", number=".", type="String",
                     description="VAFs for all tissues with short read nonzero VAF based on pileups, reads with BQ>30")
orig.header.info.add("CrossTissue", number=0, type="Flag",
                     description="Variant has VAF > 0 in another tissue")

out = pysam.VariantFile(args.out, "w", header=orig.header)

############################################################
# 5. Write output with annotations
############################################################
for rec in orig:
    key_pos = (rec.chrom, rec.pos)

    # Add current tissue VAF only if >0 for that tissue
    vaf_dict = aggregated_vaf.get(key_pos, {})

    try:
        sr_adf = rec.info.get("SR_ADF")
        sr_adr = rec.info.get("SR_ADR")
        if sr_adf and sr_adr and len(sr_adf) > 1 and len(sr_adr) > 1:
            ref = sr_adf[0] + sr_adr[0]
            alt = sr_adf[1] + sr_adr[1]
            if ref + alt > 0:
                rec.info["SR_VAF"] = float(alt / (ref + alt))
    except Exception:
        rec.info["SR_VAF"] = 0.0
        pass

    # Add TISSUE_SR_VAFS if any
    if key_pos in summary:
        rec.info["TISSUE_SR_VAFS"] = summary[key_pos]

    # Add CrossTissue flag
    if key_pos in crosstissue_flag:
        rec.info["CrossTissue"] = True

    ############################################################
    # POOLED_PB_VAF (using PB_ADF/PB_ADR, assuming Number=R [ref, alt...])
    # This still assumes index 1 is the mosaic ALT as before.
    # If needed, we can extend alt-index logic here similarly.
    ############################################################
    try:
        lr_adf = rec.info.get("PB_ADF")
        lr_adr = rec.info.get("PB_ADR")
        if lr_adf and lr_adr and len(lr_adf) > 1 and len(lr_adr) > 1:
            ref = lr_adf[0] + lr_adr[0]
            alt = lr_adf[1] + lr_adr[1]
            if ref + alt > 0:
                rec.info["POOLED_PB_VAF"] = float(alt / (ref + alt))
    except Exception:
        pass

    ############################################################
    # POOLED_ONT_VAF (using ONT_ADF/ONT_ADR, assuming Number=R [ref, alt...])
    ############################################################
    try:
        ont_adf = rec.info.get("ONT_ADF")
        ont_adr = rec.info.get("ONT_ADR")
        if ont_adf and ont_adr and len(ont_adf) > 1 and len(ont_adr) > 1:
            ref = ont_adf[0] + ont_adr[0]
            alt = ont_adf[1] + ont_adr[1]
            if ref + alt > 0:
                rec.info["POOLED_ONT_VAF"] = float(alt / (ref + alt))
    except Exception:
        pass

    out.write(rec)

out.close()

print(f"Annotated VCF written to {args.out}")

