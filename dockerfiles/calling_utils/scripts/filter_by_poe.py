#!/usr/bin/env python3
import argparse
import os
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Sequence, Tuple, Optional

import pysam

# ported to python by ChatGPT, originally BSMN bash script
# usage  ./filter_by_poe.py --vcf example_data/8_RUFUS-ST002-1D_vs_ST001_1D_BCM_150x_VEP_^Cltered.vcf.gz --fasta PON.q20q20.05.5.fa.gz --out new_test.vcf --threads 2 --failed-out failed.vcf

# IUPAC mapping for the PON base
IUPAC_CONTAINS: Dict[str, set] = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"}, "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    # N handled specially below; * handled specially below
}

_tls = threading.local()

def _get_fasta(path: str) -> pysam.FastaFile:
    fa = getattr(_tls, "fasta_handle", None)
    if fa is None:
        _tls.fasta_handle = pysam.FastaFile(path)
        fa = _tls.fasta_handle
    return fa

def decider(pon_base: str, alts: Sequence[str]) -> bool:
    """
    Return True if the variant PASSES (ALT NOT present in PON),
    False if it FAILS (ALT is present in PON or unsupported).
    Matches the original awk semantics for SNVs; non-single-base ALTs are ignored.

    Rules:
      *  -> Pass
      N  -> Fail
      A/C/G/T/R/Y/S/W/K/M/B/D/H/V -> Fail iff ANY ALT is a single base within the code's set.
      Unknown/missing PON base -> Fail
    """
    if pon_base == "*":
        return True
    if pon_base == "N" or len(pon_base) != 1:
        return False

    allowed = IUPAC_CONTAINS.get(pon_base)
    if allowed is None:
        return False

    for alt in alts or ():
        if alt is None:
            continue
        a = alt.upper()
        if len(a) == 1 and a in {"A", "C", "G", "T"} and a in allowed:
            return False  # ALT present in PON -> Fail
        # non-single-base ALTs are ignored (do not cause Fail)
    return True

def fetch_pon_base(fasta_path: str, chrom: str, pos_1based: int) -> Optional[str]:
    try:
        fa = _get_fasta(fasta_path)
        seq = fa.fetch(chrom, pos_1based - 1, pos_1based)
        if not seq:
            return None
        return seq.upper()
    except Exception:
        return None

def choose_mode_from_ext(path: str) -> str:
    return "wz" if path.endswith(".gz") else "w"

def process_record(
    idx: int,
    rec: pysam.VariantRecord,
    fasta_path: str,
    chr_regex: Optional[re.Pattern],
) -> Tuple[int, bool]:
    chrom = rec.chrom
    pos = rec.pos  # 1-based

    if chr_regex is not None and not chr_regex.match(chrom):
        return idx, False  # mimic P="?" -> default Fail

    pon_base = fetch_pon_base(fasta_path, chrom, pos)
    if pon_base is None:
        return idx, False

    alts = rec.alts or ()
    passes = decider(pon_base, alts)
    return idx, passes

def main():
    ap = argparse.ArgumentParser(description="Filter VCF by PON FASTA (IUPAC-aware) in one pass.")
    ap.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ/BCF (indexed if compressed)")
    ap.add_argument("--fasta", required=True, help="PON FASTA (bgz + .fai; .gzi if bgz)")
    ap.add_argument("--out", required=True, help="Output VCF of PASS variants (.vcf or .vcf.gz)")
    ap.add_argument("--failed-out", default=None,
                    help="Optional VCF path to also write FAIL variants.")
    ap.add_argument("--threads", type=int, default=max(os.cpu_count() or 1, 1),
                    help="Worker threads for FASTA lookups (default: CPU count)")
    ap.add_argument("--chr-regex", default=r"^chr([0-9]+|[XY])\b",
                    help="Regex contigs must match to be checked. Set '' to check all.")
    args = ap.parse_args()

    if args.failed_out and os.path.abspath(args.failed_out) == os.path.abspath(args.out):
        sys.stderr.write("[ERROR] --out and --failed-out must be different files.\n")
        sys.exit(2)

    chr_pat = re.compile(args.chr_regex) if args.chr_regex else None

    # Open input
    try:
        invcf = pysam.VariantFile(args.vcf)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open VCF '{args.vcf}': {e}\n")
        sys.exit(2)

    # Output(s)
    try:
        outvcf = pysam.VariantFile(args.out, choose_mode_from_ext(args.out), header=invcf.header)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open output '{args.out}': {e}\n")
        sys.exit(2)

    failvcf = None
    if args.failed_out:
        try:
            failvcf = pysam.VariantFile(
                args.failed_out, choose_mode_from_ext(args.failed_out), header=invcf.header
            )
        except Exception as e:
            sys.stderr.write(f"[ERROR] Failed to open failed-out '{args.failed_out}': {e}\n")
            sys.exit(2)

    # Warm up FASTA (fail fast if index missing)
    try:
        _get_fasta(args.fasta)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open FASTA '{args.fasta}': {e}\n")
        sys.exit(2)

    total = passed = failed = 0
    idx = 0
    pending = {}
    next_to_write = 0

    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futures = []
        for rec in invcf.fetch():
            my_idx = idx
            idx += 1
            fut = ex.submit(process_record, my_idx, rec, args.fasta, chr_pat)
            futures.append((my_idx, rec, fut))

        for i, rec, fut in futures:
            ok = fut.result()[1]
            pending[i] = (rec, ok)
            while next_to_write in pending:
                r, ok2 = pending.pop(next_to_write)
                total += 1
                if ok2:
                    outvcf.write(r)
                    passed += 1
                else:
                    if failvcf is not None:
                        failvcf.write(r)
                    failed += 1
                next_to_write += 1

    outvcf.close()
    if failvcf is not None:
        failvcf.close()
    invcf.close()

    sys.stderr.write(
        "[pon_filter_vcf] Done.\n"
        f"  Input variants: {total}\n"
        f"  Passed:         {passed}\n"
        f"  Failed:         {failed}\n"
        f"  PASS VCF:       {args.out}\n"
        + (f"  FAIL VCF:       {args.failed_out}\n" if args.failed_out else "")
    )

if __name__ == "__main__":
    main()

