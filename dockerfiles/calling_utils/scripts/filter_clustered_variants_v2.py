#!/usr/bin/env python3
import argparse
import pysam


def _get_cores(rec):
    """Return set of core IDs from CORE_CALLS INFO field."""
    val = rec.info.get("CORE_CALLS")
    if not val:
        return set()
    cores = set()
    for entry in str(val).split("|"):
        if ":" in entry:
            cores.add(entry.split(":", 1)[0].strip())
    return cores


def _sliding_window_singletons(records_with_idx, window):
    """
    Run sliding-window clustering on a list of (orig_idx, record) pairs
    (already in coordinate order). Returns set of orig_idx values that are
    singletons (not clustered within this core's variant set).
    """
    singleton_idxs = set()
    group = []
    group_clustered = False
    prev_rec = None

    for orig_idx, rec in records_with_idx:
        if (prev_rec is None
                or rec.chrom != prev_rec.chrom
                or rec.pos - prev_rec.pos > window):
            if group and not group_clustered:
                singleton_idxs.add(group[0])
            group = [orig_idx]
            group_clustered = False
        else:
            group.append(orig_idx)
            group_clustered = True
        prev_rec = rec

    if group and not group_clustered:
        singleton_idxs.add(group[0])

    return singleton_idxs


def filter_clustered(input_vcf, output_vcf, window=50):
    vcf_in = pysam.VariantFile(input_vcf)
    header = vcf_in.header
    all_records = list(vcf_in.fetch())
    vcf_in.close()

    # Build per-core index: core -> [(orig_idx, record), ...]
    core_index = {}
    for i, rec in enumerate(all_records):
        for core in _get_cores(rec):
            core_index.setdefault(core, []).append((i, rec))

    if not core_index:
        # No CORE_CALLS — fall back to tissue-level sliding window
        passes = set()
        group, group_clustered, prev = [], False, None
        for i, rec in enumerate(all_records):
            if prev is None or rec.chrom != prev.chrom or rec.pos - prev.pos > window:
                if group and not group_clustered:
                    passes.add(group[0])
                group, group_clustered = [i], False
            else:
                group.append(i)
                group_clustered = True
            prev = rec
        if group and not group_clustered:
            passes.add(group[0])
    else:
        # A variant passes if it is a singleton in at least one core's variant set.
        # The record is written unchanged (CORE_CALLS preserved), so downstream
        # steps treat it as passing for all cores.
        passes = set()
        for core, records_with_idx in core_index.items():
            passes |= _sliding_window_singletons(records_with_idx, window)

    vcf_out = pysam.VariantFile(output_vcf, "w", header=header)
    for i, rec in enumerate(all_records):
        if i in passes:
            vcf_out.write(rec)
    vcf_out.close()

    if output_vcf.endswith(".gz"):
        pysam.tabix_index(output_vcf, preset="vcf", force=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove clustered variants (SNVs + indels) within N bp")
    parser.add_argument("input_vcf", help="Input VCF file (bgzipped or not)")
    parser.add_argument("output_vcf", help="Output VCF file (will be bgzipped if ends with .gz)")
    parser.add_argument("--window", type=int, default=50, help="Window size in bp [default=50]")
    args = parser.parse_args()

    filter_clustered(args.input_vcf, args.output_vcf, args.window)
