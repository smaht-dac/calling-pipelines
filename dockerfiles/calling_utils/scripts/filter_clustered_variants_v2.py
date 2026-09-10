#!/usr/bin/env python3
import argparse
import pysam


# Long-read caller names, mirroring the default of --lr-callers in
# tier_filter_variants_SR_PB_ONT_v2.py and params.lr_callers in main.nf.
# Keep the three in sync if a long-read caller is ever added.
#
# Why this is needed: merge_callers_v2.py gives a core ONE core ID regardless of
# sequencing technology and unions its caller sets (merge_callers_v2.py:100), so
# a core sequenced both short- and long-read collapses into a single key. Its
# short-read and long-read calls would then cluster against each other and both
# be dropped, even though they are independent observations.
#
# This is the same reconciliation split_shared_core_calls() does downstream
# (tier_filter_variants_SR_PB_ONT_v2.py:346-372), which splits a shared core into
# {core}_SR / {core}_PB. We use a "|LR" suffix rather than "_PB" because '_' is
# legal in a core ID -- O2's main.nf already emits real cores named {core}_SR /
# {core}_PB -- whereas validate_core_id() (merge_callers_v2.py:186-196) forbids
# '|', so "<core>|LR" can never collide with a real core ID. Splitting on the
# caller list alone (rather than on a shared-core list, which this script has no
# way to obtain) is equivalent for clustering: a core with no long-read caller is
# untouched, and a long-read-only core simply has all of its calls renamed to the
# one key. The key is internal only and never written to the output VCF.
LR_CALLERS = {"longcalld"}   # compared case-insensitively


def _get_cores(rec):
    """Return set of core keys from CORE_CALLS INFO field.

    CORE_CALLS format: core1:caller1,caller2|core2:caller1

    A core's long-read calls are keyed as "<core>|LR" so they are not clustered
    against that same core's short-read calls. A core with both kinds of caller
    yields both keys.
    """
    val = rec.info.get("CORE_CALLS")
    if not val:
        return set()
    # NOTE: CORE_CALLS is declared Number=1 but contains commas, so pysam hands
    # back a TUPLE for any core listing 2+ callers and str() then yields a Python
    # repr, producing a garbage core ID. That is Open Issue #1 in
    # HANDOFF_o2_to_aws_v2.md and is deliberately left unfixed on tag v2
    # (README_long_claude.md, "Ported as-is (do not fix)"); the fix lives in
    # a27efa2 on branch core-specific. Left as-is here on purpose -- do not add a
    # rejoin guard without also revisiting that decision. Consequence for the
    # split below: multi-caller cores never reach it, so only single-caller
    # entries (the large majority) get the |LR treatment.
    cores = set()
    for entry in str(val).split("|"):
        if ":" not in entry:
            continue
        core, _, caller_list = entry.partition(":")
        core = core.strip()
        callers = {c.strip().lower() for c in caller_list.split(",") if c.strip()}
        if callers & LR_CALLERS:
            cores.add(core + "|LR")
        # An empty caller list keeps the plain core key, matching prior behavior;
        # without this the variant would get no key at all and be filtered out.
        if not callers or (callers - LR_CALLERS):
            cores.add(core)
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
