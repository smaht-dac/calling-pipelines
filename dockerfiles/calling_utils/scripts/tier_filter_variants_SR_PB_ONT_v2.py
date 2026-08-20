#!/usr/bin/env python3

################################################################################
### Libraries
################################################################################
import argparse, sys, os
import pysam
import math
import re
from collections import OrderedDict
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
### Helpers
################################################################################

def parse_core_cram_map(path):
    """Parse core_cram_map.tsv: core \\t cram_basename \\t type (SR|PB|ONT).
    Returns OrderedDict {core: [(cram_basename, type), ...]} preserving row order.
    """
    mapping = OrderedDict()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            core, basename, ctype = parts[0], parts[1], parts[2]
            mapping.setdefault(core, []).append((basename, ctype))
    return mapping


def parse_core_calls(info_value):
    """Parse CORE_CALLS INFO string: '001C1:Strelka2,RUFUS|001A3:Strelka2'.
    Returns {core: [callers]} or {} if value is None/empty.
    """
    if not info_value:
        return {}
    result = {}
    for entry in str(info_value).split('|'):
        entry = entry.strip()
        if ':' not in entry:
            continue
        core, callers_str = entry.split(':', 1)
        result[core.strip()] = [c.strip() for c in callers_str.split(',') if c.strip()]
    return result


def _cram_basename(path):
    """Strip directories and a trailing .cram/.bam, matching minipileup's sample naming."""
    base = os.path.basename(path)
    for ext in ('.cram', '.bam'):
        if base.endswith(ext):
            return base[:-len(ext)]
    return base


def _core_parts(key):
    """Constituent set of a map key, ignoring any _SR/_PB shared-core suffix."""
    base = key[:-3] if key.endswith(('_SR', '_PB')) else key
    return set(base.split('-'))


def constituent_sources(core, core_cram_map, core_ctypes):
    """Choose the map keys to sum when reconstructing a call-only merged core.

    Splitting a merged name on every hyphen and looking each atom up individually only
    works when every constituent is its own map key. It fails as soon as a *pooled*
    CRAM covers several constituents at once: for 003A1-003A2-003B1 the atoms 003A1 and
    003A2 are not keys -- the key is the composite 003A1-003A2 -- so both lookups miss
    and the pooled library contributes nothing.

    So tile the constituent set with the largest available keys instead, skipping any
    key that overlaps one already chosen (which would double-count a constituent).
    Preference order within a size matches the original per-atom fallback: the bare
    name, then {name}_SR, then {name}_PB, with the shared-core halves admitted only for
    a platform this merged core actually has.

    Returns (chosen keys, constituents left uncovered).
    """
    want = _core_parts(core)

    def rank(key):
        if key.endswith('_SR'):
            return 1
        if key.endswith('_PB'):
            return 2
        return 0

    candidates = []
    for key in core_cram_map:
        if key == core:
            continue
        parts = _core_parts(key)
        if not parts <= want:
            continue
        if key.endswith('_SR') and 'SR' not in core_ctypes:
            continue
        if key.endswith('_PB') and 'PB' not in core_ctypes:
            continue
        candidates.append((-len(parts), rank(key), key, parts))
    candidates.sort()

    chosen, covered = [], set()
    for _, _, key, parts in candidates:
        if parts & covered:
            continue
        chosen.append(key)
        covered |= parts
    return chosen, want - covered


def canonical_core(name):
    """Canonical spelling of a (possibly merged) core name.

    Merged cores are hyphen-joined constituents, and the two sides of the pipeline
    do not agree on their order: CORE_CALLS carries them in plain sorted() order
    (the convention merge_callers_v2.py uses at its CORE_CALLS write site), while
    an O2 core_cram_map.tsv carries the raw VCF-samplesheet order. Where the two
    disagree the `core not in core_calls` gate never matches and the core's column
    is silently 0/0 on every record -- observed on O2 for 003D1-003A3
    (SMHT012-3AF) and 002B3-002B1 (SMHT042-3S).

    Sorting the constituents makes both sides agree by construction. Individual
    core names contain no hyphen, so they are returned unchanged.
    """
    if '-' not in name:
        return name
    # Preserve any _SR/_PB shared-core suffix; it is not part of the core name.
    suffix = ''
    if name.endswith('_SR') or name.endswith('_PB'):
        name, suffix = name[:-3], name[-3:]
    return '-'.join(sorted(name.split('-'))) + suffix


def parse_minipileup_sample(column):
    """Split a minipileup sample column into (core, stem, platform, tissue).

    Columns are named by minipileup2-parallel.sh as

        {core}__{stem}-SR
        {core}__{stem}-PB-{tissue}
        {core}__{stem}-ONT-{tissue}

    Core IDs contain no '_', so the first '__' is always the core/stem separator;
    a column with no '__' predates the core-aware naming and yields core None.
    Returns None if the platform suffix cannot be determined.
    """
    core, stem = None, column
    if '__' in column:
        core, stem = column.split('__', 1)
    if stem.endswith('-SR'):
        return core, stem[:-3], 'SR', None
    m = re.search(r'-(PB|ONT)-(.+)$', stem)
    if m:
        return core, stem[:m.start()], m.group(1), m.group(2)
    if stem.endswith('-PB') or stem.endswith('-ONT'):
        platform = 'PB' if stem.endswith('-PB') else 'ONT'
        return core, stem[:-(len(platform) + 1)], platform, None
    return None


def build_core_cram_map_from_samples(samples, current_tissue):
    """Build the core -> column mapping from the minipileup VCF's sample columns.

    This is the AWS counterpart to main.nf's core_cram_map generation. The core ID
    is stamped into each sample column by minipileup2-parallel.sh, so the mapping
    is recovered from the data rather than from a second set of CWL arrays that
    would have to be kept positionally in sync with the CRAM arrays.

    Reproduces main.nf exactly:
      - SR columns are already tissue-specific (the SR CRAM input is), so all are
        taken. PB columns are donor-pooled across every tissue, while O2's
        core_cram_map.tsv holds only current-tissue rows -- so PB columns are
        restricted to -PB-{current_tissue}. Without that filter the map would gain
        per-core keys for other tissues, including placeholder cores ('XX') and
        core IDs that recur across tissues ('001A1').
      - ONT carries no core in main.nf and is not added to the map. Other-tissue
        PB and all ONT columns still feed the donor-level pooled counts via
        aggregate_by_group(); only this per-core map is tissue-scoped.
      - A core with both an SR and a current-tissue PB column is 'shared'
        (main.nf's SHARED_CORES, keyed on tissue::core) and is split into
        {core}_SR / {core}_PB so SR and PacBio evidence are not collapsed onto one
        core name.

    A hyphenated core is allowed to own a CRAM. `own_cram_cores` reports which cores
    have a sample column of their own, so compute_core_counts() can tell a pooled
    library (real reads, use them) from a call-only merged core (no reads of its own,
    sum the constituents) instead of assuming the latter for every hyphenated name.

    Returns (OrderedDict {core: [(column_name, type), ...]}, shared_cores set,
    own_cram_cores set).
    """
    sr_pairs, pb_pairs, own_cram_cores = [], [], set()
    for column in samples:
        parsed = parse_minipileup_sample(column)
        if parsed is None:
            continue
        core, _stem, platform, tissue = parsed
        if core is None or core == '.':
            continue
        if platform == 'SR':
            sr_pairs.append((core, column))
            own_cram_cores.add(core)
        elif platform == 'PB' and (current_tissue is None or tissue == current_tissue):
            pb_pairs.append((core, column))
            own_cram_cores.add(core)

    # SHARED_CORES: cores with both an SR and a current-tissue PacBio column.
    shared_cores = {c for c, _ in pb_pairs} & {c for c, _ in sr_pairs}

    # The shared-core split and the merged-core gate collide for MAMC. The gate tests
    # the core name exactly (is_merged = '-' in core or core == 'MAMC'), so once MAMC is
    # split neither MAMC_SR nor MAMC_PB matches and both halves are judged as ordinary
    # individual cores -- MAMC_SR on per-core short reads against the strict standalone
    # cutoff rather than on donor aggregates. Its _SR half also tends to end up an
    # always-0/0 column, since CORE_CALLS carries a bare 'MAMC' that
    # split_shared_core_calls routes to _PB whenever the caller is a long-read one.
    #
    # main.nf splits on exactly the same condition, so this is O2's behaviour too and the
    # counts are left alone; it is only surfaced. MAMC is long-read-only in the current
    # samplesheets (never an SR core), so this does not arise in practice.
    if 'MAMC' in shared_cores:
        print(
            "WARNING: core 'MAMC' has both a short-read and a current-tissue PacBio CRAM, "
            "so it is split into MAMC_SR/MAMC_PB -- and the read-support gate's "
            "`core == 'MAMC'` test then matches neither half, so both are gated as "
            "ordinary individual cores instead of as a pooled multi-core library. "
            "MAMC_SR is likely to be 0/0 throughout. This mirrors main.nf, which splits "
            "on the same condition; counts are unchanged.",
            file=sys.stderr,
        )

    mapping = OrderedDict()
    for core, column in sr_pairs:
        key = core + '_SR' if core in shared_cores else core
        mapping.setdefault(key, []).append((column, 'SR'))
        own_cram_cores.add(key)
    for core, column in pb_pairs:
        key = core + '_PB' if core in shared_cores else core
        mapping.setdefault(key, []).append((column, 'PB'))
        own_cram_cores.add(key)

    return mapping, shared_cores, own_cram_cores


def core_cram_map_to_columns(mapping, current_tissue):
    """Convert a core_cram_map.tsv mapping from CRAM basenames to column names.

    parse_core_cram_map() yields CRAM basenames, but compute_core_counts() looks
    columns up directly. Rebuild the column name minipileup would have produced --
    the same strings the old basename-reconstructing code searched for -- so the
    O2-TSV path and the column-derived path share one lookup.
    """
    converted = OrderedDict()
    for core, entries in mapping.items():
        out = []
        for basename, ctype in entries:
            if ctype == 'SR':
                out.append((f"{basename}-SR", 'SR'))
            elif ctype == 'PB':
                out.append((
                    f"{basename}-PB-{current_tissue}" if current_tissue else f"{basename}-PB",
                    'PB',
                ))
            else:
                out.append((basename, ctype))
        converted[core] = out
    return converted


def add_merged_cores_from_core_calls(mapping, original_vcf, lr_callers, canon=None):
    """Add merged (hyphenated) cores to the map, derived from CORE_CALLS.

    A merged core has no CRAM of its own, so it never appears in the --sr-core /
    --lr-core arrays. It does appear in CORE_CALLS, because it has its own VCF
    upstream and merge_callers passes it through unchanged.

    compute_core_counts() recomputes every hyphenated core by summing its
    constituents and overwriting the first pass, so the CRAM entries recorded
    here are discarded for counting. They exist only to (a) make the merged core
    a key -- which is what puts it in `cores` and therefore in the sample columns
    -- and (b) supply core_ctypes for the platform-matched constituent fallback.

    Technology is classified the way main.nf classifies MERGED_SR_CORES /
    MERGED_PB_CORES: by whether the calling caller is a long-read caller.
    """
    canon = canon or (lambda c: c)

    # Group by canonical name so one merged core cannot become two keys if its
    # constituents are spelled in different orders across records, but key the map by
    # the spelling CORE_CALLS actually used -- that spelling becomes the output sample
    # column name, and renaming it would gratuitously diverge from O2's column names.
    merged_callers, spelling = {}, {}
    for rec in original_vcf.snvs.values():
        raw = rec.info.get("CORE_CALLS")
        if raw is None:
            continue
        if isinstance(raw, (tuple, list)):
            raw = ','.join(str(v) for v in raw)
        for core, callers in parse_core_calls(raw).items():
            if '-' in core:
                key = canon(core)
                merged_callers.setdefault(key, set()).update(callers)
                spelling.setdefault(key, core)

    merged_callers = {spelling[k]: v for k, v in merged_callers.items()}

    # Compare canonically: a pooled merged core already in the map from its own CRAM may
    # be spelled differently there than in CORE_CALLS, and adding it again would give one
    # physical library two sample columns.
    already = {canon(k) for k in mapping}

    for core, callers in merged_callers.items():
        if core in mapping or canon(core) in already:
            continue
        ctype = 'PB' if (callers & lr_callers) else 'SR'
        entries = []
        for part in core.split('-'):
            for key in (part, part + '_SR', part + '_PB'):
                entries.extend(
                    (bn, ct) for bn, ct in mapping.get(key, []) if ct == ctype
                )
        if not entries:
            # No constituent CRAM available. The key and its technology are all
            # that matter here; counts are overwritten by the merged summation.
            entries = [(core, ctype)]
        mapping[core] = entries

    return mapping


def split_shared_core_calls(core_calls, shared_cores, lr_callers):
    """Reconcile CORE_CALLS keys with split map keys.

    On O2, main.nf suffixes the VCF samplesheet as well, so CORE_CALLS already
    holds {core}_SR / {core}_PB. On AWS, merge_callers writes the raw core, so a
    shared core arrives unsuffixed and would match neither half of the split map
    -- every variant on it would be dropped at the `core not in core_calls` gate.

    Split it here using the same rule main.nf uses on the VCF side: long-read
    callers become the _PB entry, everything else the _SR entry. Cores that are
    not shared, and cores already suffixed upstream, pass through untouched.
    """
    if not shared_cores:
        return core_calls
    result = {}
    for core, callers in core_calls.items():
        if core not in shared_cores:
            result[core] = callers
            continue
        sr_callers = [c for c in callers if c not in lr_callers]
        pb_callers = [c for c in callers if c in lr_callers]
        if sr_callers:
            result[core + '_SR'] = sr_callers
        if pb_callers:
            result[core + '_PB'] = pb_callers
    return result


################################################################################
### Safeguards
################################################################################
# A core name that fails to match between the core map and CORE_CALLS costs a
# variant its GT and can cost it the record entirely -- and does so silently, since
# every lookup is a plain dict miss that reads as "this core did not call it". These
# checks run before any variant is written so a mismatch stops the run instead.


def check_map_columns_exist(core_cram_map, samples):
    """Fail if the core map names a sample column the minipileup VCF does not have.

    Unreachable on the column-derived path, where the map is built from these very
    columns. It guards the --core_cram_map path, where a basename that does not line
    up with minipileup's naming yields zero counts for that core with no other sign.
    """
    have = set(samples)
    missing = OrderedDict()
    for core, entries in core_cram_map.items():
        for column, ctype in entries:
            # Merged cores with no constituent CRAM are keyed by the core name itself
            # and are recomputed by summation, so they are not expected as columns.
            if ctype in ('SR', 'PB') and column not in have and column != core:
                missing.setdefault(core, []).append(column)
    if missing:
        detail = '; '.join(f"{core}: {', '.join(cols)}" for core, cols in missing.items())
        sys.exit(
            "ERROR: the core map names minipileup sample columns that do not exist "
            f"({detail}). Those cores would silently get zero counts. Columns present: "
            f"{', '.join(sorted(have))}."
        )


def audit_core_calls_vs_map(core_cram_map, original_vcf, shared_cores, lr_callers,
                            canon, on_unmatched="error"):
    """Reconcile the cores named in CORE_CALLS against the core map keys.

    A core in CORE_CALLS with no map key gets no sample column, so it can never be
    credited with calling anything; if it is the only core calling a variant, the
    variant is dropped. Fatal by default -- silently losing variants is worse than
    refusing to run. The reverse direction (a map key never seen in CORE_CALLS) only
    costs an always-0/0 column, so it warns.
    """
    mapped = {canon(c) for c in core_cram_map}
    called_records = {}
    for rec in original_vcf.snvs.values():
        raw = rec.info.get("CORE_CALLS")
        if raw is None:
            continue
        if isinstance(raw, (tuple, list)):
            raw = ','.join(str(v) for v in raw)
        calls = split_shared_core_calls(
            parse_core_calls(raw), shared_cores or set(), lr_callers or set()
        )
        for core in calls:
            called_records[canon(core)] = called_records.get(canon(core), 0) + 1

    unmatched = OrderedDict(
        (c, n) for c, n in sorted(called_records.items()) if c not in mapped
    )
    unused = sorted(c for c in mapped if c not in called_records)

    if unmatched:
        detail = ', '.join(f"{c} ({n} input records)" for c, n in unmatched.items())
        message = (
            f"{len(unmatched)} core(s) named in CORE_CALLS have no entry in the core "
            f"map: {detail}. Map keys: {', '.join(sorted(mapped))}. These cores get no "
            "sample column, so their calls cannot be credited and variants called only "
            "by them are dropped."
        )
        if on_unmatched == "error":
            sys.exit(
                "ERROR: " + message
                + " Re-run with --on-unmatched-core warn to proceed anyway."
            )
        print("WARNING: " + message, file=sys.stderr)

    if unused:
        print(
            f"WARNING: {len(unused)} core(s) in the core map never appear in "
            f"CORE_CALLS: {', '.join(unused)}. Their columns will be 0/0 throughout.",
            file=sys.stderr,
        )

    return unmatched, unused


def report_core_map(core_cram_map, shared_cores, strict_o2_compat, own_cram_cores=None):
    """One-line summary of how cores resolved, so a silent regression is visible."""
    # 'MAMC' counts as merged the same way the read-support gate counts it
    # (is_merged = '-' in core or core == 'MAMC'): a library pooled from several cores
    # under one opaque ID rather than a hyphenated constituent list.
    own_cram_cores = own_cram_cores or set()
    merged = [c for c in core_cram_map if '-' in c or c == 'MAMC']
    pooled = [c for c in merged if c in own_cram_cores]
    print(
        f"REPORT: core map has {len(core_cram_map)} core(s) "
        f"({len(merged)} merged -- {len(pooled)} with a CRAM of their own (counted "
        f"directly), {len(merged) - len(pooled)} call-only (summed from constituents); "
        f"{len(shared_cores or set())} shared core(s) split into _SR/_PB), "
        f"merged-core matching="
        f"{'literal (--strict-o2-compat)' if strict_o2_compat else 'canonical'}.",
        file=sys.stderr,
    )


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
        self.aggregate_counts = dict() # (chrom, pos, ref, alt) -> {PB: SampleCounts, SR: SampleCounts, ONT: SampleCounts}
        self.core_counts = dict()  # (chrom,pos,ref,alt) -> {core: {SR: SampleCounts, PB: SampleCounts, ONT: SampleCounts}}

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
                # takes care of real-base context padding (minipileup2 format, e.g. CG > TG)
                if len(ALT) == 1 and alt[0] == ALT and alt[1:] == record.ref[1:]: return i
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

            self.aggregate_counts[key] = agg

    def compute_core_counts(self, core_cram_map: dict, own_cram_cores: set=None):
        """
        Compute per-core counts for SR and tissue-matched PB.
        ONT is pooled-only and not tracked per core.

        core_cram_map: OrderedDict {core: [(sample_column, type), ...]}
          type is SR or PB; ONT entries (if any) are ignored.

        Sets self.core_counts:
            {key: {core: {SR: SampleCounts, PB: SampleCounts}}}

        Entries name minipileup sample columns directly. On the column-derived path
        those columns are where the map came from; on the --core_cram_map path
        core_cram_map_to_columns() has already rebuilt them from CRAM basenames.
        Either way the PacBio tissue match (only -PB-{current_tissue} columns count
        toward a core) has been applied while building the map.

        Merged cores (hyphen in name, e.g. 001A2-001C3) are recomputed by summing
        their constituent individual cores rather than from the CRAM map directly.
        This is correct by definition and avoids bugs where the upstream pipeline
        only passes one CRAM per constituent into the core_cram_map for merged cores.

        own_cram_cores: cores that have a minipileup sample column of their own. A
        hyphenated core in this set is a *pooled library* -- one CRAM sequenced from
        several tissue cores -- so it has real reads and keeps the counts read from its
        own column; summing constituents would throw those reads away. A hyphenated
        core absent from the set is call-only (a joint call over CRAMs it does not own),
        which is the case summation exists for. Empty/None reproduces the original
        behaviour of summing every hyphenated core, which is what the --core_cram_map
        path wants: an O2 map lists constituents' basenames under a merged key, so
        ownership cannot be recovered from it.
        """
        own_cram_cores = own_cram_cores or set()
        for key, counts_ in self.counts.items():
            self.core_counts[key] = {}
            for core, cram_list in core_cram_map.items():
                core_agg = {
                    'SR': SampleCounts(core + '-SR'),
                    'PB': SampleCounts(core + '-PB'),
                }
                for column, ctype in cram_list:
                    if ctype == 'SR':
                        grp = 'SR'
                    elif ctype == 'PB':
                        grp = 'PB'
                    else:
                        continue  # ONT and unknown types are not per-core
                    sc = counts_.get(column)
                    if sc is not None:
                        core_agg[grp].REF_ADF += sc.REF_ADF
                        core_agg[grp].REF_ADR += sc.REF_ADR
                        core_agg[grp].ALT_ADF += sc.ALT_ADF
                        core_agg[grp].ALT_ADR += sc.ALT_ADR
                self.core_counts[key][core] = core_agg

        # Recompute call-only merged-core counts by summing the cores they were called
        # over. This overrides whatever was aggregated from the CRAM map above, which may
        # be incomplete if the pipeline join only provided one CRAM per constituent.
        # Cores owning a CRAM are excluded -- their first-pass counts are their own reads.
        merged_cores = [c for c in core_cram_map if '-' in c and c not in own_cram_cores]
        for core in merged_cores:
            # Technology comes from this merged core's own map entries so that shared-core
            # constituents (stored as {part}_SR / {part}_PB) are looked up with the
            # matching suffix only — prevents an SR merged core picking up PB reads.
            core_ctypes = {ctype for _, ctype in core_cram_map.get(core, [])}
            # Tile the constituents with the largest available keys, so a pooled CRAM
            # covering several of them contributes its reads instead of being missed.
            sources, uncovered = constituent_sources(core, core_cram_map, core_ctypes)
            if uncovered:
                print(
                    f"WARNING: merged core '{core}' has no counts for constituent(s) "
                    f"{', '.join(sorted(uncovered))} — no core map entry covers them, so "
                    "they contribute nothing to its per-core counts.",
                    file=sys.stderr,
                )
            for key in self.core_counts:
                merged_sr = SampleCounts(core + '-SR')
                merged_pb = SampleCounts(core + '-PB')
                for src in sources:
                    part_cc = self.core_counts[key].get(src)
                    if part_cc is None:
                        continue
                    for grp, merged in (('SR', merged_sr), ('PB', merged_pb)):
                        merged.REF_ADF += part_cc[grp].REF_ADF
                        merged.REF_ADR += part_cc[grp].REF_ADR
                        merged.ALT_ADF += part_cc[grp].ALT_ADF
                        merged.ALT_ADR += part_cc[grp].ALT_ADR
                self.core_counts[key][core] = {'SR': merged_sr, 'PB': merged_pb}


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
            '##INFO=<ID=ONT_ADR,Number=2,Type=Integer,Description="Oxford Nanopore (long-read) reverse depths (REF, ALT)">'
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

    def _binom_pval(self, alt_total: int, total: int):
        """Return binomial p-value (less) against 0.5, or None if not enough reads."""
        if alt_total >= self.min_alt_binom and total > 0:
            return binom_pvalue(alt_total, total, 0.5, alternative='less')
        return None

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

    def write_tiered_vcf(self, out_vcf_path: str, keep_info: bool=False, core_cram_map: dict=None,
                         shared_cores: set=None, lr_callers: set=None, canon=None):
        """Write tiered VCF. With core_cram_map writes multi-sample per-core VCF; otherwise single-sample."""
        if core_cram_map:
            self._write_multisample(out_vcf_path, core_cram_map, shared_cores, lr_callers, canon)
        else:
            self._write_singlesample(out_vcf_path, keep_info)

    def _write_multisample(self, out_vcf_path: str, core_cram_map: dict,
                           shared_cores: set=None, lr_callers: set=None, canon=None):
        """Write multi-sample tiered VCF with one sample column per core."""
        no_pileup_counts, fail_filters = 0, 0
        written, t1, t2 = 0, 0, 0

        # canon maps a core name to the spelling CORE_CALLS is keyed by. Identity
        # under --strict-o2-compat, which reproduces O2's literal matching.
        canon = canon or (lambda c: c)
        cores = list(core_cram_map.keys())

        # FORMAT field definitions (per-core)
        format_defs = [
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

        # INFO field definitions (tissue-level)
        info_defs = [
            '##INFO=<ID=CrossTech,Number=0,Type=Flag,Description="Alt supported in both short read and tissue-matched PacBio at or above combined thresholds">',
            '##INFO=<ID=CrossCore,Number=0,Type=Flag,Description="Variant has GT=0/1 in more than one core">',
            '##INFO=<ID=CrossCaller,Number=0,Type=Flag,Description="Alt found in 2+ unique callers across all cores">',
            '##INFO=<ID=POOLED_PB_VAF,Number=1,Type=Float,Description="PacBio VAF pooled across all of a donor\'s tissues">',
            '##INFO=<ID=POOLED_ONT_VAF,Number=1,Type=Float,Description="ONT VAF pooled across all of a donor\'s tissues">',
            '##INFO=<ID=POOLED_PB_ADF,Number=2,Type=Integer,Description="PacBio forward depths pooled across all of a donor\'s tissues (REF,ALT)">',
            '##INFO=<ID=POOLED_PB_ADR,Number=2,Type=Integer,Description="PacBio reverse depths pooled across all of a donor\'s tissues (REF,ALT)">',
            '##INFO=<ID=POOLED_ONT_ADF,Number=2,Type=Integer,Description="ONT forward depths pooled across all of a donor\'s tissues (REF,ALT)">',
            '##INFO=<ID=POOLED_ONT_ADR,Number=2,Type=Integer,Description="ONT reverse depths pooled across all of a donor\'s tissues (REF,ALT)">',
            '##INFO=<ID=REGION,Number=1,Type=String,Description="SMaHT region classification: easy, diff, or ext">',
            '##INFO=<ID=CORE_CALLS,Number=1,Type=String,Description="Per-core caller presence: core1:caller1,caller2|core2:caller1">',
            '##INFO=<ID=GERMLINE_PVAL,Number=1,Type=Float,Description="Min binomial p-value for germline deviation across pooled platforms (tissue-level SR; donor-level PB and ONT)">',
            '##INFO=<ID=GERMLINE_PVAL_SR,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled SR data (tissue-level)">',
            '##INFO=<ID=GERMLINE_PVAL_PB,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled PacBio data across all tissues (donor-level)">',
            '##INFO=<ID=GERMLINE_PVAL_ONT,Number=1,Type=Float,Description="Binomial p-value for germline deviation in pooled ONT data across all tissues (donor-level)">',
            '##INFO=<ID=TIER,Number=1,Type=String,Description="Tissue-level tier: TIER1=SR+LR cross-tech supported, TIER2=SR-only">',
            '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher exact test p-value for strand balance on pooled counts">',
            '##INFO=<ID=SB_SRC,Number=1,Type=String,Description="Platform used for pooled Fisher strand test: SR or PB">',
        ]

        with pysam.VariantFile(self.original_vcf.vcf_path) as vf_in:
            in_header = vf_in.header

            # Build new header
            new_header = pysam.VariantHeader()
            for contig_name, contig in in_header.contigs.items():
                if contig.length is not None:
                    new_header.add_line(f'##contig=<ID={contig_name},length={contig.length}>')
                else:
                    new_header.add_line(f'##contig=<ID={contig_name}>')
            for fd in format_defs:
                new_header.add_line(fd)
            for id_ in info_defs:
                new_header.add_line(id_)
            for core in cores:
                new_header.add_sample(core)

            with pysam.VariantFile(out_vcf_path, "w", header=new_header) as vf_out:
                for key in sorted(self.snvs, key=lambda k: (self.chrom_order(k[0]), k[1])):
                    if key not in self.minipileup_vcf.aggregate_counts:
                        no_pileup_counts += 1
                        continue

                    tier = self.tiers.get(key)

                    orig_rec = self.snvs[key]
                    chrom, pos, ref, alt = key

                    # Parse CORE_CALLS to determine which cores called this variant
                    core_calls_raw = orig_rec.info.get("CORE_CALLS")
                    if isinstance(core_calls_raw, (tuple, list)):
                        core_calls_raw = ','.join(str(v) for v in core_calls_raw)
                    core_calls = split_shared_core_calls(
                        parse_core_calls(core_calls_raw),
                        shared_cores or set(),
                        lr_callers or set(),
                    )
                    core_calls = {canon(c): v for c, v in core_calls.items()}
                    per_core   = self.minipileup_vcf.core_counts.get(key, {})

                    # Per-core read support gate.
                    # - Individual SR cores: per-core SR reads vs strict standalone cutoff (target_p=1e-5).
                    # - Individual PB cores: combined SR+PB gate; falls back to SR-only standalone
                    #   cutoff when pb_alt=0 (e.g. LR caller found signal but pileup disagrees).
                    # - Merged cores (hyphen in name) or MAMC: combined SR+PB gate using aggregate
                    #   PB counts when PB is available; falls back to aggregate SR vs standalone
                    #   cutoff when pb_alt=0 or the donor has no PacBio data.
                    core_alt_support = {}  # core -> bool
                    agg_counts = self.minipileup_vcf.aggregate_counts[key]
                    agg_sr   = agg_counts['SR']
                    sr_total = agg_sr.REF_ADF + agg_sr.REF_ADR + agg_sr.ALT_ADF + agg_sr.ALT_ADR
                    sr_alt   = agg_sr.ALT_ADF + agg_sr.ALT_ADR
                    for core in cores:
                        if canon(core) not in core_calls:
                            core_alt_support[core] = False
                            continue
                        core_types = {ctype for _, ctype in core_cram_map.get(core, [])}
                        is_pb_core = 'PB' in core_types and 'SR' not in core_types
                        is_merged  = '-' in core or core == 'MAMC'

                        if is_merged:
                            # Merged core: combined SR+PB if PB available, else SR-only fallback.
                            # If the combined gate fails (e.g. pb_alt=0), fall back to SR-only.
                            pb_sc    = agg_counts['PB']
                            pb_total = pb_sc.REF_ADF + pb_sc.REF_ADR + pb_sc.ALT_ADF + pb_sc.ALT_ADR
                            pb_alt   = pb_sc.ALT_ADF + pb_sc.ALT_ADR
                            if pb_total > 0:
                                thresholds = get_read_cutoffs(sr_total, pb_total)
                                combined_pass = (
                                    sr_total > 0
                                    and sr_alt >= thresholds["combined_SR"]
                                    and pb_alt >= thresholds["combined_PB"]
                                )
                                sr_only_pass = (
                                    sr_total > 0
                                    and sr_alt >= get_read_cutoffs(sr_total, 0)["SR"]
                                )
                                core_alt_support[core] = combined_pass or sr_only_pass
                            else:
                                core_alt_support[core] = (
                                    sr_total > 0
                                    and sr_alt >= get_read_cutoffs(sr_total, 0)["SR"]
                                )
                        elif is_pb_core:
                            # Individual PB core: combined SR+PB gate; falls back to SR-only
                            # if pb_alt=0 (e.g. LR caller found signal but pileup shows no LR support).
                            pb_sc    = agg_counts['PB']
                            pb_total = pb_sc.REF_ADF + pb_sc.REF_ADR + pb_sc.ALT_ADF + pb_sc.ALT_ADR
                            pb_alt   = pb_sc.ALT_ADF + pb_sc.ALT_ADR
                            sr_only_pass = (
                                sr_total > 0
                                and sr_alt >= get_read_cutoffs(sr_total, 0)["SR"]
                            )
                            if pb_total > 0:
                                thresholds = get_read_cutoffs(sr_total, pb_total)
                                combined_pass = (
                                    sr_total > 0
                                    and sr_alt >= thresholds["combined_SR"]
                                    and pb_alt >= thresholds["combined_PB"]
                                )
                                core_alt_support[core] = combined_pass or sr_only_pass
                            else:
                                core_alt_support[core] = sr_only_pass
                        else:
                            # Individual SR core: per-core SR vs standalone cutoff (unchanged)
                            cc         = per_core.get(core, {'SR': SampleCounts(core), 'PB': SampleCounts(core)})
                            sc         = cc['SR']
                            core_total = sc.REF_ADF + sc.REF_ADR + sc.ALT_ADF + sc.ALT_ADR
                            core_alt   = sc.ALT_ADF + sc.ALT_ADR
                            core_alt_support[core] = (
                                core_total > 0
                                and core_alt >= get_read_cutoffs(core_total, 0)["SR"]
                            )

                    if not any(core_alt_support.values()):
                        continue

                    fisher_result = self.tests[key]["fisher"]
                    binom_result  = self.tests[key]["binomial"]
                    fisher_pass   = fisher_result.is_pass(self.strand_alpha)

                    if tier == "TIER2":
                        binom_pass      = binom_result.is_pass("SR", self.germline_alpha_SR)
                        binom_pass_long = True
                    elif tier == "TIER1":
                        binom_pass_pb  = binom_result.is_pass("PB",  self.germline_alpha)
                        binom_pass_ont = binom_result.is_pass("ONT", self.germline_alpha)
                        binom_pass, binom_pass_long = True, True
                        if binom_pass_pb is False or binom_pass_ont is False:
                            binom_pass_long = False
                    else:
                        binom_pass, binom_pass_long = True, True

                    if fisher_pass is False or binom_pass is False or binom_pass_long is False:
                        fail_filters += 1
                        continue

                    # REGION (preserved from upstream filters; may be absent in older VCFs)
                    try:
                        region_val = orig_rec.info.get("REGION")
                    except (ValueError, KeyError):
                        region_val = None

                    # Donor-level pooled counts (all PB and ONT tissues for this donor)
                    tpb  = agg_counts['PB']
                    tont = agg_counts['ONT']

                    tpb_total  = tpb.REF_ADF  + tpb.REF_ADR  + tpb.ALT_ADF  + tpb.ALT_ADR
                    tont_total = tont.REF_ADF + tont.REF_ADR + tont.ALT_ADF + tont.ALT_ADR

                    pooled_pb_vaf  = float(tpb.ALT_ADF  + tpb.ALT_ADR)  / tpb_total  if tpb_total  > 0 else None
                    pooled_ont_vaf = float(tont.ALT_ADF + tont.ALT_ADR) / tont_total if tont_total > 0 else None

                    # Create new record
                    new_rec = new_header.new_record(
                        contig=chrom,
                        start=orig_rec.start,
                        stop=orig_rec.stop,
                        alleles=(ref, alt),
                        id=orig_rec.id,
                        qual=orig_rec.qual,
                    )

                    # Set INFO fields
                    if tier == 'TIER1':
                        new_rec.info['CrossTech'] = True
                    if pooled_pb_vaf is not None:
                        new_rec.info['POOLED_PB_VAF']  = pooled_pb_vaf
                    if pooled_ont_vaf is not None:
                        new_rec.info['POOLED_ONT_VAF'] = pooled_ont_vaf
                    new_rec.info['POOLED_PB_ADF']  = (tpb.REF_ADF,  tpb.ALT_ADF)
                    new_rec.info['POOLED_PB_ADR']  = (tpb.REF_ADR,  tpb.ALT_ADR)
                    new_rec.info['POOLED_ONT_ADF'] = (tont.REF_ADF, tont.ALT_ADF)
                    new_rec.info['POOLED_ONT_ADR'] = (tont.REF_ADR, tont.ALT_ADR)
                    if region_val is not None:
                        new_rec.info['REGION'] = region_val
                    if core_calls_raw is not None:
                        new_rec.info['CORE_CALLS'] = core_calls_raw

                    gp_sr      = binom_result.p_value_SR
                    gp_pb      = binom_result.p_value_PB
                    gp_ont     = binom_result.p_value_ONT
                    all_gp     = [p for p in [gp_sr, gp_pb, gp_ont] if p is not None]
                    gp_min     = min(all_gp) if all_gp else None

                    # Pooled germline p-values → INFO (SR: tissue-level; PB/ONT: donor-level)
                    if gp_sr is not None:
                        new_rec.info['GERMLINE_PVAL_SR'] = gp_sr
                    if gp_pb is not None:
                        new_rec.info['GERMLINE_PVAL_PB'] = gp_pb
                    if gp_ont is not None:
                        new_rec.info['GERMLINE_PVAL_ONT'] = gp_ont
                    if gp_min is not None:
                        new_rec.info['GERMLINE_PVAL'] = gp_min

                    new_rec.info['TIER']    = tier
                    new_rec.info['SB_PVAL'] = fisher_result.p_value
                    new_rec.info['SB_SRC']  = fisher_result.group

                    # Per-core FORMAT fields
                    per_core = self.minipileup_vcf.core_counts.get(key, {})

                    for core in cores:
                        called = canon(core) in core_calls
                        cc = per_core.get(core, {
                            'SR': SampleCounts(core),
                            'PB': SampleCounts(core),
                        })
                        sr = cc['SR']
                        pb = cc['PB']

                        if called and core_alt_support.get(core, False):
                            gt = (0, 1)
                        else:
                            gt = (0, 0)
                        new_rec.samples[core]['GT'] = gt

                        core_types = {ctype for _, ctype in core_cram_map.get(core, [])}
                        if core_types == {'ONT'}:
                            raise ValueError(f"ONT-only core '{core}' is not supported; TECH must be SR or PB")
                        tech = 'PB' if 'PB' in core_types else 'SR'
                        new_rec.samples[core]['TECH'] = tech

                        if tech == 'SR':
                            new_rec.samples[core]['SR_ADF'] = (sr.REF_ADF, sr.ALT_ADF)
                            new_rec.samples[core]['SR_ADR'] = (sr.REF_ADR, sr.ALT_ADR)
                            new_rec.samples[core]['PB_ADF'] = (None, None)
                            new_rec.samples[core]['PB_ADR'] = (None, None)
                            sr_total = sr.REF_ADF + sr.REF_ADR + sr.ALT_ADF + sr.ALT_ADR
                            new_rec.samples[core]['SR_VAF'] = float(sr.ALT_ADF + sr.ALT_ADR) / sr_total if sr_total > 0 else 0.0
                            new_rec.samples[core]['PB_VAF'] = None
                        else:
                            new_rec.samples[core]['SR_ADF'] = (None, None)
                            new_rec.samples[core]['SR_ADR'] = (None, None)
                            new_rec.samples[core]['PB_ADF'] = (pb.REF_ADF, pb.ALT_ADF)
                            new_rec.samples[core]['PB_ADR'] = (pb.REF_ADR, pb.ALT_ADR)
                            new_rec.samples[core]['SR_VAF'] = None
                            pb_total = pb.REF_ADF + pb.REF_ADR + pb.ALT_ADF + pb.ALT_ADR
                            new_rec.samples[core]['PB_VAF'] = float(pb.ALT_ADF + pb.ALT_ADR) / pb_total if pb_total > 0 else 0.0

                        new_rec.samples[core]['ALT_SUPPORT_PASS'] = 1 if core_alt_support.get(core, False) else 0

                        if called:
                            new_rec.samples[core]['CALLERS'] = core_calls[canon(core)]
                        else:
                            new_rec.samples[core]['CALLERS'] = ['.']

                    called_pass_cores = [c for c in cores if new_rec.samples[c]['GT'] == (0, 1)]
                    any_pb_called = any(new_rec.samples[c]['TECH'] == 'PB' for c in called_pass_cores)
                    # Shared cores are split into {name}_SR / {name}_PB columns; collapse back
                    # to base name so both count as one physical biopsy for CrossCore.
                    def _base_core(c):
                        return c[:-3] if (c.endswith('_SR') or c.endswith('_PB')) else c
                    seen_bases = set()
                    n_called_pass = 0
                    for c in called_pass_cores:
                        if (('-' in c or c == 'MAMC') and
                                new_rec.samples[c]['TECH'] == 'SR' and not any_pb_called):
                            continue
                        b = _base_core(c)
                        if b not in seen_bases:
                            seen_bases.add(b)
                            n_called_pass += 1
                    if n_called_pass > 1:
                        new_rec.info['CrossCore'] = True

                    all_callers = {
                        caller
                        for core in cores
                        if canon(core) in core_calls
                        for caller in core_calls[canon(core)]
                    }
                    if len(all_callers) >= 2:
                        new_rec.info['CrossCaller'] = True

                    vf_out.write(new_rec)
                    written += 1
                    t1 += (tier == "TIER1")
                    t2 += (tier == "TIER2")

        print(f"INFO: Wrote multi-sample tiered VCF to {out_vcf_path}.")
        print(f"REPORT: wrote {written} records (TIER1={t1}, TIER2={t2}).")
        print(f"REPORT: {no_pileup_counts} variants with no pileup counts.")
        print(f"REPORT: {fail_filters} variants failing filters.")

    def _write_singlesample(self, out_vcf_path: str, keep_info: bool=False):
        """Write tiered VCF to out_vcf_path (single-sample legacy mode).
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
                    if callers_value and len(callers_value) > 1:
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

    parser.add_argument("--current_tissue", default=None,
                    help="Tissue ID for this run (e.g. SMHT005-3AF).")
    parser.add_argument("--core_cram_map", default=None,
                    help="TSV file mapping core -> CRAM basename -> type (SR/PB/ONT), as generated by main.nf on O2. Alternative to reading cores from the minipileup sample columns; retained for O2 comparison runs.")

    parser.add_argument("--lr-callers", dest="lr_callers", default="longcallD", metavar="CALLERS",
                    help="Comma-separated long-read caller names, mirroring main.nf's params.lr_callers (default: longcallD). Used to classify merged cores and to split shared-core CORE_CALLS entries.")
    parser.add_argument("--on-unmatched-core", dest="on_unmatched_core",
                    choices=("error", "warn"), default=None,
                    help="What to do when CORE_CALLS names a core with no entry in the core map. Such a core can never be assigned a GT, so variants called only by it are dropped. 'error' (the default) refuses to run; 'warn' reports and continues. Defaults to 'warn' under --strict-o2-compat, whose purpose is to reproduce those mismatches.")
    parser.add_argument("--strict-o2-compat", dest="strict_o2_compat", action="store_true",
                    help="Match core names literally instead of canonicalising merged-core constituent order. Reproduces O2's behaviour, including the silent 0/0 column when the core map and CORE_CALLS spell a merged core differently. For bug-for-bug comparison runs only.")

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

    lr_callers = {c.strip() for c in args.lr_callers.split(',') if c.strip()}
    canon = (lambda c: c) if args.strict_o2_compat else canonical_core
    # Strict mode exists to reproduce O2's literal matching, mismatches included, so
    # an unmatched core there is the expected finding rather than a reason to stop.
    on_unmatched_core = args.on_unmatched_core or (
        "warn" if args.strict_o2_compat else "error"
    )

    with pysam.VariantFile(args.minipileup_vcf) as _mp:
        minipileup_samples = list(_mp.header.samples)

    ovcf  = OriginalVCF(args.input_vcf)
    mpvcf = MinipileupVCF(args.minipileup_vcf, ovcf, current_tissue=args.current_tissue)

    core_cram_map_data = None
    shared_cores, own_cram_cores = set(), set()
    if any('__' in s for s in minipileup_samples):
        # Cores are stamped into the minipileup sample columns by
        # minipileup2-parallel.sh. That takes precedence over --core_cram_map.
        core_cram_map_data, shared_cores, own_cram_cores = build_core_cram_map_from_samples(
            minipileup_samples, args.current_tissue
        )
        add_merged_cores_from_core_calls(core_cram_map_data, ovcf, lr_callers, canon)
    elif args.core_cram_map:
        # An O2-generated map already carries its merged cores (main.nf expands them
        # from the VCF samplesheet), so no injection here -- adding CORE_CALLS-derived
        # keys would duplicate them under a second spelling. Canonical matching is
        # what reconciles the two spellings on this path.
        core_cram_map_data = core_cram_map_to_columns(
            parse_core_cram_map(args.core_cram_map), args.current_tissue
        )

    if core_cram_map_data:
        check_map_columns_exist(core_cram_map_data, minipileup_samples)
        audit_core_calls_vs_map(
            core_cram_map_data, ovcf, shared_cores, lr_callers, canon,
            on_unmatched=on_unmatched_core,
        )
        report_core_map(core_cram_map_data, shared_cores, args.strict_o2_compat,
                        own_cram_cores)
        mpvcf.compute_core_counts(core_cram_map_data, own_cram_cores)

    tvcf  = TieredVCF(
        ovcf, mpvcf,
        strand_alpha=args.strand_alpha,
        germline_alpha=args.germline_alpha,
        germline_alpha_SR=args.germline_alpha_SR,
        min_alt_PB=args.min_alt_PB,
        min_alt_binom=args.min_alt_binom,
    )
    tvcf.write_tiered_vcf(args.output_vcf, args.keep_info, core_cram_map=core_cram_map_data,
                          shared_cores=shared_cores, lr_callers=lr_callers, canon=canon)

    if args.output_vcf.endswith(".vcf.gz"):
        pysam.tabix_index(args.output_vcf, preset="vcf", force=True)
