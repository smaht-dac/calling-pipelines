#!/usr/bin/env python3
import pysam
import csv
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from multiprocessing import cpu_count, set_start_method

# ---------------------------------------------------------------
# Parallel Step 5–7 (Sex-aware haplotype classification + phasing tag output)
#   * Parallel per-variant row
#   * Each worker opens its own BAMs (no shared state)  -> avoids race conditions
#   * Parent writes outputs (single writer)             -> avoids race conditions
# ---------------------------------------------------------------

def get_base_at_pos(aln, chrom, pos):
    """Return base in read at 1-based genome position (None if not covered)."""
    for qpos, refpos in aln.get_aligned_pairs(matches_only=True):
        if refpos == pos - 1:
            return aln.query_sequence[qpos]
    return None

def classify_row(row, bam_paths, sex):
    """
    Worker: process a single TSV row (dict), open BAMs locally, and return the enriched row.
    This function is safe to run in parallel processes.
    """
    # open BAMs per process to avoid shared handle issues
    bam_files = [pysam.AlignmentFile(b, "rb") for b in bam_paths]
    try:
        chrom = row["chrom"]
        var_pos = int(row["Var_pos"])
        germ_pos = row["Germ_pos"]
        if germ_pos == "NA":
            row.update({k: "0" for k in ["case_both","case_germ_only","case_var_only","case_none","total_reads"]})
            row["hap_classification"] = "no_germline"
            return row

        germ_pos = int(germ_pos)
        var_alt = row["Var_alt"]
        germ_alt = row["Germ_alt"]

        counts = Counter()
        fetch_start = min(var_pos, germ_pos) - 1
        fetch_end   = max(var_pos, germ_pos)

        for bam in bam_files:
            for read in bam.fetch(chrom, fetch_start, fetch_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                base_var = get_base_at_pos(read, chrom, var_pos)
                base_germ = get_base_at_pos(read, chrom, germ_pos)
                if base_var is None or base_germ is None:
                    continue

                if base_var == var_alt and base_germ == germ_alt:
                    counts["case_both"] += 1
                elif base_var == var_alt and base_germ != germ_alt:
                    counts["case_var_only"] += 1
                elif base_var != var_alt and base_germ == germ_alt:
                    counts["case_germ_only"] += 1
                else:
                    counts["case_none"] += 1

        total_reads = sum(counts.values())
        for key in ["case_both","case_germ_only","case_var_only","case_none"]:
            row[key] = counts.get(key, 0)
        row["total_reads"] = total_reads

        case_both = counts.get("case_both", 0)
        case_none = counts.get("case_none", 0)
        case_var_only = counts.get("case_var_only", 0)
        case_germ_only = counts.get("case_germ_only", 0)

        # Sex-specific chromosome logic
        if sex == "male" or sex == 'M':
            # Haploid X/Y in males
            if chrom in ["chrX", "chrY"]:
                if case_both > 1 and case_germ_only > 1 and case_none <= 1 and case_var_only <= 1: # 2 haplo
                    row["hap_classification"] = "mosaic"
                if case_both > 1 and case_germ_only <= 1 and case_none <= 1 and case_var_only <= 1: # 1 haplo
                    row["hap_classification"] = "germline"
                else: # 3+ haplo
                    row["hap_classification"] = "artifact"
                return row

        elif sex == "female" or sex == 'F':
            # Females: X diploid, Y = artifact
            if chrom == "chrY":
                row["hap_classification"] = "artifact"
                return row

        else:
            # Unknown sex: treat X/Y as artifacts
            if chrom in ["chrX", "chrY"]:
                row["hap_classification"] = "artifact"
                return row

        # Diploid classification logic
        if case_both > 1 and case_none > 1 and case_var_only <= 1 and case_germ_only <= 1: # 2 haplo
            row["hap_classification"] = "germline"
        elif case_both <= 1 and case_none <= 1 and case_var_only > 1 and case_germ_only > 1: # 2 haplo
            row["hap_classification"] = "germline"
        elif case_both > 1 and case_none > 1 and case_germ_only > 1 and case_var_only <= 1: # 3 haplo
            row["hap_classification"] = "mosaic"
        elif case_both <= 1 and case_none > 1 and case_var_only > 1 and case_germ_only > 1: # 3 haplo
            row["hap_classification"] = "mosaic"
        elif case_both <= 1 and case_var_only <= 1: # no var reads
            row["hap_classification"] = "artifact"
        else: # 1 or 4+ haplo
            row["hap_classification"] = "artifact"

        return row
    finally:
        for bam in bam_files:
            bam.close()

def write_phasing_tags(results, out_prefix):
    """Generate phasing_tags.tsv for bcftools annotate."""
    out_path = f"{out_prefix}_phasing_tags.tsv"
    with open(out_path, "w") as f:
        f.write("#CHROM\tPOS\tPHASING\n")
        for r in results:
            chrom = r["chrom"]
            pos = r["Var_pos"]
            phase = r.get("hap_classification", "").upper()
            if phase == "MOSAIC":
                tag = "MOSAIC_PHASED"
            elif phase in ["GERMLINE", "ARTIFACT"]:
                tag = phase
            elif phase == "HAPLOID":
                # If you want haploid calls to be taggable, choose one; default to UNABLE_TO_PHASE
                tag = "UNABLE_TO_PHASE"
            elif phase == "NO_GERMLINE":
                tag = "UNABLE_TO_PHASE"
            elif phase == "NO_COVERAGE":
                tag = "UNABLE_TO_PHASE"
            else:
                tag = "UNABLE_TO_PHASE"
            f.write(f"{chrom}\t{pos}\t{tag}\n")
    return out_path

def read_tsv_rows(tsv_file):
    """Read TSV into a list of (idx, row_dict) preserving order."""
    rows = []
    with open(tsv_file) as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        for idx, row in enumerate(reader):
            rows.append((idx, row))
    return rows

def main():
    # Use spawn to avoid forking issues with pysam (safer across platforms)
    try:
        set_start_method("spawn")
    except RuntimeError:
        # already set; safe to ignore
        pass

    parser = argparse.ArgumentParser(
        description="Phasing step: parallel read-level haplotype analysis, classification, and tag output."
    )
    parser.add_argument("-t", "--tsv", required=True, help="TSV from Step 4")
    parser.add_argument("-b", "--bams", nargs="+", required=True, help="PacBio BAMs")
    parser.add_argument("-s", "--sex", default="unknown",
                        choices=["male", "female", "unknown"],
                        help="Sex of the sample (for haploid X/Y logic)")
    parser.add_argument("-i", "--id", default="haplotype_summary",
                        help="Output prefix for results")
    parser.add_argument("-w", "--workers", type=int, default=cpu_count(),
                        help="Number of worker processes (default: CPU count)")
    args = parser.parse_args()

    pysam.set_verbosity(0)

    # Read TSV once in parent
    indexed_rows = read_tsv_rows(args.tsv)
    if not indexed_rows:
        print("[Warn] Input TSV has no rows.")
        # still emit empty outputs with headers
        with open(f"{args.id}.haplotyped.tsv", "w", newline="") as f_out:
            f_out.write("")  # nothing to write
        write_phasing_tags([], args.id)
        return

    # Prepare worker function (top-level, pickleable)
    worker = partial(classify_row, bam_paths=args.bams, sex=args.sex)

    # Submit in parallel; collect as (idx, enriched_row)
    results_buffer = [None] * len(indexed_rows)

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        future_to_idx = {ex.submit(worker, row): idx for idx, row in indexed_rows}
        for fut in as_completed(future_to_idx):
            idx = future_to_idx[fut]
            try:
                out_row = fut.result()
            except Exception as e:
                # Include the row index to debug data-specific failures
                raise RuntimeError(f"Worker failed on row index {idx}: {e}") from e
            results_buffer[idx] = out_row

    # Write main classification table (deterministic order)
    out_tsv = f"{args.id}.haplotyped.tsv"
    fieldnames = list(results_buffer[0].keys()) if results_buffer else []
    # Ensure consistent column order (optional: force some first)
    if fieldnames:
        # put some known columns first if present
        preferred = ["chrom","Var_pos","Var_alt","Germ_pos","Germ_alt",
                     "case_both","case_germ_only","case_var_only","case_none",
                     "total_reads","hap_classification"]
        ordered = [c for c in preferred if c in fieldnames] + [c for c in fieldnames if c not in preferred]
        fieldnames = ordered

    with open(out_tsv, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results_buffer:
            writer.writerow(row)

    # Write phasing tags
    phasing_path = write_phasing_tags(results_buffer, args.id)


        # --- Write debug intermediate file ---
    debug_tsv = f"{args.id}.read_counts.phasing.tsv"
    debug_fields = [
        "chrom", "Var_pos", "Var_alt", "hap_classification",
        "case_both", "case_none", "case_var_only", "case_germ_only"
    ]
    with open(debug_tsv, "w", newline="") as f_debug:
        writer = csv.DictWriter(f_debug, fieldnames=debug_fields, delimiter="\t")
        writer.writeheader()
        for row in results_buffer:
            # Default 0 for missing values to avoid KeyErrors
            writer.writerow({
                "chrom": row.get("chrom", ""),
                "Var_pos": row.get("Var_pos", ""),
                "Var_alt": row.get("Var_alt", ""),
                "hap_classification": row.get("hap_classification", ""),
                "case_both": row.get("case_both", 0),
                "case_none": row.get("case_none", 0),
                "case_var_only": row.get("case_var_only", 0),
                "case_germ_only": row.get("case_germ_only", 0)
            })

    print(f"[Done] Debug read count file: {debug_tsv}")


    print(f"[Done] Classification TSV: {out_tsv}")
    print(f"[Done] Phasing tags TSV: {phasing_path}")

if __name__ == "__main__":
    main()

