#!/usr/bin/env python3
"""
Step 4b: Deduplicate sequences after trimming.

Removes exact duplicate sequences from trimmed FASTAs. When multiple sequences
share the same nucleotide sequence, keeps one representative and logs the rest.
Deduplicates CBP and PCR independently so we can track what was collapsed.

Input:  data/04_trimmed/{gene}.fasta
Output: data/04b_deduped/{gene}.fasta
        data/04b_deduped/notes/dedup_log.tsv
        data/04b_deduped/notes/{gene}_collapsed.tsv

Usage:
    ./venv/bin/python pipeline/step4b_dedup.py
"""

import csv
import os
from collections import defaultdict

IN_DIR = "data/04_trimmed"
OUT_DIR = "data/04b_deduped"
NOTES_DIR = os.path.join(OUT_DIR, "notes")

OUTGROUP_IDS = {"FN597644.1", "CP026362.1"}


def read_fasta(path):
    entries = []
    header, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(parts)))
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        entries.append((header, "".join(parts)))
    return entries


def write_fasta(path, entries, line_width=80):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


def dedup_entries(entries):
    """Deduplicate a list of (header, seq) entries.

    Returns (kept, collapsed) where:
      - kept: list of (header, seq) with one representative per unique sequence
      - collapsed: dict mapping kept_header -> [list of removed headers]
    """
    seen = {}  # seq -> header (first seen)
    kept = []
    collapsed = defaultdict(list)

    for header, seq in entries:
        seq_upper = seq.upper()
        if seq_upper in seen:
            collapsed[seen[seq_upper]].append(header)
        else:
            seen[seq_upper] = header
            kept.append((header, seq))

    return kept, collapsed


def main():
    os.makedirs(NOTES_DIR, exist_ok=True)

    fasta_files = sorted(f for f in os.listdir(IN_DIR) if f.endswith(".fasta"))

    log_path = os.path.join(NOTES_DIR, "dedup_log.tsv")
    log = open(log_path, "w")
    log.write("gene\tcbp_in\tcbp_unique\tcbp_removed\t"
              "pcr_in\tpcr_unique\tpcr_removed\t"
              "outgroup\ttotal_in\ttotal_out\n")

    for fname in fasta_files:
        gene = fname.replace(".fasta", "")
        in_path = os.path.join(IN_DIR, fname)
        out_path = os.path.join(OUT_DIR, fname)

        entries = read_fasta(in_path)

        # Separate by type
        outgroup = [(h, s) for h, s in entries if h in OUTGROUP_IDS]
        cbp = [(h, s) for h, s in entries if h.startswith("CBP-")]
        pcr = [(h, s) for h, s in entries
               if h not in OUTGROUP_IDS and not h.startswith("CBP-")]

        # Dedup CBP and PCR separately
        cbp_kept, cbp_collapsed = dedup_entries(cbp)
        pcr_kept, pcr_collapsed = dedup_entries(pcr)

        # Combine: outgroup + deduped CBP + deduped PCR
        deduped = outgroup + cbp_kept + pcr_kept
        write_fasta(out_path, deduped)

        # Write collapsed details
        all_collapsed = {}
        all_collapsed.update(cbp_collapsed)
        all_collapsed.update(pcr_collapsed)

        if all_collapsed:
            detail_path = os.path.join(NOTES_DIR, f"{gene}_collapsed.tsv")
            with open(detail_path, "w") as df:
                df.write("kept_header\tremoved_header\n")
                for kept_h, removed_list in sorted(all_collapsed.items()):
                    for rem_h in removed_list:
                        df.write(f"{kept_h}\t{rem_h}\n")

        cbp_removed = sum(len(v) for v in cbp_collapsed.values())
        pcr_removed = sum(len(v) for v in pcr_collapsed.values())

        print(f"  {gene:12s}  CBP: {len(cbp):>4} -> {len(cbp_kept):>4} "
              f"(-{cbp_removed})  PCR: {len(pcr):>6} -> {len(pcr_kept):>6} "
              f"(-{pcr_removed})  total: {len(entries)} -> {len(deduped)}")

        log.write(f"{gene}\t{len(cbp)}\t{len(cbp_kept)}\t{cbp_removed}\t"
                  f"{len(pcr)}\t{len(pcr_kept)}\t{pcr_removed}\t"
                  f"{len(outgroup)}\t{len(entries)}\t{len(deduped)}\n")
        log.flush()

    log.close()
    print(f"\nDedup log: {log_path}")


if __name__ == "__main__":
    main()
