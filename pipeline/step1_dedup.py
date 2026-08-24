#!/usr/bin/env python3
"""
Step 1: Combine all environments per gene and deduplicate by exact sequence.

Input:  data/00_merged/{gene}/*.fasta  (merged reads from step 0)
Output: data/01_deduped/{gene}.fasta
        data/01_deduped/notes/dedup_log.tsv

Usage:
    ./venv/bin/python pipeline/step1_dedup.py
"""

import os
import glob

IN_DIR = "data/00_merged"
OUT_DIR = "data/01_deduped"
NOTES_DIR = os.path.join(OUT_DIR, "notes")


def read_fasta(path):
    """Yield (header, sequence) from a FASTA file."""
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def main():
    os.makedirs(NOTES_DIR, exist_ok=True)

    gene_dirs = sorted(
        d for d in os.listdir(IN_DIR)
        if os.path.isdir(os.path.join(IN_DIR, d)) and d != "notes"
    )

    log_path = os.path.join(NOTES_DIR, "dedup_log.tsv")
    log = open(log_path, "w")
    log.write("gene\ttotal_seqs\tunique_seqs\tduplicates\tdedup_rate\n")

    for gene in gene_dirs:
        gene_path = os.path.join(IN_DIR, gene)
        # Only merged files (not _unmerged)
        fasta_files = sorted(glob.glob(os.path.join(gene_path, "[0-9].fasta")))

        seen_seqs = {}  # seq -> header (keep first)
        total = 0

        for fpath in fasta_files:
            for header, seq in read_fasta(fpath):
                total += 1
                if seq not in seen_seqs:
                    seen_seqs[seq] = header

        unique = len(seen_seqs)
        dupes = total - unique
        rate = 100 * dupes / total if total > 0 else 0

        # Write deduplicated FASTA
        out_path = os.path.join(OUT_DIR, f"{gene}.fasta")
        with open(out_path, "w") as f:
            for seq, header in seen_seqs.items():
                f.write(f">{header}\n{seq}\n")

        print(f"  {gene:12s}  total={total:>10,}  unique={unique:>10,}  "
              f"dupes={dupes:>10,} ({rate:.1f}%)")
        log.write(f"{gene}\t{total}\t{unique}\t{dupes}\t{rate:.1f}\n")

    log.close()
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    main()
