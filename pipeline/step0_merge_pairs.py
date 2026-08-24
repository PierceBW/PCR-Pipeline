#!/usr/bin/env python3
"""
Step 0: Extract paired reads from BAM files and merge R1+R2 into
full-length amplicons (~500-600 bp) using the overlap region.

Input:  source_data/PCR-Primer/datasets/*.bam  (76 BAM files: 4 envs × 19 genes)
Output: data/00_merged/{gene}/{env}.fasta
        data/00_merged/{gene}/{env}_unmerged.fasta
        data/00_merged/notes/merge_log.tsv

Usage:
    ./venv/bin/python pipeline/step0_merge_pairs.py
"""

import os
import glob
import re
from collections import defaultdict

import pysam

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BAM_DIR = "source_data/PCR-Primer/datasets"
OUT_DIR = "data/00_merged"
NOTES_DIR = os.path.join(OUT_DIR, "notes")

# Standard read length expected from Illumina
EXPECTED_READ_LEN = 301

# Max mismatch fraction in overlap to accept a merge
MAX_OVERLAP_MISMATCH_FRAC = 0.10  # 10%

# Gene name normalisation: BAM filenames use "Amj" but primers.tsv uses "amj"
GENE_NAME_MAP = {"Amj": "amj", "yndE_2": "yndE_2", "acsA_2": "acsA_2"}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def reverse_complement(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def parse_bam_filename(path):
    """Return (env_id, gene_name) from a BAM filename, or (None, None)."""
    basename = os.path.basename(path)
    # Pattern: 18Primers-{env}_sorted.bam__{gene}_{hash}.bam
    m = re.match(r"18Primers-(\d+)_sorted\.bam__(.+?)_bdb11cef.*\.bam$", basename)
    if not m:
        return None, None
    env_id = m.group(1)
    gene = m.group(2)
    # Normalise gene name
    gene = GENE_NAME_MAP.get(gene, gene)
    return env_id, gene


def merge_pair(r1_seq, r1_qual, r2_seq, r2_qual, overlap_bp):
    """
    Merge R1 (forward strand) and R2 (pysam-returned, already on + strand)
    given the overlap in base pairs (computed from reference coordinates).

    Returns merged sequence string, or None if overlap quality is too poor.
    """
    if overlap_bp <= 0:
        return None

    # R1 covers the left side; its last `overlap_bp` bases overlap with
    # the first `overlap_bp` bases of R2.
    r1_unique = r1_seq[:-overlap_bp]
    r1_over = r1_seq[-overlap_bp:]
    r1_over_q = r1_qual[-overlap_bp:]

    r2_over = r2_seq[:overlap_bp]
    r2_over_q = r2_qual[:overlap_bp]
    r2_unique = r2_seq[overlap_bp:]

    # Check mismatch rate
    mismatches = sum(1 for a, b in zip(r1_over, r2_over) if a != b)
    if mismatches / overlap_bp > MAX_OVERLAP_MISMATCH_FRAC:
        return None

    # Build consensus for overlap: pick base with higher quality
    consensus = []
    for a, qa, b, qb in zip(r1_over, r1_over_q, r2_over, r2_over_q):
        if a == b:
            consensus.append(a)
        elif qa >= qb:
            consensus.append(a)
        else:
            consensus.append(b)

    return r1_unique + "".join(consensus) + r2_unique


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_bam(bam_path, env_id, gene):
    """Process one BAM file. Returns stats dict."""
    af = pysam.AlignmentFile(bam_path, "rb")

    # Group reads by query_name
    reads_by_name = defaultdict(dict)
    for read in af:
        qn = read.query_name
        if read.query_sequence is None:
            continue
        label = "R1" if read.is_read1 else "R2"
        if label not in reads_by_name[qn]:
            reads_by_name[qn][label] = read
    af.close()

    stats = {
        "gene": gene, "env": env_id,
        "total_pairs": 0, "merged": 0, "gap_filled": 0,
        "unmerged_orientation": 0,
        "unmerged_length": 0, "unmerged_no_ref": 0,
        "unmerged_mismatch": 0, "unmerged_missing_mate": 0,
        "merged_lengths": [],
    }

    # Output files
    gene_dir = os.path.join(OUT_DIR, gene)
    os.makedirs(gene_dir, exist_ok=True)
    merged_path = os.path.join(gene_dir, f"{env_id}.fasta")
    unmerged_path = os.path.join(gene_dir, f"{env_id}_unmerged.fasta")

    merged_out = open(merged_path, "w")
    unmerged_out = open(unmerged_path, "w")

    for qn, mates in reads_by_name.items():
        stats["total_pairs"] += 1

        if "R1" not in mates or "R2" not in mates:
            stats["unmerged_missing_mate"] += 1
            # Write whichever read we have
            for label, read in mates.items():
                header = f">{gene}_{env_id}_{qn}_{label}"
                unmerged_out.write(f"{header}\n{read.query_sequence}\n")
            continue

        r1 = mates["R1"]
        r2 = mates["R2"]

        # Check orientation: expect R1=FWD, R2=REV
        if r1.is_reverse or not r2.is_reverse:
            stats["unmerged_orientation"] += 1
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R1\n{r1.query_sequence}\n")
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R2\n{r2.query_sequence}\n")
            continue

        # Check length
        if len(r1.query_sequence) != EXPECTED_READ_LEN or len(r2.query_sequence) != EXPECTED_READ_LEN:
            stats["unmerged_length"] += 1
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R1\n{r1.query_sequence}\n")
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R2\n{r2.query_sequence}\n")
            continue

        # Compute overlap from reference coordinates
        if r1.reference_end is None or r2.reference_start is None:
            stats["unmerged_no_ref"] += 1
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R1\n{r1.query_sequence}\n")
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R2\n{r2.query_sequence}\n")
            continue

        overlap_bp = r1.reference_end - r2.reference_start

        # Negative overlap = gap between reads; fill with dashes
        if overlap_bp <= 0:
            gap_size = -overlap_bp
            merged_seq = r1.query_sequence + ("-" * gap_size) + r2.query_sequence
            stats["gap_filled"] += 1
            stats["merged_lengths"].append(len(merged_seq))
            header = f">{gene}_{env_id}_{qn}"
            merged_out.write(f"{header}\n{merged_seq}\n")
            continue

        # pysam returns query_sequence on + strand for REV reads already
        merged_seq = merge_pair(
            r1.query_sequence, r1.query_qualities,
            r2.query_sequence, r2.query_qualities,
            overlap_bp,
        )

        if merged_seq is None:
            stats["unmerged_mismatch"] += 1
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R1\n{r1.query_sequence}\n")
            unmerged_out.write(f">{gene}_{env_id}_{qn}_R2\n{r2.query_sequence}\n")
            continue

        stats["merged"] += 1
        stats["merged_lengths"].append(len(merged_seq))
        header = f">{gene}_{env_id}_{qn}"
        merged_out.write(f"{header}\n{merged_seq}\n")

    merged_out.close()
    unmerged_out.close()

    # Remove empty unmerged files
    if os.path.getsize(unmerged_path) == 0:
        os.remove(unmerged_path)

    return stats


def main():
    os.makedirs(NOTES_DIR, exist_ok=True)

    # Find all gene-specific BAM files (skip the unsplit parent BAMs)
    bam_files = sorted(glob.glob(os.path.join(BAM_DIR, "*.bam")))
    gene_bams = []
    for path in bam_files:
        env_id, gene = parse_bam_filename(path)
        if env_id and gene:
            gene_bams.append((path, env_id, gene))

    print(f"Found {len(gene_bams)} gene-specific BAM files")

    # Log file
    log_path = os.path.join(NOTES_DIR, "merge_log.tsv")
    log = open(log_path, "w")
    log.write("gene\tenv\ttotal_pairs\tmerged\tgap_filled\tunmerged_orientation\t"
              "unmerged_length\tunmerged_no_ref\tunmerged_mismatch\t"
              "unmerged_missing_mate\tmerge_rate\tmin_len\tmax_len\tmean_len\n")

    for path, env_id, gene in gene_bams:
        print(f"  {gene:12s} env={env_id} ... ", end="", flush=True)
        stats = process_bam(path, env_id, gene)

        lengths = stats["merged_lengths"]
        min_len = min(lengths) if lengths else 0
        max_len = max(lengths) if lengths else 0
        mean_len = sum(lengths) / len(lengths) if lengths else 0
        total_merged = stats["merged"] + stats["gap_filled"]
        merge_rate = (
            100 * total_merged / stats["total_pairs"]
            if stats["total_pairs"] > 0 else 0
        )

        gap_note = f"  (gap_filled={stats['gap_filled']})" if stats["gap_filled"] else ""
        print(f"merged={total_merged:>7,} / {stats['total_pairs']:>7,} "
              f"({merge_rate:.1f}%)  len={min_len}-{max_len} (mean {mean_len:.0f}){gap_note}")

        log.write(f"{gene}\t{env_id}\t{stats['total_pairs']}\t{stats['merged']}\t"
                  f"{stats['gap_filled']}\t"
                  f"{stats['unmerged_orientation']}\t{stats['unmerged_length']}\t"
                  f"{stats['unmerged_no_ref']}\t{stats['unmerged_mismatch']}\t"
                  f"{stats['unmerged_missing_mate']}\t{merge_rate:.1f}\t"
                  f"{min_len}\t{max_len}\t{mean_len:.1f}\n")

    log.close()
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    main()
