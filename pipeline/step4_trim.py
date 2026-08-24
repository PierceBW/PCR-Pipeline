#!/usr/bin/env python3
"""
Step 4: Primer-guided trimming — find forward primer at the start and
reverse primer RC at the end of each sequence, then trim to the
primer-to-primer region.

Phase 1 (--explore): For each gene, report where primers land in PCR reads,
        CBP, and outgroup. Show per-sequence offset distributions.
Phase 2 (--trim):    Actually trim every sequence to its primer boundaries.

Input:  data/03_with_cbp/{gene}.fasta
Output: data/04_trimmed/{gene}.fasta

Usage:
    ./venv/bin/python pipeline/step4_trim.py --explore
    ./venv/bin/python pipeline/step4_trim.py --trim
"""

import os
import csv
import argparse
import glob
from collections import Counter

IN_DIR = "data/03_with_cbp"
OUT_DIR = "data/04_trimmed"
NOTES_DIR = os.path.join(OUT_DIR, "notes")
PRIMERS_TSV = "notes/primers.tsv"

# Minimum fraction of primer bases that must match to accept a hit
MIN_PRIMER_MATCH_FRAC = 0.75

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def revcomp(seq):
    return seq.translate(COMP)[::-1]


def read_fasta(path):
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


def classify_seq(header):
    if header.startswith("CBP"):
        return "CBP"
    elif header.startswith("FN") or header.startswith("CP"):
        return "outgroup"
    else:
        return "PCR"


def find_primer(seq, primer, search_range, from_end=False):
    """
    Find the best position of a primer within a search range of a sequence.
    Returns (offset, n_matches) or (None, 0).

    If from_end=True, search from the 3' end (for reverse primer RC).
    offset is always measured from the relevant end.
    """
    plen = len(primer)
    best_pos = None
    best_matches = 0
    min_matches = int(plen * MIN_PRIMER_MATCH_FRAC)

    if from_end:
        # Search backwards from the end
        for offset in range(0, min(search_range, len(seq) - plen + 1)):
            start = len(seq) - plen - offset
            if start < 0:
                break
            n = sum(1 for a, b in zip(seq[start:start + plen], primer) if a.upper() == b.upper())
            if n > best_matches:
                best_matches = n
                best_pos = offset
    else:
        # Search forwards from the start
        for offset in range(0, min(search_range, len(seq) - plen + 1)):
            n = sum(1 for a, b in zip(seq[offset:offset + plen], primer) if a.upper() == b.upper())
            if n > best_matches:
                best_matches = n
                best_pos = offset

    if best_matches >= min_matches:
        return best_pos, best_matches
    return None, 0


def load_primers():
    primers = {}
    with open(PRIMERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            primers[row["gene"]] = (row["forward_primer"], row["reverse_primer"])
    return primers


def explore():
    """Per-gene report of primer positions in all sequence types."""
    primers = load_primers()
    gene_files = sorted(glob.glob(os.path.join(IN_DIR, "*.fasta")))

    for gene_file in gene_files:
        gene = os.path.basename(gene_file).replace(".fasta", "")
        if gene not in primers:
            continue

        fwd_primer, rev_primer = primers[gene]
        rev_rc = revcomp(rev_primer)

        # Collect primer offsets by sequence type
        results = {"CBP": {"left": [], "right": [], "total": 0, "no_left": 0, "no_right": 0},
                   "outgroup": {"left": [], "right": [], "total": 0, "no_left": 0, "no_right": 0},
                   "PCR": {"left": [], "right": [], "total": 0, "no_left": 0, "no_right": 0}}

        for header, seq in read_fasta(gene_file):
            stype = classify_seq(header)
            results[stype]["total"] += 1

            # Search for forward primer near start (within first 50 bp)
            left_off, left_match = find_primer(seq, fwd_primer, search_range=50, from_end=False)
            if left_off is not None:
                results[stype]["left"].append(left_off)
            else:
                results[stype]["no_left"] += 1

            # Search for reverse primer RC near end (within last 50 bp)
            right_off, right_match = find_primer(seq, rev_rc, search_range=50, from_end=True)
            if right_off is not None:
                results[stype]["right"].append(right_off)
            else:
                results[stype]["no_right"] += 1

        # Print results
        print(f"\n{'='*80}")
        print(f"Gene: {gene}  |  fwd_primer: {fwd_primer} ({len(fwd_primer)}bp)  |  rev_primer_RC: {rev_rc} ({len(rev_rc)}bp)")
        print(f"{'='*80}")

        for stype in ["CBP", "outgroup", "PCR"]:
            r = results[stype]
            if r["total"] == 0:
                continue

            left_ctr = Counter(r["left"])
            right_ctr = Counter(r["right"])
            left_mode = left_ctr.most_common(1)[0] if left_ctr else ("--", 0)
            right_mode = right_ctr.most_common(1)[0] if right_ctr else ("--", 0)

            print(f"\n  {stype} (n={r['total']:,}):")
            print(f"    Left offset (fwd primer from start):  mode={left_mode[0]} ({left_mode[1]:,}x)  "
                  f"no_hit={r['no_left']}  dist={dict(left_ctr.most_common(5))}")
            print(f"    Right offset (rev primer RC from end): mode={right_mode[0]} ({right_mode[1]:,}x)  "
                  f"no_hit={r['no_right']}  dist={dict(right_ctr.most_common(5))}")

        # Suggest trim
        # Use CBP offsets as ground truth (should be 0,0)
        cbp_left = Counter(results["CBP"]["left"]).most_common(1)
        cbp_right = Counter(results["CBP"]["right"]).most_common(1)
        pcr_left = Counter(results["PCR"]["left"]).most_common(1)
        pcr_right = Counter(results["PCR"]["right"]).most_common(1)

        cbp_l = cbp_left[0][0] if cbp_left else "?"
        cbp_r = cbp_right[0][0] if cbp_right else "?"
        pcr_l = pcr_left[0][0] if pcr_left else "?"
        pcr_r = pcr_right[0][0] if pcr_right else "?"

        print(f"\n  SUMMARY: CBP trim=({cbp_l},{cbp_r})  PCR trim=({pcr_l},{pcr_r})")
        if cbp_left and cbp_right:
            cbp_example_len = len(list(read_fasta(gene_file)).__next__[1]) if False else 0
            # Expected output length = CBP length (primer to primer)
            print(f"  -> Trim each sequence: cut left_offset from start, right_offset from end")


def trim():
    """Trim every sequence to its primer boundaries."""
    os.makedirs(NOTES_DIR, exist_ok=True)
    primers = load_primers()
    gene_files = sorted(glob.glob(os.path.join(IN_DIR, "*.fasta")))

    log_path = os.path.join(NOTES_DIR, "trim_log.tsv")
    log = open(log_path, "w")
    log.write("gene\ttotal\tkept\tdropped\tdrop_pct\t"
              "no_fwd\tno_rev\ttoo_short\tclipped_3prime\t"
              "cbp_left_mode\tcbp_right_mode\tpcr_left_mode\tpcr_right_mode\t"
              "uniform_len\n")

    for gene_file in gene_files:
        gene = os.path.basename(gene_file).replace(".fasta", "")
        if gene not in primers:
            continue

        fwd_primer, rev_primer = primers[gene]
        rev_rc = revcomp(rev_primer)
        fwd_len = len(fwd_primer)
        rev_len = len(rev_rc)

        # First pass: find modal offsets per type for reporting
        cbp_lefts = []
        cbp_rights = []
        pcr_lefts = []
        pcr_rights = []

        seqs = list(read_fasta(gene_file))

        for header, seq in seqs:
            stype = classify_seq(header)
            left_off, _ = find_primer(seq, fwd_primer, search_range=50, from_end=False)
            right_off, _ = find_primer(seq, rev_rc, search_range=50, from_end=True)
            if stype == "CBP":
                if left_off is not None: cbp_lefts.append(left_off)
                if right_off is not None: cbp_rights.append(right_off)
            elif stype == "PCR":
                if left_off is not None: pcr_lefts.append(left_off)
                if right_off is not None: pcr_rights.append(right_off)

        cbp_l_mode = Counter(cbp_lefts).most_common(1)[0][0] if cbp_lefts else 0
        cbp_r_mode = Counter(cbp_rights).most_common(1)[0][0] if cbp_rights else 0
        pcr_l_mode = Counter(pcr_lefts).most_common(1)[0][0] if pcr_lefts else 0
        pcr_r_mode = Counter(pcr_rights).most_common(1)[0][0] if pcr_rights else 0

        # Second pass: primer-guided trim into memory
        primer_trimmed = []
        no_fwd = 0
        no_rev = 0

        for header, seq in seqs:
            left_off, left_m = find_primer(seq, fwd_primer, search_range=50, from_end=False)
            right_off, right_m = find_primer(seq, rev_rc, search_range=50, from_end=True)

            if left_off is None:
                no_fwd += 1
                continue
            if right_off is None:
                no_rev += 1
                continue

            # Trim: keep from left_off to (len - right_off)
            trimmed_seq = seq[left_off:len(seq) - right_off] if right_off > 0 else seq[left_off:]
            primer_trimmed.append((header, trimmed_seq))

        # Third pass: enforce uniform length
        # Mode length = CBP length = the correct amplicon size
        pre_clip_lens = Counter(len(s) for _, s in primer_trimmed)
        mode_len = pre_clip_lens.most_common(1)[0][0] if pre_clip_lens else 0

        out_path = os.path.join(OUT_DIR, f"{gene}.fasta")
        trimmed_count = 0
        too_short = 0

        with open(out_path, "w") as out:
            for header, seq in primer_trimmed:
                if len(seq) < mode_len:
                    too_short += 1
                    continue
                # Clip to mode length from 3' end (keep 5' primer boundary)
                clipped = seq[:mode_len]
                out.write(f">{header}\n{clipped}\n")
                trimmed_count += 1

        clipped_count = sum(1 for _, s in primer_trimmed if len(s) > mode_len)
        total_dropped = len(seqs) - trimmed_count
        drop_pct = 100.0 * total_dropped / len(seqs) if len(seqs) > 0 else 0

        print(f"  {gene:12s}  total={len(seqs):>10,}  kept={trimmed_count:>10,}  "
              f"dropped={total_dropped:>6} ({drop_pct:.2f}%)  "
              f"[no_fwd={no_fwd} no_rev={no_rev} too_short={too_short}]  "
              f"clipped_3prime={clipped_count}  "
              f"uniform_len={mode_len}bp")

        log.write(f"{gene}\t{len(seqs)}\t{trimmed_count}\t{total_dropped}\t{drop_pct:.2f}\t"
                  f"{no_fwd}\t{no_rev}\t{too_short}\t{clipped_count}\t"
                  f"{cbp_l_mode}\t{cbp_r_mode}\t{pcr_l_mode}\t{pcr_r_mode}\t"
                  f"{mode_len}\n")

    log.close()
    print(f"\nLog written to {log_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--explore", action="store_true")
    parser.add_argument("--trim", action="store_true")
    args = parser.parse_args()

    if args.explore:
        explore()
    elif args.trim:
        trim()
    else:
        print("Specify --explore or --trim")


if __name__ == "__main__":
    main()
