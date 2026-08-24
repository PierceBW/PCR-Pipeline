#!/usr/bin/env python3
"""
Step 2: Gene-specific identity filter (>=95%) using pre-built BLAST databases.

Only processes the correct species for each gene (from primers.tsv).
BLASTs each read against the gene-specific reference amplicon and keeps
reads with >=95% identity.

Input:  data/01_deduped/{gene}.fasta
Output: data/02_filtered/{gene}.fasta
        data/02_filtered/notes/filter_log.tsv

Usage:
    ./venv/bin/python pipeline/step2_gene_filter.py
"""

import os
import subprocess
import tempfile
import csv

IN_DIR = "data/01_deduped"
OUT_DIR = "data/02_filtered"
NOTES_DIR = os.path.join(OUT_DIR, "notes")
GENE_REF_DIR = "data/13_gene_filtered/gene_refs"
PRIMERS_TSV = "notes/primers.tsv"
MIN_PIDENT = 95.0


def load_gene_species_map():
    """Load gene -> species mapping from primers.tsv."""
    mapping = {}
    with open(PRIMERS_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["gene"]] = row["species"]
    return mapping


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


def blast_filter(fasta_path, db_path, min_pident):
    """
    BLAST sequences against gene-specific DB.
    Returns dict: header -> best pident, and set of headers that pass.
    """
    # Run BLAST
    cmd = [
        "blastn",
        "-query", fasta_path,
        "-db", db_path,
        "-outfmt", "6 qseqid pident length",
        "-max_target_seqs", "1",
        "-evalue", "1e-10",
        "-num_threads", "4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse: keep best hit per query
    best_hits = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        qid = parts[0]
        pident = float(parts[1])
        if qid not in best_hits or pident > best_hits[qid]:
            best_hits[qid] = pident

    passed = {qid for qid, pid in best_hits.items() if pid >= min_pident}
    return best_hits, passed


def main():
    os.makedirs(NOTES_DIR, exist_ok=True)
    gene_species = load_gene_species_map()

    log_path = os.path.join(NOTES_DIR, "filter_log.tsv")
    log = open(log_path, "w")
    log.write("gene\tspecies\ttotal\tkept\tdropped\tno_hit\t"
              "min_pident\tmax_pident\tmean_pident\n")

    for gene, species in sorted(gene_species.items()):
        fasta_path = os.path.join(IN_DIR, f"{gene}.fasta")
        if not os.path.exists(fasta_path):
            print(f"  {gene:12s} ({species:12s})  — no input file, skipping")
            continue

        # Check if input has any sequences
        seqs = list(read_fasta(fasta_path))
        total = len(seqs)
        if total == 0:
            print(f"  {gene:12s} ({species:12s})  — 0 sequences, skipping")
            log.write(f"{gene}\t{species}\t0\t0\t0\t0\t\t\t\n")
            continue

        # Find the BLAST database
        db_path = os.path.join(GENE_REF_DIR, f"{gene}_{species}_db")
        if not os.path.exists(db_path + ".nin") and not os.path.exists(db_path + ".ndb"):
            print(f"  {gene:12s} ({species:12s})  — no BLAST DB found, skipping")
            log.write(f"{gene}\t{species}\t{total}\t0\t0\t{total}\t\t\t\n")
            continue

        # For sequences with dashes (gap-filled), write a temp file without dashes for BLAST
        # BLAST can't handle dash characters in query sequences
        has_dashes = any("-" in seq for _, seq in seqs)
        if has_dashes:
            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
            for header, seq in seqs:
                tmp.write(f">{header}\n{seq.replace('-', '')}\n")
            tmp.close()
            blast_input = tmp.name
        else:
            blast_input = fasta_path

        # Run BLAST
        best_hits, passed = blast_filter(blast_input, db_path, MIN_PIDENT)

        if has_dashes:
            os.unlink(blast_input)

        # Write filtered output (keep original sequences with dashes intact)
        out_path = os.path.join(OUT_DIR, f"{gene}.fasta")
        kept = 0
        dropped = 0
        no_hit = 0
        pidents = []

        with open(out_path, "w") as out:
            for header, seq in seqs:
                if header in passed:
                    out.write(f">{header}\n{seq}\n")
                    kept += 1
                    pidents.append(best_hits[header])
                elif header in best_hits:
                    dropped += 1
                else:
                    no_hit += 1

        min_p = min(pidents) if pidents else 0
        max_p = max(pidents) if pidents else 0
        mean_p = sum(pidents) / len(pidents) if pidents else 0

        print(f"  {gene:12s} ({species:12s})  total={total:>10,}  "
              f"kept={kept:>10,} ({100*kept/total:.1f}%)  "
              f"dropped={dropped:>7,}  no_hit={no_hit:>7,}  "
              f"pident={min_p:.1f}-{max_p:.1f} (mean {mean_p:.1f})")

        log.write(f"{gene}\t{species}\t{total}\t{kept}\t{dropped}\t{no_hit}\t"
                  f"{min_p:.2f}\t{max_p:.2f}\t{mean_p:.2f}\n")

    log.close()
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    main()
