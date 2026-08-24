#!/usr/bin/env python3
"""
Step 3: Add CBP isolates and outgroup to filtered PCR reads.

Uses full-length CBP amplicons from data/04_isolates/isolate_amplicons2/
and full-length outgroup amplicons from data/outgroup_amplicons/.

Input:  data/02_filtered/{gene}.fasta
Output: data/03_with_cbp/{gene}.fasta
        data/03_with_cbp/notes/cbp_log.tsv

Usage:
    ./venv/bin/python pipeline/step3_add_cbp.py
"""

import os
import csv

IN_DIR = "data/02_filtered"
OUT_DIR = "data/03_with_cbp"
NOTES_DIR = os.path.join(OUT_DIR, "notes")
CBP_SOURCE = "data/04_isolates/isolate_amplicons2"
OUTGROUP_SOURCE = "data/outgroup_amplicons"
PRIMERS_TSV = "notes/primers.tsv"

# Gene name mapping: isolate_amplicons2 uses "acsA" not "acsA_2"
CBP_GENE_NAME_MAP = {"acsA_2": "acsA"}


def load_gene_species_map():
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


def load_cbp(species, gene):
    """Load full-length CBP isolate amplicons."""
    cbp_gene = CBP_GENE_NAME_MAP.get(gene, gene)
    fasta_path = os.path.join(CBP_SOURCE, f"{cbp_gene}_{species}.fasta")
    if not os.path.exists(fasta_path):
        return []
    return list(read_fasta(fasta_path))


def load_outgroup(gene):
    """Load full-length outgroup amplicon."""
    fasta_path = os.path.join(OUTGROUP_SOURCE, f"{gene}_outgroup.fasta")
    if not os.path.exists(fasta_path):
        return []
    return list(read_fasta(fasta_path))


def main():
    os.makedirs(NOTES_DIR, exist_ok=True)
    gene_species = load_gene_species_map()

    log_path = os.path.join(NOTES_DIR, "cbp_log.tsv")
    log = open(log_path, "w")
    log.write("gene\tspecies\tpcr_reads\tcbp_count\tcbp_len\toutgroup_count\toutgroup_len\ttotal\n")

    for gene, species in sorted(gene_species.items()):
        filtered_path = os.path.join(IN_DIR, f"{gene}.fasta")
        if not os.path.exists(filtered_path):
            print(f"  {gene:12s} ({species:12s})  — no filtered file, skipping")
            continue

        # Count PCR reads
        pcr_seqs = list(read_fasta(filtered_path))
        pcr_count = len(pcr_seqs)

        # Get full-length CBP and outgroup
        cbp = load_cbp(species, gene)
        outgroup = load_outgroup(gene)

        # Write combined output
        out_path = os.path.join(OUT_DIR, f"{gene}.fasta")
        total = 0
        with open(out_path, "w") as out:
            # Outgroup first (for tree rooting)
            for header, seq in outgroup:
                out.write(f">{header}\n{seq}\n")
                total += 1
            # CBP isolates
            for header, seq in cbp:
                out.write(f">{header}\n{seq}\n")
                total += 1
            # PCR reads
            for header, seq in pcr_seqs:
                out.write(f">{header}\n{seq}\n")
                total += 1

        cbp_len = len(cbp[0][1]) if cbp else 0
        out_len = len(outgroup[0][1]) if outgroup else 0

        print(f"  {gene:12s} ({species:12s})  "
              f"PCR={pcr_count:>10,}  CBP={len(cbp):>4} ({cbp_len}bp)  "
              f"outgroup={len(outgroup)} ({out_len}bp)  total={total:>10,}")

        log.write(f"{gene}\t{species}\t{pcr_count}\t{len(cbp)}\t{cbp_len}\t"
                  f"{len(outgroup)}\t{out_len}\t{total}\n")

    log.close()
    print(f"\nLog written to {log_path}")


if __name__ == "__main__":
    main()
