#!/usr/bin/env python3
"""
Step 7: Prune trees to the Jocelyn in-group (CBP LCA clade).

For each gene, find the LCA of all CBP isolates in the n=200 tree,
keep only sequences within that clade, rebuild the tree, and reroot.

PCR reads outside the LCA are excluded — they may be from a different
species or taxon not represented in the culture collection.

Input:  data/05_ecosim/200/{gene}/{gene}.fasta + .nwk
Output: data/07_ingroup/{gene}/{gene}.fasta
        data/07_ingroup/{gene}/{gene}.nwk
        data/07_ingroup/{gene}/{gene}_excluded.txt

Usage:
    ./venv/bin/python pipeline/step7_prune_ingroup.py
"""

import csv
import os
import re
import subprocess
import shutil

from Bio import Phylo

IN_DIR = "data/05_ecosim/200"
OUT_DIR = "data/07_ingroup"
NOTES_DIR = os.path.join(OUT_DIR, "notes")
PRIMERS_TSV = "notes/primers.tsv"
OUTGROUP_IDS = {"FN597644.1", "CP026362.1"}

USABLE_GENES = {"acuA", "sorA", "yvqK", "albG", "thiD", "iolB", "acsA_2", "alkH"}


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


def build_tree(fasta_path, nwk_path):
    cmd = ["fasttree", "-gtr", "-nt", "-quiet", fasta_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    with open(nwk_path, "w") as f:
        f.write(result.stdout)


def reroot_tree(nwk_in, nwk_out):
    tree = Phylo.read(nwk_in, "newick")
    outgroup = None
    for og_id in OUTGROUP_IDS:
        outgroup = next(
            (c for c in tree.get_terminals() if c.name == og_id), None
        )
        if outgroup:
            break
    if outgroup is None:
        # No outgroup — just copy as-is (midpoint rooted by FastTree)
        shutil.copy(nwk_in, nwk_out)
        return False
    tree.root_with_outgroup(outgroup)
    Phylo.write(tree, nwk_out, "newick")
    return True


def get_pcr_env(name):
    m = re.search(r"PCR_E(\d+)_", name)
    return f"E{m.group(1)}" if m else "?"


def main():
    os.makedirs(NOTES_DIR, exist_ok=True)

    gene_species = {}
    with open(PRIMERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene_species[row["gene"]] = row["species"]

    log_path = os.path.join(NOTES_DIR, "prune_log.tsv")
    log = open(log_path, "w")
    log.write("gene\tspecies\ttotal_orig\tcbp\tpcr_orig\toutgroup\t"
              "pcr_inside\tpcr_outside\tpct_outside\t"
              "total_pruned\tstatus\n")

    print("Step 7: In-group pruning")
    print()

    for gene in sorted(USABLE_GENES):
        species = gene_species.get(gene, "?")
        nwk_path = os.path.join(IN_DIR, gene, f"{gene}.nwk")
        fasta_path = os.path.join(IN_DIR, gene, f"{gene}.fasta")

        if not os.path.exists(nwk_path) or not os.path.exists(fasta_path):
            print(f"  {gene:12s}  SKIP — missing files")
            continue

        # Load tree
        tree = Phylo.read(nwk_path, "newick")
        terminals = tree.get_terminals()

        cbp_leaves = [t for t in terminals if t.name.startswith("CBP-")]
        pcr_leaves = [t for t in terminals if t.name.startswith("PCR_")]
        og_leaves = [t for t in terminals if t.name in OUTGROUP_IDS]

        # Find CBP LCA
        if len(cbp_leaves) >= 2:
            lca_clade = tree.common_ancestor(cbp_leaves)
        elif len(cbp_leaves) == 1:
            lca_clade = cbp_leaves[0]
        else:
            print(f"  {gene:12s}  SKIP — no CBP leaves")
            continue

        lca_terminal_names = {t.name for t in lca_clade.get_terminals()}

        pcr_inside = [t for t in pcr_leaves if t.name in lca_terminal_names]
        pcr_outside = [t for t in pcr_leaves if t.name not in lca_terminal_names]
        pct_outside = 100 * len(pcr_outside) / len(pcr_leaves) if pcr_leaves else 0

        # Determine what to keep
        keep_names = lca_terminal_names.copy()

        # Always include outgroup if present (for rerooting the pruned tree)
        for og in og_leaves:
            keep_names.add(og.name)

        # Load FASTA and filter
        entries = read_fasta(fasta_path)
        entry_map = {h: s for h, s in entries}

        kept_entries = []
        # Outgroup first (EcoSim requirement)
        for og in og_leaves:
            if og.name in entry_map:
                kept_entries.append((og.name, entry_map[og.name]))

        # Then CBP
        for h, s in entries:
            if h.startswith("CBP-") and h in keep_names:
                kept_entries.append((h, s))

        # Then PCR inside LCA
        for h, s in entries:
            if h.startswith("PCR_") and h in lca_terminal_names:
                kept_entries.append((h, s))

        # Write outputs
        gene_dir = os.path.join(OUT_DIR, gene)
        out_fasta = os.path.join(gene_dir, f"{gene}.fasta")
        out_nwk_raw = os.path.join(gene_dir, f"{gene}_unrooted.nwk")
        out_nwk = os.path.join(gene_dir, f"{gene}.nwk")
        excluded_path = os.path.join(gene_dir, f"{gene}_excluded.txt")

        # Handle the case where ALL PCR reads are outside
        if len(pcr_inside) == 0:
            status = "ALL_PCR_OUTSIDE"
            print(f"  {gene:12s} ({species:12s})  "
                  f"CBP={len(cbp_leaves):>4}  PCR={len(pcr_leaves):>4}  "
                  f"PCR_outside={len(pcr_outside):>4} ({pct_outside:.0f}%)  "
                  f"** ALL PCR OUTSIDE — CBP-only output **")

            # Still write CBP-only + outgroup for reference
            write_fasta(out_fasta, kept_entries)
            build_tree(out_fasta, out_nwk_raw)
            reroot_tree(out_nwk_raw, out_nwk)
        else:
            status = "OK"
            write_fasta(out_fasta, kept_entries)

            # Rebuild tree on pruned set
            build_tree(out_fasta, out_nwk_raw)
            reroot_tree(out_nwk_raw, out_nwk)

            print(f"  {gene:12s} ({species:12s})  "
                  f"CBP={len(cbp_leaves):>4}  PCR={len(pcr_leaves):>4}  "
                  f"PCR_inside={len(pcr_inside):>4}  PCR_outside={len(pcr_outside):>4} ({pct_outside:.0f}%)  "
                  f"pruned_total={len(kept_entries)}")

        # Write excluded list
        os.makedirs(gene_dir, exist_ok=True)
        with open(excluded_path, "w") as f:
            f.write("header\tenv\treason\n")
            for t in pcr_outside:
                env = get_pcr_env(t.name)
                f.write(f"{t.name}\t{env}\toutside_CBP_LCA\n")

        # Write env breakdown of excluded
        if pcr_outside:
            env_counts = {}
            for t in pcr_outside:
                env = get_pcr_env(t.name)
                env_counts[env] = env_counts.get(env, 0) + 1
            env_str = ", ".join(f"{k}:{v}" for k, v in sorted(env_counts.items()))
            print(f"               excluded envs: {env_str}")

        log.write(f"{gene}\t{species}\t{len(terminals)}\t{len(cbp_leaves)}\t"
                  f"{len(pcr_leaves)}\t{len(og_leaves)}\t"
                  f"{len(pcr_inside)}\t{len(pcr_outside)}\t{pct_outside:.1f}\t"
                  f"{len(kept_entries)}\t{status}\n")
        log.flush()

    log.close()
    print(f"\nPrune log: {log_path}")


if __name__ == "__main__":
    main()
