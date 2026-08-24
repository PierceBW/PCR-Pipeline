#!/usr/bin/env python3
"""
Step 5: Subsample, build trees, and run EcoSim.

For each gene, subsample PCR reads to --n-pcr (default 200), keep all CBP +
outgroup, build tree with FastTree, reroot at outgroup, run EcoSim.

Only processes the 8 usable genes (those with substantial PCR data).

Input:  data/04b_deduped/{gene}.fasta
Output: data/05_ecosim/{n_pcr}/{gene}/{gene}.fasta
        data/05_ecosim/{n_pcr}/{gene}/{gene}.nwk
        data/05_ecosim/{n_pcr}/{gene}/{gene}.xml
        data/05_ecosim/{n_pcr}/notes/ecosim_results.tsv

Usage:
    ./venv/bin/python pipeline/step5_ecosim.py --n-pcr 200
    ./venv/bin/python pipeline/step5_ecosim.py --n-pcr 2000
    ./venv/bin/python pipeline/step5_ecosim.py --n-pcr 200 --no-ecosim
"""

import argparse
import csv
import os
import random
import re
import shutil
import subprocess
import sys
import time

from Bio import Phylo

IN_DIR = "data/04b_deduped"
OUT_BASE = "data/05_ecosim"
PRIMERS_TSV = "notes/primers.tsv"

# EcoSim lives in sibling repo
ECOSIM_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "ecosim")
)
ECOSIM_JAR = os.path.join(ECOSIM_DIR, "ecosim.jar")
JAVA = "/usr/local/opt/openjdk/bin/java"

OUTGROUP_IDS = {"FN597644.1", "CP026362.1"}

# Only run these 8 genes — the ones with substantial PCR data
USABLE_GENES = {
    "acuA", "sorA", "yvqK", "albG", "thiD",   # inaquosorum
    "iolB", "acsA_2", "alkH",                   # spizizenii
}

# Thesis expected ecotype counts (Jocelyn Wang)
THESIS_ECOTYPES = {
    "atrophaeus": 8,
    "inaquosorum": 14,
    "spizizenii": 5,
}


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
            (c for c in tree.find_clades() if c.is_terminal() and c.name == og_id),
            None,
        )
        if outgroup:
            break
    if outgroup is None:
        return False
    tree.root_with_outgroup(outgroup)
    Phylo.write(tree, nwk_out, "newick")
    return True


def run_ecosim(fasta, tree, out_xml, threads):
    """
    Run EcoSim JAR.

    EcoSim writes output relative to its own directory (cwd), so we use
    absolute paths for sequences and phylogeny, and a temp filename for output
    that we move afterward.
    """
    abs_fasta = os.path.abspath(fasta)
    abs_tree = os.path.abspath(tree)
    tmp_name = os.path.basename(out_xml)
    tmp_path = os.path.join(ECOSIM_DIR, tmp_name)

    cmd = [
        JAVA, "-Xmx4G", "-jar", ECOSIM_JAR,
        f"-s={abs_fasta}",
        f"-p={abs_tree}",
        f"-o={tmp_name}",
        "-n",               # --nogui (implies --runall)
        "-d",               # --debug (verbose output)
        f"-t={threads}",    # --threads
    ]

    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=ECOSIM_DIR)
    elapsed = time.time() - t0
    stdout = result.stdout + "\n" + result.stderr

    success = result.returncode == 0 and os.path.isfile(tmp_path)
    if success:
        os.makedirs(os.path.dirname(out_xml), exist_ok=True)
        shutil.move(tmp_path, out_xml)

    # Parse results from XML if available, else from stdout
    npop = "?"
    n_ecotypes = 0

    if success and os.path.isfile(out_xml):
        with open(out_xml) as xf:
            xml_text = xf.read()
        # npop from hillclimb result
        m = re.search(r'<hillclimb>\s*<result\s+npop="(\d+)"', xml_text)
        if m:
            npop = m.group(1)
        # ecotype count from demarcation
        m = re.search(r'<ecotypes\s+size="(\d+)"', xml_text)
        if m:
            n_ecotypes = int(m.group(1))

    if npop == "?":
        # Fallback: parse stdout
        for line in stdout.splitlines():
            if "npop" in line.lower() and "running" not in line.lower():
                m = re.search(r'npop[:\s=]+(\d+)', line, re.IGNORECASE)
                if m:
                    npop = m.group(1)
                    break

    # Save full stdout for debugging
    log_txt = os.path.join(os.path.dirname(out_xml),
                           os.path.basename(out_xml).replace(".xml", "_ecosim_log.txt"))
    with open(log_txt, "w") as f:
        f.write(stdout)

    return success, elapsed, npop, n_ecotypes, stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pcr", type=int, default=200)
    ap.add_argument("--min-seqs", type=int, default=50)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-ecosim", action="store_true")
    ap.add_argument("--genes", default=None, help="Comma-separated gene list")
    args = ap.parse_args()

    # Output directory includes n_pcr so 200 and 2000 runs don't clobber each other
    out_dir = os.path.join(OUT_BASE, str(args.n_pcr))
    notes_dir = os.path.join(out_dir, "notes")
    os.makedirs(notes_dir, exist_ok=True)

    if not args.no_ecosim and not os.path.isfile(ECOSIM_JAR):
        print(f"ERROR: ecosim.jar not found: {ECOSIM_JAR}", file=sys.stderr)
        print(f"  Expected at: {ECOSIM_JAR}", file=sys.stderr)
        sys.exit(1)

    # Load gene->species mapping
    gene_species = {}
    with open(PRIMERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene_species[row["gene"]] = row["species"]

    # Determine which genes to process
    if args.genes:
        genes = [g.strip() for g in args.genes.split(",")]
    else:
        genes = sorted(USABLE_GENES)

    rng = random.Random(args.seed)

    print(f"EcoSim run: n_pcr={args.n_pcr}  seed={args.seed}  threads={args.threads}")
    print(f"EcoSim JAR: {ECOSIM_JAR}")
    print(f"Output dir: {out_dir}")
    print(f"Genes: {', '.join(genes)}")
    print()

    log_path = os.path.join(notes_dir, f"ecosim_results.tsv")
    log = open(log_path, "w")
    log.write("gene\tspecies\tpcr_total\tpcr_sampled\tcbp\toutgroup\ttotal_seqs\t"
              "tree_time\tecosim_time\tnpop\tn_ecotypes\tthesis_expected\tstatus\n")

    for gene in genes:
        species = gene_species.get(gene)
        if not species:
            print(f"  {gene:12s}  SKIP — not in primers.tsv")
            continue

        fasta_path = os.path.join(IN_DIR, f"{gene}.fasta")
        if not os.path.exists(fasta_path):
            print(f"  {gene:12s}  SKIP — no trimmed file")
            continue

        entries = read_fasta(fasta_path)

        # Separate by type
        outgroup = [(h, s) for h, s in entries if h in OUTGROUP_IDS]
        cbp = [(h, s) for h, s in entries if h.startswith("CBP-")]
        pcr = [(h, s) for h, s in entries
               if h not in OUTGROUP_IDS and not h.startswith("CBP-")]

        # Subsample PCR
        if len(pcr) > args.n_pcr:
            sampled = rng.sample(pcr, args.n_pcr)
        else:
            sampled = pcr

        # Outgroup first (EcoSim expects outgroup as first sequence)
        sub_entries = outgroup + cbp + sampled
        if len(sub_entries) < args.min_seqs:
            print(f"  {gene:12s} ({species:12s})  SKIP — only {len(sub_entries)} seqs")
            continue

        # Output paths
        gene_dir = os.path.join(out_dir, gene)

        # Shorten PCR headers — FastTree truncates at ':' causing duplicates
        # Keep outgroup and CBP headers as-is
        # PCR headers: extract env number (gene_ENV_...) and add unique index
        # Format: PCR_E{env}_{index} so we keep env traceability
        short_entries = []
        pcr_idx = 0
        header_map = []  # (short, original) for reference
        for h, s in sub_entries:
            if h in OUTGROUP_IDS or h.startswith("CBP-"):
                short_entries.append((h, s))
            else:
                pcr_idx += 1
                # Extract env from header like "acuA_2_AmpSeq-..." where 2 is the env
                parts = h.split("_", 2)
                env = parts[1] if len(parts) > 1 else "X"
                short_h = f"PCR_E{env}_{pcr_idx:05d}"
                short_entries.append((short_h, s))
                header_map.append((short_h, h))

        # Save header mapping for traceability
        os.makedirs(gene_dir, exist_ok=True)
        map_path = os.path.join(gene_dir, f"{gene}_header_map.tsv")
        with open(map_path, "w") as mf:
            mf.write("short_header\toriginal_header\n")
            for sh, oh in header_map:
                mf.write(f"{sh}\t{oh}\n")
        out_fasta = os.path.join(gene_dir, f"{gene}.fasta")
        out_nwk_raw = os.path.join(gene_dir, f"{gene}_unrooted.nwk")
        out_nwk = os.path.join(gene_dir, f"{gene}.nwk")
        out_xml = os.path.join(gene_dir, f"{gene}.xml")

        # Write subsampled FASTA
        write_fasta(out_fasta, short_entries)

        n_total = len(sub_entries)
        print(f"  {gene:12s} ({species:12s})  "
              f"PCR={len(pcr):>10,} -> {len(sampled):>4}  "
              f"CBP={len(cbp):>4}  OG={len(outgroup)}  total={n_total}",
              end="", flush=True)

        # Build tree
        t0 = time.time()
        try:
            build_tree(out_fasta, out_nwk_raw)
        except subprocess.CalledProcessError as e:
            print(f"  TREE FAIL: {e.stderr[:200] if e.stderr else e}")
            log.write(f"{gene}\t{species}\t{len(pcr)}\t{len(sampled)}\t{len(cbp)}\t"
                      f"{len(outgroup)}\t{n_total}\t\t\t\t\t"
                      f"{THESIS_ECOTYPES.get(species, '')}\tTREE_FAIL\n")
            continue
        tree_time = time.time() - t0

        # Reroot
        if outgroup:
            rooted = reroot_tree(out_nwk_raw, out_nwk)
            if not rooted:
                print(f"  REROOT FAIL (using unrooted)", end="", flush=True)
                shutil.copy(out_nwk_raw, out_nwk)
        else:
            shutil.copy(out_nwk_raw, out_nwk)

        print(f"  tree={tree_time:.1f}s", end="", flush=True)

        # Run EcoSim
        if args.no_ecosim:
            print("  (ecosim skipped)")
            log.write(f"{gene}\t{species}\t{len(pcr)}\t{len(sampled)}\t{len(cbp)}\t"
                      f"{len(outgroup)}\t{n_total}\t{tree_time:.1f}\t\t\t\t"
                      f"{THESIS_ECOTYPES.get(species, '')}\ttree_only\n")
            continue

        print("  running ecosim...", end="", flush=True)
        success, ecosim_time, npop, n_ecotypes, stdout = run_ecosim(
            out_fasta, out_nwk, out_xml, args.threads)

        expected = THESIS_ECOTYPES.get(species, "?")

        if success:
            print(f"  {ecosim_time:.0f}s  npop={npop}  ecotypes={n_ecotypes}  "
                  f"(thesis={expected})")
            status = "OK"
        else:
            print(f"  ECOSIM FAIL ({ecosim_time:.0f}s)")
            status = "FAIL"
            n_ecotypes = "FAIL"

        log.write(f"{gene}\t{species}\t{len(pcr)}\t{len(sampled)}\t{len(cbp)}\t"
                  f"{len(outgroup)}\t{n_total}\t{tree_time:.1f}\t{ecosim_time:.1f}\t"
                  f"{npop}\t{n_ecotypes}\t{expected}\t{status}\n")
        log.flush()

    log.close()
    print(f"\nResults written to {log_path}")


if __name__ == "__main__":
    main()
