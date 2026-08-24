#!/usr/bin/env python3
"""
Step 8: Run EcoSim on pruned in-group data from step 7.

For inaquosorum genes: data is unchanged (all PCR inside LCA).
For spizizenii genes: CBP-only (all PCR were outside LCA and removed).

Uses existing trees from step 7 (already rebuilt on pruned data).

Input:  data/07_ingroup/{gene}/{gene}.fasta + .nwk
Output: data/08_ecosim_ingroup/{gene}/{gene}.xml
"""

import csv
import os
import re
import shutil
import subprocess
import time

INGROUP_DIR = "data/07_ingroup"
OUT_DIR = "data/08_ecosim_ingroup"
NOTES_DIR = os.path.join(OUT_DIR, "notes")
PRIMERS_TSV = "notes/primers.tsv"

ECOSIM_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "ecosim")
)
ECOSIM_JAR = os.path.join(ECOSIM_DIR, "ecosim.jar")
JAVA = "/usr/local/opt/openjdk/bin/java"

OUTGROUP_IDS = {"FN597644.1", "CP026362.1"}
USABLE_GENES = {"acuA", "sorA", "yvqK", "albG", "thiD", "iolB", "acsA_2", "alkH"}

THESIS_ECOTYPES = {
    "atrophaeus": 8,
    "inaquosorum": 14,
    "spizizenii": 5,
}


def count_seqs(fasta_path):
    cbp = pcr = og = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip().split()[0]
                if h in OUTGROUP_IDS:
                    og += 1
                elif h.startswith("CBP-"):
                    cbp += 1
                else:
                    pcr += 1
    return cbp, pcr, og


def run_ecosim(fasta, tree, out_xml, threads=8):
    abs_fasta = os.path.abspath(fasta)
    abs_tree = os.path.abspath(tree)
    tmp_name = os.path.basename(out_xml)
    tmp_path = os.path.join(ECOSIM_DIR, tmp_name)

    cmd = [
        JAVA, "-Xmx4G", "-jar", ECOSIM_JAR,
        f"-s={abs_fasta}",
        f"-p={abs_tree}",
        f"-o={tmp_name}",
        "-n", "-d", f"-t={threads}",
    ]

    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=ECOSIM_DIR)
    elapsed = time.time() - t0
    stdout = result.stdout + "\n" + result.stderr

    success = result.returncode == 0 and os.path.isfile(tmp_path)
    if success:
        os.makedirs(os.path.dirname(out_xml), exist_ok=True)
        shutil.move(tmp_path, out_xml)

    # Parse from XML
    npop = "?"
    n_ecotypes = 0
    if success and os.path.isfile(out_xml):
        with open(out_xml) as xf:
            xml_text = xf.read()
        m = re.search(r'<hillclimb>\s*<result\s+npop="(\d+)"', xml_text)
        if m:
            npop = m.group(1)
        m = re.search(r'<ecotypes\s+size="(\d+)"', xml_text)
        if m:
            n_ecotypes = int(m.group(1))

    # Save log
    log_txt = out_xml.replace(".xml", "_ecosim_log.txt")
    with open(log_txt, "w") as f:
        f.write(stdout)

    return success, elapsed, npop, n_ecotypes


def main():
    os.makedirs(NOTES_DIR, exist_ok=True)

    gene_species = {}
    with open(PRIMERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene_species[row["gene"]] = row["species"]

    log_path = os.path.join(NOTES_DIR, "ecosim_ingroup_results.tsv")
    log = open(log_path, "w")
    log.write("gene\tspecies\tcbp\tpcr\toutgroup\ttotal\t"
              "ecosim_time\tnpop\tn_ecotypes\tthesis_expected\tstatus\n")

    print("Step 8: EcoSim on pruned in-group")
    print(f"EcoSim JAR: {ECOSIM_JAR}")
    print()

    for gene in sorted(USABLE_GENES):
        species = gene_species.get(gene, "?")
        fasta_path = os.path.join(INGROUP_DIR, gene, f"{gene}.fasta")
        nwk_path = os.path.join(INGROUP_DIR, gene, f"{gene}.nwk")

        if not os.path.exists(fasta_path):
            print(f"  {gene:12s}  SKIP — no pruned file")
            continue

        cbp, pcr, og = count_seqs(fasta_path)
        total = cbp + pcr + og
        expected = THESIS_ECOTYPES.get(species, "?")

        # Copy fasta and tree to output dir
        gene_dir = os.path.join(OUT_DIR, gene)
        os.makedirs(gene_dir, exist_ok=True)
        out_fasta = os.path.join(gene_dir, f"{gene}.fasta")
        out_nwk = os.path.join(gene_dir, f"{gene}.nwk")
        out_xml = os.path.join(gene_dir, f"{gene}.xml")

        shutil.copy(fasta_path, out_fasta)
        shutil.copy(nwk_path, out_nwk)

        print(f"  {gene:12s} ({species:12s})  CBP={cbp:>4}  PCR={pcr:>4}  OG={og}  total={total}",
              end="", flush=True)

        print("  running ecosim...", end="", flush=True)
        success, elapsed, npop, n_ecotypes = run_ecosim(out_fasta, out_nwk, out_xml)

        if success:
            print(f"  {elapsed:.0f}s  npop={npop}  ecotypes={n_ecotypes}  (thesis={expected})")
            status = "OK"
        else:
            print(f"  ECOSIM FAIL ({elapsed:.0f}s)")
            status = "FAIL"
            n_ecotypes = "FAIL"

        log.write(f"{gene}\t{species}\t{cbp}\t{pcr}\t{og}\t{total}\t"
                  f"{elapsed:.1f}\t{npop}\t{n_ecotypes}\t{expected}\t{status}\n")
        log.flush()

    log.close()
    print(f"\nResults: {log_path}")


if __name__ == "__main__":
    main()
