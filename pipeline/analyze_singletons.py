#!/usr/bin/env python3
"""Analyze singleton composition: which are CBP vs PCR, and what thesis labels / envs."""

import xml.etree.ElementTree as ET
import os
from collections import Counter

BASE = os.path.join(os.path.dirname(__file__), "..", "data", "08_ecosim_ingroup")
GENES = ["sorA", "yvqK", "thiD", "albG", "acuA"]

for gene in GENES:
    tree = ET.parse(os.path.join(BASE, gene, f"{gene}.xml"))
    root = tree.getroot()
    ecotypes_elem = root.find("demarcation/ecotypes")

    # Collect all singletons
    cbp_singletons = []
    pcr_singletons = []

    # Collect all non-singletons
    cbp_in_nonsingle = []
    pcr_in_nonsingle = []

    for eco in ecotypes_elem.findall("ecotype"):
        size = int(eco.get("size"))
        members = [m.get("name") for m in eco.findall("member")]

        if size == 1:
            m = members[0]
            if m.startswith("CBP-"):
                # extract thesis label
                parts = m.split("_PE_")
                label = parts[1] if len(parts) == 2 else "?"
                cbp_singletons.append(label)
            else:
                # extract env
                parts = m.split("_")
                env = parts[1] if len(parts) >= 2 else "?"
                pcr_singletons.append(env)
        else:
            for m in members:
                if m.startswith("CBP-"):
                    parts = m.split("_PE_")
                    label = parts[1] if len(parts) == 2 else "?"
                    cbp_in_nonsingle.append(label)
                else:
                    parts = m.split("_")
                    env = parts[1] if len(parts) >= 2 else "?"
                    pcr_in_nonsingle.append(env)

    # Count total CBP and PCR
    phylo_size = int(root.find("phylogeny").get("size"))

    cbp_label_counts = Counter(cbp_singletons)
    pcr_env_counts = Counter(pcr_singletons)

    cbp_nonsingle_counts = Counter(cbp_in_nonsingle)
    pcr_nonsingle_counts = Counter(pcr_in_nonsingle)

    total_cbp = len(cbp_singletons) + len(cbp_in_nonsingle)
    total_pcr = len(pcr_singletons) + len(pcr_in_nonsingle)

    print(f"\n{'='*60}")
    print(f"{gene} — {phylo_size} seqs ({total_cbp} CBP, {total_pcr} PCR)")
    print(f"{'='*60}")

    print(f"\nCBP SINGLETONS: {len(cbp_singletons)} / {total_cbp} CBP ({100*len(cbp_singletons)/total_cbp:.0f}%)")
    for label in sorted(cbp_label_counts.keys()):
        s = cbp_label_counts[label]
        ns = cbp_nonsingle_counts.get(label, 0)
        tot = s + ns
        print(f"  {label}: {s} singleton / {tot} total ({100*s/tot:.0f}% singleton)")

    print(f"\nPCR SINGLETONS: {len(pcr_singletons)} / {total_pcr} PCR ({100*len(pcr_singletons)/total_pcr:.0f}%)")
    for env in sorted(pcr_env_counts.keys()):
        s = pcr_env_counts[env]
        ns = pcr_nonsingle_counts.get(env, 0)
        tot = s + ns
        print(f"  {env}: {s} singleton / {tot} total ({100*s/tot:.0f}% singleton)")
