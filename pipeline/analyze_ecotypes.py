#!/usr/bin/env python3
"""Parse EcoSim XMLs for inaquosorum genes and compute ecotype summary stats."""

import xml.etree.ElementTree as ET
import os

GENES = {
    "sorA":  {"type": "high",  "thesis_ecotypes": 41, "thesis_label": "41 in thesis (high)"},
    "thiD":  {"type": "high",  "thesis_ecotypes": 41, "thesis_label": "41 in thesis (high)"},
    "yvqK":  {"type": "high",  "thesis_ecotypes": 39, "thesis_label": "39 in thesis (high)"},
    "acuA":  {"type": "modal", "thesis_ecotypes": 14, "thesis_label": "14 in thesis (modal)"},
    "albG":  {"type": "modal", "thesis_ecotypes": 14, "thesis_label": "14 in thesis (modal)"},
}

BASE = os.path.join(os.path.dirname(__file__), "..", "data", "08_ecosim_ingroup")

for gene, info in GENES.items():
    xml_path = os.path.join(BASE, gene, f"{gene}.xml")
    tree = ET.parse(xml_path)
    root = tree.getroot()

    # Get hillclimb npop
    hc = root.find("hillclimb/result")
    npop = int(hc.get("npop"))
    omega = float(hc.get("omega"))
    sigma = float(hc.get("sigma"))
    likelihood = float(hc.get("likelihood"))

    # Get npop CI
    npop_ci_lower = int(root.find("npopCI/lower").get("value"))
    npop_ci_upper = int(root.find("npopCI/upper").get("value"))

    # Get phylogeny size
    phylo_size = int(root.find("phylogeny").get("size"))
    seq_len = int(root.find("phylogeny").get("length"))

    # Parse ecotypes
    demarcation = root.find("demarcation")
    ecotypes_elem = demarcation.find("ecotypes")
    total_ecotypes = int(ecotypes_elem.get("size"))

    singletons = 0
    doubletons = 0
    larger = 0
    mixed_ecotypes = []  # ecotypes with both CBP and PCR
    cbp_only_ecotypes = 0
    pcr_only_ecotypes = 0
    all_sizes = []

    for eco in ecotypes_elem.findall("ecotype"):
        eco_num = int(eco.get("number"))
        eco_size = int(eco.get("size"))
        all_sizes.append(eco_size)

        members = [m.get("name") for m in eco.findall("member")]
        cbp_members = [m for m in members if m.startswith("CBP-")]
        pcr_members = [m for m in members if m.startswith("PCR_")]

        has_cbp = len(cbp_members) > 0
        has_pcr = len(pcr_members) > 0

        if eco_size == 1:
            singletons += 1
        elif eco_size == 2:
            doubletons += 1
        else:
            larger += 1

        if has_cbp and has_pcr:
            # Get thesis ecotype labels from CBP members
            thesis_labels = set()
            for m in cbp_members:
                parts = m.split("_PE_")
                if len(parts) == 2:
                    thesis_labels.add(parts[1])

            # Get environments from PCR members
            envs = set()
            for m in pcr_members:
                parts = m.split("_")
                if len(parts) >= 2:
                    envs.add(parts[1])

            mixed_ecotypes.append({
                "num": eco_num,
                "size": eco_size,
                "cbp": len(cbp_members),
                "pcr": len(pcr_members),
                "thesis_labels": sorted(thesis_labels),
                "envs": sorted(envs),
            })
        elif has_cbp:
            cbp_only_ecotypes += 1
        elif has_pcr:
            pcr_only_ecotypes += 1

    # Count singleton composition
    singleton_cbp = 0
    singleton_pcr = 0
    for eco in ecotypes_elem.findall("ecotype"):
        if int(eco.get("size")) == 1:
            m = eco.find("member").get("name")
            if m.startswith("CBP-"):
                singleton_cbp += 1
            else:
                singleton_pcr += 1

    print(f"\n{'='*70}")
    print(f"GENE: {gene} ({info['type']} ecotype gene)")
    print(f"{'='*70}")
    print(f"Sequences: {phylo_size} | Alignment length: {seq_len} bp")
    print(f"HillClimb npop: {npop} (CI: {npop_ci_lower}–{npop_ci_upper})")
    print(f"Omega: {omega:.5f} | Sigma: {sigma:.5f} | Likelihood: {likelihood:.5f}")
    print(f"Thesis expects: {info['thesis_ecotypes']} ecotypes")
    print(f"Demarcated ecotypes: {total_ecotypes}")
    print(f"  Singletons (size=1): {singletons} ({singletons*100/total_ecotypes:.0f}%)")
    print(f"    - CBP singletons: {singleton_cbp}")
    print(f"    - PCR singletons: {singleton_pcr}")
    print(f"  Doubletons (size=2): {doubletons}")
    print(f"  Larger (size>2):     {larger}")
    print(f"  CBP-only ecotypes:   {cbp_only_ecotypes}")
    print(f"  PCR-only ecotypes:   {pcr_only_ecotypes}")
    print(f"  Mixed (CBP+PCR):     {len(mixed_ecotypes)}")
    print(f"\nMixed ecotypes detail:")
    for me in mixed_ecotypes:
        labels = ", ".join(me["thesis_labels"])
        envs = ", ".join(me["envs"])
        print(f"  Ecotype #{me['num']}: {me['size']} seqs "
              f"({me['cbp']} CBP [{labels}] + {me['pcr']} PCR [{envs}])")

    # Largest ecotypes
    sorted_sizes = sorted(all_sizes, reverse=True)
    print(f"\nLargest ecotypes: {sorted_sizes[:5]}")
