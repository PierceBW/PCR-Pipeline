#!/usr/bin/env python3
"""
Generate a combined PDF report with tree visualizations + analysis text
for sharing with professors.
"""

import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from Bio import Phylo

IN_DIR = "data/05_ecosim/200"
PRIMERS_TSV = "notes/primers.tsv"
OUTGROUP_IDS = {"FN597644.1", "CP026362.1"}

# Same color scheme as step6
ENV_COLORS = {"E1": "#1f77b4", "E2": "#17becf", "E4": "#2ca02c", "E5": "#9467bd"}
ENV_DEFAULT = "#7f7f7f"

INAQ_COLORS = {}
_cmap = matplotlib.colormaps["tab20"]
for i in range(14):
    INAQ_COLORS[f"I{i+1}"] = matplotlib.colors.rgb2hex(_cmap(i / 20))

SPIZ_COLORS = {}
_cmap2 = matplotlib.colormaps["Set1"]
for i in range(5):
    SPIZ_COLORS[f"S{i+1}"] = matplotlib.colors.rgb2hex(_cmap2(i / 9))

CBP_COLORS = {**INAQ_COLORS, **SPIZ_COLORS}
CBP_DEFAULT = "#d62728"
OUTSIDE_COLOR = "#ff7f0e"
OUTGROUP_COLOR = "#000000"

THESIS_ECOTYPES = {"inaquosorum": 14, "spizizenii": 5}

USABLE_GENES = ["acuA", "sorA", "yvqK", "albG", "thiD", "acsA_2", "alkH", "iolB"]

# Pruned EcoSim results (from step 8)
PRUNED_RESULTS = {
    "acsA_2": {"npop": 2, "ecotypes": 1},
    "acuA": {"npop": 115, "ecotypes": 208},
    "albG": {"npop": 37, "ecotypes": 158},
    "alkH": {"npop": 6, "ecotypes": 34},
    "iolB": {"npop": 14, "ecotypes": 86},
    "sorA": {"npop": 16, "ecotypes": 58},
    "thiD": {"npop": 43, "ecotypes": 124},
    "yvqK": {"npop": 12, "ecotypes": 56},
}


def get_leaf_type(name):
    if name in OUTGROUP_IDS:
        return "outgroup"
    if name.startswith("CBP-"):
        return "cbp"
    return "pcr"


def get_cbp_label(name):
    m = re.search(r"_PE_([A-Z]\d+)", name)
    return m.group(1) if m else "?"


def get_pcr_env(name):
    m = re.search(r"PCR_E(\d+)_", name)
    return f"E{m.group(1)}" if m else "?"


def get_leaf_color(name, outside_lca=False):
    ltype = get_leaf_type(name)
    if ltype == "outgroup":
        return OUTGROUP_COLOR
    if ltype == "cbp":
        return CBP_COLORS.get(get_cbp_label(name), CBP_DEFAULT)
    if outside_lca:
        return OUTSIDE_COLOR
    return ENV_COLORS.get(get_pcr_env(name), ENV_DEFAULT)


def compute_y_positions(clade, y_pos, y_step, positions):
    if clade.is_terminal():
        positions[clade] = y_pos[0]
        y_pos[0] += y_step
    else:
        child_ys = []
        for child in clade.clades:
            compute_y_positions(child, y_pos, y_step, positions)
            child_ys.append(positions[child])
        positions[clade] = (min(child_ys) + max(child_ys)) / 2


def compute_x_positions(clade, x_start, positions):
    positions[clade] = x_start
    if not clade.is_terminal():
        for child in clade.clades:
            bl = child.branch_length if child.branch_length else 0
            compute_x_positions(child, x_start + bl, positions)


def draw_tree(tree, outside_set, lca_clade, ax):
    terminals = tree.get_terminals()
    n_leaves = len(terminals)

    y_positions = {}
    x_positions = {}
    compute_y_positions(tree.root, [0], 1.0, y_positions)
    compute_x_positions(tree.root, 0, x_positions)
    max_x = max(x_positions.values()) if x_positions else 1

    # Branches
    for clade in tree.find_clades(order="level"):
        cx, cy = x_positions[clade], y_positions[clade]
        for child in clade.clades:
            child_x, child_y = x_positions[child], y_positions[child]
            ax.plot([cx, child_x], [child_y, child_y], color="#444444", lw=0.3, solid_capstyle="butt")
            ax.plot([cx, cx], [cy, child_y], color="#444444", lw=0.3, solid_capstyle="butt")

    # LCA highlight
    if lca_clade is not None:
        lca_ys = [y_positions[t] for t in lca_clade.get_terminals()]
        lca_x = x_positions[lca_clade]
        rect = mpatches.Rectangle(
            (lca_x, min(lca_ys) - 0.5), max_x - lca_x + max_x * 0.15,
            max(lca_ys) - min(lca_ys) + 1,
            lw=0.5, ec="#cccccc", fc="#f0f0ff", alpha=0.4, zorder=0
        )
        ax.add_patch(rect)

    # Leaves
    label_off = max_x * 0.01
    for t in terminals:
        tx, ty = x_positions[t], y_positions[t]
        outside = t.name in outside_set
        color = get_leaf_color(t.name, outside_lca=outside)
        ltype = get_leaf_type(t.name)
        marker = "o" if ltype == "pcr" else ("s" if ltype == "cbp" else "D")
        size = 4 if ltype == "outgroup" else 2
        ax.plot(tx, ty, marker=marker, color=color, markersize=size, zorder=5)

        if ltype == "cbp":
            ax.text(tx + label_off, ty, get_cbp_label(t.name), fontsize=1.2, va="center", color=color, zorder=6)
        elif ltype == "outgroup":
            ax.text(tx + label_off, ty, t.name, fontsize=2, va="center", color=color, fontweight="bold", zorder=6)

    ax.set_xlim(-max_x * 0.02, max_x * 1.2)
    ax.set_ylim(-1, n_leaves + 1)
    ax.invert_yaxis()
    ax.set_xlabel("Substitutions per site", fontsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])


def build_legend(ax, species):
    handles = [mpatches.Patch(color=OUTGROUP_COLOR, label="Outgroup")]
    if species == "inaquosorum":
        for i in range(1, 15):
            handles.append(mpatches.Patch(color=INAQ_COLORS[f"I{i}"], label=f"CBP I{i}"))
    elif species == "spizizenii":
        for i in range(1, 6):
            handles.append(mpatches.Patch(color=SPIZ_COLORS[f"S{i}"], label=f"CBP S{i}"))
    for env in ["E1", "E2", "E4", "E5"]:
        handles.append(mpatches.Patch(color=ENV_COLORS[env], label=f"PCR {env}"))
    handles.append(mpatches.Patch(color=OUTSIDE_COLOR, label="PCR outside LCA"))
    handles.append(mpatches.Patch(fc="#f0f0ff", ec="#cccccc", alpha=0.4, label="CBP LCA clade"))
    ax.legend(handles=handles, loc="upper left", fontsize=4, ncol=2, framealpha=0.8, borderpad=0.5)


def add_text_page(pdf, lines, fontsize=10):
    """Add a page of formatted text."""
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    text = "\n".join(lines)
    ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=fontsize,
            verticalalignment="top", fontfamily="monospace", wrap=True)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs("data/06_viz", exist_ok=True)
    out_path = "data/06_viz/tree_analysis_report.pdf"

    gene_species = {}
    with open(PRIMERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene_species[row["gene"]] = row["species"]

    with PdfPages(out_path) as pdf:
        # --- Title page ---
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")
        ax.text(0.5, 0.65, "PCR-Pipeline Tree Analysis", transform=ax.transAxes,
                fontsize=24, ha="center", fontweight="bold")
        ax.text(0.5, 0.55, "CBP In-Group vs Environmental PCR Reads",
                transform=ax.transAxes, fontsize=16, ha="center")
        ax.text(0.5, 0.45, "8 genes  |  n=200 PCR subsets  |  EcoSim 2.1.7",
                transform=ax.transAxes, fontsize=12, ha="center", color="#666666")
        ax.text(0.5, 0.30, "May 2026", transform=ax.transAxes,
                fontsize=12, ha="center", color="#666666")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # --- Overview page ---
        add_text_page(pdf, [
            "OVERVIEW",
            "=" * 70,
            "",
            "Question: Are the environmental PCR reads inside or outside",
            "the Jocelyn CBP in-group (least common ancestor clade)?",
            "",
            "Method:",
            "  1. Build FastTree on 200 PCR + all CBP + outgroup per gene",
            "  2. Find LCA of all CBP isolates",
            "  3. Count PCR reads inside vs outside that clade",
            "  4. Re-run EcoSim on pruned in-group",
            "",
            "RESULT: Clean species-level split",
            "-" * 40,
            "",
            "  INAQUOSORUM (acuA, albG, sorA, thiD, yvqK):",
            "    100% of PCR reads INSIDE the CBP LCA",
            "    -> Environmental reads ARE inaquosorum",
            "",
            "  SPIZIZENII (acsA_2, alkH, iolB):",
            "    100% of PCR reads OUTSIDE the CBP LCA",
            "    -> Environmental reads are a DIFFERENT TAXON",
            "",
            "This is binary: no gene shows partial overlap.",
        ], fontsize=11)

        # --- Summary table page ---
        add_text_page(pdf, [
            "SUMMARY TABLE",
            "=" * 70,
            "",
            f"{'Gene':12s} {'Species':12s} {'CBP':>5s} {'PCR':>5s} {'OG':>3s} {'In LCA':>7s} {'Out LCA':>8s} {'npop':>5s} {'Eco':>5s} {'Thesis':>7s}",
            "-" * 70,
            f"{'acuA':12s} {'inaquosorum':12s} {'141':>5s} {'200':>5s} {'1':>3s} {'200':>7s} {'0':>8s} {'115':>5s} {'208':>5s} {'14':>7s}",
            f"{'albG':12s} {'inaquosorum':12s} {'141':>5s} {'200':>5s} {'0':>3s} {'200':>7s} {'0':>8s} {'37':>5s} {'158':>5s} {'14':>7s}",
            f"{'sorA':12s} {'inaquosorum':12s} {'141':>5s} {'200':>5s} {'1':>3s} {'200':>7s} {'0':>8s} {'16':>5s} {'58':>5s} {'14':>7s}",
            f"{'thiD':12s} {'inaquosorum':12s} {'142':>5s} {'200':>5s} {'1':>3s} {'200':>7s} {'0':>8s} {'43':>5s} {'124':>5s} {'14':>7s}",
            f"{'yvqK':12s} {'inaquosorum':12s} {'141':>5s} {'200':>5s} {'1':>3s} {'200':>7s} {'0':>8s} {'12':>5s} {'56':>5s} {'14':>7s}",
            "-" * 70,
            f"{'acsA_2':12s} {'spizizenii':12s} {'209':>5s} {'200':>5s} {'1':>3s} {'0':>7s} {'200':>8s} {'2':>5s} {'1':>5s} {'5':>7s}",
            f"{'alkH':12s} {'spizizenii':12s} {'209':>5s} {'200':>5s} {'1':>3s} {'0':>7s} {'200':>8s} {'6':>5s} {'34':>5s} {'5':>7s}",
            f"{'iolB':12s} {'spizizenii':12s} {'211':>5s} {'200':>5s} {'1':>3s} {'0':>7s} {'200':>8s} {'14':>5s} {'86':>5s} {'5':>7s}",
            "",
            "npop = estimated number of populations (EcoSim hillclimb)",
            "Eco  = demarcated ecotype count",
            "In/Out LCA = PCR reads inside/outside the CBP ancestor clade",
            "",
            "Closest to thesis: alkH npop=6 (thesis=5),",
            "                   sorA npop=16 (thesis=14),",
            "                   yvqK npop=12 (thesis=14)",
        ], fontsize=10)

        # --- Per-gene: text + tree pages ---
        gene_notes = {
            "acuA": [
                "acuA (inaquosorum) - 342 sequences",
                "-" * 50,
                "CBP: 141  |  PCR: 200  |  Outgroup: 1",
                "PCR inside LCA: 200 (100%)  |  Outside: 0",
                "",
                "All 200 PCR reads fall within the CBP in-group.",
                "Environmental diversity is fully interleaved with",
                "Jocelyn's known ecotypes.",
                "",
                "Pruned EcoSim: npop=115, ecotypes=208 (thesis=14)",
            ],
            "sorA": [
                "sorA (inaquosorum) - 342 sequences  ** BEST GENE **",
                "-" * 50,
                "CBP: 141  |  PCR: 200  |  Outgroup: 1",
                "PCR inside LCA: 200 (100%)  |  Outside: 0",
                "",
                "Most CBP/PCR interleaving of any gene. 9 mixed",
                "ecotypes where CBP and PCR cluster together.",
                "Thesis ecotypes I3,I4,I5,I7,I8,I9,I10 each have",
                "environmental reads alongside them.",
                "",
                "Pruned EcoSim: npop=16, ecotypes=58 (thesis=14)",
                "npop=16 is close to thesis expectation of 14.",
            ],
            "yvqK": [
                "yvqK (inaquosorum) - 342 sequences",
                "-" * 50,
                "CBP: 141  |  PCR: 200  |  Outgroup: 1",
                "PCR inside LCA: 200 (100%)  |  Outside: 0",
                "",
                "3 mixed ecotypes. PCR reads cluster with I5, I7,",
                "and I2/I4/I6/I9/I11 groups.",
                "",
                "Pruned EcoSim: npop=12, ecotypes=56 (thesis=14)",
                "npop=12 is close to thesis expectation of 14.",
            ],
            "albG": [
                "albG (inaquosorum) - 341 sequences",
                "-" * 50,
                "CBP: 141  |  PCR: 200  |  NO outgroup",
                "PCR inside LCA: 200 (100%)  |  Outside: 0",
                "",
                "All PCR inside despite no outgroup (unrooted tree).",
                "",
                "Pruned EcoSim: npop=37, ecotypes=158 (thesis=14)",
            ],
            "thiD": [
                "thiD (inaquosorum) - 343 sequences",
                "-" * 50,
                "CBP: 142  |  PCR: 200  |  Outgroup: 1",
                "PCR inside LCA: 200 (100%)  |  Outside: 0",
                "",
                "3 mixed ecotypes at n=200 (improved to 9 at n=2000).",
                "",
                "Pruned EcoSim: npop=43, ecotypes=124 (thesis=14)",
            ],
            "acsA_2": [
                "acsA_2 (spizizenii) - 410 sequences",
                "-" * 50,
                "CBP: 209  |  PCR: 200  |  Outgroup: 1",
                "PCR inside LCA: 0 (0%)  |  Outside: 200 (100%)",
                "Excluded envs: ALL from E2",
                "",
                "** ALL PCR READS OUTSIDE CBP CLADE **",
                "CBP isolates form a monophyletic clade with zero",
                "PCR reads inside. All 200 PCR reads (all from E2)",
                "branch off separately - a different lineage/taxon.",
                "",
                "Pruned EcoSim (CBP-only): npop=2, ecotypes=1",
            ],
            "alkH": [
                "alkH (spizizenii) - 410 sequences",
                "-" * 50,
                "CBP: 209  |  PCR: 200  |  Outgroup: 1",
                "PCR inside LCA: 0 (0%)  |  Outside: 200 (100%)",
                "Excluded envs: E1:2, E2:26, E4:96, E5:76",
                "",
                "** ALL PCR READS OUTSIDE CBP CLADE **",
                "Same separation as acsA_2 but PCR from ALL envs.",
                "Rules out single-environment contamination.",
                "",
                "Pruned EcoSim (CBP-only): npop=6, ecotypes=34",
                "npop=6 is closest to thesis expectation (5)!",
            ],
            "iolB": [
                "iolB (spizizenii) - 412 sequences",
                "-" * 50,
                "CBP: 211  |  PCR: 200  |  Outgroup: 1",
                "PCR inside LCA: 0 (0%)  |  Outside: 200 (100%)",
                "Excluded envs: E1:21, E2:43, E4:129, E5:7",
                "",
                "** ALL PCR READS OUTSIDE CBP CLADE **",
                "Same pattern. PCR predominantly from E4 (129/200).",
                "",
                "Pruned EcoSim (CBP-only): npop=14, ecotypes=86",
            ],
        }

        for gene in USABLE_GENES:
            species = gene_species.get(gene, "?")
            nwk_path = os.path.join(IN_DIR, gene, f"{gene}.nwk")
            if not os.path.exists(nwk_path):
                continue

            # Text page
            notes = gene_notes.get(gene, [f"{gene} ({species})"])
            add_text_page(pdf, notes, fontsize=12)

            # Tree page
            tree = Phylo.read(nwk_path, "newick")
            terminals = tree.get_terminals()
            cbp_leaves = [t for t in terminals if t.name.startswith("CBP-")]
            pcr_leaves = [t for t in terminals if t.name.startswith("PCR_")]

            lca_clade = None
            if len(cbp_leaves) >= 2:
                lca_clade = tree.common_ancestor(cbp_leaves)

            lca_names = {t.name for t in lca_clade.get_terminals()} if lca_clade else set()
            outside_set = {t.name for t in pcr_leaves if t.name not in lca_names}

            n_leaves = len(terminals)
            fig_h = max(8, n_leaves * 0.04)
            fig, ax = plt.subplots(figsize=(11, fig_h))

            draw_tree(tree, outside_set, lca_clade, ax)
            build_legend(ax, species)

            pct_out = 100 * len(outside_set) / len(pcr_leaves) if pcr_leaves else 0
            title = (f"{gene} ({species}) — {n_leaves} seqs — "
                     f"PCR inside LCA: {len(pcr_leaves) - len(outside_set)}, "
                     f"outside: {len(outside_set)} ({pct_out:.0f}%)")
            ax.set_title(title, fontsize=8, fontweight="bold")

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # --- Conclusions page ---
        add_text_page(pdf, [
            "CONCLUSIONS",
            "=" * 70,
            "",
            "1. SPIZIZENII PCR READS ARE NOT SPIZIZENII",
            "   All 3 spizizenii genes show 100% CBP/PCR separation.",
            "   The PCR reads pass 95% BLAST identity but are",
            "   phylogenetically distinct from ALL known spizizenii",
            "   isolates. They may be a closely related unnamed species.",
            "",
            "2. INAQUOSORUM PCR READS ARE INAQUOSORUM",
            "   All 5 inaquosorum genes show 100% PCR inside CBP LCA.",
            "   Environmental amplicon diversity is genuine inaquosorum",
            "   diversity, interleaved with Jocelyn's culture collection.",
            "",
            "3. NPOP IS MORE MEANINGFUL THAN ECOTYPE COUNT",
            "   The high ecotype counts reflect EcoSim over-splitting",
            "   with added environmental diversity. But npop estimates",
            "   are in the right range:",
            "     alkH:  npop=6   (thesis=5)",
            "     yvqK:  npop=12  (thesis=14)",
            "     sorA:  npop=16  (thesis=14)",
            "",
            "4. NEXT STEPS",
            "   - For inaquosorum: the environmental reads are valid.",
            "     Consider what the extra ecotypes represent biologically.",
            "   - For spizizenii: investigate what taxon the PCR reads",
            "     actually belong to. A BLAST against NCBI nt could help.",
            "   - Consider running EcoSim on CBP-only for all genes to",
            "     reproduce Jocelyn's thesis results as a baseline.",
        ], fontsize=11)

    print(f"Report written to: {out_path}")


if __name__ == "__main__":
    main()
