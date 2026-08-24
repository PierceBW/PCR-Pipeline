#!/usr/bin/env python3
"""
Step 6: Visualize phylogenetic trees with CBP/PCR coloring and LCA marking.

For each gene, renders a rectangular phylogram as PDF/SVG:
- CBP isolates colored by thesis ecotype label (warm palette)
- PCR reads colored by environment (cool palette)
- Outgroup in black
- CBP LCA clade highlighted with a background bar
- PCR reads outside the LCA clade marked distinctly

Also produces a summary of PCR reads inside vs outside the CBP LCA per gene.

Input:  data/05_ecosim/200/{gene}/{gene}.nwk
Output: data/06_viz/{gene}_tree.pdf
        data/06_viz/lca_summary.tsv

Usage:
    ./venv/bin/python pipeline/step6_visualize_trees.py
"""

import csv
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from Bio import Phylo

IN_DIR = "data/05_ecosim/200"
OUT_DIR = "data/06_viz"
PRIMERS_TSV = "notes/primers.tsv"
OUTGROUP_IDS = {"FN597644.1", "CP026362.1"}

# Colors for PCR environments
ENV_COLORS = {
    "E1": "#1f77b4",  # blue
    "E2": "#17becf",  # cyan
    "E4": "#2ca02c",  # green
    "E5": "#9467bd",  # purple
}
ENV_DEFAULT = "#7f7f7f"  # grey for unknown env

# Colors for thesis ecotype labels (CBP)
# Inaquosorum I1-I14
INAQ_COLORS = {}
_inaq_cmap = plt.cm.get_cmap("tab20", 14)
for i in range(14):
    INAQ_COLORS[f"I{i+1}"] = matplotlib.colors.rgb2hex(_inaq_cmap(i))

# Spizizenii S1-S5
SPIZ_COLORS = {}
_spiz_cmap = plt.cm.get_cmap("Set1", 5)
for i in range(5):
    SPIZ_COLORS[f"S{i+1}"] = matplotlib.colors.rgb2hex(_spiz_cmap(i))

CBP_COLORS = {**INAQ_COLORS, **SPIZ_COLORS}
CBP_DEFAULT = "#d62728"  # red for unknown label

# Color for PCR outside LCA
OUTSIDE_COLOR = "#ff7f0e"  # orange — visually distinct "suspect" color
OUTGROUP_COLOR = "#000000"


def get_leaf_type(name):
    """Classify a leaf as 'outgroup', 'cbp', or 'pcr'."""
    if name in OUTGROUP_IDS:
        return "outgroup"
    if name.startswith("CBP-"):
        return "cbp"
    return "pcr"


def get_cbp_label(name):
    """Extract thesis ecotype label from CBP header like CBP-1234_PE_I7."""
    m = re.search(r"_PE_([A-Z]\d+)", name)
    return m.group(1) if m else "?"


def get_pcr_env(name):
    """Extract environment from PCR header like PCR_E2_00041."""
    m = re.search(r"PCR_E(\d+)_", name)
    return f"E{m.group(1)}" if m else "?"


def get_leaf_color(name, outside_lca=False):
    """Get color for a leaf node."""
    ltype = get_leaf_type(name)
    if ltype == "outgroup":
        return OUTGROUP_COLOR
    if ltype == "cbp":
        label = get_cbp_label(name)
        return CBP_COLORS.get(label, CBP_DEFAULT)
    # PCR
    if outside_lca:
        return OUTSIDE_COLOR
    env = get_pcr_env(name)
    return ENV_COLORS.get(env, ENV_DEFAULT)


def compute_y_positions(clade, y_pos, y_step, positions):
    """Recursively compute y positions for all nodes (leaves get sequential y)."""
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
    """Recursively compute x positions (cumulative branch length from root)."""
    positions[clade] = x_start
    if not clade.is_terminal():
        for child in clade.clades:
            bl = child.branch_length if child.branch_length else 0
            compute_x_positions(child, x_start + bl, positions)


def draw_tree(tree, outside_lca_set, lca_clade, gene, species, ax):
    """Draw a rectangular phylogram with colored leaves and LCA highlighting."""
    terminals = tree.get_terminals()
    n_leaves = len(terminals)

    # Compute positions
    y_positions = {}
    x_positions = {}
    compute_y_positions(tree.root, [0], 1.0, y_positions)
    compute_x_positions(tree.root, 0, x_positions)

    max_x = max(x_positions.values()) if x_positions else 1

    # Draw branches
    lines = []
    colors = []
    for clade in tree.find_clades(order="level"):
        if clade == tree.root:
            continue
        parent = tree.root  # find parent
        # Bio.Phylo doesn't expose parent directly, so use a pre-built map

    # Build parent map
    parent_map = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parent_map[child] = clade

    # Draw horizontal + vertical branch lines
    for clade in tree.find_clades(order="level"):
        cx = x_positions[clade]
        cy = y_positions[clade]

        for child in clade.clades:
            child_x = x_positions[child]
            child_y = y_positions[child]

            # Horizontal line (parent x to child x, at child y)
            ax.plot([cx, child_x], [child_y, child_y], color="#444444", linewidth=0.3, solid_capstyle="butt")
            # Vertical line (parent y to child y, at parent x)
            ax.plot([cx, cx], [cy, child_y], color="#444444", linewidth=0.3, solid_capstyle="butt")

    # Highlight LCA clade region
    if lca_clade is not None:
        lca_terminals = lca_clade.get_terminals()
        lca_ys = [y_positions[t] for t in lca_terminals]
        lca_ymin = min(lca_ys) - 0.5
        lca_ymax = max(lca_ys) + 0.5
        lca_x = x_positions[lca_clade]
        rect = mpatches.Rectangle(
            (lca_x, lca_ymin), max_x - lca_x + max_x * 0.15, lca_ymax - lca_ymin,
            linewidth=0.5, edgecolor="#cccccc", facecolor="#f0f0ff", alpha=0.4, zorder=0
        )
        ax.add_patch(rect)

    # Draw leaf markers and labels
    label_offset = max_x * 0.01
    for terminal in terminals:
        tx = x_positions[terminal]
        ty = y_positions[terminal]
        outside = terminal.name in outside_lca_set
        color = get_leaf_color(terminal.name, outside_lca=outside)
        ltype = get_leaf_type(terminal.name)

        # Marker
        marker = "o" if ltype == "pcr" else ("s" if ltype == "cbp" else "D")
        size = 4 if ltype == "outgroup" else 2
        ax.plot(tx, ty, marker=marker, color=color, markersize=size, zorder=5)

        # Label (only for CBP and outgroup — too many PCR to label)
        if ltype == "cbp":
            label = get_cbp_label(terminal.name)
            ax.text(tx + label_offset, ty, label, fontsize=1.5, va="center", color=color, zorder=6)
        elif ltype == "outgroup":
            ax.text(tx + label_offset, ty, terminal.name, fontsize=2, va="center",
                    color=color, fontweight="bold", zorder=6)

    # Axes formatting
    ax.set_xlim(-max_x * 0.02, max_x * 1.2)
    ax.set_ylim(-1, n_leaves + 1)
    ax.invert_yaxis()
    ax.set_xlabel("Substitutions per site", fontsize=6)
    ax.set_title(f"{gene} ({species}) — {n_leaves} sequences", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])


def build_legend(ax, species):
    """Add a legend for the color scheme."""
    handles = []

    # Outgroup
    handles.append(mpatches.Patch(color=OUTGROUP_COLOR, label="Outgroup"))

    # CBP by thesis label
    if species == "inaquosorum":
        for i in range(1, 15):
            label = f"I{i}"
            handles.append(mpatches.Patch(color=INAQ_COLORS[label], label=f"CBP {label}"))
    elif species == "spizizenii":
        for i in range(1, 6):
            label = f"S{i}"
            handles.append(mpatches.Patch(color=SPIZ_COLORS[label], label=f"CBP {label}"))

    # PCR envs
    for env in ["E1", "E2", "E4", "E5"]:
        handles.append(mpatches.Patch(color=ENV_COLORS[env], label=f"PCR {env}"))

    # Outside LCA
    handles.append(mpatches.Patch(color=OUTSIDE_COLOR, label="PCR outside LCA"))

    # LCA region
    handles.append(mpatches.Patch(facecolor="#f0f0ff", edgecolor="#cccccc",
                                   alpha=0.4, label="CBP LCA clade"))

    ax.legend(handles=handles, loc="upper left", fontsize=4, ncol=2,
              framealpha=0.8, borderpad=0.5)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load gene->species
    gene_species = {}
    with open(PRIMERS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene_species[row["gene"]] = row["species"]

    # Usable genes
    usable = {"acuA", "sorA", "yvqK", "albG", "thiD", "iolB", "acsA_2", "alkH"}

    log_path = os.path.join(OUT_DIR, "lca_summary.tsv")
    log = open(log_path, "w")
    log.write("gene\tspecies\ttotal\tcbp\tpcr\toutgroup\t"
              "pcr_inside_lca\tpcr_outside_lca\tpct_outside\n")

    for gene in sorted(usable):
        species = gene_species.get(gene, "?")
        nwk_path = os.path.join(IN_DIR, gene, f"{gene}.nwk")
        if not os.path.exists(nwk_path):
            print(f"  {gene}: no tree file, skipping")
            continue

        tree = Phylo.read(nwk_path, "newick")
        terminals = tree.get_terminals()

        # Classify leaves
        cbp_leaves = [t for t in terminals if get_leaf_type(t.name) == "cbp"]
        pcr_leaves = [t for t in terminals if get_leaf_type(t.name) == "pcr"]
        og_leaves = [t for t in terminals if get_leaf_type(t.name) == "outgroup"]

        # Find CBP LCA
        if len(cbp_leaves) >= 2:
            lca_clade = tree.common_ancestor(cbp_leaves)
        elif len(cbp_leaves) == 1:
            lca_clade = cbp_leaves[0]
        else:
            lca_clade = None

        # Determine which PCR reads are inside vs outside LCA
        if lca_clade is not None:
            lca_terminal_names = {t.name for t in lca_clade.get_terminals()}
        else:
            lca_terminal_names = set()

        pcr_inside = [t for t in pcr_leaves if t.name in lca_terminal_names]
        pcr_outside = [t for t in pcr_leaves if t.name not in lca_terminal_names]
        outside_set = {t.name for t in pcr_outside}

        pct_outside = 100 * len(pcr_outside) / len(pcr_leaves) if pcr_leaves else 0

        print(f"  {gene:12s} ({species:12s})  "
              f"CBP={len(cbp_leaves):>4}  PCR={len(pcr_leaves):>4}  OG={len(og_leaves)}  "
              f"PCR_inside={len(pcr_inside):>4}  PCR_outside={len(pcr_outside):>4}  "
              f"({pct_outside:.1f}% outside)")

        log.write(f"{gene}\t{species}\t{len(terminals)}\t{len(cbp_leaves)}\t"
                  f"{len(pcr_leaves)}\t{len(og_leaves)}\t"
                  f"{len(pcr_inside)}\t{len(pcr_outside)}\t{pct_outside:.1f}\n")

        # Render tree
        n_leaves = len(terminals)
        fig_height = max(6, n_leaves * 0.04)
        fig, ax = plt.subplots(figsize=(12, fig_height))

        draw_tree(tree, outside_set, lca_clade, gene, species, ax)
        build_legend(ax, species)

        pdf_path = os.path.join(OUT_DIR, f"{gene}_tree.pdf")
        fig.savefig(pdf_path, bbox_inches="tight", dpi=150)
        plt.close(fig)

        # Also save SVG
        svg_path = os.path.join(OUT_DIR, f"{gene}_tree.svg")
        fig2, ax2 = plt.subplots(figsize=(12, fig_height))
        draw_tree(tree, outside_set, lca_clade, gene, species, ax2)
        build_legend(ax2, species)
        fig2.savefig(svg_path, bbox_inches="tight")
        plt.close(fig2)

    log.close()
    print(f"\nLCA summary: {log_path}")
    print(f"Tree PDFs/SVGs: {OUT_DIR}/")


if __name__ == "__main__":
    main()
