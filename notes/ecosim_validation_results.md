# EcoSim Validation Experiment: Single-Locus vs Multi-Locus Ecotype Recovery

## Problem

EcoSim2 was producing vastly inflated ecotype counts when run on our amplicon sequences. Across all three species (*B. atrophaeus*, *B. inaquosorum*, *B. spizizenii*), we observed 10–100× more ecotypes than expected based on Jocelyn Wang's thesis, which identified 8, 14, and 5 ecotypes respectively using multi-locus concatenated analysis.

We needed to determine whether this was caused by:
1. Contamination or misclassified environmental sequences in the dataset
2. A fundamental incompatibility between single-locus amplicon trees and multi-locus ecotype assignments

## What We Did

### Step 1: Identity Filtering (95% Cutoff)

We applied a stricter identity filter to remove potentially misclassified sequences:

- BLASTed all PCR reads in the forward tree-input files against the 6-genome reference database (2 genomes per species)
- Kept only reads with ≥95% identity to their assigned species' reference genomes (up from the original 93% cutoff)
- CBP isolates and the outgroup were retained unconditionally
- Applied to all 3 species, forward strand only

**Result:** The filter removed very few sequences — 99%+ of reads passed for most genes. The biggest drop was *B. inaquosorum* comB at 18% removed; most genes lost <1%.

### Step 2: Small-Scale EcoSim Test (200 Sequences)

We subsampled filtered files to ~200 PCR reads (plus all CBP isolates and the outgroup), built FastTree phylogenies, and ran EcoSim2.

**Result:** Ecotype counts remained inflated across all species and genes (54–248 ecotypes from 200 sequences). The identity filter did not resolve the problem.

### Step 3: CBP-Only EcoSim Experiment (The Key Test)

To isolate whether the problem was the environmental PCR reads or EcoSim's behavior on single-locus trees, we ran EcoSim on **only the CBP isolates** — the same curated lab strains whose ecotype assignments are known from the thesis.

For each gene with ≥10 CBP isolates:
1. Extracted outgroup + CBP sequences only (removed all PCR reads)
2. Built a FastTree phylogeny and rerooted at the outgroup
3. Ran EcoSim2
4. Compared EcoSim's ecotype clusters to the known thesis ecotypes
5. Computed pairwise phylogenetic distances within and between thesis ecotypes

## Results

### EcoSim Ecotype Counts vs Thesis Expectations (CBP Isolates Only)

The table below includes the thesis single-gene ES2 ecotype count (from Wang Tables S2/S3), which is the number of ecotypes that gene produced when ES2 was run on the thesis's whole-genome isolates using just that one gene's alignment. The "Primer Class" column shows whether the gene was selected as a high-ecotype or modal-ecotype gene for our amplicon primer design (from `source_data/bacillus-primer/REFERENCE_GUIDE.md`).

| Species | Gene | Primer Class | Thesis Single-Gene ES2 | CBP Isolates | Thesis MLSA Ecotypes | Our CBP EcoSim | Ratio |
|---------|------|:---:|:---:|:---:|:---:|:---:|:---:|
| *B. spizizenii* | iolB | high | 18 | 211 | 5 | **5** | **1.0×** |
| *B. inaquosorum* | yvqK | high | 39 | 141 | 14 | **8** | 0.6× |
| *B. inaquosorum* | sorA | high | 41 | 141 | 14 | **17** | 1.2× |
| *B. atrophaeus* | opuAB | high | 17 | 123 | 8 | **12** | 1.5× |
| *B. inaquosorum* | thiD | high | 41 | 142 | 14 | 37 | 2.6× |
| *B. atrophaeus* | ecfA1 | modal* | 20 | 122 | 8 | 22 | 2.8× |
| *B. inaquosorum* | albG | — | — | 141 | 14 | 77 | 5.5× |
| *B. inaquosorum* | accA | modal | 14 | 142 | 14 | 95 | 6.8× |
| *B. inaquosorum* | acuA | modal | 14 | 141 | 14 | 95 | 6.8× |
| *B. atrophaeus* | acdA | modal | 8 | 123 | 8 | 64 | 8.0× |
| *B. atrophaeus* | amj | high* | 8 | 123 | 8 | 74 | 9.2× |
| *B. atrophaeus* | yndE_2 | high | 19 | 123 | 8 | 84 | 10.5× |
| *B. spizizenii* | comB | high | 13 | 210 | 5 | 67 | 13.4× |
| *B. spizizenii* | alaS | modal | 5 | 211 | 5 | 1 | 0.2× |
| *B. spizizenii* | alkH | modal | 5 | 209 | 5 | 1 | 0.2× |
| *B. spizizenii* | rmlD | high | 16 | 210 | 5 | 1 | 0.2× |

**\*** Primer classification discrepancies: **ecfA1** was classified as "modal" in our primer design but produced 20 ecotypes in the thesis (Table S2, high). **amj** was classified as "high" in our primer design but produced 8 ecotypes in the thesis (Table S3, modal). These misclassifications do not affect the pipeline but are worth noting.

### Thesis Ecotype Fragmentation Across EcoSim Clusters

The detailed comparison showed that thesis ecotypes are not monophyletic on single-gene trees. For example:

- *B. atrophaeus* **acdA**: Thesis ecotype PE_A1 (48 strains) was fragmented across **31 different** EcoSim clusters
- *B. inaquosorum* **accA**: Thesis ecotype PE_I1 (31 strains) was fragmented across **30 different** EcoSim clusters
- *B. atrophaeus* **yndE_2**: PE_A1 (48 strains) → 38 clusters; PE_A2 (36 strains) → 31 clusters

In contrast, for genes where EcoSim matched the thesis:
- *B. spizizenii* **iolB**: All 5 thesis ecotypes mapped to distinct EcoSim clusters with minimal mixing
- *B. inaquosorum* **yvqK**: Thesis ecotypes mapped to 1–4 EcoSim clusters each

### Phylogenetic Distance Analysis

Within-ecotype pairwise distances often overlapped with between-ecotype distances, indicating that single-locus trees do not cleanly separate the thesis ecotypes:

- *B. inaquosorum* accA, PE_I3: within-ecotype max distance = **0.948** (extremely high — nearly 95% divergent within a single "ecotype")
- *B. atrophaeus* acdA: within-ecotype mean distances (0.001–0.017) overlapped broadly with between-ecotype mean distance (0.016)
- *B. spizizenii* alaS: within-ecotype distances were so small (max 0.003) that EcoSim collapsed all 5 ecotypes into 1

For genes where EcoSim succeeded, within-ecotype distances were consistently lower than between-ecotype distances, providing the phylogenetic signal EcoSim needs.

## Context from Wang Thesis: Single-Gene ES2 Results

Wang's thesis (Sections 3.4, 3.8, 4.3, 4.4) directly addressed how ES2 behaves on single-gene inputs using the same CBP isolates, providing critical context for our results.

### Rarefaction Analysis (Section 3.4)

Wang ran ES2 on concatenated samples of varying numbers of genes (1, 3, 7, 20, 100, 200, 400, and 700), with 500 replicate samples at each gene number. Key findings:

- **B. atrophaeus**: Modal ecotype count was **8** at all gene numbers, even single-gene. However, with 1 gene, ~100 replicates yielded 1 ecotype and ~100 yielded 7, showing high variance.
- **B. inaquosorum**: Modal count was **14**, stabilizing quickly. With 7 genes, a bimodal distribution appeared (~100 replicates at 14, ~100 at 18–26).
- **B. spizizenii**: Modal count was **5**, but nearly 200 single-gene samples yielded **9 ecotypes** — a persistent secondary mode.

The key conclusion: ES2's ecotype count was "not affected by the level of molecular resolution" (i.e., number of genes), unlike MED (Minimum Entropy Decomposition) which required ~20 genes to stabilize. This was interpreted as evidence that "only a limited number of lineages form sequence clusters at any molecular resolution level."

### Gene Divergence Analysis (Section 3.8)

Wang tested whether genes with higher nucleotide divergence produce more ecotypes. Using Spearman's Rank correlation:

- **B. atrophaeus**: r = 0.039 (P = 0.30) — no correlation
- **B. inaquosorum**: r = -0.008 (P = 0.87) — no correlation
- **B. spizizenii**: r = -0.388 (P < 2.2e-16) — weak negative correlation (unexplained)

**Conclusion**: "Number of ecotypes demarcated by ES2 was not affected by the average pairwise nucleotide differences." Wang could not explain why certain genes consistently produced high ecotype counts.

### Thesis Single-Gene Ecotype Counts (Tables S2 & S3)

**Table S2 — Genes producing high ecotype counts** (well above modal):

| Species | Gene | Single-Gene ES2 Ecotypes | Modal |
|---------|------|:---:|:---:|
| *B. atrophaeus* | ecfA1 | 20 | 8 |
| *B. atrophaeus* | yndE_2 | 19 | 8 |
| *B. atrophaeus* | opuAB | 17 | 8 |
| *B. atrophaeus* | dapH_2 | 14 | 8 |
| *B. inaquosorum* | sorA | 41 | 14 |
| *B. inaquosorum* | thiD | 41 | 14 |
| *B. inaquosorum* | yvqK | 39 | 14 |
| *B. inaquosorum* | albE | 32 | 14 |
| *B. inaquosorum* | dapH_2 | 29 | 14 |
| *B. spizizenii* | gyra | 67 | 5 |
| *B. spizizenii* | iolB | 18 | 5 |
| *B. spizizenii* | rmlD | 16 | 5 |
| *B. spizizenii* | comB | 13 | 5 |

**Table S3 — Genes producing the modal ecotype count** (top 10 per species):

| Species | Modal | Example Genes |
|---------|:---:|---------------|
| *B. atrophaeus* | 8 | acdA, adc, alkH, alsS, amiF, **amj**, ansA, argB, argC |
| *B. inaquosorum* | 14 | arnC_1, **accA**, acoR, **acuA**, adhA, aes, albG, aldH1, amaA, amiC |
| *B. spizizenii* | 5 | acoR, acsA_2, **alaS**, **alkH**, amiF, amyE, ansA, aprX, araQ_2, araR |

Genes in **bold** are ones we used in our amplicon pipeline.

### Comparison: Thesis Single-Gene ES2 vs Our CBP-Only EcoSim

A striking discrepancy: in Wang's thesis, most single genes produced the **modal** ecotype count (8, 14, or 5) when run through ES2, even as individual genes. But our CBP-only EcoSim runs on the same isolates produced wildly different numbers for many genes.

The likely explanation is **methodological**: Wang extracted single-gene alignments from **whole-genome assemblies** (high-quality, full-length gene sequences aligned across all isolates). Our pipeline uses **short amplicon sequences** (~300–600 bp PCR products from a single primer pair), which capture only a portion of each gene. The reduced sequence length means less phylogenetic information per gene, and the tree topology may differ from what Wang obtained with full-length gene sequences.

Additionally, Wang used the core-genome alignment produced by Roary, where all genes were aligned across all genomes simultaneously. Our amplicon trees are built independently from PCR products, which may have different error profiles and alignment characteristics.

## Interpretation

The thesis ecotypes were defined using **concatenated multi-locus sequence analysis** (MLSA), which integrates phylogenetic signal across many genes simultaneously. A single housekeeping gene captures only a fraction of that signal. As a result:

1. **Some genes have too little variation** to resolve ecotypes (alaS, alkH, rmlD in *B. spizizenii* — all sequences nearly identical, EcoSim sees 1 ecotype)

2. **Some genes have discordant genealogies** due to recombination or incomplete lineage sorting, causing thesis ecotypes to appear polyphyletic on single-gene trees (acdA, amj, yndE_2 in *B. atrophaeus* — EcoSim splits each thesis ecotype into many clusters)

3. **A few genes happen to track the multi-locus ecotypes well** (iolB in *B. spizizenii*, yvqK and sorA in *B. inaquosorum*, opuAB in *B. atrophaeus*), likely because they have strong phylogenetic signal that aligns with the species-level ecotype structure

This is a well-known phenomenon in microbial population genetics: individual gene trees can differ substantially from the species tree due to horizontal gene transfer, recombination, and stochastic lineage sorting.

## Conclusion

The inflated ecotype counts in our amplicon pipeline are **not a data quality issue** — they persist even when running EcoSim on the known CBP isolates alone. The root cause is that Ecotype Simulation on single-locus amplicon phylogenies cannot reliably reproduce ecotype assignments that were originally derived from multi-locus analysis.

Importantly, Wang's thesis showed that ES2 on **full-length** single-gene alignments from whole genomes typically recovered the modal ecotype count. Our amplicon-based approach, using shorter PCR fragments, produces far more variable and often inflated results — suggesting that the reduced phylogenetic information in short amplicons exacerbates gene-tree discordance effects that are already present but more manageable with full-length gene sequences.

### Primer Classification Notes

Two genes were misclassified in our primer design reference guide relative to the thesis:
- **ecfA1** (*B. atrophaeus*): Listed as "modal" but produced 20 ecotypes in Wang's single-gene analysis (high)
- **amj** (*B. atrophaeus*): Listed as "high" but produced 8 ecotypes in Wang's single-gene analysis (modal)

These misclassifications do not affect the amplicon pipeline's operation but should be corrected in future reference documentation.

## Files

All results are in `data/12_cbp_only_test/`:
- `{species}/{gene}_forward.fasta` — CBP-only input FASTAs
- `{species}/{gene}_forward.nwk` — Rooted phylogenies
- `{species}/{gene}_forward.xml` — EcoSim2 output
- `notes/cbp_ecosim_summary.tsv` — Summary table (species, gene, CBP count, thesis vs EcoSim ecotypes)
- `notes/ecotype_comparison.tsv` — Per-thesis-ecotype fragmentation across EcoSim clusters
- `notes/distance_analysis.tsv` — Within- and between-ecotype pairwise distances per gene

Scripts used:
- `filter_by_type_strain.py` — 95% identity filter
- `test_ecosim_small.py` — Small-scale EcoSim test runner
- `cbp_ecosim_test.py` — CBP-only EcoSim experiment with comparison analysis
