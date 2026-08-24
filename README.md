# PCR-Pipeline

Pipeline for classifying environmental *Bacillus* amplicon reads to ecotypes using EcoSim2, integrating PCR environmental reads with culture-based protocol (CBP) isolate sequences from Jocelyn Wang's thesis.

## Background

This pipeline takes paired-end amplicon reads from Death Valley soil samples (sequenced on AVITI), merges, filters, and combines them with CBP isolate sequences that have known ecotype assignments. It then builds phylogenetic trees and runs EcoSim2 to demarcate ecotypes among the combined sequences.

**Species studied:**
- *B. inaquosorum* (14 thesis ecotypes: I1–I14) — 5 usable genes
- *B. spizizenii* (5 thesis ecotypes: S1–S5) — 3 usable genes
- *B. atrophaeus* (8 thesis ecotypes: A1–A8) — insufficient PCR data

**Environmental samples:** 4 soil samples from Death Valley at ~7030 ft elevation. See [notes/environment_metadata_investigation.md](notes/environment_metadata_investigation.md) for full details.

## Pipeline Steps

All scripts are in `pipeline/`. Each step reads from the previous step's output directory under `data/`.

| Step | Script | Input | Output | Description |
|------|--------|-------|--------|-------------|
| 0 | `step0_merge_pairs.py` | `ecosim_ready2/` | `data/00_merged/` | Merge paired-end reads per gene per environment |
| 1 | `step1_dedup.py` | `00_merged/` | `01_deduped/` | Deduplicate merged reads |
| 2 | `step2_gene_filter.py` | `01_deduped/` | `02_filtered/` | Gene-specific identity filter (95% to amplicon ref via BLAST) |
| 3 | `step3_add_cbp.py` | `02_filtered/` | `03_with_cbp/` | Add CBP isolate sequences with thesis ecotype labels |
| 4 | `step4_trim.py` | `03_with_cbp/` | `04_trimmed/` | Trim all sequences to uniform length per gene |
| 4b | `step4b_dedup.py` | `04_trimmed/` | `04b_deduped/` | Remove exact duplicates (CBP and PCR independently) |
| 5 | `step5_ecosim.py` | `04b_deduped/` | `05_ecosim/` | Subsample PCR reads, build FastTree, run EcoSim2 |
| 6 | `step6_visualize_trees.py` | `05_ecosim/200/` | `06_viz/` | Tree visualization with CBP/PCR coloring |
| 7 | `step7_prune_ingroup.py` | `05_ecosim/200/` | `07_ingroup/` | Find CBP LCA, keep only ingroup sequences |
| 8 | `step8_ecosim_ingroup.py` | `07_ingroup/` | `08_ecosim_ingroup/` | Re-run EcoSim on pruned ingroup |

**Helper scripts:**
- `analyze_ecotypes.py` — Parse EcoSim XMLs and compute ecotype summary stats
- `analyze_singletons.py` — Analyze singleton sequences across ecotypes
- `clip_to_uniform.py` — Clip FASTAs to uniform mode length
- `make_report_pdf.py` — Generate combined PDF report with trees

### Running the pipeline

```bash
# Set up environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Also need: fasttree, blast+, java (for EcoSim JAR)
# macOS: brew install fasttree blast

# Run steps in order
python pipeline/step0_merge_pairs.py
python pipeline/step1_dedup.py
python pipeline/step2_gene_filter.py
python pipeline/step3_add_cbp.py
python pipeline/step4_trim.py --explore   # preview trim lengths
python pipeline/step4_trim.py --trim      # apply trim
python pipeline/step4b_dedup.py
python pipeline/step5_ecosim.py --n-pcr 200
python pipeline/step6_visualize_trees.py
python pipeline/step7_prune_ingroup.py
python pipeline/step8_ecosim_ingroup.py
```

EcoSim2 JAR is expected at `../ecosim/ecosim.jar` (one level up from this repo).

## Key Results

**PCR reads split cleanly by species:**
- **Inaquosorum**: All PCR reads fall inside the CBP LCA for all 5 genes (same species as thesis isolates)
- **Spizizenii**: All PCR reads fall outside the CBP LCA for all 3 genes (different taxon)

This is a binary, all-or-nothing result across all 8 genes.

See [notes/inaquosorum_ecotype_summary.md](notes/inaquosorum_ecotype_summary.md) for per-gene ecotype breakdown and comparison to thesis results.

## Data (on lab computer only)

The following directories contain large data files and are **not in this repo**. They live on the lab computer at the same path alongside this repo:

```
PCR-Pipeline/
├── data/                  # All intermediate and final pipeline outputs (25 GB)
│   ├── 00_merged/         # Step 0 output
│   ├── 01_deduped/        # Step 1 output
│   ├── 02_filtered/       # Step 2 output
│   ├── 03_with_cbp/       # Step 3 output
│   ├── 04_trimmed/        # Step 4 output
│   ├── 04b_deduped/       # Step 4b output
│   ├── 05_ecosim/         # Step 5 output (trees, EcoSim XMLs)
│   ├── 06_viz/            # Step 6 output (tree PDFs)
│   ├── 07_ingroup/        # Step 7 output (pruned trees)
│   ├── 08_ecosim_ingroup/ # Step 8 output (ingroup EcoSim)
│   └── old_pre_dedup/     # Backup of results before step 4b was added
├── source_data/           # BAM files, thesis data, primer files (14 GB)
├── genomes/               # 487 CBP isolate genome assemblies (1.9 GB)
├── reference/             # BLAST databases, outgroup genome (44 MB)
├── archive/               # Old pipeline scripts (superseded)
└── old/                   # Old pipeline outputs (superseded, 26 GB)
```

## Dependencies

- Python 3.13+ with packages in `requirements.txt`
- [FastTree](http://www.microbesonline.org/fasttree/) (`brew install fasttree`)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`brew install blast`)
- Java (for EcoSim2 JAR)

## References

- Wang, J. et al. — Jocelyn Wang's honors thesis on Bacillus ecotype demarcation (source of CBP isolate ecotype assignments)
- Koeppel et al. 2008; Francisco et al. 2014; Wood et al. 2021 — EcoSim/ES2
