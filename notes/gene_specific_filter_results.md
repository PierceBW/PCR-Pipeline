# Gene-Specific Identity Filter Results

## Motivation

Our previous identity filter (Stage 9, `filter_by_type_strain.py`) BLASTed PCR reads against **whole reference genomes** (~4 Mb each). This meant a read could pass at ≥95% identity by matching any conserved region of the genome — not necessarily the target gene. Since our three *Bacillus* species share substantial genome-wide conservation, this approach could not detect:

- Off-target amplification products (primer mispriming to a different gene)
- Reads classified to the wrong species that happen to match a conserved locus
- Cross-species amplification where primers designed for one species amplify a homolog in another

## Approach

We built **gene-specific reference sequences** by extracting the actual amplicon region from the NCBI reference genomes:

1. For each of the 18 gene/species primer pairs (from `notes/primers.tsv`), BLAST the forward and reverse primers against the 2 reference genomes for that species
2. Extract the amplicon region (between primer hits) — typically ~500–600 bp
3. Build a per-gene BLAST database from these extracted amplicons
4. BLAST each PCR read against the gene-specific reference (not the whole genome)
5. Keep reads with ≥95% identity to the gene-specific reference; CBP isolates and outgroup retained unconditionally

**Script:** `filter_by_gene_ref.py`
**Output:** `data/13_gene_filtered/`

## Results

### Atrophaeus (very few PCR reads — all pass)

| Gene | PCR Reads | Kept | Dropped | % Kept | Mean pident |
|------|:---------:|:----:|:-------:|:------:|:-----------:|
| acdA | 2 | 2 | 0 | 100% | 98.3 |
| amj | 21 | 21 | 0 | 100% | 97.8 |
| ecfA1 | 3 | 3 | 0 | 100% | 99.9 |
| opuAB | 32 | 32 | 0 | 100% | 98.1 |
| yndE_2 | 2 | 2 | 0 | 100% | 98.6 |

All atrophaeus reads are on-target. Low PCR counts (2–32) are expected — these genes have low amplification from environmental samples.

### Inaquosorum (large read counts — reveals cross-species contamination)

| Gene | Primer Species | PCR Reads | Kept | Dropped | No Hit | % Kept | Notes |
|------|:--------------:|:---------:|:----:|:-------:|:------:|:------:|-------|
| accA | inaq | 0 | — | — | — | — | CBP only |
| **acsA_2** | **spi** | 145,498 | 7,608 | 137,807 | 83 | **5.2%** | **94.8% off-target** |
| acuA | inaq | 327,867 | 327,069 | 798 | 0 | 99.8% | On-target |
| alaS | spi | 1,897 | 1,855 | 42 | 0 | 97.8% | Cross-species but similar |
| albG | inaq | 93,521 | 93,231 | 290 | 0 | 99.7% | On-target |
| **alkH** | **spi** | 68,084 | 53,449 | 14,631 | 4 | **78.5%** | **21.5% off-target** |
| **comB** | **spi** | 414 | 72 | 341 | 1 | **17.4%** | **82.6% off-target** |
| iolB | spi | 295,159 | 287,792 | 7,364 | 3 | 97.5% | Mostly on-target |
| **rmlD** | **spi** | 196,893 | 0 | 0 | 196,893 | **0.0%** | **Reference extraction issue** |
| sorA | inaq | 301,658 | 299,871 | 1,787 | 0 | 99.4% | On-target |
| thiD | inaq | 167,643 | 166,503 | 1,140 | 0 | 99.3% | On-target |
| yvqK | inaq | 132,437 | 131,709 | 728 | 0 | 99.5% | On-target |

### Spizizenii (mostly CBP-dominated)

| Gene | Primer Species | PCR Reads | Kept | Dropped | No Hit | % Kept | Notes |
|------|:--------------:|:---------:|:----:|:-------:|:------:|:------:|-------|
| **acuA** | **inaq** | 134 | 37 | 95 | 2 | **27.6%** | **72.4% off-target** |
| alaS | spi | 11 | 9 | 0 | 2 | 81.8% | Low count |
| alkH | spi | 12,067 | 11,976 | 90 | 1 | 99.2% | On-target |
| comB | spi | 10 | 9 | 1 | 0 | 90.0% | Low count |
| iolB | spi | 580 | 577 | 3 | 0 | 99.5% | On-target |
| **rmlD** | **spi** | 251 | 0 | 0 | 251 | **0.0%** | **Reference extraction issue** |

## Key Findings

### 1. Cross-Species Gene Contamination is Real

Genes with primers designed for one species often appear in another species' directory because the original classification (Stage 4) was genome-wide, not gene-specific:

- **acsA_2** (spi primers) in inaquosorum: 95% of 145K reads are off-target
- **comB** (spi primers) in inaquosorum: 83% off-target
- **alkH** (spi primers) in inaquosorum: 21% off-target
- **acuA** (inaq primers) in spizizenii: 72% off-target

These reads matched the correct species' genome at ≥93% (passing Stage 4 classification) but are NOT from the target gene — they're from a different region of the genome that the cross-species primers happened to amplify.

### 2. Same-Species Genes are Clean

When the primers match the species (e.g., acuA/albG/sorA/thiD/yvqK in inaquosorum, alkH/iolB in spizizenii), **>97% of reads pass** the gene-specific filter. This confirms the PCR amplification is on-target for same-species genes.

### 3. rmlD Reference Extraction Problem

The rmlD amplicon extraction produced anomalously large sequences (889 bp and 1102 bp vs expected ~500 bp), suggesting the primer BLAST found hits at incorrect locations in the reference genomes. All 197K+ rmlD reads got "no hit" — they can't align well to an 889 bp reference when they're only 301 bp. This needs a manual fix to the reference extraction, not a data quality issue.

### 4. Impact on Previous Analysis

The whole-genome 95% filter (`filter_by_type_strain.py`) removed very few reads because it was checking against the whole genome, not the target gene. This gene-specific filter reveals that a substantial fraction of cross-species gene assignments contain off-target reads. For downstream EcoSim analysis, using gene-specific filtered data should produce more accurate ecotype assignments.

## Comparison: Whole-Genome vs Gene-Specific Filter

| Gene (inaq) | Whole-Genome Filter % Kept | Gene-Specific Filter % Kept | Difference |
|-------------|:--------------------------:|:---------------------------:|:----------:|
| acsA_2 | ~99% | 5.2% | **94% were off-target** |
| comB | ~82% | 17.4% | **65% were off-target** |
| alkH | ~99% | 78.5% | **21% were off-target** |
| acuA | ~99% | 99.8% | Minimal difference |
| sorA | ~99% | 99.4% | Minimal difference |
| iolB | ~99% | 97.5% | Small improvement |

## Files

- `data/13_gene_filtered/gene_refs/` — Extracted amplicon reference sequences per gene
- `data/13_gene_filtered/{species}/{gene}_forward.fasta` — Filtered output FASTAs
- `data/13_gene_filtered/notes/filter_log.tsv` — Per-gene filtering statistics
- `filter_by_gene_ref.py` — Gene-specific filter script

## Known Issues

1. **rmlD**: Reference amplicon extraction found primer hits too far apart. Needs manual investigation or a different extraction method for this gene.
2. **comB**: One reference genome (CP030925.1) extracted a 1365 bp amplicon vs 506 bp for the other. Multi-hit primer issue in that genome.
