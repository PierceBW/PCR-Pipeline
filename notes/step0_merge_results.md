# Step 0: Paired-Read Merging Results

## What We Did

The original BAM files contain paired-end Illumina reads (R1 + R2, each 301 bp) that together cover the full amplicon (~500-600 bp) with overlap in the middle. We extracted paired reads from all 72 gene-specific BAM files (18 genes x 4 environments) and merged R1+R2 into single full-length amplicon sequences.

**Script:** `pipeline/step0_merge_pairs.py`
**Input:** `source_data/PCR-Primer/datasets/*.bam` (72 gene-level BAMs)
**Output:** `data/00_merged/{gene}/{env}.fasta`

## How Merging Works

1. Parse BAM with `pysam`, group reads by `query_name` to find R1/R2 mates
2. Verify orientation: R1 maps to forward strand, R2 maps to reverse strand (98.7% of pairs)
3. Compute overlap from reference coordinates: `overlap = R1.reference_end - R2.reference_start`
4. For **positive overlap** (~61-85 bp): stitch R1 + consensus_overlap + R2, picking the higher-quality base where reads disagree
5. For **negative overlap** (gap between reads): fill with dashes (`---`) like an alignment gap. This occurs when the amplicon is longer than 602 bp (301 + 301)
6. Reject pairs with >10% mismatch in the overlap region

## Results by Gene (All Environments Combined)

| Gene | Species | Total Pairs | Merged (overlap) | Gap-Filled | Unmerged | Merge Rate | Mean Length |
|------|---------|------------:|------------------:|-----------:|---------:|-----------:|------------:|
| acuA | inaq | 1,157,977 | 1,134,029 | 1,363 | 22,585 | 98.0% | 534 bp |
| iolB | spi | 980,619 | 959,121 | 3,904 | 17,594 | 98.2% | 541 bp |
| sorA | inaq | 957,549 | 944,592 | 1,004 | 11,953 | 98.8% | 534 bp |
| rmlD | spi | 580,485 | 567,964 | 3,839 | 8,682 | 98.5% | 571 bp |
| yvqK | inaq | 456,719 | 446,974 | 838 | 8,907 | 98.0% | 520 bp |
| albG | inaq | 364,370 | 359,603 | 337 | 4,430 | 98.8% | 530 bp |
| thiD | inaq | 367,818 | 358,202 | 1,289 | 8,327 | 97.7% | 569 bp |
| acsA_2 | spi | 350,518 | 341,239 | 749 | 8,530 | 97.6% | 517 bp |
| alkH | spi | 177,536 | 162,703 | 4,373 | 10,460 | 94.1% | 596 bp |
| **alaS** | **spi** | **2,121** | **0** | **2,059** | **62** | **97.1%** | **611 bp** |
| comB | spi | 453 | 415 | 3 | 35 | 92.3% | 519 bp |
| amj | atr | 25 | 21 | 0 | 4 | 84.0% | 523 bp |
| opuAB | atr | 62 | 32 | 22 | 8 | 87.1% | 591 bp |
| ecfA1 | atr | 71 | 3 | 0 | 68 | 4.2% | 552 bp |
| acdA | atr | 5 | 2 | 0 | 3 | 40.0% | 534 bp |
| yndE_2 | atr | 2 | 2 | 0 | 0 | 100.0% | 596 bp |
| accA | inaq | 32 | 0 | 0 | 32 | 0.0% | -- |
| adc | atr | 0 | 0 | 0 | 0 | -- | -- |

**Total merged sequences: ~5.3 million across all genes**

## Key Findings

### 1. High merge rates for genes with real data
All genes with substantial read counts merge at **94-99%**. The 1-6% that fail are due to wrong orientation (adapter dimers/chimeras) or high mismatch in the overlap region.

### 2. alaS has a gap, not overlap
The alaS amplicon is ~610 bp, meaning 301 + 301 = 602 bp of read coverage leaves an ~8 bp gap in the middle. All 2,059 alaS merged sequences use dash-filling (`--------`) in that gap. This is consistent across all environments and all reads. The dashes will be treated as missing data by FastTree, which is standard practice for alignment gaps.

### 3. alkH is borderline
alkH has ~0 bp overlap on average, which means some pairs barely overlap and some have small gaps. This explains its slightly lower merge rate (94.1%) and the 4,373 gap-filled sequences. The merged sequences are still valid.

### 4. Atrophaeus genes have very few environmental reads
ecfA1 (3 merged), acdA (2), yndE_2 (2), amj (21), opuAB (32) — atrophaeus primers rarely amplify from these environmental samples. accA (inaquosorum) also has 0 merged reads (32 pairs all had wrong orientation). These genes will have very few PCR reads in downstream steps.

### 5. Merged amplicon lengths match expectations
Mean lengths per gene (517-611 bp) are consistent with the expected primer-to-primer amplicon sizes from the reference genomes (505-598 bp per `notes/primers.tsv`).

## Unmerged Read Breakdown

Reads fail to merge for these reasons:
- **Wrong orientation** (~1-2%): Both reads map to the same strand — adapter dimers or chimeric fragments
- **High mismatch** (~0.5-1%): >10% disagreement in overlap region — sequencing errors or chimeras
- **Non-standard length** (<0.1%): Reads shorter than 301 bp — partial sequences
- **Missing mate** (<0.01%): Only one read of the pair present in the BAM

## Files

- `data/00_merged/{gene}/{env}.fasta` — Merged amplicon sequences
- `data/00_merged/{gene}/{env}_unmerged.fasta` — Failed-to-merge reads (kept for reference)
- `data/00_merged/notes/merge_log.tsv` — Per-gene-per-environment statistics
