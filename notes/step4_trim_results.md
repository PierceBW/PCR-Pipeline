# Step 4: Primer-Guided Trimming Results

## What We Did

After merging paired reads (Step 0), deduplicating (Step 1), filtering by gene-specific identity (Step 2), and adding full-length CBP isolates + outgroup (Step 3), we trimmed every sequence to the primer-to-primer region.

**The problem:** Merged PCR reads include ~6 bp of adapter/reference overhang beyond each primer. CBP isolates are already primer-to-primer. Outgroup sequences were extracted by different methods with variable boundaries. Simply cutting from one end wouldn't work — we needed to find the actual primer positions in each sequence.

**The approach:** For every sequence, search for the forward primer near the 5' end and the reverse primer RC near the 3' end. Cut to those boundaries.

**Script:** `pipeline/step4_trim.py --trim`

## Per-Gene Primer Offset Analysis

The `--explore` mode confirmed a consistent pattern across all genes:

| Sequence Type | Left Offset (fwd primer from start) | Right Offset (rev primer RC from end) |
|--------------|:-----------------------------------:|:-------------------------------------:|
| **CBP** | **0** (every gene, every sequence) | **0** (every gene, every sequence) |
| **PCR** | **6** (99%+ of reads, all genes) | **6** (99%+ of reads, all genes) |
| **Outgroup** | varies by gene (primer-extracted ones: 0) | varies by gene |

CBP isolates are already perfectly primer-bounded. PCR reads consistently have 6 bp overhang on each end from the BAM mapping reference extending past the primer sites.

## Trim Results

| Gene | Species | Input | Trimmed | Dropped | Mode Length | Notes |
|------|---------|------:|--------:|--------:|------------:|-------|
| acuA | inaq | 1,127,449 | 1,123,345 | 4,104 | **522 bp** | Largest dataset |
| iolB | spi | 954,566 | 952,028 | 2,538 | **529 bp** | |
| sorA | inaq | 941,627 | 936,359 | 5,268 | **522 bp** | |
| yvqK | inaq | 445,998 | 443,266 | 2,732 | **508 bp** | |
| albG | inaq | 358,665 | 357,022 | 1,643 | **518 bp** | |
| thiD | inaq | 358,240 | 356,977 | 1,263 | **557 bp** | |
| acsA_2 | spi | 333,251 | 332,335 | 916 | **505 bp** | |
| alkH | spi | 106,638 | 106,320 | 318 | **583 bp** | |
| comB | spi | 292 | 291 | 1 | **506 bp** | Small dataset |
| alaS | spi | 222 | 221 | 1 | **598 bp** | Gap-filled reads |
| opuAB | atr | 156 | 155 | 1 | **505 bp** | |
| amj | atr | 145 | 143 | 2 | **511 bp** | |
| acdA | atr | 126 | 126 | 0 | **522 bp** | |
| ecfA1 | atr | 126 | 126 | 0 | **540 bp** | |
| yndE_2 | atr | 126 | 125 | 1 | **584 bp** | |
| rmlD | spi | 211 | 211 | 0 | **558 bp** | CBP + outgroup only |

**Total: 4,609,094 trimmed sequences** across 16 genes. Drop rate <0.5% per gene — only sequences where the primer couldn't be found.

## Outgroup Extraction

The outgroup (*B. amyloliquefaciens* FN597644.1 or *B. vallismortis* CP026362.1) was re-extracted at full amplicon length using two methods:

1. **Primer BLAST against genome** (10 genes): BLASTed forward and reverse primers against the outgroup genome, extracted the amplicon between hits. Produces exact primer-to-primer sequences.

2. **CBP BLAST against genome** (5 genes: acuA, albG, alkH, comB, thiD): When primers were too divergent to hit, BLASTed a CBP sequence against the genome and extracted the corresponding region.

3. **No reliable extraction** (3 genes: amj, opuAB, yndE_2): Used the old 600 bp extraction; primer-guided trimming cut these to approximately correct boundaries.

| Gene | Outgroup Source | Outgroup Length | Method |
|------|----------------|:-:|--------|
| accA, acdA, acsA_2, alaS, ecfA1, iolB, rmlD, sorA, yvqK | Primer BLAST | Matches CBP | Exact |
| acuA, alkH, comB, thiD | CBP BLAST vs vallismortis | Within 1 bp of CBP | Good |
| albG | CBP BLAST vs amyloliquefaciens | 20 bp shorter than CBP | Partial |
| amj, opuAB, yndE_2 | Old 600 bp, primer-trimmed | ~511-584 bp | Approximate |

## Remaining Length Variation

While the **mode** length per gene matches CBP exactly, some PCR reads have slightly different lengths due to:
- Primer found at offset 5 or 7 instead of 6 (minor variation)
- Unusual CIGARs during merging causing non-standard overlap

For most genes, >99% of sequences are at the mode length. A final uniform-length clip may be needed before tree building if FastTree requires identical lengths.

## Files

- `data/04_trimmed/{gene}.fasta` — Trimmed sequences
- `data/04_trimmed/notes/trim_log.tsv` — Per-gene trim statistics
