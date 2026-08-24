# Step 2: Gene-Specific Identity Filter Results

## What We Did

After merging paired reads (Step 0) and deduplicating (Step 1), we BLASTed every read against its gene-specific reference amplicon and kept only reads with >=95% identity. Each gene was checked against the reference amplicon for its **designed species only** (from `notes/primers.tsv`).

This is the same gene-specific filtering approach we validated earlier (`filter_by_gene_ref.py`), but now applied to the merged full-length amplicons instead of single 301 bp forward reads.

**Script:** `pipeline/step2_gene_filter.py`
**Input:** `data/01_deduped/{gene}.fasta` (merged, deduplicated reads)
**Output:** `data/02_filtered/{gene}.fasta`
**BLAST DBs:** Pre-built gene-specific databases in `data/13_gene_filtered/gene_refs/`

## Results

### Inaquosorum Genes — Nearly All Pass

| Gene | Total | Kept | Dropped | No Hit | % Kept | Mean pident |
|------|------:|-----:|--------:|-------:|-------:|------------:|
| acuA | 1,127,528 | 1,127,307 | 218 | 3 | **100.0%** | 99.7 |
| sorA | 941,919 | 941,485 | 404 | 30 | **100.0%** | 99.1 |
| yvqK | 446,429 | 445,856 | 450 | 123 | **99.9%** | 99.2 |
| thiD | 359,016 | 358,097 | 919 | 0 | **99.7%** | 98.1 |
| albG | 358,875 | 358,523 | 350 | 2 | **99.9%** | 97.7 |
| accA | 0 | — | — | — | — | — |

All five inaquosorum genes with PCR reads pass at >99.7%. These reads are on-target — the primers amplified the correct gene from the correct species. The handful of dropped reads (0.03-0.3%) are low-quality or chimeric sequences.

### Spizizenii Genes — Mixed Results

| Gene | Total | Kept | Dropped | No Hit | % Kept | Mean pident | Interpretation |
|------|------:|-----:|--------:|-------:|-------:|------------:|----------------|
| iolB | 958,152 | 954,354 | 3,790 | 8 | **99.6%** | 97.3 | Clean — on-target |
| acsA_2 | 341,453 | 333,041 | 8,404 | 8 | **97.5%** | 95.8 | Mostly clean, 2.5% off-target |
| alkH | 166,916 | 106,428 | 60,474 | 14 | **63.8%** | 95.8 | **36% off-target** |
| comB | 418 | 81 | 337 | 0 | **19.4%** | 97.9 | **81% off-target** |
| alaS | 2,059 | 10 | 2,049 | 0 | **0.5%** | 96.6 | **99.5% off-target** |
| rmlD | 570,686 | 0 | 0 | 570,686 | **0.0%** | — | **Reference issue** (see below) |

Three clear categories:
1. **iolB and acsA_2**: Mostly clean (97-100% pass)
2. **alkH, comB, alaS**: Major cross-species contamination (36-99.5% off-target)
3. **rmlD**: No reads match at all — not a data quality issue, but a reference extraction problem

### Atrophaeus Genes — Very Few PCR Reads

| Gene | Total | Kept | Dropped | No Hit | % Kept | Mean pident |
|------|------:|-----:|--------:|-------:|-------:|------------:|
| amj | 21 | 21 | 0 | 0 | **100%** | 98.0 |
| opuAB | 54 | 32 | 13 | 9 | **59%** | 98.7 |
| ecfA1 | 3 | 3 | 0 | 0 | **100%** | 99.4 |
| acdA | 2 | 2 | 0 | 0 | **100%** | 98.5 |
| yndE_2 | 2 | 2 | 0 | 0 | **100%** | 98.6 |
| adc | 0 | — | — | — | — | — |
| accA (inaq) | 0 | — | — | — | — | — |

Atrophaeus genes have extremely few environmental PCR reads (0-54 total). This is expected — atrophaeus has low environmental abundance in these samples. The reads that do exist are high quality (>98% identity).

## Understanding the Off-Target Reads

### Why alkH, comB, and alaS fail

These spizizenii genes have high drop rates because most of their reads are actually **inaquosorum reads that got classified into the wrong species directory** in the original BAM mapping. The BAM reference includes the full gene region (~810 bp), and reads from a closely related inaquosorum homolog can map to the spizizenii reference with decent alignment scores.

When we BLAST those reads against the gene-specific spizizenii amplicon reference (~500-600 bp extracted between primers), they fall below 95% because the polymorphisms between species are concentrated in the gene region.

| Gene | Reads from BAM | Pass as spizizenii | Actually inaquosorum? |
|------|---------------:|-------------------:|:---------------------:|
| alkH | 166,916 | 106,428 (64%) | ~60,000 (36%) |
| comB | 418 | 81 (19%) | ~337 (81%) |
| alaS | 2,059 | 10 (0.5%) | ~2,049 (99.5%) |

This is exactly the filtering working as intended — removing reads that are not from the target gene/species combination.

### The rmlD Problem

All 570,686 rmlD reads get "no hit" against the reference. This is because the rmlD reference amplicon was extracted incorrectly — the primer BLAST found hits too far apart in the reference genomes, producing an 889-1102 bp "amplicon" instead of the expected ~500 bp. A read of ~570 bp can't align well to an 889 bp reference.

**This is a reference extraction issue, not a data quality issue.** The rmlD reads themselves are likely fine. Fixing the rmlD reference extraction (manually locating the correct amplicon region) would recover these reads. For now, rmlD is excluded from downstream analysis.

## Impact on Dataset Size

| Category | Genes | Total PCR Reads | After Filter | Notes |
|----------|:-----:|----------------:|-------------:|-------|
| Clean inaquosorum | acuA, sorA, yvqK, thiD, albG | 3,233,767 | 3,231,268 (99.9%) | On-target, large datasets |
| Clean spizizenii | iolB, acsA_2 | 1,299,605 | 1,287,395 (99.1%) | On-target, large datasets |
| Mixed spizizenii | alkH | 166,916 | 106,428 (63.8%) | Significant cleanup |
| Mostly off-target | comB, alaS | 2,477 | 91 (3.7%) | Small datasets after filter |
| Reference issue | rmlD | 570,686 | 0 (0%) | Fixable — needs reference correction |
| Low-count atrophaeus | amj, opuAB, ecfA1, acdA, yndE_2 | 82 | 60 (73.2%) | Too few for EcoSim alone |
| No data | accA, adc | 0 | 0 | No PCR amplification |

**Usable genes for EcoSim** (>1000 filtered PCR reads):
- **inaquosorum:** acuA (1.13M), sorA (941K), yvqK (446K), albG (359K), thiD (358K)
- **spizizenii:** iolB (954K), acsA_2 (333K), alkH (106K)
- **Potentially:** rmlD (571K if reference fixed)

## Decision Points for Next Steps

### 1. Sequence Length Mismatch
Merged PCR reads are ~520-600 bp, but CBP isolates in the pipeline are only 301 bp (forward reads). Before trimming (Step 4), we need to determine if full-length CBP sequences exist and can be used. If not, we either:
- Trim PCR reads to 301 bp (losing half the merged data — defeating the purpose of merging)
- Find/generate full-length CBP sequences

### 2. rmlD Recovery
The rmlD gene has 571K reads but needs a corrected reference amplicon. This is fixable by manually identifying the correct primer binding sites in the spizizenii reference genomes.

### 3. Low-Count Genes
Atrophaeus genes (2-54 PCR reads) and comB/alaS (81/10 reads after filter) are too small for meaningful EcoSim runs on their own, but the CBP isolates (122-211 per gene) still provide ecotype structure. These genes could be used for CBP-only EcoSim validation.

## Files

- `data/02_filtered/{gene}.fasta` — Filtered PCR reads per gene
- `data/02_filtered/notes/filter_log.tsv` — Per-gene filter statistics
- `data/13_gene_filtered/gene_refs/` — Gene-specific BLAST databases used for filtering
