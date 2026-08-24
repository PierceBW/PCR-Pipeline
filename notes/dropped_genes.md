# Genes Excluded from EcoSim Analysis

## Summary

Of the 18 genes in our primer panel, **8 are usable** for EcoSim with substantial PCR data, and **10 are excluded** for various reasons.

### Usable Genes (8)

| Gene | Species | PCR Reads | CBP | Outgroup | Status |
|------|---------|----------:|----:|:--------:|--------|
| acuA | inaquosorum | 1,123,203 | 141 | yes | Clean, large dataset |
| iolB | spizizenii | 951,816 | 211 | yes | Clean, large dataset |
| sorA | inaquosorum | 936,217 | 141 | yes | Clean, large dataset |
| yvqK | inaquosorum | 443,124 | 141 | yes | Clean, large dataset |
| albG | inaquosorum | 356,881 | 141 | no | Clean, missing outgroup |
| thiD | inaquosorum | 356,834 | 142 | yes | Clean, large dataset |
| acsA_2 | spizizenii | 332,125 | 209 | yes | Clean, large dataset |
| alkH | spizizenii | 106,110 | 209 | yes | 36% dropped in filter, remaining are clean |

### Excluded Genes (10)

| Gene | Species | PCR Reads | Reason | Category |
|------|---------|----------:|--------|----------|
| accA | inaquosorum | 0 | No PCR amplification | No data |
| adc | atrophaeus | 0 | No PCR amplification | No data |
| rmlD | spizizenii | 0 | Reference extraction bug | Fixable |
| acdA | atrophaeus | 2 | Near-zero environmental amplification | Too few PCR |
| ecfA1 | atrophaeus | 3 | Near-zero environmental amplification | Too few PCR |
| yndE_2 | atrophaeus | 2 | Near-zero environmental amplification | Too few PCR |
| alaS | spizizenii | 10 | 99.5% dropped by gene filter | Cross-species contamination |
| amj | atrophaeus | 21 | Near-zero environmental amplification | Too few PCR |
| opuAB | atrophaeus | 32 | Near-zero environmental amplification | Too few PCR |
| comB | spizizenii | 80 | 81% dropped by gene filter | Cross-species contamination |

---

## Detailed Explanations

### No PCR Amplification (accA, adc)

**accA** (inaquosorum) had 32 read pairs in the BAM files, but all had wrong orientation (likely adapter dimers). Zero reads survived merging. The accA primers simply did not amplify anything usable from these environmental samples.

**adc** (atrophaeus) had zero read pairs across all 4 environments. Complete failure of PCR amplification.

### Atrophaeus Genes with Near-Zero PCR (acdA, ecfA1, yndE_2, amj, opuAB)

All 6 atrophaeus genes produced extremely few environmental PCR reads (2-32 after filtering). This is consistent across all 4 environments — *B. atrophaeus* is present at very low abundance in these environmental samples.

These genes still have their full CBP complement (~122-123 isolates per gene), so they could be used for CBP-only EcoSim validation. But they don't have enough PCR diversity to discover new ecotypes.

| Gene | PCR Reads | CBP | Outgroup | Notes |
|------|----------:|----:|:--------:|-------|
| acdA | 2 | 123 | yes | |
| ecfA1 | 3 | 122 | yes | |
| yndE_2 | 2 | 123 | no | Outgroup primer too divergent |
| amj | 21 | 122 | no | Outgroup primer too divergent |
| opuAB | 32 | 123 | no | Outgroup primer too divergent |

### Reference Extraction Bug (rmlD)

rmlD has 571K merged reads, but **zero passed the gene-specific identity filter**. This is NOT a data quality problem — the rmlD reference amplicon was extracted incorrectly.

When we BLASTed the rmlD primers against the spizizenii reference genomes, the primer hits were found too far apart:
- CP077772.1: extracted "amplicon" = 889 bp (expected ~558 bp)
- CP030925.1: extracted "amplicon" = 1102 bp

A ~570 bp PCR read can't align well to an 889 bp reference, so BLAST returns no significant hit for any read.

**Fix:** Manually identify the correct rmlD amplicon boundaries in the reference genomes. The CBP amplicons (210 isolates, 558 bp each) provide the ground truth for what the correct amplicon should look like.

**Impact:** rmlD is one of the 6 spizizenii genes and has the second-largest read count. Recovering it would add significant data.

### Cross-Species Contamination (alaS, comB)

**alaS** (spizizenii primers): 2,059 merged reads, but only 10 passed the 95% identity filter against the spizizenii alaS reference. The other 99.5% are inaquosorum reads that mapped to the alaS BAM reference because the two species share enough similarity at this locus for the reads to map, but not enough to pass the gene-specific filter. The 8 bp gap in the merged reads (alaS amplicon is ~610 bp, reads are 301+301=602) is not the cause — the filter failure is purely due to species-level divergence.

**comB** (spizizenii primers): 418 merged reads, only 80 passed (19%). Same cross-species contamination pattern as alaS, though less extreme. comB also has 210 CBP isolates, so the file isn't empty — just very few environmental PCR reads from the correct species.

Both genes' BAM files contain reads from multiple *Bacillus* species because the mapping reference (~810 bp) is long enough that reads from closely related species can map with decent alignment scores. The gene-specific filter (which checks the shorter ~500-600 bp amplicon region where species-diagnostic SNPs are concentrated) correctly removes these cross-species reads.

### Missing Outgroup (albG, amj, opuAB, yndE_2)

Four genes have no outgroup sequence because the outgroup species (*B. amyloliquefaciens* and *B. vallismortis*) are too divergent at these loci for either primer BLAST or CBP BLAST to find the homologous region.

**albG** has a partial outgroup (498 bp vs 518 bp CBP, 83% identity) which may be usable but is 20 bp short.

For the 3 atrophaeus genes (amj, opuAB, yndE_2), missing outgroup is moot since they have too few PCR reads anyway. For **albG**, the tree could be built unrooted or midpoint-rooted as an alternative.

---

## Can Any Excluded Genes Be Recovered?

| Gene | Recovery Possible? | Effort |
|------|:------------------:|--------|
| rmlD | **Yes** | Fix reference amplicon extraction — manual primer location |
| albG | **Partially** | Use partial outgroup or midpoint rooting |
| alaS | No | Reads are genuinely from the wrong species |
| comB | Marginal | Only 80 PCR reads even after recovery |
| Atrophaeus genes | No | Low species abundance, not a pipeline issue |
| accA, adc | No | No amplification from these samples |
