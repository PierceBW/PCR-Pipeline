# EcoSim n=2000 Ecotype Analysis

## Overview

Same approach as the n=200 run, but with 2000 PCR reads per gene (+ all CBP + outgroup). Total sequences per gene ~2140–2212.

## Sequence Composition per Gene

At n=2000, PCR reads dominate — CBP isolates are only 6.6–9.5% of the total. This is a very different balance compared to n=200, where CBP was 41–51%.

| Gene | Species | Total Seqs | Outgroup | CBP | PCR | CBP % | PCR % |
|------|---------|:---:|:---:|:---:|:---:|:---:|:---:|
| acsA_2 | spizizenii | 2,210 | 1 | 209 | 2,000 | 9.5% | 90.5% |
| acuA | inaquosorum | 2,142 | 1 | 141 | 2,000 | 6.6% | 93.4% |
| albG | inaquosorum | 2,141 | 0 | 141 | 2,000 | 6.6% | 93.4% |
| alkH | spizizenii | 2,210 | 1 | 209 | 2,000 | 9.5% | 90.5% |
| iolB | spizizenii | 2,212 | 1 | 211 | 2,000 | 9.5% | 90.4% |
| sorA | inaquosorum | 2,142 | 1 | 141 | 2,000 | 6.6% | 93.4% |
| thiD | inaquosorum | 2,143 | 1 | 142 | 2,000 | 6.6% | 93.3% |
| yvqK | inaquosorum | 2,142 | 1 | 141 | 2,000 | 6.6% | 93.4% |

At this ratio, CBP isolates are a small minority. This may contribute to the massive over-splitting — the environmental reads dominate the tree topology and EcoSim has less "anchor" structure from the CBP isolates.

## Summary Table

| Gene | Species | Total | Ecotypes | CBP-only | PCR-only (NEW) | Mixed | Singletons | npop | Thesis |
|------|---------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| acsA_2 | spizizenii | 2,210 | 1,017 | 122 | 894 | 1 | 944 | 69 | 5 |
| acuA | inaquosorum | 2,142 | 1,255 | 56 | 1,194 | 5 | 1,166 | 160 | 14 |
| albG | inaquosorum | 2,141 | 1,352 | 30 | 1,319 | 3 | 1,272 | 132 | 14 |
| alkH | spizizenii | 2,210 | 1,226 | 65 | 1,161 | 0 | 1,152 | 90 | 5 |
| iolB | spizizenii | 2,212 | 1,188 | 94 | 1,094 | 0 | 1,090 | 170 | 5 |
| sorA | inaquosorum | 2,142 | 1,084 | 42 | 1,032 | **10** | 999 | 78 | 14 |
| thiD | inaquosorum | 2,143 | 1,134 | 29 | 1,096 | **9** | 1,023 | 196 | 14 |
| yvqK | inaquosorum | 2,142 | 1,347 | 42 | 1,303 | 2 | 1,274 | 120 | 14 |

## Comparison: n=200 vs n=2000

| Gene | Ecotypes (200) | Ecotypes (2000) | Mixed (200) | Mixed (2000) | Singletons (200) | Singletons (2000) |
|------|:-:|:-:|:-:|:-:|:-:|:-:|
| acsA_2 | 149 | 1,017 | 0 | 1 | 134 | 944 |
| acuA | 197 | 1,255 | 2 | 5 | 165 | 1,166 |
| albG | 158 | 1,352 | 2 | 3 | 136 | 1,272 |
| alkH | 26 | 1,226 | 0 | 0 | 16 | 1,152 |
| iolB | 116 | 1,188 | 0 | 0 | 88 | 1,090 |
| sorA | 73 | 1,084 | **9** | **10** | 51 | 999 |
| thiD | 81 | 1,134 | 3 | **9** | 58 | 1,023 |
| yvqK | 66 | 1,347 | 3 | 2 | 50 | 1,274 |

**Key trend:** More PCR reads = massively more ecotypes (mostly singletons). The mixed ecotype count stays relatively stable or increases modestly — sorA remains the champion.

---

## Per-Gene Detailed Breakdown

### acsA_2 (spizizenii) — 1,017 ecotypes

- **209 CBP + 2,000 PCR** across 1,017 ecotypes
- **1 mixed ecotype** (up from 0 at n=200)
- **122 CBP-only** (194 isolates), **894 PCR-only** (1,999 reads, 825 singletons)

**Mixed ecotype:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #897 | 16 | 15 | 1 | S1, S2, S4 | E2:1 |

Just one PCR read merged with CBP — essentially still separated.

**Largest CBP-only:** Eco#1017 with 64 CBP (S1, S2, S3)

**PCR-only pattern:** All large PCR ecotypes come from E2 (environment 2).

---

### acuA (inaquosorum) — 1,255 ecotypes

- **141 CBP + 2,000 PCR** across 1,255 ecotypes
- **5 mixed ecotypes** (up from 2 at n=200)
- **56 CBP-only** (119 isolates), **1,194 PCR-only** (1,954 reads, 1,117 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #1138 | 23 | 3 | 20 | I7 | E4:3 E5:17 |
| #666 | 17 | 8 | 9 | I8 | E4:9 |
| #1070 | 15 | 9 | 6 | I6, I7, I9, I10 | E1:1 E2:4 E4:1 |
| #1061 | 11 | 1 | 10 | I9 | E1:1 E2:4 E4:3 E5:2 |
| #1255 | 2 | 1 | 1 | I14 | E2:1 |

Good diversity of thesis labels in mixed ecotypes — I7, I8, I9, I10, I14 all have environmental matches.

---

### albG (inaquosorum) — 1,352 ecotypes

- **140 CBP + 2,000 PCR** across 1,352 ecotypes
- **No outgroup** — unrooted tree
- **3 mixed ecotypes**
- **30 CBP-only** (100 isolates), **1,319 PCR-only** (1,986 reads, 1,251 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #1341 | 32 | 31 | 1 | I1, I2, I3, I5, I6, I13 | E4:1 |
| #1343 | 20 | 8 | 12 | I7 | E2:1 E4:5 E5:6 |
| #1342 | 2 | 1 | 1 | I5 | E4:1 |

---

### alkH (spizizenii) — 1,226 ecotypes

- **209 CBP + 2,000 PCR** across 1,226 ecotypes
- **0 mixed ecotypes** — complete separation (same as n=200)
- **65 CBP-only** (209 isolates), **1,161 PCR-only** (2,000 reads, 1,096 singletons)

alkH remains completely separated between CBP and PCR at both sample sizes.

**Largest CBP-only:** Eco#1175 and #1221 each with 49 CBP (S1, S2)

**Largest PCR-only:** Eco#519 with 68 reads (E4:65) — strong E4 dominance continues.

---

### iolB (spizizenii) — 1,188 ecotypes

- **211 CBP + 2,000 PCR** across 1,188 ecotypes
- **0 mixed ecotypes** — complete separation (same as n=200)
- **94 CBP-only** (211 isolates), **1,094 PCR-only** (2,000 reads, 1,007 singletons)

Like alkH, iolB shows zero overlap between culture collection and environmental reads.

---

### sorA (inaquosorum) — 1,084 ecotypes ★ BEST GENE

- **141 CBP + 2,000 PCR** across 1,084 ecotypes
- **10 mixed ecotypes** (up from 9 at n=200) — still the most overlap
- 41 CBP + 137 PCR in mixed ecotypes
- **42 CBP-only** (100 isolates), **1,032 PCR-only** (1,863 reads, 961 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #814 | 61 | 1 | 60 | I4 | E1:1 E2:1 E4:5 E5:53 |
| #811 | 51 | 6 | 45 | I9 | E1:2 E4:3 E5:40 |
| #5 | 14 | 11 | 3 | I4, I5, I6, I7, I9 | E4:3 |
| #1077 | 13 | 2 | 11 | I7 | E4:11 |
| #223 | 9 | 8 | 1 | I8 | E4:1 |
| #1070 | 7 | 1 | 6 | I10 | E2:1 E4:1 E5:4 |
| #222 | 6 | 3 | 3 | I8 | E1:1 E5:2 |
| #806 | 6 | 1 | 5 | I6 | E4:5 |
| #1074 | 6 | 4 | 2 | I3 | E5:2 |
| #387 | 5 | 4 | 1 | I1, I5, I12 | E1:1 |

**10 different thesis ecotype labels** (I1, I3, I4, I5, I6, I7, I8, I9, I10, I12) have environmental reads matching them. This is extraordinary coverage — nearly all of Jocelyn's inaquosorum ecotypes (14 expected) have environmental representation in sorA.

---

### thiD (inaquosorum) — 1,134 ecotypes

- **142 CBP + 2,000 PCR** across 1,134 ecotypes
- **9 mixed ecotypes** (up from 3 at n=200) — big improvement
- 76 CBP + 37 PCR in mixed ecotypes
- **29 CBP-only** (66 isolates), **1,096 PCR-only** (1,963 reads, 998 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #673 | 24 | 17 | 7 | I1, I2, I3, I5, I7, I10 | E4:7 |
| #170 | 21 | 20 | 1 | I5, I7, I9, I11, I12 | E2:1 |
| #193 | 20 | 1 | 19 | I7 | E2:2 E5:17 |
| #674 | 19 | 15 | 4 | I7, I8 | E4:4 |
| #696 | 18 | 17 | 1 | I1, I4, I6, I13 | E4:1 |
| #168 | 4 | 2 | 2 | I3, I5 | E1:1 E5:1 |
| #679 | 3 | 2 | 1 | I5 | E4:1 |
| #14 | 2 | 1 | 1 | I5 | E2:1 |
| #666 | 2 | 1 | 1 | I10 | E4:1 |

thiD improved dramatically from n=200 (3 mixed) to n=2000 (9 mixed). More PCR reads allowed detection of matches to thesis labels I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12, I13 — 13 of 14 expected inaquosorum ecotypes represented!

---

### yvqK (inaquosorum) — 1,347 ecotypes

- **141 CBP + 2,000 PCR** across 1,347 ecotypes
- **2 mixed ecotypes** (down from 3 at n=200)
- **42 CBP-only** (118 isolates), **1,303 PCR-only** (1,998 reads, 1,239 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #969 | 18 | 17 | 1 | I2, I4, I6 | E5:1 |
| #4 | 7 | 6 | 1 | I7, I8 | E2:1 |

Minimal overlap at n=2000 — environmental reads are largely novel relative to CBP for this gene.

---

## Key Takeaways

### 1. Massive over-splitting at n=2000

Every gene has 1,000+ ecotypes with the majority being singletons (~90%). EcoSim's demarcation algorithm is being overwhelmed by the diversity of 2,000 environmental reads. This is a known issue — the algorithm performs best with moderate-sized datasets.

### 2. sorA and thiD are the best-performing genes

- **sorA**: 10 mixed ecotypes covering thesis labels I1, I3, I4, I5, I6, I7, I8, I9, I10, I12 (10 of 14)
- **thiD**: 9 mixed ecotypes covering thesis labels I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12, I13 (13 of 14)

These two genes provide the strongest evidence that environmental amplicon reads can be linked back to known ecotypes from the culture collection.

### 3. Spizizenii genes show zero or near-zero overlap

alkH and iolB have zero mixed ecotypes at both sample sizes. acsA_2 has just 1 (with only 1 PCR read). This suggests the environmental spizizenii diversity is substantially different from the culture collection — or that spizizenii is less abundant in these environments.

### 4. The n=200 subsets may be more informative

Fewer ecotypes, fewer singletons, and the mixed ecotype counts are similar. For practical EcoSim analysis, smaller subsets with balanced CBP/PCR ratios may give cleaner demarcation.

### 5. Environment-specific clustering persists

Large PCR-only ecotypes continue to draw predominantly from single environments, even at n=2000. This ecological signal is robust.

---

## Files

- Results TSV: `data/05_ecosim/2000/notes/ecosim_results.tsv`
- Per-gene data: `data/05_ecosim/2000/{gene}/`
  - `{gene}.fasta` — subsampled sequences
  - `{gene}.nwk` — rerooted tree
  - `{gene}.xml` — EcoSim XML results
  - `{gene}_ecosim_log.txt` — full EcoSim output
  - `{gene}_header_map.tsv` — PCR header mapping (short → original)
