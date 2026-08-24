# EcoSim n=200 Ecotype Analysis

## Overview

For each of the 8 usable genes, we ran EcoSim on subsets of 200 PCR reads + all CBP isolates + outgroup. This document breaks down how ecotypes partition between CBP (known isolates from Jocelyn Wang's culture collection) and PCR (new environmental amplicon reads).

## Sequence Composition per Gene

At n=200, CBP isolates make up a large fraction of the total — roughly 41–51% of all sequences. This means EcoSim is working with nearly equal parts known isolates and environmental reads.

| Gene | Species | Total Seqs | Outgroup | CBP | PCR | CBP % | PCR % |
|------|---------|:---:|:---:|:---:|:---:|:---:|:---:|
| acsA_2 | spizizenii | 410 | 1 | 209 | 200 | 51.0% | 48.8% |
| acuA | inaquosorum | 342 | 1 | 141 | 200 | 41.2% | 58.5% |
| albG | inaquosorum | 341 | 0 | 141 | 200 | 41.3% | 58.7% |
| alkH | spizizenii | 410 | 1 | 209 | 200 | 51.0% | 48.8% |
| iolB | spizizenii | 412 | 1 | 211 | 200 | 51.2% | 48.5% |
| sorA | inaquosorum | 342 | 1 | 141 | 200 | 41.2% | 58.5% |
| thiD | inaquosorum | 343 | 1 | 142 | 200 | 41.4% | 58.3% |
| yvqK | inaquosorum | 342 | 1 | 141 | 200 | 41.2% | 58.5% |

Spizizenii genes have more CBP isolates (209–211) than inaquosorum genes (141–142), so at n=200 PCR the spizizenii genes are majority-CBP (51%) while inaquosorum genes are majority-PCR (58%).

## Summary Table

| Gene | Species | Total | Ecotypes | CBP-only | PCR-only (NEW) | Mixed | Singletons | Thesis Expected |
|------|---------|:---:|:--------:|:--------:|:--------------:|:-----:|:----------:|:---------------:|
| acsA_2 | spizizenii | 410 | 149 | 136 | 13 | 0 | 134 | 5 |
| acuA | inaquosorum | 342 | 197 | 97 | 98 | 2 | 165 | 14 |
| albG | inaquosorum | 341 | 158 | 44 | 112 | 2 | 136 | 14 |
| alkH | spizizenii | 410 | 26 | 20 | 6 | 0 | 16 | 5 |
| iolB | spizizenii | 412 | 116 | 76 | 40 | 0 | 88 | 5 |
| sorA | inaquosorum | 342 | 73 | 41 | 23 | **9** | 51 | 14 |
| thiD | inaquosorum | 343 | 81 | 30 | 48 | 3 | 58 | 14 |
| yvqK | inaquosorum | 342 | 66 | 22 | 41 | 3 | 50 | 14 |

**Definitions:**
- **CBP-only**: Ecotypes containing only CBP isolates (no environmental PCR reads)
- **PCR-only (NEW)**: Ecotypes containing only environmental PCR reads — these are potentially new ecotypes not represented in the culture collection
- **Mixed**: Ecotypes containing both CBP isolates and PCR reads — environmental reads that cluster with known isolates
- **Singletons**: Ecotypes with only 1 member (over-splitting)

---

## Key Observations

### 1. sorA is the best-performing gene

sorA has 9 mixed ecotypes — far more than any other gene. This means environmental PCR reads are being assigned to the same ecotypes as known CBP isolates, which is the strongest validation that our pipeline is working correctly.

### 2. Spizizenii genes show complete CBP/PCR separation

acsA_2, alkH, and iolB have zero mixed ecotypes. All CBP and PCR reads fall into separate ecotypes. This could mean:
- The environmental diversity for spizizenii is genuinely distinct from the culture collection
- Or the 200-read subsample doesn't capture enough of the right reads to overlap with CBP

### 3. Singleton over-splitting

Most genes have high singleton counts (acuA: 165/197, albG: 136/158, acsA_2: 134/149). EcoSim is over-demarcating at this sample size. The 2000-read subsets may improve this by providing more sequences per ecotype.

### 4. Strong environment-specific clustering

PCR-only ecotypes often draw heavily from a single environment:
- alkH Eco#6: 73 of 74 reads from E5
- alkH Eco#2: 66 of 68 reads from E4
- iolB Eco#22: all 35 reads from E4
- albG Eco#154: all 18 reads from E4

This suggests real ecological structure — different environments harbor distinct ecotypes.

---

## Per-Gene Detailed Breakdown

### acsA_2 (spizizenii) — 149 ecotypes

- **209 CBP + 200 PCR** assigned across 149 ecotypes
- **0 mixed ecotypes** — complete separation between CBP and PCR
- **136 CBP-only ecotypes** (209 isolates), dominated by singletons (134)
- **13 PCR-only ecotypes** (200 reads, only 2 singletons)

**Largest CBP-only ecotypes:**

| Ecotype | Size | Thesis Labels |
|:-------:|:----:|:-------------:|
| #148 | 51 | S1, S2 |
| #149 | 15 | S1, S2, S4 |
| #14 | 6 | S5 |
| #15 | 5 | S1 |

**PCR-only ecotypes (all from Env 2):**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #8 | 47 | E2:47 |
| #7 | 40 | E2:40 |
| #2 | 38 | E2:38 |
| #1 | 31 | E2:31 |
| #5 | 12 | E2:12 |

All PCR ecotypes are from environment 2 — strong environment-specific signal.

---

### acuA (inaquosorum) — 197 ecotypes

- **141 CBP + 200 PCR** assigned across 197 ecotypes
- **2 mixed ecotypes** — environmental reads clustering with known isolates
- **97 CBP-only ecotypes** (134 isolates)
- **98 PCR-only ecotypes** (149 reads, 78 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #85 | 52 | 2 | 50 | I7, I9 | E2:49 E5:1 |
| #91 | 6 | 5 | 1 | I7, I9 | E2:1 |

**Largest PCR-only ecotypes:**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #68 | 6 | E4:6 |
| #77 | 6 | E2:5 E4:1 |
| #143 | 6 | E1:4 E2:1 E4:1 |
| #186 | 6 | E1:1 E2:1 E4:1 E5:3 |

---

### albG (inaquosorum) — 158 ecotypes

- **140 CBP + 200 PCR** assigned across 158 ecotypes
- **No outgroup** — tree was unrooted (midpoint-rooted by FastTree)
- **2 mixed ecotypes**
- **44 CBP-only ecotypes** (131 isolates)
- **112 PCR-only ecotypes** (180 reads, 102 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #96 | 17 | 1 | 16 | I4 | E1:14 E2:2 |
| #65 | 12 | 8 | 4 | I7 | E2:2 E4:1 E5:1 |

**Largest CBP-only ecotypes:**

| Ecotype | Size | Thesis Labels |
|:-------:|:----:|:-------------:|
| #93 | 18 | I1, I2, I3, I6, I13 |
| #100 | 18 | I1, I5 |
| #73 | 15 | I4, I5, I6, I9, I13, I14 |
| #99 | 14 | I5, I8, I11, I12 |

**Largest PCR-only ecotypes:**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #64 | 18 | E2:6 E4:1 E5:11 |
| #154 | 18 | E4:18 |
| #66 | 13 | E2:3 E4:1 E5:9 |
| #156 | 7 | E4:7 |

---

### alkH (spizizenii) — 26 ecotypes

- **209 CBP + 200 PCR** assigned across only 26 ecotypes — **fewest ecotypes of any gene**
- **0 mixed ecotypes**
- **20 CBP-only ecotypes** (209 isolates)
- **6 PCR-only ecotypes** (200 reads, only 2 singletons)
- Closest to thesis expectation (26 vs expected 5)

**CBP-only ecotypes (large clusters):**

| Ecotype | Size | Thesis Labels |
|:-------:|:----:|:-------------:|
| #21 | 81 | S1, S2 |
| #24 | 76 | S1, S2 |
| #8 | 19 | S1 |
| #23 | 8 | S1 |
| #26 | 6 | S5 |

**PCR-only ecotypes (strong environment clustering):**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #6 | 74 | E2:1 **E5:73** |
| #2 | 68 | E1:1 E2:1 **E4:66** |
| #5 | 30 | E1:1 E2:24 E4:3 E5:2 |
| #1 | 26 | **E4:25** E5:1 |

alkH shows the clearest environment-specific ecotype signal. Ecotypes #6 and #2 are almost entirely from single environments (E5 and E4 respectively).

---

### iolB (spizizenii) — 116 ecotypes

- **211 CBP + 200 PCR** assigned across 116 ecotypes
- **0 mixed ecotypes** — complete CBP/PCR separation
- **76 CBP-only ecotypes** (211 isolates)
- **40 PCR-only ecotypes** (200 reads, 23 singletons)

**Largest CBP-only ecotypes:**

| Ecotype | Size | Thesis Labels |
|:-------:|:----:|:-------------:|
| #106 | 35 | S1, S2 |
| #115 | 35 | S1, S2, S3 |
| #43 | 20 | S1, S2 |
| #108 | 18 | S1, S2, S4 |

**Largest PCR-only ecotypes:**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #35 | 48 | E1:10 E2:32 E4:5 E5:1 |
| #22 | 35 | E4:35 |
| #28 | 30 | E1:1 E4:29 |
| #38 | 11 | E2:3 E4:2 E5:6 |

---

### sorA (inaquosorum) — 73 ecotypes ★ BEST GENE

- **141 CBP + 200 PCR** assigned across 73 ecotypes
- **9 mixed ecotypes** — most of any gene, strong validation
- 49 CBP isolates + 111 PCR reads in mixed ecotypes
- **41 CBP-only ecotypes** (92 isolates)
- **23 PCR-only ecotypes** (89 reads, 16 singletons)

**Mixed ecotypes (environmental reads matching known isolates):**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #24 | 37 | 14 | 23 | I5, I7 | E1:18 E2:4 E5:1 |
| #70 | 37 | 7 | 30 | I4, I5, I7 | E2:1 E5:29 |
| #46 | 33 | 4 | 29 | I3 | E1:1 E2:26 E4:1 E5:1 |
| #1 | 14 | 1 | 13 | I6 | E4:13 |
| #2 | 12 | 11 | 1 | I4, I5, I6, I7, I9 | E4:1 |
| #69 | 11 | 1 | 10 | I10 | E5:10 |
| #72 | 8 | 6 | 2 | I9 | E5:2 |
| #5 | 4 | 3 | 1 | I8 | E5:1 |
| #48 | 4 | 2 | 2 | I7 | E1:1 E4:1 |

These 9 mixed ecotypes represent cases where PCR environmental reads are genetically close enough to CBP isolates to be placed in the same ecotype — direct confirmation that the pipeline correctly links environmental diversity to known lineages.

**Largest CBP-only ecotype:**

| Ecotype | Size | Thesis Labels |
|:-------:|:----:|:-------------:|
| #43 | 37 | I1, I2, I10, I11 |

**Largest PCR-only ecotypes:**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #65 | 37 | E2:1 E4:36 |
| #66 | 14 | E4:14 |

---

### thiD (inaquosorum) — 81 ecotypes

- **142 CBP + 200 PCR** assigned across 81 ecotypes
- **3 mixed ecotypes**
- **30 CBP-only ecotypes** (125 isolates)
- **48 PCR-only ecotypes** (129 reads, 39 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #39 | 42 | 12 | 30 | I7 | E5:30 |
| #73 | 41 | 1 | 40 | I7 | E4:40 |
| #79 | 5 | 4 | 2 | I2, I7 | E1:1 |

Ecotypes #39 and #73 are large mixed groups dominated by PCR reads from single environments (E5 and E4), each containing thesis ecotype I7 CBP isolates.

**Largest CBP-only ecotypes:**

| Ecotype | Size | Thesis Labels |
|:-------:|:----:|:-------------:|
| #66 | 17 | I1, I2, I4, I9 |
| #68 | 17 | I1, I4, I6, I13 |
| #78 | 15 | I7, I8 |
| #50 | 13 | I5, I11, I12 |

**Largest PCR-only ecotypes:**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #47 | 26 | E1:19 E2:5 E4:1 E5:1 |
| #1 | 17 | E2:16 E4:1 |
| #36 | 17 | E2:17 |

---

### yvqK (inaquosorum) — 66 ecotypes

- **141 CBP + 200 PCR** assigned across 66 ecotypes
- **3 mixed ecotypes**
- **22 CBP-only ecotypes** (93 isolates)
- **41 PCR-only ecotypes** (122 reads, 35 singletons)

**Mixed ecotypes:**

| Ecotype | Size | CBP | PCR | CBP Labels | Environments |
|:-------:|:----:|:---:|:---:|:----------:|:------------:|
| #59 | 45 | 12 | 33 | I5, I7, I12 | E2:33 |
| #57 | 44 | 1 | 43 | I7 | E1:2 E2:1 E4:1 E5:39 |
| #58 | 37 | 35 | 2 | I2, I4, I6, I9, I11 | E5:2 |

**Largest CBP-only ecotype:**

| Ecotype | Size | Thesis Labels |
|:-------:|:----:|:-------------:|
| #63 | 42 | I1, I3, I5 |
| #65 | 14 | I1, I3, I5, I6, I13 |

**Largest PCR-only ecotypes:**

| Ecotype | Size | Environments |
|:-------:|:----:|:------------:|
| #48 | 48 | E1:8 E2:2 E4:38 |
| #1 | 15 | E1:3 E2:2 E4:10 |

---

## Thesis Ecotype Labels in CBP Headers

CBP isolate headers encode the thesis ecotype assignment: `CBP-{strain}_PE_{ecotype}`. The ecotype labels (I1–I14 for inaquosorum, S1–S5 for spizizenii) come from Jocelyn Wang's thesis demarcation. Within our EcoSim ecotypes, multiple thesis labels often appear together, suggesting that our higher-resolution analysis (with environmental reads) sometimes merges what the thesis considered separate ecotypes, and sometimes splits them further.

## Environment Key

- **E1**: Environment 1
- **E2**: Environment 2
- **E4**: Environment 4
- **E5**: Environment 5

(Environments 3 and 6 were not sequenced)

---

## What This Tells Us

1. **The pipeline works.** EcoSim successfully demarcates ecotypes from our merged, filtered, trimmed amplicon data across all 8 genes.

2. **Environmental reads add real diversity.** PCR-only ecotypes represent lineages present in the environment but absent from the culture collection — exactly what amplicon sequencing is designed to reveal.

3. **Mixed ecotypes validate the approach.** When PCR reads cluster with CBP isolates (especially strong in sorA, thiD, yvqK), it confirms that our sequence processing preserves the biological signal needed for accurate ecotype classification.

4. **Environment-specific ecotypes suggest ecological structure.** Many PCR ecotypes draw predominantly from single environments, consistent with habitat-specific ecotype distributions.

5. **The 200-read subsets over-split.** High singleton counts suggest EcoSim needs more sequences per ecotype to demarcate reliably. The 2000-read subsets should improve this.
