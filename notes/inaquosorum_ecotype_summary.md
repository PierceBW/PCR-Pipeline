# Inaquosorum Ecotype Summary — Post-Deduplication Results

## What This Document Shows

This summarizes the EcoSim ecotype demarcation results for all 5 usable **B. inaquosorum** genes, using the **pruned in-group dataset** (step 8) after **deduplication of identical sequences** (step 4b). These are the results where all PCR environmental reads fall within the CBP (culture-based phylogenetics) LCA clade — meaning the environmental reads are confirmed to be the same species as Jocelyn Wang's thesis isolates.

Each gene's tree contains deduplicated CBP isolates (17–26 unique sequences from ~141 originals) + 200 randomly subsampled PCR amplicon reads + 1 outgroup. EcoSim was run on these combined trees to demarcate ecotypes.

---

## Key Definitions

- **CBP isolates**: Culture-based phylogenetics isolates from Jocelyn Wang's thesis, each assigned to a thesis ecotype (I1–I14)
- **PCR reads**: Environmental amplicon sequences, labeled by source environment (E1, E2, E4, E5)
- **npop**: Estimated number of populations from EcoSim's hill-climbing algorithm
- **npop CI**: 95% confidence interval for npop
- **Demarcated ecotypes**: Number of ecotypes identified by EcoSim's monophyly-based demarcation
- **Singleton**: An ecotype containing only 1 sequence
- **Doubleton**: An ecotype containing exactly 2 sequences
- **Mixed ecotype**: An ecotype containing both CBP isolates AND PCR environmental reads
- **High ecotype gene**: A gene that produced significantly MORE ecotypes than expected in Jocelyn's thesis CBP-only analysis (e.g., sorA found 41 ecotypes vs the modal 14)
- **Modal ecotype gene**: A gene that produced closer to the average ~14 ecotypes in Jocelyn's thesis CBP-only analysis

---

## CBP Deduplication

Trimming Jocelyn's full-length isolate amplicons to match PCR read length collapses most CBP isolates into identical sequences. These duplicates were removed before tree building and EcoSim (step 4b):

| Gene | CBP Isolates | Unique Sequences | Duplicates Removed | % Duplicate |
|------|:------------:|:----------------:|:------------------:|:-----------:|
| sorA | 141 | 17 | 124 | 88% |
| yvqK | 141 | 21 | 120 | 85% |
| thiD | 142 | 26 | 116 | 82% |
| albG | 141 | 25 | 116 | 82% |
| acuA | 141 | 18 | 123 | 87% |

PCR reads were also deduplicated. The dedup log is at `data/04b_deduped/notes/dedup_log.tsv`.

---

## Summary Table

| Gene | Type | Seqs | Len (bp) | npop | npop CI | Ecotypes | Singletons (CBP / PCR) | Doubletons (CBP / PCR / Mixed) | Larger (CBP / PCR / Mixed) | Chao1 | Wang thesis ecotypes (ref) |
|------|------|:----:|:--------:|:----:|:-------:|:--------:|:----------------------:|:-----------------------------:|:--------------------------:|:-----:|:------:|
| **sorA** | High | 218 | 522 | 77 | 44–89 | 76 | 21 (5 / 16) | 21 (0 / 18 / **3**) | 34 (1 / 29 / **4**) | 86 | 41 |
| **yvqK** | High | 222 | 508 | 80 | 2–95 | 36 | 10 (3 / 7) | 6 (1 / 5 / 0) | 20 (0 / 15 / **5**) | 44 | 39 |
| **thiD** | High | 227 | 557 | 156 | 48–225 | 33 | 11 (6 / 5) | 5 (1 / 4 / 0) | 17 (2 / 10 / **5**) | 45 | 41 |
| **albG** | Modal | 225 | 518 | 16 | 4–25 | 34 | 8 (1 / 7) | 7 (0 / 7 / 0) | 19 (0 / 16 / **3**) | 39 | 14 |
| **acuA** | Modal | 219 | 522 | **4** | 4–4 | 12 | 6 (5 / 1) | 1 (0 / 0 / **1**) | 5 (0 / 1 / **4**) | 30 | 14 |

Singletons are split as (CBP-only / PCR-only). Doubletons and larger ecotypes are split as (CBP-only / PCR-only / Mixed). Mixed ecotypes contain both CBP and PCR sequences.

**Chao1** = S_obs + f1²/(2·f2), where f1 = singletons, f2 = doubletons. Estimates total ecotype richness including unobserved ecotypes.

**"Type"** = Whether Jocelyn's thesis classified this gene as producing high or modal ecotype counts using CBP data alone (Thesis Table S2/S3). sorA (41), thiD (41), and yvqK (39) are "high ecotype" genes (≥1.5× the modal 14). albG (14) and acuA (14) are "modal" genes.

---

## Singleton Breakdown

| Gene | Total Singletons | CBP Singletons | PCR Singletons | % Singletons |
|------|:----------------:|:--------------:|:--------------:|:------------:|
| sorA | 21 | 5 (29% of CBP) | 16 (8% of PCR) | 28% |
| yvqK | 10 | 3 (14% of CBP) | 7 (4% of PCR) | 28% |
| thiD | 11 | 6 (23% of CBP) | 5 (2% of PCR) | 33% |
| albG | 8 | 1 (4% of CBP) | 7 (4% of PCR) | 24% |
| acuA | 6 | 5 (28% of CBP) | 1 (0.5% of PCR) | 50% |

### CBP singletons by thesis ecotype

| Thesis Ecotype | sorA | yvqK | thiD | albG | acuA |
|:--------------:|:----:|:----:|:----:|:----:|:----:|
| I1 | 1/1 (100%) | — | — | — | 2/2 (100%) |
| I2 | — | — | 1/2 (50%) | — | 1/1 (100%) |
| I3 | 1/2 (50%) | — | — | — | — |
| I4 | — | — | 1/1 (100%) | — | — |
| I5 | 1/3 (33%) | — | 2/5 (40%) | — | 1/4 (25%) |
| I7 | — | 1/4 (25%) | — | 1/2 (50%) | — |
| I8 | 2/2 (100%) | — | — | — | — |
| I10 | — | 1/2 (50%) | — | — | — |
| I14 | — | 1/1 (100%) | 1/1 (100%) | — | 1/1 (100%) |
| I9 | — | — | 1/2 (50%) | — | — |

"—" = 0 singletons for that ecotype/gene. Fraction is singletons/total unique CBP sequences of that ecotype. Only ecotypes with at least one singleton shown.

### PCR singletons by environment

| Environment | sorA | yvqK | thiD | albG | acuA |
|:-----------:|:----:|:----:|:----:|:----:|:----:|
| E1 | 3/51 (6%) | 4/36 (11%) | — | 1/25 (4%) | — |
| E2 | 4/41 (10%) | 1/43 (2%) | 3/69 (4%) | — | — |
| E4 | 1/60 (2%) | 2/74 (3%) | 2/42 (5%) | 6/80 (8%) | 1/98 (1%) |
| E5 | 8/48 (17%) | — | — | — | — |

"—" = 0 singletons. Fraction is singletons/total PCR reads from that environment.

---

## Ecotype Composition

| Gene | CBP-Only Ecotypes | PCR-Only Ecotypes | Mixed Ecotypes |
|------|:-----------------:|:-----------------:|:--------------:|
| sorA | 6 | 63 | **7** |
| yvqK | 4 | 27 | **5** |
| thiD | 9 | 19 | **5** |
| albG | 1 | 30 | **3** |
| acuA | 5 | 2 | **5** |

---

## Mixed Ecotype Details

### sorA — 7 mixed ecotypes

| Ecotype | Size | CBP | PCR | Thesis Labels | Environments |
|:-------:|:----:|:---:|:---:|:-------------:|:------------:|
| #1 | 2 | 1 | 1 | I6 | E4 |
| #5 | 3 | 2 | 1 | I12, I5 | E2 |
| #9 | 5 | 1 | 4 | I7 | E4 |
| #21 | 2 | 1 | 1 | I6 | E1 |
| #44 | 7 | 1 | 6 | I7 | E4 |
| #65 | 5 | 2 | 3 | I4, I9 | E5 |
| #67 | 2 | 1 | 1 | I10 | E5 |

7 of 14 thesis ecotypes (I4, I5, I6, I7, I9, I10, I12) appear in mixed ecotypes for sorA.

### yvqK — 5 mixed ecotypes

| Ecotype | Size | CBP | PCR | Thesis Labels | Environments |
|:-------:|:----:|:---:|:---:|:-------------:|:------------:|
| #3 | 18 | 1 | 17 | I7 | E1, E4 |
| #22 | 38 | 1 | 37 | I5 | E1, E2, E5 |
| #29 | 10 | 1 | 9 | I7 | E1, E4, E5 |
| #33 | 3 | 2 | 1 | I10, I7 | E2 |
| #34 | 12 | 11 | 1 | I1, I13, I3, I6, I8 | E4 |

6 of 14 thesis ecotypes (I1, I3, I5, I6, I7, I8, I10, I13) appear in mixed ecotypes for yvqK.

### thiD — 5 mixed ecotypes

| Ecotype | Size | CBP | PCR | Thesis Labels | Environments |
|:-------:|:----:|:---:|:---:|:-------------:|:------------:|
| #1 | 56 | 2 | 54 | I7 | E1, E2, E4, E5 |
| #13 | 45 | 4 | 41 | I11, I5, I9 | E1, E2 |
| #20 | 8 | 2 | 6 | I10, I7 | E2, E4 |
| #27 | 5 | 1 | 4 | I8 | E4 |
| #30 | 3 | 2 | 1 | I7, I8 | E4 |

6 of 14 thesis ecotypes (I5, I7, I8, I9, I10, I11) appear in mixed ecotypes for thiD.

### albG — 3 mixed ecotypes

| Ecotype | Size | CBP | PCR | Thesis Labels | Environments |
|:-------:|:----:|:---:|:---:|:-------------:|:------------:|
| #16 | 16 | 15 | 1 | I1, I10, I11, I12, I2, I3, I5, I8, I9 | E4 |
| #18 | 28 | 7 | 21 | I1, I3, I4, I5, I6, I8, I9 | E1, E2, E5 |
| #29 | 9 | 1 | 8 | I7 | E4 |

11 of 14 thesis ecotypes (I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12) appear in mixed ecotypes for albG.

### acuA — 5 mixed ecotypes

| Ecotype | Size | CBP | PCR | Thesis Labels | Environments |
|:-------:|:----:|:---:|:---:|:-------------:|:------------:|
| #1 | 68 | 3 | 65 | I13, I6, I7 | E1, E2, E4, E5 |
| #4 | 2 | 1 | 1 | I8 | E4 |
| #8 | 21 | 2 | 19 | I11, I9 | E1, E2, E4 |
| #9 | 3 | 2 | 1 | I5 | E1 |
| #11 | 24 | 5 | 19 | I5, I6, I7 | E2, E4, E5 |

8 of 14 thesis ecotypes (I5, I6, I7, I8, I9, I11, I13) appear in mixed ecotypes for acuA.

---

## EcoSim Parameters

| Gene | npop | Omega | Sigma | Likelihood |
|------|:----:|:-----:|:-----:|:----------:|
| sorA | 77 | 0.003 | 3.611 | 1.000 |
| yvqK | 80 | 0.036 | 0.203 | 0.940 |
| thiD | 156 | 0.056 | 0.079 | 0.330 |
| albG | 16 | 0.329 | 0.020 | 0.003 |
| acuA | 4 | 0.023 | 0.446 | 0.011 |

- **Omega**: Rate of periodic selection (ecotype-formation events)
- **Sigma**: Rate of drift
- **Likelihood**: How well the model fits the data

---

## Comparison to Thesis Expectations

Jocelyn's thesis found a modal 14 ecotypes for B. inaquosorum using CBP core-genome data, but individual genes varied widely (Table S2/S3). Our combined CBP+PCR results after deduplication:

| Gene | Type | Thesis ecotypes (single gene) | Our npop | npop CI | Our ecotypes |
|------|------|:----------------------------:|:--------:|:-------:|:------------:|
| **acuA** | Modal | 14 | **4** | 4–4 | 12 |
| **albG** | Modal | 14 | **16** | 4–25 | 34 |
| **sorA** | High | 41 | 77 | 44–89 | 76 |
| **yvqK** | High | 39 | 80 | 2–95 | 36 |
| **thiD** | High | 41 | 156 | 48–225 | 33 |

---

## Comparison: Pre-Dedup vs Post-Dedup

| Gene | Pre-Dedup Seqs | Post-Dedup Seqs | Pre npop | Post npop | Pre Ecotypes | Post Ecotypes | Pre Singletons (%) | Post Singletons (%) | Pre Mixed | Post Mixed |
|------|:--------------:|:---------------:|:--------:|:---------:|:------------:|:-------------:|:------------------:|:-------------------:|:---------:|:----------:|
| sorA | 342 | 218 | 16 | 77 | 58 | 76 | 36 (62%) | 21 (28%) | 9 | 7 |
| yvqK | 342 | 222 | 12 | 80 | 56 | 36 | 42 (75%) | 10 (28%) | 3 | 5 |
| thiD | 343 | 227 | 43 | 156 | 124 | 33 | 89 (72%) | 11 (33%) | 2 | 5 |
| albG | 341 | 225 | 37 | 16 | 158 | 34 | 136 (86%) | 8 (24%) | 2 | 3 |
| acuA | 342 | 219 | 115 | 4 | 208 | 12 | 176 (85%) | 6 (50%) | 2 | 5 |

Pre-dedup results are in `data/old_pre_dedup/`.

---

## Notes

- Deduplication dramatically reduced ecotype counts (e.g., acuA: 208 → 12, albG: 158 → 34)
- Singleton rates dropped from 62–86% to 24–50%
- Mixed ecotype counts increased for yvqK (3→5), thiD (2→5), albG (2→3), acuA (2→5); sorA decreased (9→7)
- acuA npop=4 (CI: 4–4) is now well below thesis expectation of 14; albG npop=16 (CI: 4–25) is closest
- albG has the most thesis ecotypes represented in mixed ecotypes (11 of 14)
- Thesis ecotypes I14 does not appear in any mixed ecotype for any gene
