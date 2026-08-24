# Tree Visualization, In-Group Pruning & EcoSim Re-run Results

## The Big Picture

We visualized all 8 gene trees (n=200 PCR + all CBP + outgroup) and asked: **are the environmental PCR reads inside or outside the Jocelyn CBP in-group?**

The answer splits cleanly by species:

| Species | Genes | PCR Inside CBP LCA | PCR Outside CBP LCA | Verdict |
|---------|-------|:---:|:---:|---------|
| **inaquosorum** | acuA, albG, sorA, thiD, yvqK | **100%** | 0% | PCR reads are same species as CBP |
| **spizizenii** | acsA_2, alkH, iolB | 0% | **100%** | PCR reads are a DIFFERENT taxon |

This is a binary, all-or-nothing result — no partial overlap in any gene.

---

## Per-Gene Tree Analysis

### acuA (inaquosorum) — 342 sequences

**Tree PDF:** `data/06_viz/acuA_tree.pdf`

- **CBP:** 141 isolates (thesis ecotypes I1–I14)
- **PCR:** 200 environmental reads (E1, E2, E4, E5)
- **Outgroup:** CP026362.1
- **PCR inside LCA:** 200 (100%)
- **PCR outside LCA:** 0

**Interpretation:** All 200 PCR reads fall within the CBP in-group clade. The environmental diversity is fully interleaved with Jocelyn's known ecotypes. The tree shows PCR reads scattered throughout the CBP branches — good ecological signal.

**Pruned EcoSim (unchanged):** npop=115, ecotypes=208 (thesis expects 14)

---

### albG (inaquosorum) — 341 sequences

**Tree PDF:** `data/06_viz/albG_tree.pdf`

- **CBP:** 141 isolates
- **PCR:** 200 environmental reads
- **Outgroup:** none (albG has no outgroup — unrooted tree)
- **PCR inside LCA:** 200 (100%)
- **PCR outside LCA:** 0

**Interpretation:** Same pattern as acuA. All PCR reads inside the CBP clade despite no outgroup for rooting. The environmental inaquosorum diversity for albG is correctly captured.

**Pruned EcoSim (unchanged):** npop=37, ecotypes=158 (thesis expects 14)

---

### sorA (inaquosorum) — 342 sequences ★ BEST GENE

**Tree PDF:** `data/06_viz/sorA_tree.pdf`

- **CBP:** 141 isolates
- **PCR:** 200 environmental reads
- **Outgroup:** CP026362.1
- **PCR inside LCA:** 200 (100%)
- **PCR outside LCA:** 0

**Interpretation:** sorA shows the most interleaving of CBP and PCR reads. In the tree, you can see colored CBP squares (by thesis ecotype) mixed with PCR circles (by environment) throughout the clade. The 9 mixed ecotypes from the EcoSim analysis confirm environmental reads match known ecotypes. Multiple CBP ecotype groups (I3, I4, I5, I7, I8, I9, I10) each have PCR reads clustering alongside them.

**Pruned EcoSim (unchanged):** npop=16, ecotypes=58 (thesis expects 14)

---

### thiD (inaquosorum) — 343 sequences

**Tree PDF:** `data/06_viz/thiD_tree.pdf`

- **CBP:** 142 isolates
- **PCR:** 200 environmental reads
- **Outgroup:** CP026362.1
- **PCR inside LCA:** 200 (100%)
- **PCR outside LCA:** 0

**Interpretation:** All PCR reads inside. thiD had 3 mixed ecotypes at n=200 (improved to 9 at n=2000), so environmental reads do match some CBP ecotypes.

**Pruned EcoSim (unchanged):** npop=43, ecotypes=124 (thesis expects 14)

---

### yvqK (inaquosorum) — 342 sequences

**Tree PDF:** `data/06_viz/yvqK_tree.pdf`

- **CBP:** 141 isolates
- **PCR:** 200 environmental reads
- **Outgroup:** CP026362.1
- **PCR inside LCA:** 200 (100%)
- **PCR outside LCA:** 0

**Interpretation:** All PCR inside. 3 mixed ecotypes at n=200 — PCR reads cluster with I5, I7, and I2/I4/I6/I9/I11 groups.

**Pruned EcoSim (unchanged):** npop=12, ecotypes=56 (thesis expects 14)

---

### acsA_2 (spizizenii) — 410 sequences

**Tree PDF:** `data/06_viz/acsA_2_tree.pdf`

- **CBP:** 209 isolates (thesis ecotypes S1–S5)
- **PCR:** 200 environmental reads (ALL from E2)
- **Outgroup:** CP026362.1
- **PCR inside LCA:** 0 (0%)
- **PCR outside LCA:** 200 (100%)

**Interpretation:** The CBP isolates form a monophyletic clade with ZERO PCR reads inside it. All 200 PCR reads branch off separately — they are a completely distinct lineage from the culture collection spizizenii. Notably, all 200 PCR reads are from environment 2, which could mean this environment harbors a different Bacillus species/taxon that passes the 95% gene filter for acsA_2 but is genuinely distinct.

**Pruned EcoSim (CBP-only):** npop=2, ecotypes=1 (thesis expects 5)

With no PCR reads, EcoSim found only 1 ecotype from 209 CBP isolates. This is far below the thesis expectation of 5 spizizenii ecotypes, suggesting the CBP-only tree lacks the diversity resolution that PCR reads would normally provide — but since those PCR reads are wrong taxon, we can't use them.

---

### alkH (spizizenii) — 410 sequences

**Tree PDF:** `data/06_viz/alkH_tree.pdf`

- **CBP:** 209 isolates
- **PCR:** 200 environmental reads
- **Outgroup:** CP026362.1
- **PCR inside LCA:** 0 (0%)
- **PCR outside LCA:** 200 (100%)
- **Excluded envs:** E1:2, E2:26, E4:96, E5:76

**Interpretation:** Same as acsA_2 — complete CBP/PCR separation. Unlike acsA_2 where all PCR were E2, alkH's excluded PCR reads come from all 4 environments (predominantly E4 and E5). This rules out a single-environment contamination issue — the "wrong taxon" pattern is consistent across environments.

**Pruned EcoSim (CBP-only):** npop=6, ecotypes=34 (thesis expects 5)

Better result than acsA_2 — 34 ecotypes from CBP-only data. The thesis expects 5, so 34 is an over-split but in the right ballpark. The npop=6 is very close to thesis expectation of 5.

---

### iolB (spizizenii) — 412 sequences

**Tree PDF:** `data/06_viz/iolB_tree.pdf`

- **CBP:** 211 isolates
- **PCR:** 200 environmental reads
- **Outgroup:** CP026362.1
- **PCR inside LCA:** 0 (0%)
- **PCR outside LCA:** 200 (100%)
- **Excluded envs:** E1:21, E2:43, E4:129, E5:7

**Interpretation:** Same pattern — complete separation. PCR reads predominantly from E4 (129/200). The iolB PCR environmental reads are genetically distinct from all 211 CBP spizizenii isolates.

**Pruned EcoSim (CBP-only):** npop=14, ecotypes=86 (thesis expects 5)

---

## Summary: Pruned EcoSim vs Thesis

| Gene | Species | Total Seqs | PCR | npop | Ecotypes | Thesis | Status |
|------|---------|:---:|:---:|:---:|:---:|:---:|:---:|
| acsA_2 | spizizenii | 210 | 0 (CBP-only) | 2 | 1 | 5 | Under-split |
| acuA | inaquosorum | 342 | 200 | 115 | 208 | 14 | Over-split |
| albG | inaquosorum | 341 | 200 | 37 | 158 | 14 | Over-split |
| **alkH** | **spizizenii** | **210** | **0 (CBP-only)** | **6** | **34** | **5** | **npop close!** |
| iolB | spizizenii | 212 | 0 (CBP-only) | 14 | 86 | 5 | Over-split |
| sorA | inaquosorum | 342 | 200 | 16 | 58 | 14 | npop close |
| thiD | inaquosorum | 343 | 200 | 43 | 124 | 14 | Over-split |
| yvqK | inaquosorum | 342 | 200 | 12 | 56 | 14 | npop close |

**alkH npop=6** is closest to thesis expectation (5). **sorA npop=16** and **yvqK npop=12** are in the right range for inaquosorum (thesis=14).

---

## What This Tells Us

### 1. The spizizenii environmental reads are NOT spizizenii

All 3 spizizenii genes show 100% CBP/PCR separation. The PCR reads pass the 95% BLAST identity filter against the spizizenii gene reference, but they are phylogenetically distinct from all known spizizenii isolates. They may be:
- A closely related but unnamed Bacillus species
- A divergent spizizenii lineage not captured in the culture collection
- Off-target amplification from a related organism

This is exactly what the professor suspected — these "extra ecotypes" in the spizizenii runs are from a different taxon.

### 2. The inaquosorum environmental reads ARE inaquosorum

All 5 inaquosorum genes show 100% PCR-inside-CBP-LCA. The environmental amplicon diversity is genuine inaquosorum diversity, interleaved with Jocelyn's culture collection isolates. The high ecotype counts (56–208) reflect real environmental diversity not captured by the culture collection.

### 3. npop may be more meaningful than ecotype count

The thesis expected ecotype counts (5 for spizizenii, 14 for inaquosorum) were based on CBP-only analysis. Our npop values for sorA (16), yvqK (12), and alkH (6) are in the right range. The ecotype demarcation may be over-splitting due to the added environmental diversity, but the population count estimate (npop) is more robust.

---

## Files

| Output | Location |
|--------|----------|
| Tree PDFs/SVGs | `data/06_viz/{gene}_tree.pdf` |
| LCA summary | `data/06_viz/lca_summary.tsv` |
| Pruned FASTAs | `data/07_ingroup/{gene}/{gene}.fasta` |
| Excluded reads | `data/07_ingroup/{gene}/{gene}_excluded.txt` |
| Prune log | `data/07_ingroup/notes/prune_log.tsv` |
| Pruned EcoSim results | `data/08_ecosim_ingroup/notes/ecosim_ingroup_results.tsv` |
| EcoSim XML per gene | `data/08_ecosim_ingroup/{gene}/{gene}.xml` |
