# Environment Metadata Investigation

## Answer

The "E" labels (E1, E2, E4, E5) come from **plate well positions** (row E on Plate 1) in the AmpSeq submission form. They map to Death Valley soil samples D30–D34. All four are from the **7000-ft Death Valley elevation station** (~7030 ft GPS, ~36°14'38"N, 117°4'26"W).

| Label | Well | Sample | Elevation (ft) | Soil Type | Replicate | Amplicon Data? |
|:-----:|:----:|:------:|:--------------:|:---------:|:---------:|:--------------:|
| **E1** | Plate 1 - E1 | **D30** | 7030 | **Bulk** | B5 | Yes |
| **E2** | Plate 1 - E2 | **D31** | 7030 | **Rhizosphere** | R1 | Yes |
| E3 | Plate 1 - E3 | D32 | 7030 | Rhizosphere | R2 | **No** |
| **E4** | Plate 1 - E4 | **D33** | 7030 | **Rhizosphere** | R3 | Yes |
| **E5** | Plate 1 - E5 | **D34** | 7030 | **Rhizosphere** | R4 | Yes |

**Key facts:**
- E1 (D30) is the **only bulk soil** sample; E2/E4/E5 are all **rhizosphere** replicates
- E3 (D32, rhizosphere R2) produced no usable amplicon data
- All samples are from the same location and elevation (~7030 ft, Death Valley)
- 64 total samples were submitted (S1–S18 Spring Mountains + D1–D45 Death Valley); only these 4 yielded PCR amplicons

---

## Source Data

### AmpSeq Submission Form (`ampseq_files/AMPSEQ-General Sample Submission Form TemplateNov9.xlsx`)

- **Customer**: Cohan lab (Wesleyan University)
- **Contact**: folabemiwo@wesleyan.edu
- **Project**: Death Valley Project
- **Submission date**: 2024-11-09
- **Total samples**: 64 (18 Spring Mountains S1–S18, 45 Death Valley D1–D45, plus S11A)
- **Plate layout**: Single plate, wells A1–F4 filled

The "E" in E1/E2/E4/E5 refers to **row E of the 96-well plate**, not "Environment":

| Well | # | Sample |
|------|---|--------|
| Plate 1 - E1 | 49 | D30 |
| Plate 1 - E2 | 50 | D31 |
| Plate 1 - E3 | 51 | D32 |
| Plate 1 - E4 | 52 | D33 |
| Plate 1 - E5 | 53 | D34 |
| Plate 1 - E6 | 54 | D35 |
| ... | ... | ... |

### Death Valley Collection Sheet (`ampseq_files/Death Valley collection-August2024.xlsx`)

Full metadata for the 4 samples with amplicon data:

| Sample | Elev Station | Lat | Lon | GPS Elev (ft) | Soil | Replicate |
|--------|:-----------:|-----|-----|:-------------:|:----:|:---------:|
| D30 | 7000 | 36°14'38"N | 117°4'26"W | 7030 | Bulk | B5 |
| D31 | 7000 | 36°14'39"N | 117°4'27"W | 7030 | Rhizosphere | R1 |
| D32 | 7000 | 36°14'39"N | 117°4'27"W | 7030 | Rhizosphere | R2 |
| D33 | 7000 | 36°14'39"N | 117°4'27"W | 7030 | Rhizosphere | R3 |
| D34 | 7000 | 36°14'38"N | 117°4'26"W | 7030 | Rhizosphere | R4 |

### Full sampling design (64 samples, 2 mountain ranges, 5 elevations)

**Spring Mountains (S1–S18):**

| Elevation | Bulk Samples | Rhizosphere Samples |
|:---------:|:------------:|:-------------------:|
| 3000 ft (~3390 ft GPS) | S1, S2, S3 | S4, S5, S6 |
| 5000 ft (~5390 ft GPS) | S7, S8, S9 | S10, S11, S11A, S12 |
| 7000 ft (~7170 ft GPS) | S13, S14, S15 | S16, S17, S18 |

**Death Valley (D1–D45):**

| Elevation | Bulk Samples | Rhizosphere Samples |
|:---------:|:------------:|:-------------------:|
| 0 ft (sea level) | D1, D2, D3, D4, D5 | — |
| 3000 ft (~2560 ft GPS) | D6, D7, D8, D9, D10 | D11, D12, D13, D14, D15 |
| 5000 ft (~4500 ft GPS) | D16, D17, D18, D19, D20 | D21, D22, D23, D24, D25 |
| 7000 ft (~7030 ft GPS) | D26, D27, D28, D29, D30 | D31, D32, D33, D34, D35 |
| 10000 ft (~9590 ft GPS) | D36, D37, D38, D39, D40 | D41, D42, D43, D44, D45 |

Of these 64 samples, **only D30–D34 (wells E1–E5) produced PCR amplicon data**, and D32 (E3) yielded nothing usable. This means all our PCR data comes from a single elevation (7000 ft) at a single site in Death Valley.

---

## Implications for Analysis

1. **No elevation or geographic comparison is possible** — all 4 samples are from the same ~7030 ft site
2. **Bulk vs rhizosphere is the only ecological contrast**: E1 (D30) = bulk soil vs E2/E4/E5 (D31/D33/D34) = rhizosphere
3. **E2, E4, E5 are biological replicates** (rhizosphere R1, R3, R4 from the same location)
4. **E3 (D32 = rhizosphere R2) failed** — no amplicons recovered

---

## BAM File Mapping

| BAM File | Well | Sample | Soil | Replicate |
|----------|:----:|:------:|:----:|:---------:|
| `18Primers-1_sorted.bam` | E1 | D30 | Bulk | B5 |
| `18Primers-2_sorted.bam` | E2 | D31 | Rhizosphere | R1 |
| `18Primers-4_sorted.bam` | E4 | D33 | Rhizosphere | R3 |
| `18Primers-5_sorted.bam` | E5 | D34 | Rhizosphere | R4 |

---

## Sequencing Platform

- **Platform**: AVITI (Element Biosciences)
- **Run ID**: `AVITI-Run1-012225`
- **Flow cell**: `2421528783`
- **Read format**: Paired-end, 301 bp each, merged to ~520–570 bp amplicons

---

## Read Count Distribution by Sample

### Raw read pairs (before filtering, from step 0 merge):

| Gene | E1/D30 (Bulk) | E2/D31 (Rhizo) | E4/D33 (Rhizo) | E5/D34 (Rhizo) | Total |
|------|-------:|-------:|-------:|-------:|--------:|
| sorA | 187,409 | 178,454 | 361,018 | 230,668 | 957,549 |
| yvqK | 61,826 | 98,498 | 171,065 | 125,330 | 456,719 |
| thiD | 48,069 | 145,970 | 76,784 | 96,995 | 367,818 |
| albG | 35,838 | 72,842 | 157,291 | 98,399 | 364,370 |
| acuA | 55,678 | 350,572 | 656,140 | 95,587 | 1,157,977 |

### After subsampling to 200 PCR reads per gene (used in EcoSim):

| Gene | E1/D30 (Bulk) | E2/D31 (Rhizo) | E4/D33 (Rhizo) | E5/D34 (Rhizo) |
|------|---:|---:|---:|---:|
| sorA | 51 | 41 | 60 | 48 |
| yvqK | 36 | 43 | 74 | 47 |
| thiD | 32 | 69 | 42 | 57 |
| albG | 25 | 43 | 80 | 52 |
| acuA | 17 | 70 | 98 | 15 |

---

*Updated 2026-06-15 with AmpSeq submission form and Death Valley collection metadata.*
