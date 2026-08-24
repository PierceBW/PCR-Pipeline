# PCR Amplicon Pipeline — Original Pipeline Documentation

> **Note:** This documents the original (pre-processing) pipeline that produces the input for the current pipeline in `pipeline/`. The original scripts are archived in `archive/` on the lab computer. The current pipeline starts from `ecosim_ready2/` — see the main [README](../README.md) for those steps.

This document describes the steps from raw sequencing data to EcoSim2-ready inputs, including all design decisions and their rationale. Intermediate outputs from each stage are preserved on the lab computer so any decision can be revisited.

---

## Overview

**Goal:** Classify environmental *Bacillus* amplicon reads to three species, combine with isolate genome amplicons from Jocelyn Wang's thesis (known ecotypes), and run EcoSim2 to determine which ecotypes the field sequences fall into.

**Three species:** *B. atrophaeus*, *B. spizizenii*, *B. inaquosorum*

**18 genes** (from `bacillus primers.xlsx` — authoritative source):
| Species | Genes |
|---------|-------|
| spizizenii | acsA_2, alaS, alkH, iolB, rmlD, comB |
| inaquosorum | accA, acuA, albG, thiD, yvqK, sorA |
| atrophaeus | ecfA1, acdA, adc, amj, yndE_2, opuAB |

> **Dropped:** `arnC_1` — not in the Excel primer file, and no BAM reads exist for it in the Galaxy history. Primers remain in git history if needed later.

---

## Stage 1 — Reference Setup

### 1a. 6-species reference genomes

**Script:** `fetch_6ref_genomes_ncbi.py`
**Output:** `reference/refs_6genomes.fasta`, `reference/refs_6genomes_manifest.tsv`
**BLAST DB:** built with `reference/build_blast_db_6refs.sh` → `reference/refs_6genomes_db.*`

6 full genomes from NCBI (2 per species):
- *B. atrophaeus*: NZ_CP007640.1, NZ_CP195110.1
- *B. spizizenii*: CP077772.1, CP030925.1
- *B. inaquosorum*: NZ_CP162598.1, NZ_CP080644.1

Used for species classification (Stage 3).

### 1b. Primer table

**Script:** `parse_primer_htmls.py` (rewritten April 2026 to use Excel as sole source)
**Output:** `notes/primers.tsv`
**Source:** `bacillus primers.xlsx` — columns: Gene, Species, Reason for Picking, Left Primer, Right Primer

**Important primer direction note:**
- HTML files stored the reverse primer as it appears on the **+ strand of the amplicon** (i.e. the revcomp of the actual binding sequence)
- The Excel file stores the reverse primer as the **actual binding sequence** (what you'd order from a company)
- `notes/primers.tsv` stores the Excel version (actual binding sequence)
- `extract_amplicons_by_primers.py` was updated to search `revcomp(rev_primer)` on the + strand of the genome (not `rev_primer` directly)

**accA discrepancy:** The HTML files had completely wrong primers for accA. The Excel primers are correct.

### 1c. Strain table (Jocelyn Wang thesis isolates)

**Script:** `build_strain_table.py`
**Output:** `notes/full_strain_table.tsv`
**Source:** `matching_seq/Supplementary Table 3...xlsx`, sheet `wcloneid`

479 CBP strains with species and ecotype assignments. 3 strains in the table have no genome file (CBP-2516, CBP-3066, CBP-3068).

---

## Stage 2 — Extract Environmental Reads from BAMs

**Script:** `build_final_fastas.py`
**Source:** `PCR-Primer/` (Galaxy history export)
**Output:** `final_fastas2/` — 36 files (18 genes × forward + reverse)

The `PCR-Primer/` Galaxy history contains per-gene per-environment BAMs demuxed with 18 primer pairs:
- `18Primers-1_sorted.bam` = env1
- `18Primers-2_sorted.bam` = env2
- `18Primers-4_sorted.bam` = env4
- `18Primers-5_sorted.bam` = env5

Per gene, reads from all 4 envs are combined. Forward reads (not reverse-strand mapped) go to `{gene}_forward.fasta`; reverse-strand reads are reverse-complemented and go to `{gene}_reverse.fasta`.

**Headers:** `>AmpSeq-MD2-env{N}` (first read per env), `>AmpSeq-MD2-env{N}_v{i}` (subsequent)

**Outgroup placeholder:** Each file starts with `>FN597644.1` copied from the old pipeline as a placeholder. This is replaced with the proper per-gene amplicon in Stage 6.

> **Backup:** `final_fastas/` and `deduped_fastas/` are the old 11-gene pipeline outputs. Kept as backup; not used in the current pipeline.

**Read counts by gene** (total across all envs):

| Gene | Forward | Reverse | Note |
|------|---------|---------|------|
| acuA | 1,150,344 | 1,157,279 | |
| iolB | 975,346 | 979,689 | |
| sorA | 950,611 | 957,183 | |
| rmlD | 577,226 | 580,129 | |
| yvqK | 455,722 | 456,580 | |
| thiD | 365,453 | 367,606 | |
| albG | 362,520 | 364,244 | |
| acsA_2 | 348,085 | 350,359 | |
| alkH | 176,763 | 177,364 | |
| alaS | 2,084 | 2,096 | |
| comB | 437 | 442 | |
| ecfA1 | 62 | 64 | very few |
| opuAB | 62 | 56 | very few |
| amj | 24 | 24 | very few |
| accA | 20 | 12 | very few |
| acdA | 4 | 3 | almost none |
| yndE_2 | 2 | 2 | almost none |
| adc | 0 | 0 | no reads |

---

## Stage 3 — Deduplication

**Script:** `deduplicate_fastas.py`
**Input:** `final_fastas2/`
**Output:** `deduped_fastas2/`, `dedup_counts/`

Exact sequence deduplication (uppercase, gaps stripped). The outgroup (first sequence) is always preserved unchanged and excluded from deduplication.

**Result:** 10.8M total reads → 3.2M unique (70% were duplicates)

Count tables for prof review:
- `dedup_counts/{gene}_{strand}_counts.tsv` — per unique sequence: representative ID, count, sequence
- `dedup_counts/dedup_summary.tsv` — per file summary: total, unique, duplicated, % unique

---

## Stage 4 — Species Classification

**Script:** `classify_6ref_blast.py`
**Input:** `deduped_fastas2/`
**Output:** `classified_for_ecosim2/`
**Reference DB:** `reference/refs_6genomes_db` (2 genomes per species × 3 species = 6 total)

Each sequence is BLASTed against the 6 reference genomes. Ranking is by **bitscore** (continuous score, more discriminating than rounded pident). Assignment rules:
- Best bitscore hit(s) ≥ 93% pident → assigned to that species
- Best hit < 93% pident → dropped (logged in `classified_for_ecosim2/notes/dropped_below93.tsv`)
- **Tied best bitscore across multiple species → assigned to ALL tied species** (sequence appears in each tied species' FASTA)

**Why bitscore not pident:** BLAST's `pident` output is rounded to 3 decimal places. Two hits can show identical pident but have different actual alignment scores. Bitscore is continuous and resolves most apparent ties. Genuine cross-species ties (same bitscore, different species) are assigned to all tied species rather than dropped — the sequence is real data and the tree will show where it sits phylogenetically.

**Why 2 refs per species:** Having 2 reference genomes per species means a read can tie between the two refs of the *same* species — this is not a real tie and is correctly handled (both refs have the same species label, so the species set is size 1 and the sequence is assigned normally).

Tied sequences are logged in `classified_for_ecosim2/notes/tied_unassigned.tsv` with which species tied.

Reverse reads are reverse-complemented before BLAST. Per-species sub-FASTAs are written under `classified_for_ecosim2/{species}/`.

---

## Stage 5 — Isolate Amplicon Extraction (Jocelyn Wang genomes)

**Script:** `extract_amplicons_by_primers.py`
**Input:** `genomes/` (487 CBP isolate genomes), `notes/primers.tsv`, `notes/full_strain_table.tsv`
**Output:** `isolate_amplicons2/`

For each of the 487 CBP isolate genomes, the primer pairs are used to extract the amplicon region (same method as PCR: find fwd primer, find revcomp(rev) downstream, extract between them). Headers are `CBP-XXXX_ecotype` (e.g. `CBP-1646_PE_A1`).

**Result:** 2,842 amplicons across 18 gene/species combos; 14 not-found (contig-boundary edge cases).

**Search strategy (fast → slow):**
1. Exact match on concatenated contigs (null-byte separated, C-level str.find)
2. Hamming sliding window (≤2 mismatches) per contig
3. BLAST fallback (blastn-short) for indel cases

**Output sequences are full amplicons** (~500–600 bp, primer-to-primer). These get trimmed in Stage 7.

---

## Stage 6 — Merge PCR Reads + Isolate Amplicons

**Script:** `merge_for_ecosim.py`
**Input:** `classified_for_ecosim2/` + `isolate_amplicons2/`
**Output:** `ecosim_ready2/`

For each `{species}/{gene}_{strand}.fasta` in the classified output, the matching isolate amplicons (`isolate_amplicons2/{gene}_{species}.fasta`) are appended. Isolates are species-independent — they are added to both forward and reverse files for each gene/species combination.

**Result:** 78 FASTAs, 10,737,704 total sequences (10,732,678 PCR reads + 5,026 isolate sequences). The higher PCR count vs Stage 3 reflects the tie-to-all assignment change (Stage 4): reads that previously tied cross-species and were dropped are now assigned to each tied species, adding sequences (notably ~22k alkH reads recovered for spizizenii).

**Per-species FASTA counts from merge:**

| Species | PCR reads | Isolates | Total |
|---------|-----------|----------|-------|
| atrophaeus | small | 123 | ~1.4k |
| spizizenii | ~35k | 209–211 | ~35.5k |
| inaquosorum | ~10.7M | 141–142 | ~10.7M |

---

## Stage 7 — Prepare for Tree Building (IN PROGRESS)

**Scripts:** `fetch_outgroup_genome.py`, `extract_outgroup_amplicons.py`, `prepare_for_trees.py`
**Input:** `ecosim_ready2/`
**Output:** `tree_inputs/`

### 7a. Sequence length problem

Inspecting `ecosim_ready2/` files revealed three populations:

| Population | Length | Header pattern | Action |
|------------|--------|---------------|--------|
| PCR reads | ~301 bp | `AmpSeq-MD2-env*` | Keep as-is |
| Isolate amplicons | ~505–598 bp | `CBP-XXXX_PE_*` | Trim to 301 bp |
| Trash reads | <200 bp | `AmpSeq-MD2-env*` | Drop |
| Outgroup placeholder | 287 bp | `FN597644.1` | Replace |

The 500+ bp isolate amplicons are full primer-to-primer sequences including both primer binding regions. The PCR reads start at the forward primer position and read ~301 bp into the amplicon. They cover the same gene region, just with different endpoints.

### 7b. Trimming decision (forward files)

**`isolate[:301]`** — take the first 301 bp of the isolate amplicon.

The isolate amplicon structure is: `[fwd_primer (~20bp)][internal region (~480bp)][revcomp(rev_primer) (~21bp)]`. The PCR forward reads also start near the fwd primer position and read ~301 bp. Both populations cover the same gene window.

PCR reads start a few bases (~0–7 bp, variable per read) before the primer binding site — this is normal sequencing overhang and varies between reads. The overlap region is ~275–295 bp of shared sequence, which is what drives tree topology.

**Decision rationale:** A ~6 bp end-effect at the primer boundary is negligible for phylogenetics. Using `[:301]` is simple, consistent, and avoids per-read primer searching. If this is ever revisited, the raw data is in `ecosim_ready2/` unchanged.

### 7c. Trimming decision (reverse files)

**`revcomp(isolate)[:301]`** — reverse-complement the full isolate amplicon, then take first 301 bp.

The isolate was extracted on the + strand. PCR reverse reads were reverse-complemented in `build_final_fastas.py` (reverse-strand BAM reads → revcomp → stored). So both PCR reverse reads and trimmed isolate reverse sequences read 5'→3' from the reverse primer position.

### 7d. Outgroup: FN597644.1

**FN597644.1 = *B. amyloliquefaciens* DSM 7 = ATCC 23350** (confirmed via NCBI lookup). This is the correct outgroup for all three species — prof confirmed *B. amyloliquefaciens* is appropriate.

The current placeholder in `ecosim_ready2/` files is a 287 bp snippet from the old pipeline (wrong). The proper approach:
1. Fetch the complete FN597644.1 genome (~3.9 Mbp) from NCBI → `reference/amyloliquefaciens_FN597644.fasta`
2. For each gene: take 20 representative PCR reads from `deduped_fastas2/`, BLAST each against FN597644.1 (`blastn -word_size 7`), take the **modal hit position** (most queries agree), extract 301 bp from that locus
3. Write per-gene outgroup FASTA to `reference/outgroup_amplicons/{gene}_outgroup.fasta`

**Why BLAST with representative reads, not with primers:** The primers (20 bp) were designed specifically for the three target species and produce hundreds of spurious near-matches in the amyloliquefaciens genome. Using actual ~301 bp amplicon sequences as BLAST queries gives a specific, unique hit at the correct homologous locus. The modal position (most common hit across 20 independently sampled reads) is used to reject spurious outliers.

**Identity to amyloliquefaciens:** Most genes show 74–82% pident — typical for the *B. subtilis* group across species. ecfA1 is higher (~88%) and amj lower (~74%).

**Per-gene outgroup status (validated with 20 queries each):**

| Gene | Outgroup locus | Status |
|------|---------------|--------|
| acuA | FN597644.1:2822423 | ✓ 20/20 consistent |
| albG | FN597644.1:3626480 | ✓ 20/20 consistent |
| iolB | FN597644.1:3843909 | ✓ 20/20 consistent |
| acsA_2 | FN597644.1:2799399 | ✓ modal (19/20) |
| alaS | FN597644.1:2594457 | ✓ modal (19/20) |
| alkH | FN597644.1:3857845 | ✓ modal (19/20) |
| ecfA1 | FN597644.1:151812 | ✓ modal (19/20) |
| rmlD | FN597644.1:3670891 | ✓ modal (19/20) |
| thiD | FN597644.1:1268538 | ✓ modal (19/20) |
| opuAB | FN597644.1:293258 | ✓ modal (15/20) |
| amj | FN597644.1:430723 | ✓ modal (18/20) |
| yvqK | FN597644.1:3212841 | ✓ modal (17/20) |
| accA | FN597644.1:~2774200 | ✓ modal (majority) |
| yndE_2 | FN597644.1:375483 | ✓ 2/2 consistent (sparse gene, no tree) |
| acdA | scattered | ✗ not reliable (4 reads, no tree anyway) |
| **comB** | **15–21 bp partial matches only** | **⚠ gene absent/too diverged — no outgroup** |
| **sorA** | **16–26 bp partial matches only** | **⚠ gene absent/too diverged — no outgroup** |
| adc | no sequences | ✗ 0 reads, no tree |

**comB and sorA:** These two genes have no detectable homolog in FN597644.1 by BLAST (all matches are 15–26 bp fragments, not full-region hits). Both genes have ≥100 sequences in inaquosorum and/or spizizenii and would otherwise make trees. **Action needed:** Ask prof whether to use a different outgroup genome for these two genes, or proceed without an outgroup for those trees.

> **Revisit from:** `check_outgroup_hits.py --genes comB,sorA` to re-examine hit quality. `reference/amyloliquefaciens_FN597644.fasta` preserved. An alternative outgroup genome (e.g. *B. subtilis* 168, already in `refs_6genomes.fasta`) could be tried.

Only **one outgroup per gene** is used — the same FN597644.1 sequence for all three species within a gene.

### 7e. PCR read length enforcement

FastTree and VeryFastTree require **all sequences in a file to be exactly the same length**. PCR reads nominally target 301 bp, but some reads are shorter (200–300 bp) due to read-through or early termination. These shorter reads pass the 200 bp trash filter but would crash tree building.

**Rule applied in `prepare_for_trees.py`:**
- PCR read ≥ 301 bp → kept, trimmed to exactly 301 bp with `seq[:301]`
- PCR read 200–300 bp → **dropped** (too short; would break alignment)
- PCR read < 200 bp → dropped as trash

This means the "valid" count in `filter_log.tsv` only includes 301 bp reads. Reads between 200–300 bp are counted under `trash_dropped`.

### 7f. Minimum sequence filter

Gene/species/strand combos with fewer than 100 sequences (PCR + isolates, excluding outgroup) are not sufficient for EcoSim2 and are excluded. These are logged in `tree_inputs/filter_log.tsv`.

### 7g. No alignment needed

Because all sequences end up at 301 bp covering the same amplicon region, a multiple sequence aligner (MAFFT, MUSCLE) is not required. FastTree can run directly on the unaligned-but-same-length files.

> **If you want to revisit this:** All intermediate outputs are preserved. `ecosim_ready2/` has the un-trimmed sequences. Running MAFFT would produce a gapped MSA where the 301 bp reads and 580 bp isolates would be aligned properly — this is the "correct" bioinformatics approach but adds complexity and requires the aligner to reconcile two very different length populations. The 301 bp trim is a pragmatic choice that produces equivalent results for tree building purposes.

---

## Stage 8 — Identity Filtering & Validation

### 8a. Whole-Genome 95% Filter (Initial Attempt)

**Script:** `filter_by_type_strain.py`
**Input:** `tree_inputs/{species}/{gene}_forward.fasta`
**Output:** `data/09_filtered_ecosim/`

BLASTed all PCR reads against the 6-genome reference database. Kept reads with ≥95% identity to their assigned species' reference genomes. CBP isolates and outgroup retained unconditionally.

**Result:** Removed very few sequences (99%+ passed for most genes). Did not resolve ecotype inflation.

### 8b. EcoSim Validation Experiment

**Scripts:** `test_ecosim_small.py`, `cbp_ecosim_test.py`
**Output:** `data/10_test_ecosim_200/`, `data/12_cbp_only_test/`
**Report:** `notes/ecosim_validation_results.md`

Tested whether EcoSim2's inflated ecotype counts were caused by data quality issues or a fundamental single-locus limitation. Key experiments:

1. **Small-scale test** (200 sequences): Ecotype counts remained inflated (54–248 from 200 sequences)
2. **CBP-only test**: Ran EcoSim on ONLY the known CBP isolates. Results varied by gene — some genes recovered near-thesis ecotype counts (iolB: 5 vs expected 5), others were massively inflated (acdA: 64 vs expected 8)

**Conclusion:** The problem is not data quality. Single-locus amplicon trees cannot reliably reproduce multi-locus ecotype assignments due to gene-tree discordance. See `notes/ecosim_validation_results.md` for full analysis including thesis context.

### 8c. Gene-Specific Identity Filter

**Script:** `filter_by_gene_ref.py`
**Input:** `tree_inputs/{species}/{gene}_forward.fasta`
**Output:** `data/13_gene_filtered/`
**Report:** `notes/gene_specific_filter_results.md`

The whole-genome filter (8a) had a fundamental flaw: it checked identity against the entire genome, not the target gene. This meant off-target amplification products could pass if they matched any conserved region.

The gene-specific filter extracts the actual amplicon region from each reference genome using the primer sequences, then filters reads against only that region.

**Key finding:** Cross-species gene contamination is real and substantial:
- acsA_2 (spi primers) in inaquosorum: **95% of reads off-target**
- comB (spi primers) in inaquosorum: **83% off-target**
- acuA (inaq primers) in spizizenii: **72% off-target**
- Same-species genes are clean (>97% pass)

See `notes/gene_specific_filter_results.md` for full results.

---

## Stage 9 — Tree Building (IN PROGRESS)

**Script:** `build_trees.py`
**Tool:** FastTree 2 (< 10k seqs) or VeryFastTree (≥ 10k seqs)
**Input:** `tree_inputs/{species}/{gene}_{strand}.fasta`
**Output:** `tree_inputs/{species}/{gene}_{strand}.nwk`

Each FASTA gets its own tree. Forward and reverse are run separately (separate EcoSim runs). Model: GTR + CAT nucleotide (`-gtr -nt`). `build_trees.py` runs all files in parallel using `--threads` workers.

**Tool selection logic:**
- FastTree: files with < 10,000 sequences — uses `OMP_NUM_THREADS` env var for threading
- VeryFastTree: files with ≥ 10,000 sequences — uses `-threads N` flag (FastTree does NOT accept `-threads`)

**Result:** 40 input FASTAs, **40 trees built** (complete).

**Trees per species:**

| Species | Gene/strand combos | Notes |
|---------|-------------------|-------|
| atrophaeus | 10 | 5 genes × 2 strands |
| inaquosorum | 22 | 11 genes × 2 strands |
| spizizenii | 8 | 4 genes × 2 strands |

**Outgroup rooting:** FastTree/VeryFastTree produce unrooted trees. Each tree must be rerooted at the `FN597644.1` leaf before EcoSim2. **Pending:** `reroot_trees.py` (in `old/`) needs updating for new paths.

**Install:** `brew install fasttree` or download VeryFastTree binary from GitHub.

---

## Stage 9 — EcoSim2 (PENDING)

**Tool:** EcoSim2 JAR
**Input:** matching `.fasta` + `.nwk` pairs from `tree_inputs/`
**Output:** per-gene per-species EcoSim2 results

EcoSim2 takes the aligned FASTA and the rooted Newick tree and assigns each sequence to an ecotype. The CBP isolate sequences (labeled `CBP-XXXX_ecotype`) will appear in the tree alongside the environmental reads, allowing comparison of amplicon-derived clusters with the thesis ecotype assignments.

---

## Directory Reference

| Directory | Contents | Status |
|-----------|----------|--------|
| `final_fastas2/` | All 18 genes × 2 strands, raw reads from PCR-Primer BAMs | Complete |
| `deduped_fastas2/` | Deduplicated version of final_fastas2/ | Complete |
| `dedup_counts/` | Per-sequence count tables + summary | Complete |
| `classified_for_ecosim2/` | Species-classified reads, per-species subdirs | Complete |
| `isolate_amplicons2/` | Per-gene CBP isolate amplicons (corrected primers) | Complete |
| `ecosim_ready2/` | Merged PCR + isolates, pre-trim | Complete |
| `tree_inputs/` | 40 trimmed FASTAs + 40 NWK trees; filter_log.tsv | Complete |
| `reference/` | 6-ref genomes, BLAST DB, outgroup genome + per-gene amplicons | Complete |
| `genomes/` | 487 CBP isolate genomes | Complete |
| `PCR-Primer/` | Galaxy history export with per-gene per-env BAMs | Source data |
| `bam_files/` | env1/2/4/5 BAM files | Source data |
| `final_fastas/` | OLD 11-gene pipeline (backup) | Backup |
| `deduped_fastas/` | OLD 11-gene deduped (backup) | Backup |
| `old/` | Superseded scripts and analysis outputs | Archive |
| `data/09_filtered_ecosim/` | Whole-genome 95% filtered FASTAs | Complete |
| `data/10_test_ecosim_200/` | 200-seq EcoSim test results | Complete |
| `data/12_cbp_only_test/` | CBP-only EcoSim validation experiment | Complete |
| `data/13_gene_filtered/` | Gene-specific filtered FASTAs + gene refs | Complete |
| `notes/` | primers.tsv, full_strain_table.tsv, procedure docs, result reports | Reference |
| `bacillus-primer/` | HTML primer reference files | Reference |
| `matching_seq/` | Thesis supplementary table + assembly matching | Reference |

---

## Script Reference

| Script | Purpose |
|--------|---------|
| `parse_primer_htmls.py` | Build `notes/primers.tsv` from `bacillus primers.xlsx` |
| `build_strain_table.py` | Build `notes/full_strain_table.tsv` from thesis Excel |
| `collect_genomes.py` | Extract CBP genome FASTAs from Galaxy history exports |
| `fetch_6ref_genomes_ncbi.py` | Download 6 NCBI reference genomes |
| `build_final_fastas.py` | Extract per-gene FASTAs from PCR-Primer BAMs → `final_fastas2/` |
| `deduplicate_fastas.py` | Deduplicate FASTAs with count tables → `deduped_fastas2/` |
| `classify_6ref_blast.py` | BLAST-classify reads to species → `classified_for_ecosim2/` |
| `extract_amplicons_by_primers.py` | Extract per-gene amplicons from CBP isolate genomes |
| `merge_for_ecosim.py` | Merge classified PCR reads + isolate amplicons → `ecosim_ready2/` |
| `fetch_outgroup_genome.py` | Download FN597644.1 (*B. amyloliquefaciens*) genome |
| `extract_outgroup_amplicons.py` | Extract per-gene amplicons from outgroup genome |
| `prepare_for_trees.py` | Filter ≥100 seqs, trim to 301 bp, replace outgroup → `tree_inputs/` |
| `build_trees.py` | Run FastTree/VeryFastTree on all `tree_inputs/` FASTAs in parallel |
| `check_outgroup_hits.py` | Diagnostic: test BLAST hit consistency for outgroup extraction |
| `filter_by_type_strain.py` | Whole-genome 95% identity filter (Stage 8a) |
| `test_ecosim_small.py` | Small-scale EcoSim test runner (Stage 8b) |
| `cbp_ecosim_test.py` | CBP-only EcoSim validation experiment (Stage 8b) |
| `filter_by_gene_ref.py` | Gene-specific identity filter using amplicon references (Stage 8c) |

---

## Key Design Decisions Log

| Decision | Choice | Rationale | Revisit from |
|----------|--------|-----------|-------------|
| Primer source | Excel (`bacillus primers.xlsx`) | HTML had wrong accA primers; Excel has all 18 genes | `notes/primers.tsv` |
| Rev primer direction | Store as actual binding seq; search `revcomp(rev)` on genome | Excel convention; HTML stored + strand representation | `extract_amplicons_by_primers.py` lines 210, 221 |
| arnC_1 | Dropped | No BAM reads; not in Excel | `notes/primers.tsv` git history |
| Deduplication | Exact match (uppercase, gaps stripped) | Remove PCR duplicates before classification | `deduped_fastas2/`, `dedup_counts/` |
| Species classification cutoff | 93% identity | Below this is likely cross-species noise | `classify_6ref_blast.py` CUTOFF_PCT |
| Outgroup | FN597644.1 (*B. amyloliquefaciens* DSM 7) | Confirmed by prof; same for all genes and species | `reference/amyloliquefaciens_FN597644.fasta` |
| Outgroup extraction method | BLAST with representative amplicons (modal position), not primers | Primers (20 bp) give hundreds of spurious hits in amyloliquefaciens; 301 bp amplicon queries give unique locus-specific hits | `check_outgroup_hits.py` output, `deduped_fastas2/` |
| comB/sorA outgroup | No outgroup found in FN597644.1 | All BLAST hits are 15–26 bp fragments; gene absent or too diverged | `check_outgroup_hits.py --genes comB,sorA` — ask prof about alternative outgroup |
| Isolate trim length | `[:301]` (forward); `revcomp(seq)[:301]` (reverse) | Matches PCR read length; same gene region; ~6 bp end offset is phylogenetically negligible | `ecosim_ready2/` has untrimmed originals |
| Minimum sequences | 100 per gene/species/strand combo | EcoSim2 needs enough sequences for meaningful ecotype clustering | `tree_inputs/filter_log.tsv` |
| Alignment | None (not required) | All sequences are 301 bp covering same amplicon region | `ecosim_ready2/` + MAFFT available if needed |
| Forward vs reverse separation | Kept separate through all stages | Forward and reverse reads may have different error profiles; separate EcoSim runs | Can be merged if needed |
| PCR read length enforcement | Drop PCR reads < 301 bp (even if ≥ 200 bp); truncate ≥ 302 bp to 301 | FastTree/VeryFastTree require all seqs exactly same length; 200 bp trash cutoff was not sufficient | `prepare_for_trees.py` min_read_len + trim_len args |
| Classification scoring metric | BLAST bitscore (not pident) | `pident` is rounded to 3 decimal places; bitscore is continuous and resolves most apparent ties | `classify_6ref_blast.py` BLAST outfmt column 5 |
| Cross-species ties | Assign to ALL tied species | A genuinely ambiguous read is real data — the tree will show where it clusters; dropping it loses information | `classified_for_ecosim2/notes/tied_unassigned.tsv` |
