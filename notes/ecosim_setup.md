# EcoSim Setup & Running Guide

## Location

EcoSim (Ecotype Simulation) JAR is at:
```
/Users/frederickcohan/Desktop/Github/ecosim/ecosim.jar
```
This is a sibling repo to PCR-Pipeline: `~/Desktop/Github/ecosim/`

Version: **Ecotype Simulation 2.1.7**

## Dependencies

- **Java**: `/usr/local/opt/openjdk/bin/java` (OpenJDK 25.0.2 via Homebrew)
- **FastTree**: `/usr/local/bin/fasttree` (for tree building)
- **Biopython**: `Bio.Phylo` (for rerooting trees, in `./venv/`)

## How EcoSim Works

EcoSim takes:
1. A **FASTA file** — aligned sequences (all same length), outgroup listed first
2. A **Newick tree** — same leaf names as FASTA headers

It outputs an **XML file** with demarcation results including ecotype assignments.

## Running EcoSim Manually

```bash
/usr/local/opt/openjdk/bin/java -Xmx4G -jar /Users/frederickcohan/Desktop/Github/ecosim/ecosim.jar \
    -s=/absolute/path/to/sequences.fasta \
    -p=/absolute/path/to/tree.nwk \
    -o=output.xml \
    -n -d -t=8
```

### Flags
| Flag | Long form | Meaning |
|------|-----------|---------|
| `-s` | `--sequences` | Input FASTA file (absolute path) |
| `-p` | `--phylogeny` | Input Newick tree (absolute path) |
| `-o` | `--output` | Output XML filename (written to cwd) |
| `-n` | `--nogui` | No GUI, implies `--runall` (runs full pipeline including demarcation) |
| `-r` | `--runall` | Run everything including demarcation |
| `-d` | `--debug` | Verbose output |
| `-t` | `--threads` | Number of threads |

### Important Notes
- EcoSim writes output **relative to its own directory** (cwd), not the path in `-o`
- Always run with `cwd` set to the ecosim directory
- Use **absolute paths** for `-s` and `-p`
- After run, move the output XML from ecosim dir to your desired location

## Running via Pipeline Script

```bash
# Small subset (200 PCR reads + all CBP + outgroup per gene)
./venv/bin/python pipeline/step5_ecosim.py --n-pcr 200

# Larger subset
./venv/bin/python pipeline/step5_ecosim.py --n-pcr 2000

# Tree building only (no EcoSim)
./venv/bin/python pipeline/step5_ecosim.py --n-pcr 200 --no-ecosim

# Specific genes only
./venv/bin/python pipeline/step5_ecosim.py --n-pcr 200 --genes acuA,iolB
```

### What the Script Does
1. Reads trimmed uniform-length FASTAs from `data/04_trimmed/`
2. Separates outgroup, CBP isolates, and PCR reads
3. Subsamples PCR reads to `--n-pcr` (keeps all CBP + outgroup)
4. Writes subsampled FASTA with outgroup first
5. Builds tree with `fasttree -gtr -nt`
6. Reroots tree at outgroup (FN597644.1 or CP026362.1)
7. Runs EcoSim
8. Parses ecotype count from EcoSim stdout
9. Logs results to TSV

### Output Structure
```
data/05_ecosim/
├── 200/                          # n_pcr=200 run
│   ├── acuA/
│   │   ├── acuA.fasta           # Subsampled sequences
│   │   ├── acuA_unrooted.nwk   # Raw FastTree output
│   │   ├── acuA.nwk            # Rerooted tree
│   │   ├── acuA.xml            # EcoSim results
│   │   └── acuA_ecosim_log.txt # Full EcoSim stdout
│   └── notes/
│       └── ecosim_results.tsv   # Summary table
├── 2000/                         # n_pcr=2000 run
│   └── ...
```

## Thesis Expected Ecotype Counts (Jocelyn Wang)

| Species | Expected Ecotypes |
|---------|:-:|
| atrophaeus | 8 |
| inaquosorum | 14 |
| spizizenii | 5 |

## Usable Genes (8)

| Gene | Species | PCR Reads | CBP | Outgroup |
|------|---------|----------:|----:|:--------:|
| acuA | inaquosorum | 1,120,479 | 141 | yes |
| sorA | inaquosorum | 934,890 | 141 | yes |
| iolB | spizizenii | 949,671 | 211 | yes |
| yvqK | inaquosorum | 443,112 | 141 | yes |
| albG | inaquosorum | 355,719 | 141 | no |
| thiD | inaquosorum | 353,533 | 142 | yes |
| acsA_2 | spizizenii | 331,833 | 209 | yes |
| alkH | spizizenii | 105,936 | 209 | yes |

Note: albG has no outgroup, so tree will be unrooted (midpoint-rooted by FastTree).
