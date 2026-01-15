# TFTag Pipeline User Manual

**TFTag** is a Python pipeline for genome-wide or targeted design of CRISPR–Cas9 gRNAs and homology arms for **N- and C-terminal gene tagging**.  
It produces a **queryable SQLite database** (and Parquet snapshot) that users can search for their gene of interest and select optimal tagging guides.

---

## Contents

1. Overview  
2. Features  
3. Installation  
4. Inputs  
5. Running the pipeline  
6. Pipeline stages  
7. Output database and schema  
8. Querying the SQLite database  
9. Guide selection recommendations  
10. Performance considerations  
11. Troubleshooting  
12. Reproducibility and best practices  

---

## 1. Overview

TFTag automates the full design workflow required for endogenous tagging:

- Identification of start/stop codons from genome annotations
- Discovery of nearby CRISPR gRNAs
- Designability checks for homology arms and primers
- Silent-edit design to prevent re-cutting
- On-target efficiency scoring (RS3)
- Off-target enumeration (Cas-OFFinder)
- Off-target risk scoring (CCLMoff)
- Validation primer design
- Persistent storage in SQLite + Parquet formats

The output is intended to serve both as:
- a **design table** for experimental use, and
- a **reusable resource** that can be queried programmatically.

---

## 2. Features

- Supports **genome-wide** or **gene-list–restricted** design
- Handles **N-terminal (start codon)** and **C-terminal (stop codon)** tagging
- Produces **240 bp homology arms** in gene orientation
- Applies **synonymous edits** to disrupt PAM or protospacer
- Optional **RS3 efficiency scoring**
- Optional **Cas-OFFinder off-target enumeration**
- Optional **CCLMoff off-target risk scoring**
- Designs **PCR validation primers**
- Writes results to:
  - SQLite database (indexed, queryable)
  - Parquet file (analysis-friendly)

---

## 3. Installation

### 3.1 Clone the repository

```bash
git clone https://github.com/YOURNAME/tftag-pipeline.git
cd tftag-pipeline
```

### 3.2 Create an environment (recommended)

```bash
conda create -n tftag python=3.10
conda activate tftag
```

### 3.3 Install TFTag
```bash
pip install -e .
```

Verify installation:

```bash
tftag --help
```

## 4. Inputs

### 4.1 Genome annotation (GTF/GFF)

- Provided via `--gtf`  
- Must contain `start_codon` and `stop_codon` features  
- A `gffutils` database will be created automatically (or reused)  

### 4.2 Reference genome (FASTA)

- Provided via `--fasta`  
- Chromosome/contig names must match the annotation  

### 4.3 Gene selection (`--genes`)

Optional. If omitted, **all genes** in the GTF database are processed.

Accepted formats:
- Path to a file with one gene ID per line  
- Comma-separated list of gene IDs  

Examples:

```bash
--genes genes.txt
--genes "FBgn0000008,FBgn0000017"
Gene identifiers must match those in the annotation database.
```

---

## 5. Running the pipeline

### 5.1 Minimal run (no off-target scoring)

```bash
tftag \
  --gtf genome.gtf \
  --fasta genome.fa \
  --outdb out/tftag_guides.sqlite
```

### 5.2 Targeted gene set

```bash
tftag \
  --gtf genome.gtf \
  --fasta genome.fa \
  --genes genes.txt \
  --outdb out/tftag_guides.sqlite

5.3 Enable RS3 efficiency scoring

RS3 runs automatically if the Python package is installed.

```bash
tftag \
  --gtf genome.gtf \
  --fasta genome.fa \
  --tracr Hsu2013
```

### 5.4 Enable Cas-OFFinder off-target enumeration

```bash
tftag \
  --gtf genome.gtf \
  --fasta genome.fa \
  --specificity \
  --cas cas-offinder \
  --device C
```

### 5.5 Enable CCLMoff off-target risk scoring

```bash
tftag \
  --gtf genome.gtf \
  --fasta genome.fa \
  --specificity \
  --cclmoff_cmd "python run_cclmoff.py --pairs {pairs} --out {output}"
```

---

## 6. Pipeline stages

### Stage A: Annotation parsing
- Builds or loads a `gffutils` database  
- Extracts start/stop codons  
- Assigns N/C tag identifiers (N1, C1, …)  

### Stage B: PAM scanning
- Scans ±window around each codon  
- Identifies NGG PAMs on both strands  
- Computes cut site and distance to codon  

### Stage C: Designability prefilter
Rejects guides if:
- Homology arms extend beyond contig boundaries  
- Primer windows cannot be accommodated  

### Stage D: Efficiency scoring (optional)
- Computes RS3 score using 30-mer context  
- Adds warnings if unavailable  

### Stage E: Off-target analysis (optional)
- Runs Cas-OFFinder  
- Summarizes mismatch counts  

### Stage F: CCLMoff scoring (optional)
- Builds (on, off) sequence pairs  
- Runs external predictor  
- Aggregates risk per spacer  

### Stage G: Homology arm construction
- Builds HAL and HAR (240 bp each)  
- Sequences are returned in **gene orientation**  

### Stage H: Silent-edit design
- Determines which arm must be edited  
- Applies synonymous mutations  
- Records mutation locations  

### Stage I: Validation primer design
- Designs upstream and downstream PCR primers  
- Uses relaxation rounds if strict criteria fail  

### Stage J: Persistence
- Appends results to SQLite database  
- Writes Parquet snapshot  

---

## 7. Output database

### 7.1 Files

- `out/tftag_guides.sqlite` – primary queryable database  
- `out/tftag_guides.parquet` – analysis snapshot  

### 7.2 Key columns

| Category | Examples |
|----------|----------|
| Gene | `FB-id`, `name`, `feature`, `tag` |
| Guide | `spacer`, `gRNA_start`, `gRNA_end`, `cut_distance` |
| Arms | `HAL_seq_gene`, `HAR_seq_gene`, `_mut` |
| Scoring | `rs3_score`, `n_hits`, `cclmoff_max` |
| QC | `designable`, `rejected`, `warnings` |
| Primers | `PCR1_forward_primer_seq`, `PCR2_reverse_primer_seq` |

---

## 8. Querying the SQLite database

### Example: list guides for a gene

```sql
SELECT
  name, tag, feature,
  spacer, cut_distance,
  rs3_score, n_hits, cclmoff_max
FROM guides
WHERE name = 'pax2'
ORDER BY cut_distance ASC, rs3_score DESC;
Example: best N-terminal guide per gene
SELECT *
FROM guides
WHERE feature = 'start_codon'
  AND designable = 1
  AND rejected = 0
ORDER BY name, cut_distance, rs3_score DESC;
```

---

## 9. Guide selection recommendations

A practical ranking heuristic:

1. `designable = TRUE`  
2. `rejected = FALSE`  
3. Minimal `cut_distance`  
4. High `rs3_score`  
5. Low `n_hits`  
6. Low `cclmoff_max`  
7. Few or no warnings  

---

## 10. Performance considerations

- Whole-genome FASTA is loaded into memory  
- Suitable for Drosophila-sized genomes  
- For larger genomes, future versions may use indexed FASTA access  
- Expensive steps:
  - Cas-OFFinder  
  - CCLMoff  
- Consider caching spacer-level results for repeated runs  

---

## 11. Troubleshooting

### No guides found
- No NGG sites in scan window  
- Chromosome name mismatch between GTF and FASTA  
- Missing start/stop codon annotations  

### Many guides marked `designable = FALSE`
- Genes near contig ends  
- Primer window constraints too strict  

### RS3 scores missing
- RS3 package not installed  
- 30-mer context truncated near contig end  

### Cas-OFFinder errors
- Incorrect binary path  
- Invalid device specification  

Check the `warnings` column for details.

## 12. Reproducibility and best practices

For public database releases, record:

- Genome build and annotation version  
- Pipeline version / git commit  
- RS3 tracrRNA  
- Cas-OFFinder mismatch settings  
- CCLMoff model/version  

This information should be stored in a metadata table or release notes.

---

## Citation

If you use TFTag, please cite the associated publication (when available) and/or the GitHub repository.