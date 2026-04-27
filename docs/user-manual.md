# TFTag Pipeline User Manual

## 1. Overview

**TFTag** is a Python pipeline for genome-wide or targeted design of CRISPR–Cas9 gRNAs and homology arms for **N- and C-terminal gene tagging**.

It produces a **queryable SQLite database** (plus Parquet/CSV) that supports:

* experimental guide selection
* large-scale resource generation
* programmatic querying

### Contents

1. Overview  
2. Key Features  
3. Installation  
4. Inputs  
5. Pipeline stages  
6. Output database and schema  
8. Querying the SQLite database  
9. Guide selection recommendations  
10. Performance considerations  
11. Troubleshooting  
12. Reproducibility and best practices  

---

## 2. Key Features

### Core design

* N- and C-terminal tagging (start/stop codons)
* Strand-aware gRNA discovery
* 240 bp homology arms (gene orientation)
* Synonymous mutation design to prevent re-cutting
* Validation primer design

### Annotation improvements

* Terminus-level grouping (`N1`, `C1`, etc.)
* Mapping of:

  * transcripts (`FBtr`)
  * isoforms (`-RA`, `-RB`, etc.)
* Detection of **potential stop-codon readthrough**

### Scoring and selection (new)

* Composite **guide selection score**
* Tiered selection system:

  * Tier 0: standard guides
  * Tier 1: multi-perfect-match guides
  * Tier 2: non-coding edit required (last resort)
* Chromosome-aware off-target scoring

### Off-target analysis (enhanced)

* Cas-OFFinder integration
* Per-mismatch counts (`n_mm*`)
* Chromosome-aware counts:

  * `n_mm*_same_chr`
  * `n_mm*_other_chr`

### Output formats

* SQLite (indexed, queryable)
* Parquet (analysis)
* CSV (human-readable)

---

## 3. Installation

### 3.1 Clone the repository

```bash
git clone https://github.com/kjihens/tftag-pipeline.git
cd tftag-pipeline
```

### 3.2 Create Conda environment

```bash
conda env create -f environment.yml
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

---

## 4. Inputs

### 4.1 Genome annotation (GTF/GFF)

- Provided via `--gtf`  
- Must include:

    * `start_codon`
    * `stop_codon`
    * `transcript_id`
    * `transcript_symbol`

- A `gffutils` database will be created automatically (or reused)  

---

### 4.2 Reference genome (FASTA)

- Provided via `--fasta`  
- Chromosome/contig names must match the annotation  

---

### 4.3 Gene selection

Optional. If omitted, **all genes** in the GTF database are processed.

Accepted formats:
- Path to a file with one gene ID per line  
- Comma-separated list of gene IDs 

Gene identifiers must match those in the annotation database.

Examples:

```bash
--genes genes.txt
--genes "FBgn0000008,FBgn0000017"
```

---

### 4.4 Selection modes

Optional. If omitted, **all guides** that were identified will be included.

```bash
--select-mode all
--select-mode closest
--select-mode rs3
--select-mode score
```

- all (default): Show all.
- closest: Chooses the guide cutting nearest to the start/stop codon.
- rs3: Chooses the guide with the best efficiency as predicted by rs3.
- score: Chooses the guide with the score

---

## 5. Pipeline stages

### 5.1 Stage A: Annotation parsing

- Builds or loads a `gffutils` database  
- Extracts start/stop codons  
- Assigns N/C tag identifiers (N1, C1, …)  

Adds:

* `terminus_transcripts`
* `terminus_isoforms`
* `potential_readthrough`

---

### 5.2 Stage B: PAM scanning

- Scans ±window around each codon  
- Identifies NGG PAMs on both strands  
- Computes cut site and distance to codon 

---

### 5.3 Stage C: Designability prefilter

Rejects guides if:
- Homology arms extend beyond contig boundaries  
- Primer windows cannot be accommodated  

---

### 5.4 Stage D: RS3 scoring

- Computes RS3 score using 30-mer context  
- Adds warnings if unavailable  

---

### 5.5 Stage E: Off-target analysis

- Runs Cas-OFFinder  
- Summarizes mismatch counts  

Adds:

| Column          | Meaning                      |
| --------------- | ---------------------------- |
| n_mm0           | total 0-mismatch, should be at least 1    |
| n_mm1           | total 1-mismatch off-targets |
| n_mm1_same_chr  | same chromosome              |
| n_mm1_other_chr | other chromosomes            |

---

### 5.6 Stage F: Homology arms

* 240 bp HAL/HAR
* gene-oriented sequences

---

### 5.7 Stage G: Silent edits

* Applied only if:

  * PAM + ≥13 bp spacer overlap contained in HAL or HAR
* Coding edits preferred
* Non-coding edits rejected

---

### 5.8 Stage H: Guide scoring

See section 8 for in-depth explanation

Components:

#### Distance

```
exp(-cut_distance / 10)
```

#### Off-target

```
1 / (1 + weighted penalty)
```

#### RS3

(logistic-transformed)

#### Edit score

| Case            | Score |
| --------------- | ----- |
| no edit         | 1.0   |
| coding edit     | 0.7   |
| non-coding edit | 0.1   |

---

### Stage I: Tier assignment

| Tier | Meaning                  |
| ---- | ------------------------ |
| 0    | preferred                |
| 1    | multiple perfect matches |
| 2    | rejected                 |

---

### Stage J: Missing guides

If no guide exists:

* row added
* `no_guide = TRUE`

---

## 6. Output database

### Key columns

| Category    | Columns                                 |
| ----------- | --------------------------------------- |
| Gene        | gene_id, gene_symbol, tag               |
| Isoforms    | terminus_transcripts, terminus_isoforms |
| Readthrough | potential_readthrough                   |
| Guide       | grna_seq_23                             |
| Off-target  | n_mm*, n_mm*_same_chr                   |
| Score       | selection_score, selection_tier         |
| QC          | warnings, selection_warning             |
| Primers     | PCR1_forward_primer_seq, PCR2_reverse_primer_seq |

---

## 7. Query examples

### Best guides

```sql
SELECT *
FROM guides
WHERE selection_tier = 0
ORDER BY gene_id, tag, selection_score DESC;
```

---

### All guides for gene

```sql
SELECT *
FROM guides
WHERE gene_symbol = 'pax2'
ORDER BY selection_tier ASC, selection_score DESC;
```

### All guides for gene, filtered and sorted

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

## 8. Guide Scoring

TFTag assigns each guide a **composite selection score** (`selection_score`) and a **selection tier** (`selection_tier`) to prioritise guides per terminus.

The system reflects practical experimental constraints:
- proximity to the tagging site  
- off-target risk (especially on the same chromosome)  
- predicted efficiency  
- need for sequence editing  

---

### 8.1 Overview of scoring workflow

For each guide:

1. Compute component scores:
   - distance score  
   - off-target score  
   - efficiency score (RS3)  
   - edit score  

2. Combine into a weighted composite score

3. Assign a **selection tier** (hard constraints)

4. Rank guides:
   - first by `selection_tier`
   - then by `selection_score`

---

### 8.2 Distance component

Guides closer to the codon are preferred.

```python
distance_score = exp(-cut_distance / d0)
```

- `cut_distance`: distance (bp) from Cas9 cut site to codon  
- `d0`: decay constant (default ≈ 10 bp)

Interpretation:
- 0 bp → 1.0  
- 10 bp → ~0.37  
- >30 bp → strongly penalised  

---

### 8.3 Off-target component

Off-targets are penalised with:
- strong weighting for low-mismatch hits  
- additional penalty for same-chromosome hits  

#### Weighted penalty

```python
penalty = (
    w1 * n_mm1 +
    w2 * n_mm2 +
    w3 * n_mm3 +
    w1_chr * n_mm1_same_chr
)
```

Typical behaviour:
- 1-mismatch hits dominate the penalty  
- same-chromosome hits are especially undesirable  

#### Final score

```python
offtarget_score = 1 / (1 + penalty)
```

Interpretation:
- 1.0 → no off-targets  
- ~0.5 → moderate risk  
- → 0 → severe off-target burden  

---

### 8.4 RS3 efficiency component

RS3 values are **log-odds**, not probabilities.

They are transformed as:

```python
rs3_score_norm = 1 / (1 + exp(-rs3_score))
```

Interpretation:
- ~0.5 → average efficiency  
- → 1 → high efficiency  
- → 0 → low efficiency  

---

### 8.5 Edit component

Guides requiring sequence edits are penalised.

| Case | Description | Score |
|------|-------------|------|
| no edit | no mutation required | 1.0 |
| coding edit | synonymous edit in CDS | 0.7 |
| non-coding edit | edit outside coding region | 0.1 |

---

### 8.6 Composite score

The final score is a weighted combination:

```python
selection_score = (
    w_dist * distance_score +
    w_off  * offtarget_score +
    w_rs3  * rs3_score_norm +
    w_edit * edit_score
)
```

Weights are configurable but typically:
- distance: high importance  
- off-targets: high importance  
- RS3: moderate  
- edits: moderate  

---

### 8.7 Selection tiers (hard constraints)

Guides are first classified into tiers:

| Tier | Condition | Meaning |
|------|----------|--------|
| 0 | n_mm0 == 1 AND not rejected | preferred |
| 1 | n_mm0 > 1 | multiple perfect matches |
| 2 | rejected = TRUE | requires non-coding edits |

Ranking is:
1. by `selection_tier`  
2. then by `selection_score`  

---

### 8.8 Warnings

Warnings are recorded in:

```
selection_warning
```

Triggered by:
- `n_mm0 > 1`
- non-coding edit required  

---

### 8.9 Practical interpretation

Best guides typically:
- cut within ~10 bp of the codon  
- have no 1-mismatch off-targets  
- have no same-chromosome off-targets  
- do not require edits  
- have high RS3 scores  

---

### 8.10 Recommended usage

Example query:

```sql
SELECT *
FROM guides
WHERE selection_tier = 0
ORDER BY selection_score DESC;
```

---

### 8.11 Notes

- Off-target penalties are asymmetric (low mismatches dominate)  
- Same-chromosome hits are penalised more heavily as more difficult to detect in the fly  
- RS3 is treated probabilistically after transformation  
- Tiering ensures biologically undesirable guides are never selected first  

## 9. Troubleshooting

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

---

## 10. Best practices

Record:

* genome build
* annotation version
* mismatch threshold
* RS3 model
* scoring parameters

## 11. Citation

If you use TFTag, please cite the associated publication (when available) and/or the GitHub repository.
