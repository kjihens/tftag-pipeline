# TFTag Pipeline User Manual

## 1. Overview

**TFTag** is a Python pipeline for large-scale or targeted CRISPR/Cas9 donor and guide RNA design for endogenous protein tagging in *Drosophila melanogaster*.

The pipeline was designed primarily for the TFTag resource, which aims to systematically tag all *Drosophila* transcription factors at their endogenous loci, but it is broadly applicable to any gene set.

TFTag integrates:

* transcript-aware annotation parsing,
* strand-aware gRNA discovery,
* strain-specific sequence validation,
* off-target analysis,
* donor engineering,
* splice-aware synonymous blocking mutation design,
* guide ranking and prioritisation,
* large-scale reproducible output generation.

The pipeline produces:

* an indexed SQLite database,
* Parquet files,
* CSV tables,
* structured JSON provenance logs,
* human-readable run logs.

The system is explicitly designed for:

* genome-scale resource generation,
* high-throughput CRISPR engineering,
* reproducibility,
* explainable filtering,
* downstream programmatic querying.

All genomic coordinates are represented internally as **1-based inclusive coordinates**.

---

# 2. Design Philosophy

TFTag was developed around several core principles.

## 2.1 Explicit Coordinate Handling

CRISPR design pipelines frequently obscure strand orientation and coordinate conventions.

TFTag instead:

* stores all genomic intervals explicitly,
* tracks PAM and protospacer coordinates separately,
* stores gRNA sequences in guide orientation,
* preserves strand-aware interpretation throughout the pipeline.

This substantially reduces ambiguity during downstream analysis.

---

## 2.2 Transcript-Aware Biology

The pipeline operates at the level of transcript-supported termini rather than simplistic gene-level models.

TFTag therefore:

* groups compatible start/stop codons into termini,
* tracks transcript and isoform support,
* evaluates splice-aware edit safety,
* reconstructs CDS structures for readthrough detection.

This is particularly important for:

* alternatively spliced genes,
* compact insect genomes,
* exon-intron structures near tagging sites.

---

## 2.3 Reproducibility

Every major design decision is recorded explicitly.

The pipeline preserves:

* rejected guides,
* rejection reasons,
* warning states,
* scoring weights,
* run parameters,
* provenance metadata.

This allows complete post hoc reconstruction of guide selection decisions.

---

## 2.4 Conservative Donor Engineering

TFTag prioritises biologically safe donor design.

The pipeline therefore:

* avoids unnecessary donor mutations,
* rejects edits in non-coding regions by default,
* checks splice-site overlap,
* prefers minimal synonymous sequence perturbation,
* supports codon-usage-aware silent edits.

---

# 3. Key Features

## Core CRISPR Design

* N-terminal and C-terminal tagging
* Strand-aware SpCas9 NGG scanning
* Explicit PAM/protospacer coordinate handling
* HDR donor design
* Validation primer design
* Composite guide ranking

---

## Transcript-Aware Annotation

* Terminus grouping (`N1`, `C1`, etc.)
* Transcript-aware codon grouping
* Isoform tracking
* Translation-based readthrough detection
* Splice-aware edit safety annotation

---

## Off-Target Analysis

* Cas-OFFinder integration
* Per-mismatch off-target counts
* Same-chromosome vs other-chromosome separation
* Strict uniqueness filtering
* Composite specificity scoring

---

## Injection strain-Aware Design

* VCF compatibility checking
* Injection strain-specific guide reconstruction
* PAM mutation detection
* Chromosome-specific injection strain assignment

---

## Large-Scale Output Support

* SQLite database generation
* Indexed query support
* CSV export
* Parquet export
* JSON provenance logs

---

# 4. Installation

## 4.1 Clone Repository

```bash
git clone https://github.com/kjihens/tftag-pipeline.git
cd tftag-pipeline
```

---

## 4.2 Create Conda Environment

```bash
conda env create -f environment.yml
conda activate tftag
```

---

## 4.3 Install TFTag

```bash
pip install -e .
```

Verify installation:

```bash
tftag --help
```

---

# 5. Inputs

## 5.1 Genome Annotation (GTF/GFF)

Provided using:

```bash
--gtf annotation.gtf
```

The annotation must contain:

* `start_codon`
* `stop_codon`
* transcript identifiers
* exon structures

TFTag automatically creates and caches a `gffutils` database.

### Important assumptions

The pipeline assumes:

* CDS structures are valid,
* transcript exon ordering is correct,
* chromosome names match the FASTA.

---

## 5.2 Reference Genome FASTA

Provided using:

```bash
--fasta genome.fasta
```

Requirements:

* chromosome names must match the annotation,
* FASTA should represent the same genome build as the GTF.

---

## 5.3 Gene Selection

Optional.

If omitted, all genes are processed.

Accepted formats:

```bash
--genes genes.txt
```

or:

```bash
--genes "FBgn0000008,FBgn0000017"
```

---

## 5.4 Injection strain VCF Input

TFTag supports strain-aware design using VCF files.

Example:

```bash
--strain-vcf attP40=Cas9_on_2.vcf.gz
--strain-vcf attP2=Cas9_on_3.vcf.gz
```

Chromosome assignment:

```bash
--strain-group attP40=3L,3R
--strain-group attP2=2L,2R,X
```

This allows the pipeline to:

* reconstruct Injection strain-specific guide sequences,
* identify PAM-disrupting polymorphisms,
* reject incompatible guides.

---

## 5.5 Selection Modes

```bash
--selection all
--selection closest
--selection rs3
--selection score
```

### Modes

| Mode      | Behaviour                            |
| --------- | ------------------------------------ |
| `all`     | Keep all guides                      |
| `closest` | Keep guide nearest to codon          |
| `rs3`     | Keep guide with best RS3 score       |
| `score`   | Keep guide with best composite score |

---

# 6. Pipeline Architecture

The pipeline is divided into distinct stages.

---

# 6.1 Stage A — Annotation Parsing

The annotation module:

1. loads the GTF database,
2. extracts start/stop codons,
3. groups compatible termini,
4. associates transcripts and isoforms.

## Terminus grouping

Multiple transcripts may share the same physical terminus.

TFTag therefore groups termini into labels such as:

| Label | Meaning                 |
| ----- | ----------------------- |
| `N1`  | first unique N-terminus |
| `C1`  | first unique C-terminus |

---

## Translation-Based Readthrough Detection

The previous heuristic readthrough detection system has been replaced.

The pipeline now:

1. reconstructs transcript CDS structures,
2. appends the annotated stop codon,
3. translates the CDS,
4. evaluates stop codon behaviour directly.

### Output states

| Status                                    | Meaning                        |
| ----------------------------------------- | ------------------------------ |
| `clean_terminal_stop`                     | Single terminal stop codon     |
| `internal_stop_plus_terminal_stop`        | Premature stop codons detected |
| `no_stop_in_translation`                  | Translation lacks valid stop   |
| `queried_stop_does_not_translate_as_stop` | Annotated stop not recognised  |

### Output columns

| Column                      | Meaning                   |
| --------------------------- | ------------------------- |
| `potential_readthrough`     | TRUE/FALSE                |
| `stop_translation_status`   | Translation QC result     |
| `stop_translation_evidence` | Transcript-level evidence |
| `readthrough_transcripts`   | Affected transcript IDs   |

---

# 6.2 Stage B — PAM Scanning

The scanner searches both strands for SpCas9 NGG sites.

## Coordinate representation

Each guide stores:

| Column                  | Meaning                    |
| ----------------------- | -------------------------- |
| `protospacer_start/end` | 20 nt spacer coordinates   |
| `pam_start/end`         | PAM coordinates            |
| `grna_seq_23`           | Full 23-mer guide sequence |
| `spacer`                | 20 nt spacer               |
| `pam_seq`               | PAM sequence               |
| `grna_strand`           | Target strand              |
| `cut_pos`               | Cas9 cleavage position     |
| `cut_distance`          | Distance to codon          |

The 23-mer is always stored in guide orientation.

---

## Search Window Design

The scanner expands beyond the nominal codon window.

This prevents guides from being missed when:

* the PAM lies inside the search window,
* but part of the spacer extends outside it.

---

## Maximum Cut Distance Filtering

Optional:

```bash
--max-cut-distance 15
```

This removes guides whose cut site is too far from the insertion point.

Motivation:

HDR efficiency generally decreases with increasing cut distance.

---

# 6.3 Stage C — Injection strain-Aware Sequence Checking

The strain-check module reconstructs guide sequences after applying VCF variants.

## Behaviour

For each guide:

1. overlapping variants are fetched,
2. Injection strain-specific sequence is reconstructed,
3. PAM integrity is evaluated.

---

## Rejection conditions

| Reason                     | Meaning                      |
| -------------------------- | ---------------------------- |
| `strain_pam_gg_mutated`     | PAM GG disrupted             |
| `strain_nonidentical_23mer` | Guide differs from reference |

---

# 6.4 Stage D — Designability Filtering

The donor prefilter evaluates whether a guide can realistically support donor construction.

Guides are rejected if:

* homology arms extend beyond contig boundaries,
* validation primer windows cannot be constructed,
* required intervals exceed chromosome limits.

## Designability outputs

| Column        | Meaning     |
| ------------- | ----------- |
| `designable`  | TRUE/FALSE  |
| `skip_reason` | Explanation |

---

# 6.5 Stage E — RS3 Efficiency Scoring

TFTag uses Rule Set 3 (RS3) scoring.

## 30-mer construction

The RS3 sequence context is:

```text
4 nt upstream
20 nt spacer
3 nt PAM
3 nt downstream
```

The sequence is always generated in guide orientation.

---

## Output columns

| Column         | Meaning            |
| -------------- | ------------------ |
| `rs3_30mer`    | RS3 input sequence |
| `rs3_score`    | Raw RS3 score      |
| `rs3_tracrRNA` | tracrRNA model     |

---

# 6.6 Stage F — Off-Target Analysis

Off-target enumeration is performed using Cas-OFFinder.

## Mismatch counting

The pipeline stores exact mismatch counts:

| Column  | Meaning         |
| ------- | --------------- |
| `n_mm0` | Perfect matches |
| `n_mm1` | 1-mismatch hits |
| `n_mm2` | 2-mismatch hits |
| `n_mm3` | 3-mismatch hits |
| `n_mm4` | 4-mismatch hits |

Additional chromosome-aware summaries:

| Column            | Meaning           |
| ----------------- | ----------------- |
| `n_mm1_same_chr`  | Same chromosome   |
| `n_mm1_other_chr` | Other chromosomes |

---

## Strict uniqueness filtering

```bash
--min-offtarget-mismatch -1
```

Requires:

* exactly one perfect genomic hit,
* no off-targets within searched mismatch space.

---

# 6.7 Stage G — Homology Arm Construction

TFTag constructs:

| Arm   | Meaning            |
| ----- | ------------------ |
| `HAL` | Left homology arm  |
| `HAR` | Right homology arm |

Default length:

```text
240 bp
```

Both genomic and gene-oriented sequences are stored.

---

# 6.8 Stage H — Edit-Arm Determination

The donor is evaluated for potential Cas9 re-cutting.

## Core rule

An arm requires editing only if it contains:

1. the complete PAM,
2. at least N PAM-proximal spacer bases.

Default:

```text
13 PAM-proximal bases
```

This substantially reduces unnecessary donor editing.

---

## Coding-only restriction

By default:

* coding edits are allowed,
* non-coding edits are rejected.

This avoids introducing uncertain regulatory mutations.

---

# 6.9 Stage I — Splice-Aware Edit Safety

The splice-check module evaluates whether proposed blocking mutations overlap biologically sensitive regions.

## Important design decision

The system only flags edits if the mutation itself overlaps:

* intronic sequence,
* splice donor sites,
* splice acceptor sites.

It does NOT penalise:

* exons overlapping introns of alternative isoforms,
* exon overlap with unrelated genes,
* transcript structures unaffected by the mutation.

This greatly reduces false positive splice warnings.

---

## Splice-related output columns

| Column                      | Meaning                    |
| --------------------------- | -------------------------- |
| `edit_overlaps_splice_site` | TRUE/FALSE                 |
| `edit_overlaps_intron`      | TRUE/FALSE                 |
| `splice_risk`               | Combined splice annotation |
| `splice_risk_details`       | Human-readable explanation |

---

# 6.10 Stage J — Silent Blocking Mutation Design

TFTag attempts to prevent donor re-cutting using synonymous edits.

---

## Biological rationale

Experimental evidence suggests blocking mutations are most effective in:

* the PAM,
* the four spacer bases closest to the PAM.

TFTag therefore prioritises mutations in these regions.

---

## Mutation hierarchy

### Tier 1 — PAM GG mutation

Attempt one synonymous mutation in:

* PAM position 22,
* PAM position 23.

PAM position 21 is not prioritised because it usually preserves NGG recognition.

---

### Tier 2 — PAM-proximal spacer mutation

Attempt one synonymous mutation in the four PAM-proximal spacer positions.

---

### Tier 3 — Dual distal mutations

If no high-efficiency blocking mutation exists:

* introduce two synonymous mutations elsewhere in the protospacer.

---

## Mutation reporting

| Column                 | Meaning                    |
| ---------------------- | -------------------------- |
| `edit_priority`        | Blocking strategy used     |
| `n_blocking_mutations` | Number of edits introduced |

Example values:

| Value                 | Meaning                      |
| --------------------- | ---------------------------- |
| `pam_gg_single`       | Single PAM mutation          |
| `pam_proximal_single` | Single PAM-proximal mutation |
| `double_distal`       | Two distal mutations         |
| `none_possible`       | No synonymous edit possible  |

---

# 6.11 Stage K — Codon Usage-Aware Synonymous Editing

TFTag supports codon-aware synonymous mutation selection.

## Motivation

Synonymous codons are not biologically equivalent.

Codon substitutions may influence:

* translation efficiency,
* mRNA stability,
* co-translational folding.

---

## Codon usage table generation

The pipeline can:

1. reconstruct CDS sequences,
2. calculate genome-wide codon frequencies,
3. cache codon usage tables.

If a cached table already exists:

* recalculation is skipped.

---

## Silent mutation optimisation modes

### Codon-usage preserving mode

Choose synonymous codons whose usage frequency most closely matches the original codon.

### GC-preserving mode

Fallback/default mode.

Chooses synonymous codons that preserve GC content most closely.

This minimises sequence perturbation.

---

# 6.12 Stage L — Composite Guide Scoring

TFTag assigns:

* a composite `selection_score`,
* a `selection_tier`.

The ranking system integrates:

* HDR geometry,
* specificity,
* efficiency,
* donor engineering complexity,
* splice safety.

---

# 7. Guide Scoring System

## 7.1 Scoring Philosophy

The scoring system combines:

* hard biological constraints,
* soft prioritisation metrics.

Hard constraints determine the tier.

Soft metrics determine ranking within the tier.

---

## 7.2 Distance Component

Guides closer to the insertion site are preferred.

```python
distance_score = exp(-cut_distance / 10)
```

Interpretation:

| Distance | Approximate score  |
| -------- | ------------------ |
| 0 bp     | 1.0                |
| 10 bp    | 0.37               |
| 30 bp    | strongly penalised |

---

## 7.3 Off-Target Component

Low-mismatch off-targets dominate the specificity penalty.

Same-chromosome hits receive additional weighting because they are more difficult to segregate genetically.

### Penalty model

```python
penalty = (
    w1 * n_mm1 +
    w2 * n_mm2 +
    w3 * n_mm3 +
    w_chr * n_mm1_same_chr
)
```

### Final score

```python
offtarget_score = 1 / (1 + penalty)
```

---

## 7.4 RS3 Component

RS3 outputs log-odds rather than probabilities.

TFTag converts them using a logistic transform:

```python
rs3_score_norm = 1 / (1 + exp(-rs3_score))
```

---

## 7.5 Edit Complexity Component

Guides requiring donor edits are penalised.

Additional penalties are applied for:

* multiple edits,
* distal edits,
* splice-risk overlap.

---

## 7.6 Composite Score

```python
selection_score = (
    w_dist   * distance_score +
    w_off    * offtarget_score +
    w_rs3    * rs3_score_norm +
    w_edit   * edit_score +
    w_splice * splice_score
)
```

Scoring weights are written into:

* run logs,
* JSON provenance files.

---

## 7.7 Selection Tiers

| Tier | Meaning                          |
| ---- | -------------------------------- |
| 0    | Preferred guides                 |
| 1    | Multiple perfect genomic matches |
| 2    | Rejected guides                  |

Guides are ranked:

1. by tier,
2. by composite score.

---

# 8. Output Structure

## Main Outputs

| File       | Description                 |
| ---------- | --------------------------- |
| `.sqlite`  | Indexed relational database |
| `.parquet` | Efficient analytical format |
| `.csv`     | Human-readable table        |
| `.log`     | Human-readable run log      |
| `.json`    | Structured provenance log   |

---

## Important Output Columns

| Category      | Columns                                            |
| ------------- | -------------------------------------------------- |
| Gene          | `gene_id`, `gene_symbol`, `tag`                    |
| Transcript    | `terminus_transcripts`, `terminus_isoforms`        |
| Readthrough   | `potential_readthrough`, `stop_translation_status` |
| Guide         | `grna_seq_23`, `cut_distance`                      |
| Specificity   | `n_mm*`, `n_mm*_same_chr`                          |
| Editing       | `requires_edit_arm`, `edit_priority`               |
| Splice safety | `splice_risk`, `splice_risk_details`               |
| Ranking       | `selection_score`, `selection_tier`                |
| QC            | `warnings`, `reject_reason`                        |

---

# 9. Query Examples

## Best guides

```sql
SELECT *
FROM guides
WHERE selection_tier = 0
ORDER BY gene_id, tag, selection_score DESC;
```

---

## All guides for a gene

```sql
SELECT *
FROM guides
WHERE gene_symbol = 'pax2'
ORDER BY selection_tier ASC, selection_score DESC;
```

---

## Guides requiring splice-risk review

```sql
SELECT *
FROM guides
WHERE splice_risk != 'none';
```

---

## Guides rejected during filtering

```sql
SELECT gene_symbol, tag, spacer, reject_reason
FROM guides
WHERE rejected = 1;
```

---

# 10. Performance Considerations

The most computationally expensive stages are:

1. Cas-OFFinder enumeration,
2. RS3 scoring,
3. Primer3 execution.

TFTag minimises unnecessary work by:

* filtering early,
* only scoring retained guides,
* only designing primers after final selection.

---

# 11. Troubleshooting

## No guides found

Potential causes:

* no NGG PAMs near the terminus,
* chromosome naming mismatch,
* aggressive cut-distance filtering,
* Injection strain-specific PAM mutations.

---

## Many guides rejected

Potential causes:

* excessive off-target burden,
* non-coding edit requirement,
* splice-risk overlap,
* donor design constraints.

---

## Missing RS3 scores

Potential causes:

* truncated 30-mer context,
* RS3 package unavailable,
* chromosome-edge effects.

---

## Cas-OFFinder failures

Check:

* binary path,
* CUDA/device specification,
* input FASTA indexing.

---

# 12. Reproducibility Recommendations

Always record:

* genome build,
* annotation version,
* mismatch threshold,
* scoring weights,
* RS3 model,
* Injection strain VCF versions,
* TFTag version.

---

# 13. Citation

If you use TFTag, please cite:

* the associated publication (when available),
* the TFTag repository,
* Cas-OFFinder,
* Rule Set 3.
