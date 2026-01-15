# TFTag Pipeline

TFTag is a Python pipeline for designing CRISPR gRNAs and homology arms for
N- and C-terminal tagging of genes.

## Features
- Start/stop codon–targeted gRNA discovery
- 240 bp homology arm design
- Silent PAM/protospacer edits
- RS3 efficiency scoring (optional)
- Cas-OFFinder off-target enumeration (optional)
- CCLMoff off-target risk scoring (optional)
- Validation primer design
- SQLite + Parquet output database

## Installation

```bash
git clone https://github.com/kjihens/tftag-pipeline.git
cd tftag-pipeline
pip install -e .
Usage
tftag --gtf genome.gtf --fasta genome.fa --genes genes.txt
