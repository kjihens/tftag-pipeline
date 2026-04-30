"""
Genome I/O utilities for TFTag.

Provides:
- GTF → gffutils database creation/loading
- FASTA → in-memory dictionary

Design notes
------------
- Optimised for Drosophila-scale genomes (full FASTA in memory).
- Fails early with clear errors if inputs are invalid.
"""

from __future__ import annotations

import os
from typing import Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils


# ---------------------------------------------------------------------
# GTF database
# ---------------------------------------------------------------------

def create_gtf_db(gtf_file: str, db_path: str) -> gffutils.FeatureDB:
    """
    Create or load a gffutils database from a GTF/GFF file.

    Behaviour
    ---------
    - If db_path does not exist → create database
    - If db_path exists → load existing database

    Notes
    -----
    - Does NOT verify that db_path corresponds to the same GTF file.
      If you change GTF, you should delete the DB manually.
    """

    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")

    os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)

    if not os.path.exists(db_path):
        print(f"Creating GTF database: {db_path}")

        db = gffutils.create_db(
            gtf_file,
            dbfn=db_path,
            force=False,
            keep_order=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )

        print("GTF database created.")

    else:
        print(f"Loading existing GTF database: {db_path}")
        db = gffutils.FeatureDB(db_path, keep_order=True)

    return db


# ---------------------------------------------------------------------
# FASTA loading
# ---------------------------------------------------------------------

def load_fasta_dict(genome_fasta_path: str) -> Dict[str, SeqRecord]:
    """
    Load genome FASTA into a dictionary keyed by contig name.

    Returns
    -------
    dict[str, SeqRecord]

    Raises
    ------
    FileNotFoundError:
        if FASTA does not exist

    ValueError:
        if FASTA is empty or malformed

    Notes
    -----
    - Entire genome is loaded into memory.
    - Suitable for Drosophila-scale genomes.
    """

    if not os.path.exists(genome_fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {genome_fasta_path}")

    records = SeqIO.parse(genome_fasta_path, "fasta")
    fasta_dict = SeqIO.to_dict(records)

    if not fasta_dict:
        raise ValueError(f"FASTA file appears empty or invalid: {genome_fasta_path}")

    return fasta_dict