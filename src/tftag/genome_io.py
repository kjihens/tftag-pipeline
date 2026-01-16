"""
I/O helpers: load GTF/FASTA and persist results.
"""
import os
from Bio import SeqIO
import gffutils

def create_GTF_db(gtf_file: str, db_path: str):
    os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)

    if not os.path.exists(db_path):
        print(f"Creating GFF/GTF database at {db_path}...")
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
        print("GFF/GTF database created.")
    else:
        db = gffutils.FeatureDB(db_path, keep_order=True)
        print(f"Loaded existing GFF/GTF database from {db_path}.")
    return db


def load_fasta_dict(genome_fasta_path: str):
    """Load genome FASTA into a dict keyed by chromosome name (Drosophila-scale)."""
    return SeqIO.to_dict(SeqIO.parse(genome_fasta_path, "fasta"))

