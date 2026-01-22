"""
Extract start/stop codons and assign N/C tags.

Output table (`attribute_df`) is the contract for downstream modules:
- gene_id, gene_symbol
- terminus: 'N' (start codon) or 'C' (stop codon)
- feature: 'start_codon'|'stop_codon' (kept for clarity / provenance)
- codon_start, codon_end: genomic 1-based inclusive
- gene_strand: '+'|'-'
- chromosome
- tag: e.g. N1, N2, C1, ...
"""
from __future__ import annotations

from typing import Iterable

import pandas as pd
from tqdm.auto import tqdm
import gffutils


def build_attribute_table(genes: Iterable[str], db, *, verbose: bool = True) -> pd.DataFrame:
    """
    Build an attribute table for codon features for a set of genes.

    Parameters
    ----------
    genes:
      Iterable of gene IDs to look up in the gffutils database.
    db:
      gffutils DB handle (as returned by your genome_io.create_GTF_db).
    verbose:
      If True, print a short warning summary for missing genes.

    Returns
    -------
    pd.DataFrame
      One row per start/stop codon feature, with stable tag numbering per gene+terminus.
    """
    rows: list[dict] = []
    missing: list[str] = []

    for gene_id in tqdm(genes, desc="Fetch gene features", leave=False):
        try:
            gene = db[gene_id]
        except gffutils.exceptions.FeatureNotFoundError:
            missing.append(gene_id)
            continue

        # Prefer gene symbol/name from the gene feature, fall back to ID.
        gene_symbol = (
            (gene.attributes.get("gene_symbol") or gene.attributes.get("gene_name") or [None])[0]
        ) or gene_id

        for feat in db.children(gene, featuretype=["start_codon", "stop_codon"]):
            terminus = "N" if feat.featuretype == "start_codon" else "C"
            rows.append(
                dict(
                    gene_id=gene_id,
                    gene_symbol=gene_symbol,
                    terminus=terminus,
                    feature=feat.featuretype,
                    codon_start=int(feat.start),
                    codon_end=int(feat.end),
                    gene_strand=feat.strand,
                    chromosome=feat.seqid,
                )
            )

    df = pd.DataFrame(rows)
    if df.empty:
        if verbose and missing:
            print(f"Warning: {len(missing)} genes not found in gffutils DB (showing first 10): {missing[:10]}")
        return df

    # De-duplicate and enforce deterministic ordering prior to per-gene numbering.
    df = df.drop_duplicates(["gene_id", "feature", "codon_start", "codon_end"])
    df = df.sort_values(
        ["gene_id", "terminus", "chromosome", "codon_start", "codon_end"]
    ).reset_index(drop=True)

    # Number tags per gene and terminus (N1, N2, ... ; C1, C2, ...).
    df["tag_number"] = df.groupby(["gene_id", "terminus"]).cumcount() + 1
    df["tag"] = df["terminus"] + df["tag_number"].astype(str)
    df = df.drop(columns=["tag_number"])

    if verbose and missing:
        print(f"Warning: {len(missing)} genes not found in gffutils DB (showing first 10): {missing[:10]}")

    return df