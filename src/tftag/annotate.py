"""
Extract start/stop codons and assign N/C tags.
"""
from __future__ import annotations
from typing import Iterable
import pandas as pd
from tqdm.auto import tqdm
import gffutils


def build_attribute_table(genes: Iterable[str], db) -> pd.DataFrame:
    rows = []
    missing = []

    for gene_id in tqdm(genes, desc="Fetch gene features", leave=False):
        try:
            gene = db[gene_id]
        except gffutils.exceptions.FeatureNotFoundError:
            missing.append(gene_id)
            continue

        # Prefer symbol/name from gene feature
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
        return df

    df = df.drop_duplicates(["gene_id", "feature", "codon_start", "codon_end"])
    df = df.sort_values(["gene_id", "feature", "chromosome", "codon_start", "codon_end"]).reset_index(drop=True)

    df["tag_number"] = df.groupby(["gene_id", "feature"]).cumcount() + 1
    df["tag"] = df["terminus"] + df["tag_number"].astype(str)
    df = df.drop(columns=["tag_number"])

    if missing:
        print(f"Warning: {len(missing)} genes not found in gffutils DB (showing first 10): {missing[:10]}")

    return df