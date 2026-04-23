"""
Annotation helpers for TFTag.

What this module does
--------------------
This module turns a set of gene IDs into a tidy "attribute table" with one row
per annotated start_codon / stop_codon feature. That table is the biological
starting point for the rest of the pipeline.

For each terminus row, we record:
- the gene ID and gene symbol
- whether it is an N- or C-terminus
- genomic location and strand
- a stable per-gene terminus tag (N1, N2, ..., C1, C2, ...)
- which transcript isoforms contain that exact terminus
- whether a C-terminus is a potential readthrough-like candidate, based on GTF structure

Important interpretation notes
------------------------------
1. "potential_readthrough" is a *structural annotation heuristic* based on the GTF.
   It does NOT prove translational stop-codon readthrough.
2. A terminus may be shared by multiple isoforms. We therefore explicitly annotate
   transcript / isoform membership for each exact codon feature.
3. Column names in this module are chosen to match the newer TFTag pipeline schema:
     gene_id, gene_symbol, terminus, feature,
     codon_start, codon_end, gene_strand, chromosome, tag, ...

Expected input
--------------
- `genes`: iterable of gene IDs resolvable in the gffutils database
- `db`: a gffutils FeatureDB created from the annotation GTF/GFF

Expected output
---------------
A pandas DataFrame with one row per start/stop codon feature.
"""

from __future__ import annotations

from typing import Iterable
import pandas as pd
from tqdm.auto import tqdm
import gffutils


# ---------------------------------------------------------------------
# Transcript discovery helpers
# ---------------------------------------------------------------------

def _transcript_like_children(gene, db):
    """
    Return transcript-like children of a gene.

    Why this exists
    ---------------
    Different annotations do not always use the same feature type for transcripts.
    Some use "transcript", others "mRNA", and some include ncRNA-style transcript
    containers as well. We therefore try a broad but still sensible set first.

    Fallback behaviour
    ------------------
    If no recognised transcript feature types are found, we fall back to direct
    gene children that themselves have CDS/start_codon/stop_codon descendants.
    That makes the code more robust to annotation idiosyncrasies.
    """
    tx_types = [
        "transcript",
        "mRNA",
        "ncRNA",
        "tRNA",
        "rRNA",
        "lnc_RNA",
        "snRNA",
        "snoRNA",
    ]

    txs = list(db.children(gene, featuretype=tx_types, order_by="start"))
    if txs:
        return txs

    # Fallback for less standard annotations:
    # take direct children that behave like transcript containers.
    out = []
    for child in db.children(gene, level=1, order_by="start"):
        has_cds = any(True for _ in db.children(child, featuretype="CDS"))
        has_start = any(True for _ in db.children(child, featuretype="start_codon"))
        has_stop = any(True for _ in db.children(child, featuretype="stop_codon"))
        if has_cds or has_start or has_stop:
            out.append(child)
    return out


def _get_transcript_label(tx) -> str:
    """
    Return a human-readable transcript label if possible.

    Preference order
    ----------------
    1. transcript_name
    2. Name
    3. transcript_id
    4. feature ID

    For FlyBase-style annotations this often yields names like:
      dsx-RA, dsx-RB, ovo-RC
    """
    attrs = tx.attributes
    return (
        (attrs.get("transcript_name") or attrs.get("Name") or attrs.get("transcript_id") or [None])[0]
        or tx.id
    )


def _extract_isoform_suffix(tx_label: str) -> str:
    """
    Extract a compact isoform suffix from a transcript label when possible.

    Typical Drosophila / FlyBase examples:
      dsx-RA -> A
      dsx-RB -> B
      ovo-RC -> C

    If the expected '-R...' pattern is not present, return the original label.
    This is deliberate: we do not want to silently invent or lose information.
    """
    if "-R" in tx_label and len(tx_label.split("-R")[-1]) >= 1:
        return tx_label.split("-R")[-1]
    return tx_label


# ---------------------------------------------------------------------
# Terminus-to-isoform mapping
# ---------------------------------------------------------------------

def _transcripts_with_exact_terminus(
    gene_id: str,
    feature: str,
    chrom: str,
    strand: str,
    codon_start: int,
    codon_end: int,
    db,
) -> list[str]:
    """
    Return all transcript labels for transcripts of a gene that contain
    this exact start_codon or stop_codon feature.

    Exact matching criteria
    -----------------------
    We require equality of:
    - feature type (start_codon or stop_codon)
    - chromosome
    - strand
    - start
    - end

    This is important because a single gene can have multiple distinct N- or C-termini,
    and multiple transcripts can share one exact terminus.
    """
    try:
        gene = db[gene_id]
    except gffutils.exceptions.FeatureNotFoundError:
        return []

    txs = _transcript_like_children(gene, db)
    matches = []

    for tx in txs:
        for feat in db.children(tx, featuretype=feature, order_by="start"):
            if (
                feat.seqid == chrom
                and feat.strand == strand
                and int(feat.start) == int(codon_start)
                and int(feat.end) == int(codon_end)
            ):
                matches.append(_get_transcript_label(tx))
                break

    return sorted(set(matches))


def add_terminus_isoform_columns(
    attribute_df: pd.DataFrame,
    db,
    *,
    show_progress: bool = True,
    transcript_col: str = "terminus_transcripts",
    isoform_col: str = "terminus_isoforms",
) -> pd.DataFrame:
    """
    Add columns listing all transcripts / isoforms that contain each exact terminus.

    Adds
    ----
    transcript_col:
      Semicolon-separated transcript labels, e.g. "dsx-RA;dsx-RB"

    isoform_col:
      Semicolon-separated compact isoform suffixes, e.g. "A;B"

    Why this matters
    ----------------
    A start or stop codon can be shared by multiple isoforms. For tagging and
    interpretation, it is very useful to know whether a terminus is unique or shared.
    """
    df = attribute_df.copy()
    df[transcript_col] = ""
    df[isoform_col] = ""

    required = ["gene_id", "feature", "chromosome", "gene_strand", "codon_start", "codon_end"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"add_terminus_isoform_columns missing required columns: {missing}")

    # Cache repeated lookups because many rows can refer to the same exact terminus.
    cache = {}

    iterator = df.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Annotating terminus isoforms", leave=False)

    for row in iterator:
        idx = row.Index

        key = (
            row.gene_id,
            row.feature,
            row.chromosome,
            row.gene_strand,
            int(row.codon_start),
            int(row.codon_end),
        )

        if key not in cache:
            txs = _transcripts_with_exact_terminus(
                gene_id=row.gene_id,
                feature=row.feature,
                chrom=row.chromosome,
                strand=row.gene_strand,
                codon_start=int(row.codon_start),
                codon_end=int(row.codon_end),
                db=db,
            )
            isoforms = sorted(set(_extract_isoform_suffix(tx) for tx in txs))
            cache[key] = (txs, isoforms)

        txs, isoforms = cache[key]
        df.at[idx, transcript_col] = ";".join(txs)
        df.at[idx, isoform_col] = ";".join(isoforms)

    return df


# ---------------------------------------------------------------------
# Readthrough-like structural heuristic
# ---------------------------------------------------------------------

def _collect_transcript_cds_and_stops(tx, db):
    """
    Collect basic CDS and stop-codon information for one transcript-like feature.

    Returns
    -------
    dict or None
      None if the transcript has no CDS.
      Otherwise a dictionary containing:
      - tx_id
      - chromosome
      - strand
      - cds_start
      - cds_end
      - stop_sites (set of exact stop-codon intervals)

    Why this is useful
    ------------------
    To identify C-termini that might be readthrough-like, we compare transcripts
    of the same gene and ask whether one transcript terminates at a given stop
    while another transcript has CDS continuation beyond it in gene orientation.
    """
    cds_feats = list(db.children(tx, featuretype="CDS", order_by="start"))
    if not cds_feats:
        return None

    stop_feats = list(db.children(tx, featuretype="stop_codon", order_by="start"))

    cds_start = min(f.start for f in cds_feats)
    cds_end = max(f.end for f in cds_feats)

    return {
        "tx_id": tx.id,
        "chromosome": tx.seqid,
        "strand": tx.strand,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "stop_sites": {(f.start, f.end) for f in stop_feats},
    }


def _is_potential_readthrough_for_stop(
    gene_id: str,
    chrom: str,
    strand: str,
    stop_start: int,
    stop_end: int,
    db,
) -> bool:
    """
    Heuristic GTF-only test for whether a given stop codon could represent
    a readthrough-like C-terminus.

    We flag True if:
      1) at least one transcript of the gene uses this exact stop_codon
      2) another transcript has CDS extending beyond this stop in gene orientation
      3) the extension is codon-compatible (multiple of 3 bp)

    Important limitation
    --------------------
    This is NOT proof of translational stop-codon readthrough.
    It is a structural annotation heuristic intended to flag potentially tricky
    C-termini for donor / tagging design and downstream interpretation.
    """
    try:
        gene = db[gene_id]
    except gffutils.exceptions.FeatureNotFoundError:
        return False

    txs = _transcript_like_children(gene, db)
    if not txs:
        return False

    tx_info = []
    for tx in txs:
        info = _collect_transcript_cds_and_stops(tx, db)
        if info is None:
            continue
        if info["chromosome"] != chrom or info["strand"] != strand:
            continue
        tx_info.append(info)

    if not tx_info:
        return False

    # At least one transcript must explicitly terminate at this exact stop.
    stopping_txs = [t for t in tx_info if (stop_start, stop_end) in t["stop_sites"]]
    if not stopping_txs:
        return False

    # Look for another transcript with CDS continuation beyond this stop.
    for t in tx_info:
        if (stop_start, stop_end) in t["stop_sites"]:
            continue

        if strand == "+":
            # On + strand, "beyond stop" means downstream, i.e. larger coordinate.
            if t["cds_end"] <= stop_end:
                continue
            extension = t["cds_end"] - stop_end
        else:
            # On - strand, "beyond stop" in gene orientation means smaller coordinate.
            if t["cds_start"] >= stop_start:
                continue
            extension = stop_start - t["cds_start"]

        # Require codon-compatible extension.
        if extension > 0 and extension % 3 == 0:
            return True

    return False


def add_potential_readthrough_column(
    attribute_df: pd.DataFrame,
    db,
    *,
    show_progress: bool = True,
    column_name: str = "potential_readthrough",
) -> pd.DataFrame:
    """
    Add a boolean column marking C-termini that are potential readthrough candidates.

    Behaviour
    ---------
    - stop_codon rows are evaluated
    - start_codon rows are always False

    Practical use
    -------------
    This warns that a given C-terminus may be shared with, or structurally related to,
    a longer coding isoform, which can matter for C-terminal tagging.
    """
    df = attribute_df.copy()
    df[column_name] = False

    required = ["gene_id", "chromosome", "gene_strand", "codon_start", "codon_end", "feature"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"add_potential_readthrough_column missing required columns: {missing}")

    cache = {}

    iterator = df.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Annotating readthrough candidates", leave=False)

    for row in iterator:
        idx = row.Index

        if row.feature != "stop_codon":
            df.at[idx, column_name] = False
            continue

        key = (
            row.gene_id,
            row.chromosome,
            row.gene_strand,
            int(row.codon_start),
            int(row.codon_end),
        )

        if key not in cache:
            cache[key] = _is_potential_readthrough_for_stop(
                gene_id=row.gene_id,
                chrom=row.chromosome,
                strand=row.gene_strand,
                stop_start=int(row.codon_start),
                stop_end=int(row.codon_end),
                db=db,
            )

        df.at[idx, column_name] = cache[key]

    return df


# ---------------------------------------------------------------------
# Main entry point used by the pipeline
# ---------------------------------------------------------------------

def build_attribute_table(genes: Iterable[str], db, *, verbose: bool = True) -> pd.DataFrame:
    """
    Build an attribute table for codon features for a set of genes.

    Parameters
    ----------
    genes:
      Iterable of gene IDs to look up in the gffutils database.

    db:
      gffutils DB handle (as returned by genome_io.create_GTF_db).

    verbose:
      If True, print a short warning summary for missing genes and enable progress bars
      for the more expensive annotation steps.

    Returns
    -------
    pd.DataFrame
      One row per start/stop codon feature, with:
      - stable tag numbering per gene+terminus
      - transcript / isoform membership for that exact terminus
      - potential readthrough flag for C-termini
    """
    rows: list[dict] = []
    missing: list[str] = []

    # -----------------------------------------------------------------
    # Step 1: collect raw start/stop codon features gene by gene
    # -----------------------------------------------------------------
    for gene_id in tqdm(genes, desc="Fetch gene features", leave=False):
        try:
            gene = db[gene_id]
        except gffutils.exceptions.FeatureNotFoundError:
            missing.append(gene_id)
            continue

        # Prefer a human-readable gene symbol/name if available.
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

    # If nothing was found, return early after reporting missing genes if requested.
    if df.empty:
        if verbose and missing:
            print(f"Warning: {len(missing)} genes not found in gffutils DB (showing first 10): {missing[:10]}")
        return df

    # -----------------------------------------------------------------
    # Step 2: de-duplicate and sort deterministically
    # -----------------------------------------------------------------
    # This is important because later per-gene numbering (N1, N2, C1, C2, ...)
    # should be stable across runs.
    df = df.drop_duplicates(["gene_id", "feature", "codon_start", "codon_end"])
    df = df.sort_values(
        ["gene_id", "terminus", "chromosome", "codon_start", "codon_end"]
    ).reset_index(drop=True)

    # -----------------------------------------------------------------
    # Step 3: assign stable per-gene terminus tags
    # -----------------------------------------------------------------
    # Example:
    #   gene X with two start codons -> N1, N2
    #   gene X with three stop codons -> C1, C2, C3
    df["tag_number"] = df.groupby(["gene_id", "terminus"]).cumcount() + 1
    df["tag"] = df["terminus"] + df["tag_number"].astype(str)
    df = df.drop(columns=["tag_number"])

    # -----------------------------------------------------------------
    # Step 4: annotate which isoforms share each exact terminus
    # -----------------------------------------------------------------
    df = add_terminus_isoform_columns(df, db, show_progress=verbose)

    # -----------------------------------------------------------------
    # Step 5: annotate potential readthrough-like C-termini
    # -----------------------------------------------------------------
    df = add_potential_readthrough_column(df, db, show_progress=verbose)

    # -----------------------------------------------------------------
    # Step 6: report missing gene IDs, if any
    # -----------------------------------------------------------------
    if verbose and missing:
        print(f"Warning: {len(missing)} genes not found in gffutils DB (showing first 10): {missing[:10]}")

    return df