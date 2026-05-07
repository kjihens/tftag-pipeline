"""
Annotation helpers for TFTag.

Purpose
-------
This module converts a set of gene IDs into a tidy "attribute table" with one row
per annotated start_codon / stop_codon feature. That table is the biological
starting point for the rest of the pipeline.

For each terminus row, we record:
- gene ID and gene symbol
- whether it is an N- or C-terminus
- genomic location and strand
- a stable per-gene terminus tag (N1, N2, ..., C1, C2, ...)
- which transcript isoforms contain that exact terminus
- whether a C-terminus is a potential readthrough-like candidate, based on GTF structure

Important interpretation notes
------------------------------
1. "potential_readthrough" is a structural annotation heuristic based on the GTF.
   It does NOT prove translational stop-codon readthrough.
2. A terminus may be shared by multiple isoforms, so transcript / isoform
   membership is explicitly annotated for each exact codon feature.
3. This implementation does NOT assume that start_codon / stop_codon / CDS are
   direct gffutils children of transcript features. Instead, it primarily matches
   transcript membership via feature attributes such as transcript_id. This is
   much more robust for standard GTFs.

Expected output columns
-----------------------
- gene_id
- gene_symbol
- terminus
- feature
- codon_start
- codon_end
- gene_strand
- chromosome
- tag
- terminus_transcripts
- terminus_isoforms
- potential_readthrough
"""

from __future__ import annotations

from typing import Iterable
import pandas as pd
from tqdm.auto import tqdm
import gffutils
from Bio.Seq import Seq
from .utils import get_sequence


# ---------------------------------------------------------------------
# Transcript / attribute helpers
# ---------------------------------------------------------------------

def _transcript_like_children(gene, db):
    """
    Return transcript-like children of a gene.

    Why this helper exists
    ----------------------
    Different annotations do not always use the same feature type for transcripts.
    Some use 'transcript', others 'mRNA', and some include other transcript-like
    containers. We therefore try a broad set first.

    Fallback behaviour
    ------------------
    If none of those types are found, we fall back to direct gene children that
    themselves look transcript-like because they own CDS / start_codon / stop_codon
    descendants.
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
    Return the most human-readable transcript label available.

    For FlyBase GTFs, transcript_symbol is the preferred label, e.g.:
      Polr1B-RA
    """
    attrs = tx.attributes

    for key in ("transcript_symbol", "transcript_name", "Name", "symbol", "Alias"):
        vals = attrs.get(key, [])
        if vals and vals[0]:
            return str(vals[0])

    vals = attrs.get("transcript_id", [])
    if vals and vals[0]:
        return str(vals[0])

    return str(tx.id)


def _extract_isoform_suffix(tx_label: str) -> str:
    """
    Extract a compact FlyBase-style isoform suffix.

    Examples:
      Polr1B-RA -> A
      dsx-RB    -> B
      ovo-RC    -> C
    """
    if "-R" in tx_label:
        suffix = tx_label.split("-R")[-1]
        if suffix:
            return suffix
    return tx_label


def _feature_transcript_ids(feat) -> set[str]:
    """
    Extract transcript IDs from a feature's attributes.

    Why this is important
    ---------------------
    Standard GTFs usually encode transcript membership through attributes such as
    transcript_id, rather than through a strict parent-child tree that gffutils
    can always traverse in the way we want.

    We therefore use attribute-based matching whenever possible.

    Supported keys
    --------------
    - transcript_id   : normal GTF
    - Parent / parent : useful in some GFF-like files

    Returns
    -------
    set[str]
      Possibly empty set of transcript identifiers associated with the feature.
    """
    attrs = feat.attributes
    ids = set()

    for key in ("transcript_id", "Parent", "parent"):
        vals = attrs.get(key, [])
        for v in vals:
            if v:
                ids.add(str(v))

    return ids


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
    Return transcript labels for transcripts of a gene that contain
    this exact start_codon or stop_codon feature.

    Matching strategy
    -----------------
    We do NOT assume that codon features are direct children of transcript
    features in the gffutils hierarchy. Instead, we:

    1. collect transcript-like children of the gene
    2. map transcript IDs -> human-readable transcript labels
    3. scan all matching codon features under the gene
    4. use attribute transcript IDs to assign those codon features to transcripts

    Exact matching criteria for the terminus itself
    -----------------------------------------------
    - same feature type (start_codon or stop_codon)
    - same chromosome
    - same strand
    - same start
    - same end
    """
    try:
        gene = db[gene_id]
    except gffutils.exceptions.FeatureNotFoundError:
        return []

    txs = _transcript_like_children(gene, db)
    if not txs:
        return []

    # Map any transcript-like identifier we can find to a readable label.
    tx_id_to_label = {}
    for tx in txs:
        label = _get_transcript_label(tx)
        tx_ids = _feature_transcript_ids(tx)
        tx_ids.add(tx.id)  # fallback
        for tid in tx_ids:
            tx_id_to_label[tid] = label

    matches = set()

    # Search codon features under the gene and map them back to transcripts via attributes.
    for feat in db.children(gene, featuretype=feature, order_by="start"):
        if (
            feat.seqid == chrom
            and feat.strand == strand
            and int(feat.start) == int(codon_start)
            and int(feat.end) == int(codon_end)
        ):
            feat_tx_ids = _feature_transcript_ids(feat)
            for tid in feat_tx_ids:
                if tid in tx_id_to_label:
                    matches.add(tx_id_to_label[tid])

    return sorted(matches)


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
# Translation-based readthrough / stop-codon QC
# ---------------------------------------------------------------------

from Bio.Seq import Seq
from .utils import get_sequence


def _transcript_label_to_id_map_for_gene(gene_id: str, db) -> dict[str, str]:
    """
    Map readable transcript labels back to transcript IDs.

    This allows us to use values stored in terminus_transcripts, e.g. dsx-RA,
    and recover the transcript_id used on CDS/start/stop_codon features.
    """
    try:
        gene = db[gene_id]
    except gffutils.exceptions.FeatureNotFoundError:
        return {}

    mapping: dict[str, str] = {}

    for tx in _transcript_like_children(gene, db):
        attrs = tx.attributes

        tx_ids = set()
        tx_ids.add(str(tx.id))

        for key in ("transcript_id", "Parent", "parent"):
            for value in attrs.get(key, []):
                if value:
                    tx_ids.add(str(value))

        labels = set(tx_ids)

        for key in ("transcript_symbol", "transcript_name", "Name", "symbol"):
            for value in attrs.get(key, []):
                if value:
                    labels.add(str(value))

        canonical_id = None
        tx_id_values = attrs.get("transcript_id", [])
        if tx_id_values:
            canonical_id = str(tx_id_values[0])
        else:
            canonical_id = str(tx.id)

        for label in labels:
            mapping[label] = canonical_id

    return mapping


def _parse_semicolon_list(value) -> list[str]:
    """Parse semicolon-separated annotation values."""
    if value is None or pd.isna(value):
        return []

    text = str(value).strip()
    if not text or text == "none":
        return []

    return [item.strip() for item in text.split(";") if item.strip()]


def _features_for_transcript_under_gene(
    gene,
    db,
    *,
    featuretype: str,
    transcript_id: str,
) -> list:
    """
    Return features of a given type assigned to transcript_id.

    Uses transcript_id / Parent attributes rather than assuming strict
    gffutils parent-child structure.
    """
    out = []

    for feat in db.children(gene, featuretype=featuretype, order_by="start"):
        feat_tx_ids = _feature_transcript_ids(feat)

        if transcript_id in feat_tx_ids:
            out.append(feat)

    return out


def _stop_codon_matches_terminus(
    stop_feat,
    *,
    chrom: str,
    strand: str,
    stop_start: int,
    stop_end: int,
) -> bool:
    """Return True if a stop_codon feature is the exact terminus being assessed."""
    return (
        stop_feat.seqid == chrom
        and stop_feat.strand == strand
        and int(stop_feat.start) == int(stop_start)
        and int(stop_feat.end) == int(stop_end)
    )


def _cds_sequence_for_transcript(
    gene,
    db,
    fasta_dict,
    *,
    transcript_id: str,
) -> str:
    """
    Reconstruct transcript CDS sequence in coding orientation.

    Notes
    -----
    Many GTFs store stop_codon separately from CDS. Therefore this function
    reconstructs only CDS features. Stop codon assessment is handled separately.
    """
    cds_feats = _features_for_transcript_under_gene(
        gene,
        db,
        featuretype="CDS",
        transcript_id=transcript_id,
    )

    if not cds_feats:
        return ""

    strand = cds_feats[0].strand

    if strand == "+":
        cds_feats = sorted(cds_feats, key=lambda feat: int(feat.start))
    elif strand == "-":
        cds_feats = sorted(cds_feats, key=lambda feat: int(feat.start), reverse=True)
    else:
        return ""

    pieces: list[str] = []

    for feat in cds_feats:
        pieces.append(
            get_sequence(
                fasta_dict,
                feat.seqid,
                int(feat.start),
                int(feat.end),
                strand,
            )
        )

    seq = "".join(pieces).upper()
    usable_len = (len(seq) // 3) * 3
    return seq[:usable_len]


def _stop_codon_sequence_for_terminus(
    fasta_dict,
    *,
    chrom: str,
    strand: str,
    stop_start: int,
    stop_end: int,
) -> str:
    """Return the annotated stop codon sequence in coding orientation."""
    seq = get_sequence(
        fasta_dict,
        chrom,
        int(stop_start),
        int(stop_end),
        strand,
    ).upper()

    return seq if len(seq) == 3 else ""


def _translation_based_stop_status_for_transcript(
    gene,
    db,
    fasta_dict,
    *,
    transcript_id: str,
    chrom: str,
    strand: str,
    stop_start: int,
    stop_end: int,
) -> dict:
    """
    Assess whether one transcript cleanly terminates at the queried stop codon.

    Logic
    -----
    For a transcript assigned to this terminus:
    1. It should contain the exact stop_codon feature.
    2. The terminus stop sequence should translate as '*'.
    3. CDS + stop should translate with exactly one terminal stop.
    4. Internal stop codons before the terminal stop are suspicious.

    Returns
    -------
    dict with:
      transcript_id
      status
      evidence
    """
    stop_feats = _features_for_transcript_under_gene(
        gene,
        db,
        featuretype="stop_codon",
        transcript_id=transcript_id,
    )

    has_exact_stop = any(
        _stop_codon_matches_terminus(
            stop_feat,
            chrom=chrom,
            strand=strand,
            stop_start=stop_start,
            stop_end=stop_end,
        )
        for stop_feat in stop_feats
    )

    if not has_exact_stop:
        return {
            "transcript_id": transcript_id,
            "status": "not_assigned_to_exact_stop",
            "evidence": "transcript lacks exact queried stop_codon feature",
        }

    cds_seq = _cds_sequence_for_transcript(
        gene,
        db,
        fasta_dict,
        transcript_id=transcript_id,
    )

    if not cds_seq:
        return {
            "transcript_id": transcript_id,
            "status": "cds_missing",
            "evidence": "no CDS sequence reconstructed",
        }

    stop_seq = _stop_codon_sequence_for_terminus(
        fasta_dict,
        chrom=chrom,
        strand=strand,
        stop_start=stop_start,
        stop_end=stop_end,
    )

    if len(stop_seq) != 3:
        return {
            "transcript_id": transcript_id,
            "status": "invalid_stop_length",
            "evidence": f"stop_codon sequence length is {len(stop_seq)}",
        }

    stop_aa = str(Seq(stop_seq).translate(table=1))

    if stop_aa != "*":
        return {
            "transcript_id": transcript_id,
            "status": "queried_stop_does_not_translate_as_stop",
            "evidence": f"queried stop sequence {stop_seq} translates as {stop_aa}",
        }

    coding_plus_stop = cds_seq + stop_seq
    usable_len = (len(coding_plus_stop) // 3) * 3
    coding_plus_stop = coding_plus_stop[:usable_len]

    protein = str(Seq(coding_plus_stop).translate(table=1))

    if "*" not in protein:
        return {
            "transcript_id": transcript_id,
            "status": "no_stop_in_translation",
            "evidence": "CDS+queried_stop translation contains no stop codon",
        }

    stop_positions = [i for i, aa in enumerate(protein) if aa == "*"]

    if stop_positions == [len(protein) - 1]:
        return {
            "transcript_id": transcript_id,
            "status": "clean_terminal_stop",
            "evidence": "CDS+queried_stop contains single terminal stop",
        }

    if stop_positions and stop_positions[-1] == len(protein) - 1:
        return {
            "transcript_id": transcript_id,
            "status": "internal_stop_plus_terminal_stop",
            "evidence": f"internal stop codon(s) at amino-acid positions {stop_positions[:-1]}",
        }

    return {
        "transcript_id": transcript_id,
        "status": "internal_stop_no_terminal_stop",
        "evidence": f"stop codon(s) at amino-acid positions {stop_positions}, none terminal",
    }


def add_translation_stop_status_columns(
    attribute_df: pd.DataFrame,
    db,
    fasta_dict,
    *,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Add translation-based stop/readthrough status columns.

    For start_codon rows:
        potential_readthrough = False
        stop_translation_status = "not_applicable"

    For stop_codon rows:
        reconstruct CDS for transcript(s) covered by that terminus, append the
        queried stop codon, translate, and assess whether termination is clean.

    Added columns
    -------------
    - potential_readthrough
    - stop_translation_status
    - stop_translation_evidence
    - readthrough_transcripts
    """
    df = attribute_df.copy()

    for col, default in [
        ("potential_readthrough", False),
        ("stop_translation_status", "not_applicable"),
        ("stop_translation_evidence", "none"),
        ("readthrough_transcripts", "none"),
    ]:
        if col not in df.columns:
            df[col] = default

    if df.empty:
        return df

    required = [
        "gene_id",
        "feature",
        "chromosome",
        "gene_strand",
        "codon_start",
        "codon_end",
        "terminus_transcripts",
    ]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"add_translation_stop_status_columns missing required columns: {missing}")

    cache: dict[tuple[str, str, str, int, int, tuple[str, ...]], list[dict]] = {}

    iterator = df.itertuples(index=True)
    if show_progress:
        iterator = tqdm(
            iterator,
            total=len(df),
            desc="Translation-based stop/readthrough annotation",
            leave=False,
        )

    for row in iterator:
        idx = row.Index

        if row.feature != "stop_codon":
            df.at[idx, "potential_readthrough"] = False
            df.at[idx, "stop_translation_status"] = "not_applicable"
            df.at[idx, "stop_translation_evidence"] = "not a stop_codon terminus"
            df.at[idx, "readthrough_transcripts"] = "none"
            continue

        transcript_labels = _parse_semicolon_list(row.terminus_transcripts)

        if not transcript_labels:
            df.at[idx, "potential_readthrough"] = False
            df.at[idx, "stop_translation_status"] = "unknown_no_terminus_transcripts"
            df.at[idx, "stop_translation_evidence"] = "no transcript labels assigned to terminus"
            df.at[idx, "readthrough_transcripts"] = "none"
            continue

        try:
            gene = db[row.gene_id]
        except gffutils.exceptions.FeatureNotFoundError:
            df.at[idx, "potential_readthrough"] = False
            df.at[idx, "stop_translation_status"] = "unknown_gene_missing"
            df.at[idx, "stop_translation_evidence"] = "gene not found in database"
            df.at[idx, "readthrough_transcripts"] = "none"
            continue

        label_to_id = _transcript_label_to_id_map_for_gene(row.gene_id, db)
        transcript_ids = [
            label_to_id.get(label, label)
            for label in transcript_labels
        ]

        key = (
            row.gene_id,
            row.chromosome,
            row.gene_strand,
            int(row.codon_start),
            int(row.codon_end),
            tuple(sorted(transcript_ids)),
        )

        if key not in cache:
            cache[key] = [
                _translation_based_stop_status_for_transcript(
                    gene,
                    db,
                    fasta_dict,
                    transcript_id=tx_id,
                    chrom=row.chromosome,
                    strand=row.gene_strand,
                    stop_start=int(row.codon_start),
                    stop_end=int(row.codon_end),
                )
                for tx_id in transcript_ids
            ]

        results = cache[key]

        statuses = [res["status"] for res in results]
        evidence = [
            f"{res['transcript_id']}:{res['status']}({res['evidence']})"
            for res in results
        ]

        suspicious_statuses = {
            "queried_stop_does_not_translate_as_stop",
            "no_stop_in_translation",
            "internal_stop_no_terminal_stop",
            "internal_stop_plus_terminal_stop",
        }

        suspicious = [
            res["transcript_id"]
            for res in results
            if res["status"] in suspicious_statuses
        ]

        if suspicious:
            df.at[idx, "potential_readthrough"] = True
            df.at[idx, "stop_translation_status"] = "suspicious_stop_translation"
            df.at[idx, "readthrough_transcripts"] = ";".join(suspicious)
        elif all(status == "clean_terminal_stop" for status in statuses):
            df.at[idx, "potential_readthrough"] = False
            df.at[idx, "stop_translation_status"] = "clean_terminal_stop"
            df.at[idx, "readthrough_transcripts"] = "none"
        else:
            df.at[idx, "potential_readthrough"] = False
            df.at[idx, "stop_translation_status"] = "uncertain"
            df.at[idx, "readthrough_transcripts"] = "none"

        df.at[idx, "stop_translation_evidence"] = "; ".join(evidence)

    return df


# ---------------------------------------------------------------------
# Main entry point used by the pipeline
# ---------------------------------------------------------------------

def build_attribute_table(
    genes: Iterable[str],
    db,
    *,
    fasta_dict=None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Build an attribute table for codon features for a set of genes.

    Parameters
    ----------
    genes:
      Iterable of gene IDs to look up in the gffutils database.

    db:
      gffutils DB handle.

    verbose:
      If True, print a short warning summary for missing genes and enable progress
      bars for the more expensive annotation steps.

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

    # -------------------------------------------------------------
    # Step 1: collect raw start/stop codon features gene by gene
    # -------------------------------------------------------------
    iterator = tqdm(genes, desc="Fetch gene features", leave=False) if verbose else genes

    for gene_id in iterator:
        try:
            gene = db[gene_id]
        except gffutils.exceptions.FeatureNotFoundError:
            missing.append(gene_id)
            continue

        # Prefer a readable gene symbol if available.
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

    # Early exit if nothing was found.
    if df.empty:
        if verbose and missing:
            print(f"Warning: {len(missing)} genes not found in gffutils DB (showing first 10): {missing[:10]}")
        return df

    # -------------------------------------------------------------
    # Step 2: de-duplicate and sort deterministically
    # -------------------------------------------------------------
    # Stable ordering is important because per-gene numbering below should be reproducible.
    df = df.drop_duplicates(["gene_id", "feature", "codon_start", "codon_end"])
    df = df.sort_values(
        ["gene_id", "terminus", "chromosome", "codon_start", "codon_end"]
    ).reset_index(drop=True)

    # -------------------------------------------------------------
    # Step 3: assign stable per-gene terminus tags
    # -------------------------------------------------------------
    # Example:
    #   gene X with two start codons -> N1, N2
    #   gene X with three stop codons -> C1, C2, C3
    df["tag_number"] = df.groupby(["gene_id", "terminus"]).cumcount() + 1
    df["tag"] = df["terminus"] + df["tag_number"].astype(str)
    df = df.drop(columns=["tag_number"])

    # -------------------------------------------------------------
    # Step 4: annotate which isoforms share each exact terminus
    # -------------------------------------------------------------
    df = add_terminus_isoform_columns(df, db, show_progress=verbose)

    # -------------------------------------------------------------
    # Step 5: annotate potential readthrough-like C-termini
    # -------------------------------------------------------------
    if fasta_dict is not None:
        df = add_translation_stop_status_columns(
            df,
            db,
            fasta_dict,
            show_progress=verbose,
        )
    else:
        df["potential_readthrough"] = False
        df["stop_translation_status"] = "not_assessed_no_fasta"
        df["stop_translation_evidence"] = "fasta_dict was not provided"
        df["readthrough_transcripts"] = "none"

    # -------------------------------------------------------------
    # Step 6: report missing genes, if any
    # -------------------------------------------------------------
    if verbose and missing:
        print(f"Warning: {len(missing)} genes not found in gffutils DB (showing first 10): {missing[:10]}")

    # Optional lightweight sanity prints during development
    # if verbose:
    #     print("Readthrough candidates:", int(df["potential_readthrough"].sum()))
    #     print("Termini with transcript assignments:", int((df["terminus_transcripts"] != "").sum()))

    return df

