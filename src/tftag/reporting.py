from __future__ import annotations
import pandas as pd
import numpy as np
import os
from .runlog import TFTagRunLogger
from .storage import to_sqlite


def add_no_guide_rows(
    candidates: pd.DataFrame,
    attribute: pd.DataFrame,
) -> pd.DataFrame:
    """
    Ensure every terminus from the attribute table appears in the final output.

    If no guide was found for a terminus, add one placeholder row with:
      - guide_found = False
      - no_guide_reason = "no guide RNA found in PAM search window"
      - selection_tier = 99
      - selection_warning = "no guide RNA found"
    """
    if attribute.empty:
        return candidates

    target_cols = [
        "gene_id",
        "gene_symbol",
        "terminus",
        "feature",
        "tag",
        "chromosome",
        "gene_strand",
        "codon_start",
        "codon_end",
        "terminus_transcripts",
        "terminus_isoforms",
        "potential_readthrough",
    ]
    target_cols = [c for c in target_cols if c in attribute.columns]

    termini = attribute[target_cols].drop_duplicates(["gene_id", "tag"]).copy()

    if candidates.empty:
        found_keys = pd.DataFrame(columns=["gene_id", "tag"])
    else:
        found_keys = candidates[["gene_id", "tag"]].drop_duplicates()

    missing = termini.merge(
        found_keys,
        on=["gene_id", "tag"],
        how="left",
        indicator=True,
    )
    missing = missing[missing["_merge"] == "left_only"].drop(columns=["_merge"])

    if missing.empty:
        out = candidates.copy()
        out["guide_found"] = True
        if "no_guide_reason" not in out.columns:
            out["no_guide_reason"] = "none"
        return out

    # Add all columns present in candidates to the missing-placeholder rows.
    for col in candidates.columns:
        if col not in missing.columns:
            missing[col] = pd.NA

    missing["guide_found"] = False
    missing["no_guide_reason"] = "no guide RNA found in PAM search window"
    missing["selection_tier"] = 99
    missing["selection_score"] = np.nan
    missing["selection_warning"] = "no guide RNA found"
    missing["warnings"] = "no guide RNA found"

    out = candidates.copy()
    out["guide_found"] = True
    if "no_guide_reason" not in out.columns:
        out["no_guide_reason"] = "none"

    # Match column order.
    missing = missing[out.columns]
    return pd.concat([out, missing], ignore_index=True)

def final_output_summary(candidates: pd.DataFrame) -> dict:
    """
    Summarise final output, including explicit no-guide placeholder rows.

    `guide_found=False` rows represent termini that were analysed but produced
    no usable guide row.
    """
    if candidates.empty:
        return {
            "output_rows": 0,
            "guide_rows": 0,
            "no_guide_termini": 0,
            "genes_with_guides": 0,
            "termini_with_guides": 0,
        }

    if "guide_found" in candidates.columns:
        guide_mask = candidates["guide_found"].fillna(True).astype(bool)
    else:
        guide_mask = pd.Series(True, index=candidates.index)

    return {
        "output_rows": int(len(candidates)),
        "guide_rows": int(guide_mask.sum()),
        "no_guide_termini": int((~guide_mask).sum()),
        "genes_with_guides": int(candidates.loc[guide_mask, "gene_id"].nunique())
        if "gene_id" in candidates.columns
        else 0,
        "termini_with_guides": int(
            candidates.loc[guide_mask, ["gene_id", "tag"]].drop_duplicates().shape[0]
        )
        if {"gene_id", "tag"}.issubset(candidates.columns)
        else 0,
    }
def drop_internal_columns(candidates: pd.DataFrame) -> pd.DataFrame:
    """
    Remove columns useful internally but not normally needed in final user output.

    Important: do not drop `guide_found`, `no_guide_reason`,
    `selection_score`, `selection_tier`, or `selection_warning`.
    """
    drop_columns = [
        "spacer",
        "pam_seq",
        "protospacer_start",
        "protospacer_end",
        "pam_start",
        "pam_end",
        "designable",
        "skip_reason",
        "stock_pam_gg_mutated",
        "HALs",
        "HARs",
        "HALe",
        "HARe",
        "HAL_pam_in_arm",
        "HAR_pam_in_arm",
        "HAL_pam_proximal_overlap",
        "HAR_pam_proximal_overlap",
    ]

    present = [c for c in drop_columns if c in candidates.columns]
    return candidates.drop(columns=present)

def terminus_summary(attribute: pd.DataFrame) -> dict:
    """Summarise the requested/analysed terminus table."""
    if attribute.empty:
        return {
            "genes_analysed": 0,
            "total_termini": 0,
            "termini_counts": {},
        }

    termini_counts = (
        attribute["terminus"].value_counts().to_dict()
        if "terminus" in attribute.columns
        else {}
    )

    return {
        "genes_analysed": int(attribute["gene_id"].nunique()),
        "total_termini": int(len(attribute)),
        "termini_counts": termini_counts,
    }


def candidate_coverage_summary(candidates: pd.DataFrame) -> dict:
    """Summarise how many genes/termini currently have candidate guide rows."""
    if candidates.empty:
        return {
            "candidate_guides": 0,
            "genes_with_candidates": 0,
            "termini_with_candidates": 0,
        }

    return {
        "candidate_guides": int(len(candidates)),
        "genes_with_candidates": int(candidates["gene_id"].nunique())
        if "gene_id" in candidates.columns
        else 0,
        "termini_with_candidates": int(
            candidates[["gene_id", "tag"]].drop_duplicates().shape[0]
        )
        if {"gene_id", "tag"}.issubset(candidates.columns)
        else 0,
    }


def print_summary_block(title: str, summary: dict) -> None:
    """Print a simple human-readable summary block."""
    print(f"\n{title}:")
    for key, value in summary.items():
        if key == "termini_counts" and isinstance(value, dict):
            for term, count in sorted(value.items()):
                print(f"  {term}-termini: {count}")
        else:
            print(f"  {key}: {value}")


def add_provenance(
    candidates: pd.DataFrame,
    *,
    run_id: str,
    gtf_db_path: str,
    genome_fasta_path: str,
) -> pd.DataFrame:
    """Add run-level provenance columns to final output rows."""
    out = candidates.copy()
    out["run_id"] = run_id
    out["gtf_db_path"] = gtf_db_path
    out["genome_fasta_path"] = genome_fasta_path
    return out

def write_outputs(
    candidates: pd.DataFrame,
    *,
    basename: str,
    outdir: str,
    output_table: str,
    write_csv: bool,
    write_parquet: bool,
) -> None:
    """Write SQLite plus optional CSV/Parquet snapshots."""
    to_sqlite(
        candidates,
        basename,
        outdir,
        output_table,
        if_exists="append",
        index=False,
        create_indices=True,
    )

    if write_parquet:
        candidates.to_parquet(os.path.join(outdir, basename + ".parquet"), index=False)

    if write_csv:
        candidates.to_csv(os.path.join(outdir, basename + ".csv"), index=False)


def write_no_guide_only_output(
    *,
    attribute: pd.DataFrame,
    run_id: str,
    gtf_db_path: str,
    genome_fasta_path: str,
    basename: str,
    outdir: str,
    output_table: str,
    write_csv: bool,
    write_parquet: bool,
    logger: TFTagRunLogger,
    reason: str,
) -> None:
    """
    Write a valid output file containing only no-guide placeholder rows.

    This handles the edge case where PAM scanning produced no candidate guides
    anywhere. It is still useful to emit a complete terminus-level report.
    """
    candidates = annotate.add_no_guide_rows(pd.DataFrame(), attribute)
    candidates = add_provenance(
        candidates,
        run_id=run_id,
        gtf_db_path=gtf_db_path,
        genome_fasta_path=genome_fasta_path,
    )
    candidates = drop_internal_columns(candidates)

    write_outputs(
        candidates,
        basename=basename,
        outdir=outdir,
        output_table=output_table,
        write_csv=write_csv,
        write_parquet=write_parquet,
    )

    final_summary = final_output_summary(candidates)

    print(reason)
    print(f"Finished. Wrote {len(candidates)} no-guide rows to:")
    print(f"   - {outdir + '/' + basename + '.sqlite'} (table: {output_table})")
    if write_parquet:
        print(f"   - {outdir + '/' + basename + '.parquet'}")
    if write_csv:
        print(f"   - {outdir + '/' + basename + '.csv'}")

    print_summary_block("Final output summary", final_summary)

    logger.add_summary("Early no-guide output", {"reason": reason})
    logger.add_dataframe_summary("Final output summary", candidates)
    logger.success()