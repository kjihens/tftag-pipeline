"""
End-to-end orchestrator for TFTag.

High-level flow
---------------
1) Build codon/terminus table from genome annotation
2) Scan PAMs around each terminus
3) Optionally check stock VCF compatibility and guide/stock compatibility
4) Prefilter designable guides
5) Score RS3 efficiency
6) Enumerate and summarise off-targets
7) Build homology arms and determine edit requirements
8) Score/rank guides
9) Optionally select one guide per terminus
10) Design validation primers
11) Add no-guide placeholder rows
12) Write SQLite/CSV/Parquet and run logs

Conventions
-----------
- Genomic coordinates are 1-based inclusive.
- `mismatches` defines the maximum Cas-OFFinder mismatch class searched.
"""

from __future__ import annotations

import os
import uuid

import numpy as np

from . import annotate, efficiency, scan, stockcheck
from .genome_io import create_gtf_db, load_fasta_dict
from .helpers import filter_by_offtarget_mismatch, parse_genes_arg
from .off_targeting import enumerate_offtargets_cas_offinder, merge_specificity
from .runlog import TFTagRunLogger
from .select import add_guide_selection_score, select_one_per_tag
from .reporting import (
    add_no_guide_rows,
    add_provenance,
    candidate_coverage_summary,
    drop_internal_columns,
    final_output_summary,
    print_summary_block,
    terminus_summary,
    write_no_guide_only_output,
    write_outputs,
)
from .editing import choose_arm_for_mutation, apply_silent_edits
from .primers import design_validation_primers

def make_run_parameters(
    *,
    gtf_file: str,
    gtf_db_path: str,
    genome_fasta_path: str,
    genes: str | None,
    basename: str,
    outdir: str,
    output_table: str,
    pam_window_up: int,
    pam_window_down: int,
    tracrRNA: str,
    batch_size_rs3: int,
    protospacer_overlap_len: int,
    run_offtargets: bool,
    cas_offinder_bin: str,
    device_spec: str,
    mismatches: int,
    pam_pattern: str,
    min_offtarget_mismatch: int | None,
    offtarget_batch_size: int,
    selection: str,
    write_csv: bool,
    write_parquet: bool,
    check_stock_vcf_compatibility: bool,
    check_stock_variants: bool,
    stock_identical_only: bool,
    stock_vcfs: dict[str, str],
    chrom_to_stock: dict[str, str],
) -> dict:
    """Collect run parameters for reproducible logging."""
    return {
        "gtf_file": gtf_file,
        "gtf_db_path": gtf_db_path,
        "genome_fasta_path": genome_fasta_path,
        "genes": genes,
        "basename": basename,
        "outdir": outdir,
        "output_table": output_table,
        "pam_window_up": pam_window_up,
        "pam_window_down": pam_window_down,
        "tracrRNA": tracrRNA,
        "batch_size_rs3": batch_size_rs3,
        "protospacer_overlap_len": protospacer_overlap_len,
        "run_offtargets": run_offtargets,
        "cas_offinder_bin": cas_offinder_bin,
        "device_spec": device_spec,
        "mismatches": mismatches,
        "pam_pattern": pam_pattern,
        "min_offtarget_mismatch": min_offtarget_mismatch,
        "offtarget_batch_size": offtarget_batch_size,
        "selection": selection,
        "write_csv": write_csv,
        "write_parquet": write_parquet,
        "check_stock_vcf_compatibility": check_stock_vcf_compatibility,
        "check_stock_variants": check_stock_variants,
        "stock_identical_only": stock_identical_only,
        "stock_vcfs": stock_vcfs,
        "chrom_to_stock": chrom_to_stock,
    }


def run_pipeline(
    gtf_file: str,
    gtf_db_path: str,
    genome_fasta_path: str,
    genes: str | None = None,
    basename: str = "tftag",
    outdir: str = "out",
    output_table: str = "guides",
    pam_window_up: int = 30,
    pam_window_down: int = 30,
    tracrRNA: str = "Hsu2013",
    batch_size_rs3: int = 2048,
    protospacer_overlap_len: int = 13,
    run_offtargets: bool = True,
    cas_offinder_bin: str = "cas-offinder",
    device_spec: str = "C",
    mismatches: int = 4,
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG",
    min_offtarget_mismatch: int | None = 0,
    offtarget_batch_size: int = 500,
    selection: str = "all",
    write_csv: bool = True,
    write_parquet: bool = True,
    stock_vcfs: dict[str, str] | None = None,
    chrom_to_stock: dict[str, str] | None = None,
    check_stock_vcf_compatibility: bool = False,
    check_stock_variants: bool = False,
    stock_identical_only: bool = False,
) -> None:
    """
    Run the TFTag design pipeline.

    Parameters
    ----------
    selection:
        "all" keeps all candidate guides.
        "closest", "rs3", and "score" reduce to one guide per (gene_id, tag).

    min_offtarget_mismatch:
        None or 0 disables off-target filtering.
        Positive values require off-target enumeration and filter guides according
        to the helper logic in `helpers.filter_by_offtarget_mismatch`.
        -1 requests strict uniqueness within the searched mismatch space.
    """
    if selection not in ("all", "closest", "rs3", "score"):
        raise ValueError("selection must be one of: all, closest, rs3, score")

    if mismatches < 0:
        raise ValueError("mismatches must be >= 0")

    if min_offtarget_mismatch is not None:
        min_offtarget_mismatch = int(min_offtarget_mismatch)

        if min_offtarget_mismatch < -1:
            raise ValueError(
                "min_offtarget_mismatch must be -1, 0, a positive integer, or None"
            )

        if min_offtarget_mismatch > 0 and (min_offtarget_mismatch - 1) > mismatches:
            raise ValueError(
                f"min_offtarget_mismatch={min_offtarget_mismatch} requires "
                f"mismatches >= {min_offtarget_mismatch - 1}, but mismatches={mismatches}."
            )

    if not run_offtargets and min_offtarget_mismatch not in (None, 0):
        raise ValueError(
            "--min-offtarget-mismatch requires off-target enumeration. "
            "Remove --no-offtargets or set --min-offtarget-mismatch 0."
        )

    os.makedirs(outdir, exist_ok=True)

    run_id = uuid.uuid4().hex[:12]
    stock_vcfs = stock_vcfs or {}
    chrom_to_stock = chrom_to_stock or {}

    logger = TFTagRunLogger(outdir=outdir, basename=basename, run_id=run_id)

    logger.start(
        parameters=make_run_parameters(
            gtf_file=gtf_file,
            gtf_db_path=gtf_db_path,
            genome_fasta_path=genome_fasta_path,
            genes=genes,
            basename=basename,
            outdir=outdir,
            output_table=output_table,
            pam_window_up=pam_window_up,
            pam_window_down=pam_window_down,
            tracrRNA=tracrRNA,
            batch_size_rs3=batch_size_rs3,
            protospacer_overlap_len=protospacer_overlap_len,
            run_offtargets=run_offtargets,
            cas_offinder_bin=cas_offinder_bin,
            device_spec=device_spec,
            mismatches=mismatches,
            pam_pattern=pam_pattern,
            min_offtarget_mismatch=min_offtarget_mismatch,
            offtarget_batch_size=offtarget_batch_size,
            selection=selection,
            write_csv=write_csv,
            write_parquet=write_parquet,
            check_stock_vcf_compatibility=check_stock_vcf_compatibility,
            check_stock_variants=check_stock_variants,
            stock_identical_only=stock_identical_only,
            stock_vcfs=stock_vcfs,
            chrom_to_stock=chrom_to_stock,
        )
    )

    try:
        # ------------------------------------------------------------
        # Validate stock VCF arguments and optionally check compatibility
        # ------------------------------------------------------------
        if check_stock_vcf_compatibility and not stock_vcfs:
            raise ValueError("No stock VCFs supplied. Use --stock_vcf NAME=PATH.")

        if check_stock_variants and not stock_vcfs:
            raise ValueError(
                "Stock variant checking requested, but no stock VCFs were supplied."
            )

        if check_stock_variants and not chrom_to_stock:
            raise ValueError(
                "Stock variant checking requested, but no chromosome-to-stock mapping was supplied."
            )

        if check_stock_variants:
            for stock_name, vcf_path in stock_vcfs.items():
                stockcheck.ensure_fetchable_vcf(vcf_path)

        if check_stock_vcf_compatibility:
            logger.log("Checking stock VCF compatibility with reference FASTA...")

            for stock_name, vcf_path in stock_vcfs.items():
                summary = stockcheck.check_vcf_reference_compatibility(
                    vcf_path=vcf_path,
                    fasta_path=genome_fasta_path,
                    max_records=10000,
                )

                logger.log(f"\nStock: {stock_name}")
                logger.log(f"  VCF: {summary['vcf_path']}")
                logger.log(f"  Records checked: {summary['n_records_checked']}")
                logger.log(f"  REF mismatches: {summary['n_ref_mismatches']}")

                if summary["contigs_only_in_vcf"]:
                    extra = summary["contigs_only_in_vcf"]
                    logger.log(
                        f"  Contigs only in VCF: {len(extra)} "
                        f"(first 10: {extra[:10]})"
                    )

                if summary["contigs_only_in_fasta"]:
                    extra = summary["contigs_only_in_fasta"]
                    logger.log(
                        f"  Contigs only in FASTA: {len(extra)} "
                        f"(first 10: {extra[:10]})"
                    )

                if summary["length_mismatches"]:
                    logger.log(f"  Length mismatches: {summary['length_mismatches']}")

                if summary["mismatch_examples"]:
                    logger.log("  REF mismatch examples:")
                    for ex in summary["mismatch_examples"]:
                        logger.log(f"    {ex}")

        if check_stock_vcf_compatibility and not check_stock_variants:
            logger.add_summary(
                "Stock VCF compatibility",
                {
                    "status": "completed",
                    "n_stock_vcfs": len(stock_vcfs),
                },
            )
            logger.log("Stock VCF compatibility check completed.")
            logger.success()
            return

        # ------------------------------------------------------------
        # Load annotation/database and genome FASTA
        # ------------------------------------------------------------
        db = create_gtf_db(gtf_file, gtf_db_path)
        fasta_dict = load_fasta_dict(genome_fasta_path)

        genes_list = parse_genes_arg(genes)
        if genes_list is None:
            genes_list = [g.id for g in db.features_of_type("gene")]

        attribute = annotate.build_attribute_table(genes_list, db)
        input_summary = terminus_summary(attribute)

        logger.add_summary("Input summary", input_summary)

        if attribute.empty:
            logger.finish_early("No start_codon or stop_codon features found for requested genes.")
            return

        # ------------------------------------------------------------
        # PAM scanning
        # ------------------------------------------------------------
        candidates = scan.scan_for_guides(
            attribute,
            fasta_dict,
            window_up=pam_window_up,
            window_down=pam_window_down,
        )

        logger.add_summary("PAM scanning", candidate_coverage_summary(candidates))

        if candidates.empty:
            write_no_guide_only_output(
                attribute=attribute,
                run_id=run_id,
                gtf_db_path=gtf_db_path,
                genome_fasta_path=genome_fasta_path,
                basename=basename,
                outdir=outdir,
                output_table=output_table,
                write_csv=write_csv,
                write_parquet=write_parquet,
                logger=logger,
                reason="No candidate guides found in PAM search windows.",
            )
            return

        candidates = scan.add_grna_23_coordinates(candidates)

        # ------------------------------------------------------------
        # Optional stock-sequence compatibility annotation/filtering
        # ------------------------------------------------------------
        if check_stock_variants:
            logger.log("Annotating candidate guides against stock VCFs...")

            candidates = stockcheck.annotate_stock_variants(
                guides_df=candidates,
                stock_vcfs=stock_vcfs,
                chrom_to_stock=chrom_to_stock,
                reference_fasta_path=genome_fasta_path,
                show_progress=True,
            )

            removed_pam = 0
            removed_nonidentical = 0

            before = len(candidates)
            candidates = candidates[~candidates["stock_pam_gg_mutated"]].copy()
            removed_pam = before - len(candidates)

            if removed_pam:
                logger.log(
                    f"Removed {removed_pam} guides with PAM GG mutation in assigned stock."
                )

            if stock_identical_only:
                before = len(candidates)
                candidates = candidates[candidates["stock_seq_matches_ref"]].copy()
                removed_nonidentical = before - len(candidates)

                if removed_nonidentical:
                    logger.log(
                        f"Removed {removed_nonidentical} non-identical guides due to "
                        "--stock_identical_only."
                    )

            logger.add_summary(
                "Stock filtering",
                {
                    "removed_pam_gg_mutated": removed_pam,
                    "removed_nonidentical": removed_nonidentical,
                    "remaining_guides": len(candidates),
                },
            )

            if candidates.empty:
                write_no_guide_only_output(
                    attribute=attribute,
                    run_id=run_id,
                    gtf_db_path=gtf_db_path,
                    genome_fasta_path=genome_fasta_path,
                    basename=basename,
                    outdir=outdir,
                    output_table=output_table,
                    write_csv=write_csv,
                    write_parquet=write_parquet,
                    logger=logger,
                    reason="No guides remain after stock sequence filtering.",
                )
                return

        # ------------------------------------------------------------
        # Designability filtering
        # ------------------------------------------------------------
        candidates = prefilter_designable(
            candidates,
            fasta_dict,
            show_progress=True,
        )

        before = len(candidates)
        candidates = candidates[candidates["designable"]].copy()
        removed_non_designable = before - len(candidates)

        if removed_non_designable:
            logger.log(f"Removed {removed_non_designable} non-designable guides.")

        logger.add_summary(
            "Designability filtering",
            {
                "removed_non_designable": removed_non_designable,
                "remaining_guides": len(candidates),
            },
        )

        if candidates.empty:
            write_no_guide_only_output(
                attribute=attribute,
                run_id=run_id,
                gtf_db_path=gtf_db_path,
                genome_fasta_path=genome_fasta_path,
                basename=basename,
                outdir=outdir,
                output_table=output_table,
                write_csv=write_csv,
                write_parquet=write_parquet,
                logger=logger,
                reason="No designable guides remain after prefilter.",
            )
            return

        # ------------------------------------------------------------
        # RS3 efficiency scoring
        # ------------------------------------------------------------
        candidates = efficiency.score_rs3(
            candidates,
            fasta_dict,
            tracrRNA=tracrRNA,
            batch_size=batch_size_rs3,
        )

        # ------------------------------------------------------------
        # Off-target enumeration and optional filtering
        # ------------------------------------------------------------
        if run_offtargets:
            _hits, spec = enumerate_offtargets_cas_offinder(
                candidates,
                genome_fasta_path,
                outdir=outdir,
                run_id=run_id,
                cas_offinder_bin=cas_offinder_bin,
                device_spec=device_spec,
                mismatches=mismatches,
                pam_pattern=pam_pattern,
                batch_size=offtarget_batch_size,
                show_progress=True,
            )
            candidates = merge_specificity(candidates, spec, mismatches=mismatches)

            if min_offtarget_mismatch is not None and min_offtarget_mismatch != 0:
                before = len(candidates)
                candidates = filter_by_offtarget_mismatch(
                    candidates,
                    min_offtarget_mismatch,
                )
                removed_by_offtarget_filter = before - len(candidates)

                logger.add_summary(
                    "Off-target filtering",
                    {
                        "removed_by_offtarget_filter": removed_by_offtarget_filter,
                        "remaining_guides": len(candidates),
                    },
                )

                if candidates.empty:
                    write_no_guide_only_output(
                        attribute=attribute,
                        run_id=run_id,
                        gtf_db_path=gtf_db_path,
                        genome_fasta_path=genome_fasta_path,
                        basename=basename,
                        outdir=outdir,
                        output_table=output_table,
                        write_csv=write_csv,
                        write_parquet=write_parquet,
                        logger=logger,
                        reason="No guides remain after off-target filtering.",
                    )
                    return

        else:
            candidates["n_hits"] = np.nan
            for k in range(0, mismatches + 1):
                candidates[f"n_mm{k}"] = np.nan

        logger.add_dataframe_summary("Post off-target summary", candidates)

        # ------------------------------------------------------------
        # Homology arms, edit requirements, and silent edits
        # ------------------------------------------------------------
        candidates = add_homology_arms(
            candidates,
            fasta_dict,
            show_progress=True,
        )

        candidates = choose_arm_for_mutation(
            candidates,
            protospacer_overlap_len=protospacer_overlap_len,
            coding_only=True,
            show_progress=True,
        )

        candidates = apply_silent_edits(
            candidates,
            show_progress=True,
        )

        # ------------------------------------------------------------
        # Score all guides, optionally reduce to one guide per terminus
        # ------------------------------------------------------------
        candidates = add_guide_selection_score(candidates)
        logger.add_dataframe_summary("Post scoring summary", candidates)

        if selection != "all":
            before = len(candidates)
            candidates = select_one_per_tag(candidates, mode=selection)
            removed_by_selection = before - len(candidates)

            logger.add_summary(
                "Guide selection",
                {
                    "selection_mode": selection,
                    "removed_by_selection": removed_by_selection,
                    "remaining_guides": len(candidates),
                },
            )

            if candidates.empty:
                write_no_guide_only_output(
                    attribute=attribute,
                    run_id=run_id,
                    gtf_db_path=gtf_db_path,
                    genome_fasta_path=genome_fasta_path,
                    basename=basename,
                    outdir=outdir,
                    output_table=output_table,
                    write_csv=write_csv,
                    write_parquet=write_parquet,
                    logger=logger,
                    reason="No guides remain after guide selection.",
                )
                return

        # ------------------------------------------------------------
        # Validation primers are designed after optional selection so expensive
        # primer design is only run for retained guide rows.
        # ------------------------------------------------------------
        candidates = design_validation_primers(
            candidates,
            fasta_dict,
            show_progress=True,
        )

        # ------------------------------------------------------------
        # Add explicit rows for termini that have no retained guide.
        # ------------------------------------------------------------
        candidates = add_no_guide_rows(candidates, attribute)

        candidates = add_provenance(
            candidates,
            run_id=run_id,
            gtf_db_path=gtf_db_path,
            genome_fasta_path=genome_fasta_path,
        )

        candidates = drop_internal_columns(candidates)

        # ------------------------------------------------------------
        # Persist outputs
        # ------------------------------------------------------------
        write_outputs(
            candidates,
            basename=basename,
            outdir=outdir,
            output_table=output_table,
            write_csv=write_csv,
            write_parquet=write_parquet,
        )

        print(f"Finished. Wrote {len(candidates)} rows to:")
        print(f"   - {outdir + '/' + basename + '.sqlite'} (table: {output_table})")
        if write_parquet:
            print(f"   - {outdir + '/' + basename + '.parquet'}")
        if write_csv:
            print(f"   - {outdir + '/' + basename + '.csv'}")

        final_summary = final_output_summary(candidates)

        print_summary_block("Input summary", input_summary)
        print_summary_block("Final output summary", final_summary)

        logger.add_dataframe_summary("Final output summary", candidates)
        logger.success()

    except Exception as e:
        logger.failure(e)
        raise

    finally:
        txt_log, json_log = logger.write()
        print(f"Log written to: {txt_log}")
        print(f"JSON log written to: {json_log}")