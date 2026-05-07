"""
End-to-end orchestrator for TFTag.

Pipeline overview
-----------------
1. Build a terminus table from genome annotation.
2. Scan each terminus for nearby NGG PAMs.
3. Optionally annotate stock-specific guide compatibility.
4. Mark non-designable guides as rejected.
5. Score retained guides with RS3.
6. Run Cas-OFFinder on retained guides and optionally reject poor-specificity guides.
7. Build homology arms and determine whether donor edits are required.
8. Mark guides requiring non-coding-arm edits as rejected.
9. Apply synonymous donor edits to retained guides.
10. Score all tested guides.
11. Optionally select one retained guide per terminus, marking others as non-selected.
12. Design validation primers for retained guides.
13. Add explicit no-PAM/no-guide rows for termini with no scanned gRNAs.
14. Write SQLite/CSV/Parquet outputs and run logs.

Output model
------------
The pipeline does not silently drop tested guides. Every scanned gRNA remains in
the final output table with:

- guide_status
- guide_reject_reason

If no PAM was found for a terminus, reporting.add_no_guide_rows() adds one
placeholder row for that terminus.

Coordinate convention
---------------------
Genomic coordinates are 1-based inclusive throughout.
"""

from __future__ import annotations

import os
import uuid

import numpy as np

from . import annotate, efficiency, scan, stockcheck, splicecheck
from .editing import apply_silent_edits, choose_arm_for_mutation
from .genome_io import create_gtf_db, load_fasta_dict
from .helpers import filter_by_offtarget_mismatch, parse_genes_arg
from .homology import add_homology_arms, prefilter_designable
from .off_targeting import enumerate_offtargets_cas_offinder, merge_specificity
from .primers import design_validation_primers
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
from .runlog import TFTagRunLogger
from .select import (
    DEFAULT_SCORING_WEIGHTS,
    add_guide_selection_score,
    select_one_per_tag,
)
from .status import initialise_guide_status, mark_rejected, retained_mask


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
    codon_choice: str = "gc",
    codon_usage_table: str | None = None,
    force_recompute_codon_usage: bool = False,
) -> dict:
    """
    Collect run parameters for reproducible logging.

    Keeping this in one place ensures the text and JSON logs describe the exact
    settings used for each run.
    """
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
        "guide_scoring_weights": DEFAULT_SCORING_WEIGHTS,
        "write_csv": write_csv,
        "write_parquet": write_parquet,
        "check_stock_vcf_compatibility": check_stock_vcf_compatibility,
        "check_stock_variants": check_stock_variants,
        "stock_identical_only": stock_identical_only,
        "stock_vcfs": stock_vcfs,
        "chrom_to_stock": chrom_to_stock,
        "codon_choice": codon_choice,
        "codon_usage_table": codon_usage_table,
        "force_recompute_codon_usage": force_recompute_codon_usage,
    }


def _merge_work_back(candidates, work):
    """
    Merge a processed retained-guide subset back into the full candidate table.

    This is the central pattern used by the pipeline: process only retained rows,
    then copy updated columns back while preserving rejected rows.
    """
    out = candidates.copy()

    for col in work.columns:
        out.loc[work.index, col] = work[col]

    return out


def _mark_guides_not_in_filtered_result(work, filtered, reason: str):
    """
    Mark retained guides rejected if they are absent from a filtered result.

    This assumes `filtered` preserves original row indices. Therefore
    select_one_per_tag() should not reset the index.
    """
    kept_index = set(filtered.index)
    failed = ~work.index.isin(kept_index)

    return mark_rejected(work, failed, reason)


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
    codon_choice: str = "gc",
    codon_usage_table: str | None = None,
    force_recompute_codon_usage: bool = False,
    
) -> None:
    """
    Run the TFTag design pipeline.

    Selection modes
    ---------------
    - all:
        Keep all retained and rejected tested guides.
    - closest, rs3, score:
        Select one retained guide per (gene_id, tag), and mark other retained
        guides as rejected with reason "non_selected".

    Off-target filtering
    --------------------
    min_offtarget_mismatch:
    - None or 0: no filtering
    - positive integer: reject guides with off-targets below this mismatch class
    - -1: strict uniqueness within searched mismatch space
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
            for _stock_name, vcf_path in stock_vcfs.items():
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
        # Load annotation/database, genome FASTA and codon usage table (if needed)
        # ------------------------------------------------------------
        db = create_gtf_db(gtf_file, gtf_db_path)
        fasta_dict = load_fasta_dict(genome_fasta_path)

        codon_usage_lookup = {}

        if codon_choice == "usage":
            from .codon_usage import codon_usage_lookup as make_lookup
            from .codon_usage import get_or_build_codon_usage_table

            cache_path = codon_usage_table or os.path.join(outdir, "codon_usage.tsv")

            codon_usage_df = get_or_build_codon_usage_table(
                db=db,
                fasta_dict=fasta_dict,
                cache_path=cache_path,
                genes=None,
                force_recompute=force_recompute_codon_usage,
                metadata={
                    "gtf_file": gtf_file,
                    "genome_fasta_path": genome_fasta_path,
                },
                show_progress=True,
            )

            codon_usage_lookup = make_lookup(codon_usage_df)

        genes_list = parse_genes_arg(genes)
        if genes_list is None:
            genes_list = [gene.id for gene in db.features_of_type("gene")]

        attribute = annotate.build_attribute_table(genes_list, db)
        input_summary = terminus_summary(attribute)

        logger.add_summary("Input summary", input_summary)

        if attribute.empty:
            logger.finish_early(
                "No start_codon or stop_codon features found for requested genes."
            )
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

        # From here on, every scanned guide has explicit lifecycle state.
        candidates = initialise_guide_status(candidates)
        candidates = scan.add_grna_23_coordinates(candidates)

        # ------------------------------------------------------------
        # Optional stock-specific guide compatibility annotation
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

            pam_mutated = candidates["stock_pam_gg_mutated"].fillna(True).astype(bool)
            rejected_pam = int(pam_mutated.sum())

            candidates = mark_rejected(
                candidates,
                pam_mutated,
                "stock_pam_gg_mutated",
            )

            if rejected_pam:
                logger.log(
                    f"Marked {rejected_pam} guides rejected due to PAM GG mutation in assigned stock."
                )

            rejected_nonidentical = 0

            if stock_identical_only:
                nonidentical = ~candidates["stock_seq_matches_ref"].astype("boolean").fillna(False)
                rejected_nonidentical = int(nonidentical.sum())

                candidates = mark_rejected(
                    candidates,
                    nonidentical,
                    "stock_nonidentical_23mer",
                )

                if rejected_nonidentical:
                    logger.log(
                        f"Marked {rejected_nonidentical} guides rejected due to "
                        "--stock_identical_only."
                    )

            logger.add_summary(
                "Stock filtering",
                {
                    "rejected_pam_gg_mutated": rejected_pam,
                    "rejected_nonidentical": rejected_nonidentical,
                    "currently_retained_guides": int(retained_mask(candidates).sum()),
                    "total_tested_guides": len(candidates),
                },
            )

        # ------------------------------------------------------------
        # Designability prefilter
        # ------------------------------------------------------------
        # Only retained guides are tested further. Previously rejected guides stay
        # in the table but are not processed by expensive downstream steps.
        work = candidates.loc[retained_mask(candidates)].copy()

        if not work.empty:
            work = prefilter_designable(
                work,
                fasta_dict,
                show_progress=True,
            )

            non_designable = ~work["designable"].astype("boolean").fillna(False)
            rejected_non_designable = int(non_designable.sum())

            work = mark_rejected(
                work,
                non_designable,
                "non_designable",
            )

            candidates = _merge_work_back(candidates, work)
        else:
            rejected_non_designable = 0

        logger.add_summary(
            "Designability filtering",
            {
                "rejected_non_designable": rejected_non_designable,
                "currently_retained_guides": int(retained_mask(candidates).sum()),
                "total_tested_guides": len(candidates),
            },
        )

        # ------------------------------------------------------------
        # RS3 efficiency scoring
        # ------------------------------------------------------------
        work = candidates.loc[retained_mask(candidates)].copy()

        if not work.empty:
            work = efficiency.score_rs3(
                work,
                fasta_dict,
                tracrRNA=tracrRNA,
                batch_size=batch_size_rs3,
            )
            candidates = _merge_work_back(candidates, work)

        # ------------------------------------------------------------
        # Off-target enumeration and optional off-target rejection
        # ------------------------------------------------------------
        work = candidates.loc[retained_mask(candidates)].copy()

        if run_offtargets and not work.empty:
            _hits, spec = enumerate_offtargets_cas_offinder(
                work,
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

            work = merge_specificity(work, spec, mismatches=mismatches)

            if min_offtarget_mismatch is not None and min_offtarget_mismatch != 0:
                filtered = filter_by_offtarget_mismatch(
                    work,
                    min_offtarget_mismatch,
                    max_mismatches=mismatches,
                )

                rejected_by_offtarget = len(work) - len(filtered)

                work = _mark_guides_not_in_filtered_result(
                    work,
                    filtered,
                    "failed_offtarget_filter",
                )
            else:
                rejected_by_offtarget = 0

            candidates = _merge_work_back(candidates, work)

            logger.add_summary(
                "Off-target filtering",
                {
                    "rejected_by_offtarget_filter": rejected_by_offtarget,
                    "currently_retained_guides": int(retained_mask(candidates).sum()),
                    "total_tested_guides": len(candidates),
                },
            )

        elif not run_offtargets:
            # NaN means "not computed", not "zero off-targets".
            candidates["n_hits"] = np.nan

            for k in range(mismatches + 1):
                candidates[f"n_mm{k}"] = np.nan
                candidates[f"n_mm{k}_same_chr"] = np.nan
                candidates[f"n_mm{k}_other_chr"] = np.nan

        logger.add_dataframe_summary("Post off-target summary", candidates)

        # ------------------------------------------------------------
        # Homology arms and donor edit requirements
        # ------------------------------------------------------------
        work = candidates.loc[retained_mask(candidates)].copy()

        if not work.empty:
            work = add_homology_arms(
                work,
                fasta_dict,
                show_progress=True,
            )

            work = choose_arm_for_mutation(
                work,
                protospacer_overlap_len=protospacer_overlap_len,
                coding_only=True,
                show_progress=True,
            )

            # choose_arm_for_mutation uses `rejected=True` for guides requiring
            # edits in the non-coding arm. Convert that internal flag into the
            # general guide-status model.
            noncoding_edit = work["rejected"].astype("boolean").fillna(False)
            rejected_noncoding_edit = int(noncoding_edit.sum())

            work = mark_rejected(
                work,
                noncoding_edit,
                "requires_noncoding_homology_arm_edit",
            )

            candidates = _merge_work_back(candidates, work)
        else:
            rejected_noncoding_edit = 0

        logger.add_summary(
            "Edit-arm filtering",
            {
                "rejected_requires_noncoding_homology_arm_edit": rejected_noncoding_edit,
                "currently_retained_guides": int(retained_mask(candidates).sum()),
                "total_tested_guides": len(candidates),
            },
        )

        # ------------------------------------------------------------
        # Silent edits
        # ------------------------------------------------------------
        # Recompute retained rows after edit-arm filtering so rejected guides are
        # not silently edited.
        work = candidates.loc[retained_mask(candidates)].copy()

        if not work.empty:
            work = apply_silent_edits(
                work,
                codon_choice=codon_choice,
                codon_usage_lookup=codon_usage_lookup,
                show_progress=True,
            )

            # Annotate potential splice-disrupting risk of any edits applied in the previous step.
            work = splicecheck.annotate_edit_splice_risk(
                work,
                db,
                show_progress=True,
            )

            candidates = _merge_work_back(candidates, work)

        # ------------------------------------------------------------
        # Guide scoring
        # ------------------------------------------------------------
        # Score all tested guides so rejected guides still have interpretable
        # ranking information where available.
        candidates = add_guide_selection_score(candidates)

        logger.add_summary("Guide scoring weights", DEFAULT_SCORING_WEIGHTS)
        logger.add_dataframe_summary("Post scoring summary", candidates)

        # ------------------------------------------------------------
        # Optional one-guide-per-terminus selection
        # ------------------------------------------------------------
        if selection != "all":
            work = candidates.loc[retained_mask(candidates)].copy()

            if not work.empty:
                selected = select_one_per_tag(work, mode=selection)
                selected_index = set(selected.index)

                non_selected = ~work.index.isin(selected_index)
                rejected_non_selected = int(non_selected.sum())

                work = mark_rejected(
                    work,
                    non_selected,
                    "non_selected",
                )

                candidates = _merge_work_back(candidates, work)
            else:
                rejected_non_selected = 0

            logger.add_summary(
                "Guide selection",
                {
                    "selection_mode": selection,
                    "rejected_non_selected": rejected_non_selected,
                    "currently_retained_guides": int(retained_mask(candidates).sum()),
                    "total_tested_guides": len(candidates),
                },
            )

        # ------------------------------------------------------------
        # Validation primer design
        # ------------------------------------------------------------
        # Primer3 is relatively expensive, so design primers only for guides that
        # remain retained after all rejection/selection stages.
        work = candidates.loc[retained_mask(candidates)].copy()

        if not work.empty:
            work = design_validation_primers(
                work,
                fasta_dict,
                show_progress=True,
            )
            candidates = _merge_work_back(candidates, work)

        # ------------------------------------------------------------
        # No-PAM/no-guide placeholder rows
        # ------------------------------------------------------------
        # This adds rows only for termini that produced no scanned guide at all.
        # It should not be used to represent guides rejected later by filters.
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

    except Exception as exc:
        logger.failure(exc)
        raise

    finally:
        txt_log, json_log = logger.write()
        print(f"Log written to: {txt_log}")
        print(f"JSON log written to: {json_log}")