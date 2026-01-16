"""
End-to-end orchestrator.
"""
from __future__ import annotations

import os
import uuid
import numpy as np
import pandas as pd

from . import annotate, scan, efficiency, design
from .genome_io import create_GTF_db, load_fasta_dict
from .storage import to_sqlite
from .select import select_one_per_tag
from .helpers import parse_genes_arg, filter_by_offtarget_mismatch
from .off_targeting import enumerate_offtargets_cas_offinder, merge_specificity


def run_pipeline(
    gtf_file: str,
    gtf_db_path: str,
    genome_fasta_path: str,
    genes: str | None = None,
    output_db_path: str = "out/tftag_guides.sqlite",
    output_table: str = "guides",
    pam_window_up: int = 30,
    pam_window_down: int = 30,
    tracrRNA: str = "Hsu2013",
    batch_size_rs3: int = 2048,
    protospacer_overlap_len: int = 13,
    # Off-targets (Cas-OFFinder only)
    run_offtargets: bool = True,
    cas_offinder_bin: str = "cas-offinder",
    device_spec: str = "C",
    mismatches: int = 4,
    # For SpCas9 NGG: 20N + GG (Cas-OFFinder uses pattern line length to define target length)
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG",
    # Filtering: 0/None = no filter; 2 => require n_mm0==1 and n_mm1==0; -1 => strict uniqueness
    min_offtarget_mismatch: int | None = 0,
    # Selection / outputs
    per_tag: str = "all",  # all|closest|rs3
    write_csv: bool = True,
    write_parquet: bool = True,
) -> None:
    """
    Run the TFTag pipeline.

    per_tag:
      - "all": keep all candidate guides for each (gene_id, tag)
      - "closest": keep one guide per (gene_id, tag) with smallest cut_distance
      - "rs3": keep one guide per (gene_id, tag) with highest rs3_score

    min_offtarget_mismatch:
      - None or 0: no off-target filtering
      - 2: keep only guides with n_mm0 == 1 and n_mm1 == 0
      - 3: additionally require n_mm2 == 0
      - -1: strict uniqueness (only mm0==1, all other n_mm* == 0)
    """
    if per_tag not in ("all", "closest", "rs3"):
        raise ValueError("per_tag must be one of: all, closest, rs3")
    if mismatches < 0:
        raise ValueError("mismatches must be >= 0")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_db_path) or ".", exist_ok=True)
    os.makedirs("out", exist_ok=True)

    run_id = uuid.uuid4().hex[:12]

    # Create/load GTF db
    db = create_GTF_db(gtf_file, gtf_db_path)

    # Load genome FASTA
    fasta_dict = load_fasta_dict(genome_fasta_path)

    # Resolve genes list
    genes_list = parse_genes_arg(genes)
    if genes_list is None:
        genes_list = [g.id for g in db.features_of_type("gene")]

    # Build codon attribute table
    attribute = annotate.build_attribute_table(genes_list, db)

    # Scan for candidate guides
    candidates = scan.scan_for_guides(
        attribute, fasta_dict, window_up=pam_window_up, window_down=pam_window_down
    )
    if candidates.empty:
        print("No candidate guides found.")
        return

    # Prefilter feasibility (design-aware)
    candidates = design.prefilter_designable(candidates, fasta_dict, show_progress=True)
    if "designable" in candidates.columns and not candidates["designable"].any():
        print("No designable guides after prefilter.")
        return

    # RS3 scoring
    candidates = efficiency.score_rs3(
        candidates, fasta_dict, tracrRNA=tracrRNA, batch_size=batch_size_rs3
    )

    # Off-target enumeration (Cas-OFFinder)
    if run_offtargets:
        _hits, spec = enumerate_offtargets_cas_offinder(
            candidates,
            genome_fasta_path,
            out_dir="out",
            run_id=run_id,
            cas_offinder_bin=cas_offinder_bin,
            device_spec=device_spec,
            mismatches=mismatches,
            pam_pattern=pam_pattern,
        )
        candidates = merge_specificity(candidates, spec, mismatches=mismatches)

        # Apply mismatch-based filter (if requested)
        if min_offtarget_mismatch is not None and int(min_offtarget_mismatch) != 0:
            candidates = filter_by_offtarget_mismatch(candidates, min_offtarget_mismatch)
            if candidates.empty:
                print("No guides remain after applying --min-offtarget-mismatch filter.")
                return

    else:
        # Populate expected columns so downstream code and selection can run
        candidates["n_hits"] = np.nan
        for k in range(0, mismatches + 1):
            candidates[f"n_mm{k}"] = np.nan

        # If user requested an off-target filter but we didn't compute off-targets, stop early
        if min_offtarget_mismatch is not None and int(min_offtarget_mismatch) != 0:
            print("Error: --min-offtarget-mismatch requires --run-offtargets (Cas-OFFinder).")
            return

    # Design steps
    candidates = design.add_homology_arms(candidates, fasta_dict, show_progress=True)

    candidates = design.choose_arm_for_mutation(
        candidates,
        protospacer_overlap_len=protospacer_overlap_len,
        coding_only=True,
        show_progress=True,
    )

    candidates = design.apply_silent_edits(candidates, show_progress=True)

    candidates = design.validation_primers(candidates, fasta_dict, show_progress=True)

    # Reduce to one guide per (gene_id, tag) if requested
    candidates = select_one_per_tag(candidates, mode=per_tag)

    # Add provenance
    candidates["run_id"] = run_id
    candidates["gtf_db_path"] = gtf_db_path
    candidates["genome_fasta_path"] = genome_fasta_path

    # Write SQLite
    to_sqlite(
        candidates,
        output_db_path,
        output_table,
        if_exists="append",
        index=False,
        create_indices=True,
    )

    # Optional snapshots
    base = os.path.splitext(output_db_path)[0]
    if write_parquet:
        parquet_path = base + ".parquet"
        os.makedirs(os.path.dirname(parquet_path) or ".", exist_ok=True)
        candidates.to_parquet(parquet_path, index=False)
    if write_csv:
        csv_path = base + ".csv"
        os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
        candidates.to_csv(csv_path, index=False)

    # Console summary
    print(f"Finished. Wrote {len(candidates)} guides to:")
    print(f"   - {output_db_path}")
    if write_parquet:
        print(f"   - {base}.parquet")
    if write_csv:
        print(f"   - {base}.csv")