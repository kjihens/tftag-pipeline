"""
End-to-end orchestrator for TFTag.

High-level flow
--------------
1) Build/start from codon features (start/stop codons) via annotate.build_attribute_table
2) Scan PAMs around each codon via scan.scan_for_guides
3) Prefilter loci that cannot support full arms/primer windows via design.prefilter_designable
4) Compute RS3 on-target scores via efficiency.score_rs3
5) Enumerate off-targets (Cas-OFFinder) and merge summary; optional filtering
6) Optionally select one guide per (gene_id, tag) before expensive design steps
7) Design steps: homology arms, choose edit arm, apply silent edits, validation primers
8) Write SQLite (+ optional CSV/Parquet snapshots)

Conventions
-----------
- Genomic coordinates: 1-based inclusive
- `mismatches` defines Cas-OFFinder maximum mismatches searched (and thus the n_mm* columns produced)
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
    basename: str = "tftag",
    outdir: str = "out",
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
    # For SpCas9 NGG: 20N + GG (pattern length must match query length)
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG",
    # Filtering: None/0 = no filter; 2 => require n_mm0==1 and n_mm1==0; -1 => strict uniqueness (within searched mm space)
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
      - -1: strict uniqueness (mm0==1 and all other n_mm* == 0)
        Note: “strict” is only as strict as the `mismatches` you searched (e.g. mismatches=4).
    """
    if per_tag not in ("all", "closest", "rs3"):
        raise ValueError("per_tag must be one of: all, closest, rs3")
    if mismatches < 0:
        raise ValueError("mismatches must be >= 0")

    if min_offtarget_mismatch is not None:
        min_offtarget_mismatch = int(min_offtarget_mismatch)
        if min_offtarget_mismatch < -1:
            raise ValueError("min_offtarget_mismatch must be -1, 0, a positive integer, or None")
        # If the user requests a threshold beyond what Cas-OFFinder computed, fail fast.
        # Example: min_offtarget_mismatch=6 requires at least n_mm1..n_mm5 columns.
        if min_offtarget_mismatch > 0 and (min_offtarget_mismatch - 1) > mismatches:
            raise ValueError(
                f"min_offtarget_mismatch={min_offtarget_mismatch} requires mismatches >= {min_offtarget_mismatch - 1}, "
                f"but mismatches={mismatches}."
            )

    # Ensure output directories exist
    os.makedirs(outdir, exist_ok=True)

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
            out_dir=out_dir,
            run_id=run_id,
            cas_offinder_bin=cas_offinder_bin,
            device_spec=device_spec,
            mismatches=mismatches,
            pam_pattern=pam_pattern,
        )
        candidates = merge_specificity(candidates, spec, mismatches=mismatches)

        # Apply mismatch-based filter (if requested)
        if min_offtarget_mismatch is not None and min_offtarget_mismatch != 0:
            candidates = filter_by_offtarget_mismatch(candidates, min_offtarget_mismatch)
            if candidates.empty:
                print("No guides remain after applying --min-offtarget-mismatch filter.")
                return

    else:
        # Populate expected columns so downstream code and selection can run.
        # (Selection should tolerate NaNs; if not, set conservative values here instead.)
        candidates["n_hits"] = np.nan
        for k in range(0, mismatches + 1):
            candidates[f"n_mm{k}"] = np.nan

        # If user requested an off-target filter but we didn't compute off-targets, stop early.
        if min_offtarget_mismatch is not None and min_offtarget_mismatch != 0:
            print("Error: --min-offtarget-mismatch requires off-target enumeration (run_offtargets=True).")
            return

    # If selecting one per tag, do it BEFORE expensive primer3 work.
    if per_tag != "all":
        candidates = select_one_per_tag(candidates, mode=per_tag)
        if candidates.empty:
            print("No guides remain after per-tag selection.")
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

    # If per_tag == "all", selection happens here (no-op) to keep output deterministic.
    # This also supports the case where per_tag == "all" but you still want the function to exist here.
    candidates = select_one_per_tag(candidates, mode=per_tag)

    # Add provenance
    candidates["run_id"] = run_id
    candidates["gtf_db_path"] = gtf_db_path
    candidates["genome_fasta_path"] = genome_fasta_path

    # Write SQLite
    to_sqlite(
        candidates,
        basename,
        outdir,
        output_table,
        if_exists="append",
        index=False,
        create_indices=True,
    )

    # Optional snapshots
    if write_parquet:
        candidates.to_parquet(os.path.join(outdir, basename + '.parquet'), index=False)
    if write_csv:
        candidates.to_csv(os.path.join(outdir, basename + '.csv'), index=False)

    # Console summary
    print(f"Finished. Wrote {len(candidates)} guides to:")
    print(f"   - {outdir + '/' + basename + '.sqlite'} (table: {output_table})")
    if write_parquet:
        print(f"   - {outdir + '/' + basename + '.parquet'}")
    if write_csv:
        print(f"   - {outdir + '/' + basename + '.csv'}")