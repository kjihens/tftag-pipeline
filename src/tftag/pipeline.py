"""
End-to-end orchestrator.
"""
from __future__ import annotations
import os, uuid
import numpy as np
import pandas as pd

from . import annotate, storage as tio, scan, efficiency, offtarget, cclmoff, design


def _parse_genes_arg(genes: str | None) -> list[str] | None:
    if genes is None:
        return None
    genes = str(genes).strip()
    if not genes:
        return None
    if os.path.exists(genes):
        with open(genes) as fh:
            return [ln.strip() for ln in fh if ln.strip()]
    # comma-separated list
    return [g.strip() for g in genes.split(",") if g.strip()]


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
    do_specificity: bool = False,
    cas_offinder_bin: str = "cas-offinder",
    device_spec: str = "C",
    batch_size_rs3: int = 2048,
    cclmoff_cmd: str | None = None,
    cclmoff_pairs_path: str = "out/cclmoff_pairs.tsv",
    cclmoff_preds_path: str = "out/cclmoff_preds.tsv",
    cclmoff_agg_method: str = "max",
    protospacer_overlap_len: int = 13,
) -> None:

    # Ensure output directories exist
    for p in [output_db_path, cclmoff_pairs_path, cclmoff_preds_path]:
        os.makedirs(os.path.dirname(p) or ".", exist_ok=True)

    run_id = uuid.uuid4().hex[:12]

    # Create/load GTF db
    db = tio.createGTFdb(gtf_file, gtf_db_path)

    # Load genome FASTA (dict backend)
    fasta_dict = tio.load_fasta_dict(genome_fasta_path)

    # Resolve genes list
    genes_list = _parse_genes_arg(genes)
    if genes_list is None:
        genes_list = [g.id for g in db.features_of_type("gene")]

    # Build codon attribute table
    attribute = annotate.build_attribute_table(genes_list, db)

    # Scan for candidate guides
    candidates = scan.scan_for_guides(attribute, fasta_dict, window_up=pam_window_up, window_down=pam_window_down)
    if candidates.empty:
        print("No candidate guides found.")
        return

    # Prefilter feasibility (design-aware)
    candidates = design.prefilter_designable(candidates, fasta_dict, show_progress=True)

    if not candidates["designable"].any():
        print("No designable guides after prefilter.")
        return

    # RS3 scoring
    candidates = efficiency.score_rs3(candidates, fasta_dict, tracrRNA=tracrRNA, batch_size=batch_size_rs3)

    # Off-target analysis
    hits = pd.DataFrame()
    if do_specificity:
        inp = offtarget.write_cas_offinder_input(candidates, genome_fasta_path, outfile=f"out/{run_id}_cas_input.txt")
        out = offtarget.run_cas_offinder(inp, cas_offinder_bin=cas_offinder_bin, device_spec=device_spec, output_file=f"out/{run_id}_cas_hits.txt")
        hits = offtarget.parse_cas_offinder_output(out)
        spec = offtarget.summarize_specificity(hits)

        # Uniqueness safety
        if not spec["spacer"].is_unique:
            raise RuntimeError("Specificity summary is not unique per spacer; check parser/aggregator.")

        candidates = candidates.merge(spec, on="spacer", how="left")
    else:
        for col in ["n_hits", "n_mm0", "n_mm1", "n_mm2", "n_mm3", "n_mm4"]:
            candidates[col] = np.nan

    # CCLMoff (optional)
    if cclmoff_cmd:
        if hits.empty:
            inp = offtarget.write_cas_offinder_input(candidates, genome_fasta_path, outfile=f"out/{run_id}_cas_input.txt")
            out = offtarget.run_cas_offinder(inp, cas_offinder_bin=cas_offinder_bin, device_spec=device_spec, output_file=f"out/{run_id}_cas_hits.txt")
            hits = offtarget.parse_cas_offinder_output(out)

        pairs_map = cclmoff.build_pairs_from_hits(candidates, hits, fasta_dict, outfile=cclmoff_pairs_path)
        preds_file = cclmoff.run_cclmoff_from_template(cclmoff_pairs_path, cclmoff_preds_path, cclmoff_cmd)
        preds = cclmoff.parse_cclmoff_output(preds_file)
        agg = cclmoff.aggregate(pairs_map, preds)

        if not agg["spacer"].is_unique:
            raise RuntimeError("CCLMoff aggregate is not unique per spacer; check pairing/aggregation.")

        candidates = candidates.merge(agg, on="spacer", how="left")

        # choose primary score if you want a single column
        if cclmoff_agg_method == "max" and "cclmoff_max" in candidates.columns:
            candidates["cclmoff_primary"] = candidates["cclmoff_max"]
        elif cclmoff_agg_method == "sum" and "cclmoff_sum" in candidates.columns:
            candidates["cclmoff_primary"] = candidates["cclmoff_sum"]

    # Homology arms
    candidates = design.add_homology_arms(candidates, fasta_dict, show_progress=True)

    # Decide which arm to mutate
    candidates = design.choose_arm_for_mutation(
        candidates, protospacer_overlap_len=protospacer_overlap_len, coding_only=True, show_progress=True
    )

    # Apply silent edits
    candidates = design.apply_silent_edits(candidates, show_progress=True)

    # Validation primers
    candidates = design.validation_primers(candidates, fasta_dict, show_progress=True)

    # Add run provenance
    candidates["run_id"] = run_id
    candidates["gtf_db_path"] = gtf_db_path
    candidates["genome_fasta_path"] = genome_fasta_path

    # Write outputs
    tio.to_sqlite(candidates, output_db_path, output_table, if_exists="append", index=False, create_indices=False)
    parquet_path = os.path.splitext(output_db_path)[0] + f".{run_id}.parquet"
    os.makedirs(os.path.dirname(parquet_path) or ".", exist_ok=True)
    candidates.to_parquet(parquet_path, index=False)

    print(f"Finished. Wrote {len(candidates)} guides to {output_db_path} and {parquet_path}. run_id={run_id}")