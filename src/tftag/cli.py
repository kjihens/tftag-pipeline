"""
CLI for TFTag pipeline.
"""
from __future__ import annotations
import argparse
from .pipeline import run_pipeline


def main():
    ap = argparse.ArgumentParser(description="TFTag genome-wide CRISPR design")

    ap.add_argument("--gtf", required=True, help="GTF/GFF annotation file.")
    ap.add_argument("--db", default="gtf.db", help="Path to gffutils sqlite DB (created if absent).")
    ap.add_argument("--fasta", required=True, help="Genome FASTA (used for sequence extraction).")
    ap.add_argument("--genes", default=None, help="Either (a) path to a file with one gene ID per line, or (b) comma-separated gene IDs.")
    ap.add_argument("--outdb", default="out/tftag_guides.sqlite", help="Output SQLite database path.")
    ap.add_argument("--table", default="guides", help="Output table name in SQLite DB.")
    ap.add_argument("--pam_up", type=int, default=30, help="bp upstream of codon to scan for PAMs.")
    ap.add_argument("--pam_down", type=int, default=30, help="bp downstream of codon to scan for PAMs.")
    ap.add_argument("--tracr", default="Hsu2013", help="RS3 tracrRNA model identifier.")
    ap.add_argument("--rs3_batch", type=int, default=2048, help="Batch size for RS3 inference.")
    ap.add_argument("--specificity", action="store_true", help="Run Cas-OFFinder and add mismatch hit counts.")
    ap.add_argument("--cas", default="cas-offinder", help="Cas-OFFinder binary.")
    ap.add_argument("--device", default="C", help="Cas-OFFinder device spec (build-dependent; e.g. C).")
    ap.add_argument("--cclmoff_cmd", default=None, help="Command template to run CCLMoff; uses {pairs} and {output}.")
    ap.add_argument("--cclmoff_pairs", default="out/cclmoff_pairs.tsv", help="CCLMoff pairs TSV path.")
    ap.add_argument("--cclmoff_preds", default="out/cclmoff_preds.tsv", help="CCLMoff predictions TSV path.")
    ap.add_argument("--cclmoff_agg", default="max", choices=["max", "sum"], help="Primary CCLMoff aggregation.")
    ap.add_argument("--protospacer_overlap_len", type=int, default=13, help="Overlap length between PAM-proximal protospacer and arm to decide which arm to mutate.")

    args = ap.parse_args()

    run_pipeline(
        gtf_file=args.gtf,
        gtf_db_path=args.db,
        genome_fasta_path=args.fasta,
        genes=args.genes,
        output_db_path=args.outdb,
        output_table=args.table,
        pam_window_up=args.pam_up,
        pam_window_down=args.pam_down,
        tracrRNA=args.tracr,
        do_specificity=args.specificity,
        cas_offinder_bin=args.cas,
        device_spec=args.device,
        batch_size_rs3=args.rs3_batch,
        cclmoff_cmd=args.cclmoff_cmd,
        cclmoff_pairs_path=args.cclmoff_pairs,
        cclmoff_preds_path=args.cclmoff_preds,
        cclmoff_agg_method=args.cclmoff_agg,
        protospacer_overlap_len=args.protospacer_overlap_len,
    )


if __name__ == "__main__":
    main()