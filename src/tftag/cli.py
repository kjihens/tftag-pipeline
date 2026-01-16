"""
CLI for TFTag pipeline.
"""
from __future__ import annotations
import argparse
from .pipeline import run_pipeline


def main():
    ap = argparse.ArgumentParser(description="TFTag genome-wide CRISPR design")

    # Inputs
    ap.add_argument("--gtf", required=True, help="GTF/GFF annotation file.")
    ap.add_argument("--db", default="gtf.db", help="Path to gffutils sqlite DB (created if absent).")
    ap.add_argument("--fasta", required=True, help="Genome FASTA (used for sequence extraction).")
    ap.add_argument("--genes", default=None, help="Either (a) path to a file with one gene ID per line, or (b) comma-separated gene IDs.")
    
    # Outputs
    ap.add_argument("--outdb", default="out/tftag_guides.sqlite", help="Output SQLite database path.")
    ap.add_argument("--table", default="guides", help="Output table name in SQLite DB.")

    # Guide scanning
    ap.add_argument("--pam_up", type=int, default=30, help="bp upstream of codon to scan for PAMs.")
    ap.add_argument("--pam_down", type=int, default=30, help="bp downstream of codon to scan for PAMs.")
    ap.add_argument("--protospacer_overlap_len", type=int, default=0, help="Length of overlap (in bp) allowed between protospacer and coding sequence.")
    
    # RS3 efficiency scoring
    ap.add_argument("--tracr", default="Hsu2013", help="RS3 tracrRNA model identifier.")
    ap.add_argument("--rs3_batch", type=int, default=2048, help="Batch size for RS3 inference.")

    # Off-target enumeration (Cas-OFFinder)
    ap.add_argument("--offtargets", dest="offtargets", action="store_true", help="Run Cas-OFFinder (default)")
    ap.add_argument("--no-offtargets", dest="offtargets", action="store_false", help="Skip Cas-OFFinder")
    ap.set_defaults(offtargets=True)
    ap.add_argument("--cas-offinder", default="cas-offinder")
    ap.add_argument("--device", default="C")
    ap.add_argument("--mismatches", type=int, default=4)
    ap.add_argument("--pam-pattern", default="NNNNNNNNNNNNNNNNNNNNNGG")
    ap.add_argument("--no-offtarget-sites", action="store_true", help="Keep only guides with no off-target hits (mm1..mmN == 0 and mm0==1)")

    # Selection
    ap.add_argument("--per-tag", choices=["all", "closest", "rs3"], default="all", help="Return all guides, or select one guide per gene per terminus (start/stop)")

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
        batch_size_rs3=args.rs3_batch,
        run_offtargets=args.offtargets,
        cas_offinder_bin=args.cas_offinder,
        device_spec=args.device,
        mismatches=args.mismatches,
        pam_pattern=args.pam_pattern,
        no_offtarget_sites=args.no_offtarget_sites,
        per_terminus=args.per_tag,
        protospacer_overlap_len=args.protospacer_overlap_len,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())