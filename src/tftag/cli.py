"""
CLI for TFTag pipeline.
"""
from __future__ import annotations
import argparse
from .pipeline import run_pipeline
from .cli_parsers import (
    parse_name_path,
    parse_name_chroms,
    build_stock_vcf_dict,
    build_chrom_to_stock_map,
)


def main():
    ap = argparse.ArgumentParser(description="TFTag genome-wide CRISPR design",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                 )

    # Inputs
    ap.add_argument("--gtf", required=True, help="GTF/GFF annotation file.")
    ap.add_argument("--db", default="gtf.db", help="Path to gffutils sqlite DB (created if absent).")
    ap.add_argument("--fasta", required=True, help="Genome FASTA (used for sequence extraction).")
    ap.add_argument("--genes", default=None, help="Either (a) path to a file with one gene ID per line, or (b) comma-separated gene IDs.")
    ap.add_argument("--stock-vcf", action="append", type=parse_name_path, default=[], metavar="NAME=PATH", help=("Named stock VCF. Repeatable. Example: "
            "--stock_vcf attP40=Cas9_on_3.mapping.vcf"
        )
    )
    
    # Outputs
    ap.add_argument("--basename", default="tftag", help="Base name for output files (DB, CSV, Parquet).")
    ap.add_argument("--outdir", default="out", help="Output directory.")
    ap.add_argument("--table", default="guides", help="Output table name in SQLite DB.")
    ap.add_argument("--write-csv", dest="write_csv", action="store_true", help="Also write output table as CSV.")
    ap.add_argument("--no-write-csv", dest="write_csv", action="store_false", help="Do not write output table as CSV.")
    ap.set_defaults(write_csv=True)
    ap.add_argument("--write-parquet", dest="write_parquet", action="store_true", help="Also write output table as Parquet.")
    ap.add_argument("--no-write-parquet", dest="write_parquet", action="store_false", help="Do not write output table as Parquet.")
    ap.set_defaults(write_parquet=False)


    # Guide scanning
    ap.add_argument("--pam-up", type=int, default=30, help="bp upstream of codon to scan for PAMs.")
    ap.add_argument("--pam-down", type=int, default=30, help="bp downstream of codon to scan for PAMs.")
    ap.add_argument("--protospacer-overlap-len", type=int, default=13, help="Minimum number of PAM-proximal protospacer bases that must overlap a Homology arm, together with the PAM, before donor mutation is required.")
    
    # RS3 efficiency scoring
    ap.add_argument("--tracr", default="Hsu2013", help="RS3 tracrRNA model identifier.")
    ap.add_argument("--rs3-batch", type=int, default=2048, help="Batch size for RS3 inference.")

    # Off-target enumeration (Cas-OFFinder)
    ap.add_argument("--offtargets", dest="offtargets", action="store_true", help="Run Cas-OFFinder (default)")
    ap.add_argument("--no-offtargets", dest="offtargets", action="store_false", help="Skip Cas-OFFinder")
    ap.set_defaults(offtargets=True)
    ap.add_argument("--cas-offinder", default="cas-offinder")
    ap.add_argument("--device", default="C")
    ap.add_argument("--mismatches", type=int, default=4)
    ap.add_argument("--pam-pattern", default="NNNNNNNNNNNNNNNNNNNNNGG")
    ap.add_argument("--min-offtarget-mismatch", type=int, default=0, help=("Optional off-target filter. 0 disables filtering. "
                                                                           "2 keeps guides with no 1-mismatch off-targets; "
                                                                           "3 also excludes 2-mismatch off-targets; "
                                                                           "-1 requires strict uniqueness within searched mismatch space.")
    )
                                                                                                                                                
    ap.add_argument("--offtarget-batch-size", type=int, default=500, help="Number of unique spacers per Cas-OFFinder batch.")

    # variant checking
    ap.add_argument("--stock-group", action="append", type=parse_name_chroms, default=[], metavar="NAME=CHR1,CHR2,...", help=("Assign chromosomes to a named stock. Repeatable. Example: "
            "--stock-group attP40=3L,3R --stock-group attP2=2L,2R,X,4,Y"
        )
    )
    ap.add_argument("--check-stock-vcf-compatibility", action="store_true", help=(
            "Validate that each supplied stock VCF matches the reference FASTA "
            "(contigs, lengths, REF alleles at sampled/all checked variant sites)."
        )
    )
    ap.add_argument("--check-stock-variants", action="store_true", help=(
            "Annotate candidate guides against the relevant stock VCF for the chromosome."
        )
    )

    # Codon Table usage
    ap.add_argument(
        "--codon-choice",
        choices=["gc", "usage"],
        default="gc",
        help="Synonymous codon choice strategy for blocking edits.",
    )

    ap.add_argument(
        "--codon-usage-table",
        default=None,
        help="Path to cached codon-usage TSV. Used when --codon-choice usage.",
    )

    ap.add_argument(
        "--force-recompute-codon-usage",
        action="store_true",
        help="Recompute codon-usage table even if cache exists.",
    )

    # Selection
    ap.add_argument("--selection", choices=["all", "closest", "rs3", "score"], default="all", help="Return all guides, or select one guide per gene per terminus (start/stop)")
    ap.add_argument("--stock-identical-only", action="store_true", help=("Guide selection mode. all keeps all guides; closest, rs3, and score "
                                                                         "select one guide per gene terminus. score uses composite ranking."
                                                                         )
    )

    args = ap.parse_args()

    stock_vcfs = build_stock_vcf_dict(args.stock_vcf)
    chrom_to_stock = build_chrom_to_stock_map(
        args.stock_group,
        known_stock_names=set(stock_vcfs.keys())
    )

    run_pipeline(
        gtf_file=args.gtf,
        gtf_db_path=args.db,
        genome_fasta_path=args.fasta,
        genes=args.genes,
        basename=args.basename,
        outdir=args.outdir,
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
        min_offtarget_mismatch=args.min_offtarget_mismatch,
        offtarget_batch_size=args.offtarget_batch_size,
        selection=args.selection,
        protospacer_overlap_len=args.protospacer_overlap_len,
        write_csv=args.write_csv,
        write_parquet=args.write_parquet,
        stock_vcfs=stock_vcfs,
        chrom_to_stock=chrom_to_stock,
        check_stock_vcf_compatibility=args.check_stock_vcf_compatibility,
        check_stock_variants=args.check_stock_variants,
        stock_identical_only=args.stock_identical_only,
        codon_choice=args.codon_choice,
        codon_usage_table=args.codon_usage_table,
        force_recompute_codon_usage=args.force_recompute_codon_usage,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())