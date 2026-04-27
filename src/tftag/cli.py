"""
CLI for TFTag pipeline.
"""
from __future__ import annotations
import argparse
from .pipeline import run_pipeline

def _parse_name_path(spec: str) -> tuple[str, str]:
    """
    Parse NAME=PATH, e.g.:
      attP40=data/Cas9_on_3.mapping.vcf
    """
    if "=" not in spec:
        raise argparse.ArgumentTypeError(
            f"Expected NAME=PATH, got {spec!r}"
        )
    name, path = spec.split("=", 1)
    name = name.strip()
    path = path.strip()

    if not name:
        raise argparse.ArgumentTypeError(
            f"Invalid NAME in {spec!r}"
        )
    if not path:
        raise argparse.ArgumentTypeError(
            f"Invalid PATH in {spec!r}"
        )
    return name, path

def _parse_name_chroms(spec: str) -> tuple[str, list[str]]:
    """
    Parse NAME=CHR1,CHR2,..., e.g.:
      attP40=3L,3R
    """
    if "=" not in spec:
        raise argparse.ArgumentTypeError(
            f"Expected NAME=CHR1,CHR2,..., got {spec!r}"
        )
    name, chroms = spec.split("=", 1)
    name = name.strip()
    chrom_list = [c.strip() for c in chroms.split(",") if c.strip()]

    if not name:
        raise argparse.ArgumentTypeError(
            f"Invalid NAME in {spec!r}"
        )
    if not chrom_list:
        raise argparse.ArgumentTypeError(
            f"No chromosomes provided in {spec!r}"
        )
    return name, chrom_list

def _build_stock_vcf_dict(stock_vcf_specs: list[tuple[str, str]]) -> dict[str, str]:
    stock_vcfs: dict[str, str] = {}
    for name, path in stock_vcf_specs:
        if name in stock_vcfs:
            raise argparse.ArgumentTypeError(
                f"Duplicate stock name in --stock_vcf: {name!r}"
            )
        stock_vcfs[name] = path
    return stock_vcfs


def _build_chrom_to_stock_map(stock_group_specs: list[tuple[str, list[str]]],
                              known_stock_names: set[str]) -> dict[str, str]:
    chrom_to_stock: dict[str, str] = {}

    for stock_name, chroms in stock_group_specs:
        if stock_name not in known_stock_names:
            raise argparse.ArgumentTypeError(
                f"--stock_group refers to unknown stock {stock_name!r}. "
                f"Define it first with --stock_vcf."
            )
        for chrom in chroms:
            if chrom in chrom_to_stock and chrom_to_stock[chrom] != stock_name:
                raise argparse.ArgumentTypeError(
                    f"Chromosome {chrom!r} assigned to multiple stocks: "
                    f"{chrom_to_stock[chrom]!r} and {stock_name!r}"
                )
            chrom_to_stock[chrom] = stock_name

    return chrom_to_stock

def main():
    ap = argparse.ArgumentParser(description="TFTag genome-wide CRISPR design")

    # Inputs
    ap.add_argument("--gtf", required=True, help="GTF/GFF annotation file.")
    ap.add_argument("--db", default="gtf.db", help="Path to gffutils sqlite DB (created if absent).")
    ap.add_argument("--fasta", required=True, help="Genome FASTA (used for sequence extraction).")
    ap.add_argument("--genes", default=None, help="Either (a) path to a file with one gene ID per line, or (b) comma-separated gene IDs.")
    ap.add_argument("--stock_vcf", action="append", type=_parse_name_path, default=[], metavar="NAME=PATH", help=("Named stock VCF. Repeatable. Example: "
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
    ap.add_argument("--pam_up", type=int, default=30, help="bp upstream of codon to scan for PAMs.")
    ap.add_argument("--pam_down", type=int, default=30, help="bp downstream of codon to scan for PAMs.")
    ap.add_argument("--protospacer_overlap_len", type=int, default=13, help="Length of overlap (in bp) allowed between protospacer and coding sequence.")
    
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
    ap.add_argument("--min-offtarget-mismatch", type=int, default=0, help="Minimum mismatch count allowed for off-target sites (e.g. 2 = allow ≥2 mismatches only).")
    ap.add_argument("--offtarget_batch_size", type=int, default=500, help="Number of unique spacers per Cas-OFFinder batch.")

    # variant checking
    ap.add_argument("--stock_group", action="append", type=_parse_name_chroms, default=[], metavar="NAME=CHR1,CHR2,...", help=("Assign chromosomes to a named stock. Repeatable. Example: "
            "--stock_group attP40=3L,3R --stock_group attP2=2L,2R,X,4,Y"
        )
    )
    ap.add_argument("--check_stock_vcf_compatibility", action="store_true", help=(
            "Validate that each supplied stock VCF matches the reference FASTA "
            "(contigs, lengths, REF alleles at sampled/all checked variant sites)."
        )
    )
    ap.add_argument("--check_stock_variants", action="store_true", help=(
            "Annotate candidate guides against the relevant stock VCF for the chromosome."
        )
    )
    # Selection
    ap.add_argument("--selection", choices=["all", "closest", "rs3", "score"], default="all", help="Return all guides, or select one guide per gene per terminus (start/stop)")
    ap.add_argument("--stock_identical_only", action="store_true", help="Keep only guides whose full 23-mer target sequence is identical between reference and assigned stock.")

    args = ap.parse_args()

    stock_vcfs = _build_stock_vcf_dict(args.stock_vcf)
    chrom_to_stock = _build_chrom_to_stock_map(
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
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())