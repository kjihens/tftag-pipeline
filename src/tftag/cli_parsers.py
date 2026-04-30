import argparse

def parse_name_path(spec: str) -> tuple[str, str]:
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

def parse_name_chroms(spec: str) -> tuple[str, list[str]]:
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

def build_stock_vcf_dict(stock_vcf_specs: list[tuple[str, str]]) -> dict[str, str]:
    stock_vcfs: dict[str, str] = {}
    for name, path in stock_vcf_specs:
        if name in stock_vcfs:
            raise argparse.ArgumentTypeError(
                f"Duplicate stock name in --stock_vcf: {name!r}"
            )
        stock_vcfs[name] = path
    return stock_vcfs


def build_chrom_to_stock_map(stock_group_specs: list[tuple[str, list[str]]],
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
