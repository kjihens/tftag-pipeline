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

def build_strain_vcf_dict(strain_vcf_specs: list[tuple[str, str]]) -> dict[str, str]:
    strain_vcfs: dict[str, str] = {}
    for name, path in strain_vcf_specs:
        if name in strain_vcfs:
            raise argparse.ArgumentTypeError(
                f"Duplicate strain name in --strain_vcf: {name!r}"
            )
        strain_vcfs[name] = path
    return strain_vcfs


def build_chrom_to_strain_map(strain_group_specs: list[tuple[str, list[str]]],
                              known_strain_names: set[str]) -> dict[str, str]:
    chrom_to_strain: dict[str, str] = {}

    for strain_name, chroms in strain_group_specs:
        if strain_name not in known_strain_names:
            raise argparse.ArgumentTypeError(
                f"--strain_group refers to unknown strain {strain_name!r}. "
                f"Define it first with --strain_vcf."
            )
        for chrom in chroms:
            if chrom in chrom_to_strain and chrom_to_strain[chrom] != strain_name:
                raise argparse.ArgumentTypeError(
                    f"Chromosome {chrom!r} assigned to multiple strains: "
                    f"{chrom_to_strain[chrom]!r} and {strain_name!r}"
                )
            chrom_to_strain[chrom] = strain_name

    return chrom_to_strain
