import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
import os

def get_gene_coordinates(gene_list_file, gff_file):
    """
    Get the coordinates of start and stop codons for genes in the gene list.
    :param gene_list_file: Path to the file containing gene names (one per line).
    :param gff_file: Path to the GFF file.
    :return: List of dictionaries with gene name, feature type, chromosome, start, end, and strand.
    """
    # Load the GFF database
    if not os.path.exists(f"{gff_file}.db"):
        db = gffutils.create_db(gff_file, dbfn=f"{gff_file}.db", force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)
    else:
        db = gffutils.FeatureDB(f"{gff_file}.db", keep_order=True)

    print("GFF database loaded.")

    # Read the gene list
    with open(gene_list_file, "r") as f:
        gene_names = [line.strip() for line in f]

    results = []
    for gene_name in gene_names:
        try:
            gene = db[gene_name]
            for feature in db.children(gene, featuretype=["start_codon", "stop_codon"]):
                results.append({
                    "gene_name": gene_name,
                    "feature": feature.featuretype,
                    "chromosome": feature.seqid,
                    "start": feature.start,
                    "end": feature.end,
                    "strand": feature.strand
                })
        except gffutils.exceptions.FeatureNotFoundError:
            print(f"Gene {gene_name} not found in the GFF file.")
    return results

def get_sequences_with_context(coordinates, fasta_file, upstream=30, downstream=30):
    """
    Get sequences 30 bp upstream and downstream of start/stop codons.
    :param coordinates: List of dictionaries with genomic coordinates.
    :param fasta_file: Path to the FASTA file.
    :param upstream: Number of bases upstream to include.
    :param downstream: Number of bases downstream to include.
    :return: List of dictionaries with sequence information.
    """
    # Load the FASTA file
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    results = []
    for coord in coordinates:
        chromosome = coord["chromosome"]
        strand = coord["strand"]
        start = coord["start"]
        end = coord["end"]

        # Ensure the chromosome exists in the FASTA file
        if chromosome not in fasta_sequences:
            print(f"Chromosome {chromosome} not found in the FASTA file.")
            continue

        # Get the sequence with upstream and downstream context
        chromosome_seq = fasta_sequences[chromosome].seq
        if strand == "+":
            seq_start = max(0, start - upstream - 1)  # Convert to 0-based index
            seq_end = end + downstream
        else:
            seq_start = max(0, start - downstream - 1)  # Convert to 0-based index
            seq_end = end + upstream

        sequence = chromosome_seq[seq_start:seq_end]
        if strand == "-":
            sequence = sequence.reverse_complement()

        results.append({
            "gene_name": coord["gene_name"],
            "feature": coord["feature"],
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "strand": strand,
            "sequence": str(sequence)
        })
    return results

def find_23bp_sequences(sequences):
    """
    Identify 23 bp subsequences ending with 'GG' and mark as + strand,
    or starting with 'CC' and mark as - strand.
    :param sequences: List of dictionaries with sequence information.
    :return: List of dictionaries with genomic information of identified subsequences.
    """
    results = []
    for seq_info in sequences:
        sequence = seq_info["sequence"]
        chromosome = seq_info["chromosome"]
        start = seq_info["start"]
        strand = seq_info["strand"]

        # Identify subsequences ending with 'GG' (mark as + strand)
        for i in range(len(sequence) - 22):
            subseq = sequence[i:i + 23]
            if subseq.endswith("GG"):
                results.append({
                    "chromosome": chromosome,
                    "start": start + i,
                    "end": start + i + 22,
                    "strand": "+",
                    "sequence": subseq
                })

        # Identify subsequences starting with 'CC' (mark as - strand)
        for i in range(len(sequence) - 22):
            subseq = sequence[i:i + 23]
            if subseq.startswith("CC"):
                results.append({
                    "chromosome": chromosome,
                    "start": start + i,
                    "end": start + i + 22,
                    "strand": "-",
                    "sequence": subseq
                })
    return results

def write_gff(output_file, subsequences):
    """
    Write the identified subsequences to a GFF file.
    :param output_file: Path to the output GFF file.
    :param subsequences: List of dictionaries with genomic information of identified subsequences.
    """
    with open(output_file, "w") as f:
        for subseq in subsequences:
            f.write(f"{subseq['chromosome']}\tcustom\t23bp_sequence\t{subseq['start']}\t{subseq['end']}\t.\t{subseq['strand']}\t.\tsequence={subseq['sequence']}\n")

# Main script
if __name__ == "__main__":
    gene_list_file = "testlist.txt"
    gff_file = "genome_files/dmel-all-r6.63.gff"
    fasta_file = "genome_files/dmel-all-chromosome-r6.63.fasta"
    output_file = "termini_gRNAs.gff"

    # Step 1: Get gene coordinates
    coordinates = get_gene_coordinates(gene_list_file, gff_file)
    print(coordinates.head())  # Print the first few coordinates for verification

    # Step 2: Get sequences with upstream and downstream context
    sequences = get_sequences_with_context(coordinates, fasta_file)

    # Step 3: Find 23 bp subsequences
    subsequences = find_23bp_sequences(sequences)

    # Step 4: Write to GFF file
    write_gff(output_file, subsequences)

    print(f"Output written to {output_file}")