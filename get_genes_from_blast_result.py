import csv
import os
import sys

from Bio import SeqIO


def load_blast_results(path: str) -> dict:
    pass


def get_gene_names(path: str) -> list[str]:
    gene_names = []

    for each_file in os.listdir(path):
        if each_file.endswith((".fasta", ".fa", ".fna", ".fas")):
            gene = os.path.splitext(os.path.basename(each_file))[0]
            gene_names.append(gene)

    return gene_names


def main() -> None:
    assert len(sys.argv) == 4, (
        "Run this script with: python3 get_genes_from_blast_result.py "
        "[path_to_blast_files] [path_to_fasta_files] [output_directory]"
    )

    blast_files_path = sys.argv[1]
    fasta_files_path = sys.argv[2]
    # outdir = sys.argv[3]

    gene_names = get_gene_names(fasta_files_path)
    blast_results = load_blast_results(blast_files_path)


if __name__ == "__main__":
    main()
