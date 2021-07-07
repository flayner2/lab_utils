from collections import defaultdict
import csv
import os
import sys

from Bio import SeqIO


def load_blast_results(path: str, genes: list) -> dict[str, list]:
    blast_results = defaultdict(list)

    for each_file in os.listdir(path):
        if each_file.endswith(".tsv"):
            file_path = os.path.join(path, each_file)
            genome_name = each_file.split(".gbff")[0]
            gene_name = "not_found"

            for gene in genes:
                if gene in each_file:
                    gene_name = gene
                    break

            with open(file_path, "r") as tsv_file:
                tsv_reader = csv.reader(tsv_file, delimiter="\t")

                for row in tsv_reader:
                    hits = defaultdict(list)

                    each_hit = {}
                    each_hit["query"] = row[0]
                    each_hit["subject"] = row[1]
                    each_hit["identity"] = row[2]
                    each_hit["align_len"] = row[3]
                    each_hit["evalue"] = row[10]

                    hits[genome_name].append(each_hit)
                    blast_results[gene_name].append(hits)

    return blast_results


def get_gene_names(path: str) -> list[str]:
    gene_names = []

    for each_file in os.listdir(path):
        if each_file.endswith((".fasta", ".fa", ".fna", ".fas")):
            gene = os.path.splitext(each_file)[0]
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
    blast_results = load_blast_results(blast_files_path, gene_names)


if __name__ == "__main__":
    main()
