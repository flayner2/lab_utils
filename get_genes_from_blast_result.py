from collections import defaultdict
import csv
import os
import sys

from Bio import SeqIO


def _shorten_str(string: str, perc: float = 0.2) -> str:
    assert 0 < perc <= 1.0, "Provide a percentage between 0 and 1"

    max_len = round(len(string) * perc)
    return string[:max_len]


def get_user_choices(blast_hits: dict[str, list]) -> dict[str, list]:
    user_choices = defaultdict(list)

    for gene, all_hits in blast_hits.items():
        print(f"Gene: {gene}")

        for hit_collection in all_hits:
            for genome, hits in hit_collection.items():
                print(f"Genome: {genome}\n")
                print("\tquery\tsubject\tidentity\talign_length\tevalue\n")

                for pos, hit in enumerate(hits):
                    query = _shorten_str(hit["query"], 0.3)
                    subject = _shorten_str(hit["subject"], 1)
                    print(
                        (
                            f"{pos}\t{query}...\t{subject}\t"
                            f"{hit['identity']}\t{hit['align_len']}"
                            f"\t{hit['evalue']}\n"
                        )
                    )

                while True:
                    choice = input(
                        "Choose one by providing a valid index or type 'quit' to exit: "
                    )

                    if choice.lower() == "quit":
                        sys.exit(0)
                    elif choice.isnumeric() and int(choice) in range(0, len(hits)):
                        choice = int(choice)

                        try:
                            user_choices[gene].append(hits[choice])
                        except IndexError:
                            print("Invalid index, try again...\n")
                            continue
                    else:
                        print(
                            (
                                "Invalid option. Provide a valid index or type "
                                "'quit' to exit.\n"
                            )
                        )
                        continue

    return user_choices


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
                hits = defaultdict(list)

                for row in tsv_reader:
                    each_hit = {}
                    each_hit["query"] = row[0]
                    each_hit["subject"] = row[1]
                    each_hit["identity"] = float(row[2])
                    each_hit["align_len"] = int(row[3])
                    each_hit["evalue"] = float(row[10])

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
    choices = get_user_choices(blast_results)


if __name__ == "__main__":
    main()
