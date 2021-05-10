import csv
import math
import os
import sys

from Bio import SeqIO


def main() -> None:
    assert len(sys.argv) > 1, (
        "Run this script with: python3 get_genes_from_blast_result.py "
        "[path_to_blast_files] [path_to_fasta_files] [output_directory]"
    )

    blast_files_path = sys.argv[1]
    fasta_files_path = sys.argv[2]
    outdir = sys.argv[3]

    for fasta_file in os.listdir(fasta_files_path):
        species_id = fasta_file.split(".fna")[0]
        sequences = SeqIO.to_dict(
            SeqIO.parse(os.path.join(fasta_files_path, fasta_file), "fasta")
        )

        for blast_file in os.listdir(blast_files_path):
            if species_id in blast_file:
                with open(os.path.join(blast_files_path, blast_file), "r") as csv_file:
                    blast_data = csv.reader(csv_file, delimiter="\t")

                    gene_name = blast_file.split(".fna")[-1].split(".blastn")[0]

                    gene_file_name = f"{species_id}{gene_name}.fasta"

                    for each_hit in blast_data:
                        if math.floor(float(each_hit[2])) == 100:
                            new_record = sequences[each_hit[0]]
                            start, end = int(each_hit[8]), int(each_hit[9])
                            new_record.seq = new_record.seq[start - 1 : end]
                            new_record.description = (
                                f"{new_record.description} {start}:{end}"
                            )

                            with open(
                                os.path.join(outdir, gene_file_name), "w+"
                            ) as outfile:
                                SeqIO.write(
                                    sequences=new_record, handle=outfile, format="fasta"
                                )

                            break


if __name__ == "__main__":
    main()
