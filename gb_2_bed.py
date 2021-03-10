from Bio import SeqIO
import os
import sys
from collections import defaultdict


def main() -> None:
    if len(sys.argv) < 3:
        sys.exit(
            "Run this script with `python gb_2_bed.py {input_directory} {output_directory}`"
        )

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    try:
        for gb_file in os.listdir(input_folder):
            if gb_file.endswith(".gbff"):
                path_to_file = os.path.join(input_folder, gb_file)
                organism = ""
                contigs = defaultdict(dict)

                for seq_record in SeqIO.parse(path_to_file, "genbank"):
                    if not organism:
                        organism = seq_record.annotations["organism"].replace(" ", "_")

                    contig_id = seq_record.id

                    for feature in seq_record.features:
                        if feature.type == "gene":
                            gene_id = feature.qualifiers["gene"][0]
                            gene_start, gene_end = (
                                feature.location.nofuzzy_start,
                                feature.location.nofuzzy_end,
                            )
                            contigs[contig_id][gene_id] = gene_start, gene_end

            outfile_name = f"{organism}_genome.bed"
            output_path = os.path.join(output_folder, outfile_name)

            with open(output_path, "a+") as outfile:
                for contig_name, gene_dict in contigs.items():
                    for gene_name, positions in gene_dict.items():
                        start, end = positions
                        to_write = f"{contig_name}\t{start}\t{end}\t{gene_name}\n"

                        outfile.write(to_write)

    except FileNotFoundError:
        sys.exit("The provided path is not valid.")


if __name__ == "__main__":
    main()
