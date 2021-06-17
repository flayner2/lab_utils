import os
import sys

from Bio import SeqIO

"""
Converts a set of .gbff files to .fna files using Biopython.
Example run inside LGM: `python3 gbff_to_fna.py`
"""


def main() -> None:
    assert (
        len(sys.argv) >= 3
    ), "Run this with: python3 gbff_to_fna.py {path_to_inputs} {path_to_outputs}"

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    for filename in os.listdir(input_path):
        if filename.endswith(".gbff"):
            in_name = os.path.join(input_path, filename)
            outname = os.path.join(output_path, filename.replace(".gbff", ".fna"))
            print(f"Converting the file {str(filename)}")
            SeqIO.convert(in_name, "genbank", outname, "fasta")


if __name__ == "__main__":
    main()
