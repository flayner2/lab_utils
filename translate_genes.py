import os
import sys
from typing import Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_seqs(fasta_file: str) -> Iterator[SeqRecord]:
    return SeqIO.parse(fasta_file, format="fasta")


def translate(fasta_file: str, output_dir: str) -> None:
    sequences = load_seqs(fasta_file)

    output_name = f"{os.path.basename(fasta_file).split('.fas')[0]}_translated.fas"
    output_file = os.path.join(output_dir, output_name)

    with open(output_file, "w+") as outfile:
        for record in sequences:
            if len(record.seq) % 3 == 0:
                record.seq = record.seq.translate()
                outfile.write(record.format("fasta"))


def main():
    assert (
        len(sys.argv) >= 3
    ), "Run this with: python3 translate_genes.py {input_file} {output_path}"

    fasta_file = sys.argv[1]
    output_dir = sys.argv[2]

    translate(fasta_file, output_dir)


if __name__ == "__main__":
    main()
