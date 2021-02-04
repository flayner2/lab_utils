import sys
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO


def perform_blast(
    query_path: str, min_e: float, subject_path: str, out_name: str
) -> None:
    """Runs a local version of Blast on a input fasta file containing a set of sequences,
    only returning those with expect < min_e. Produces a Blast alignment xml file.

    Args:
        query_path (str): path to the query fasta file (input file)
        min_e (float): the minimun e-value thresshold
        subject_path (str): path to the subject fasta file (to compare against)
        out_name (str): name of the Blast alignment output file
    """
    cline = NcbiblastnCommandline(
        query=query_path, subject=subject_path, outfmt=5, out=out_name, evalue=min_e
    )
    _, _ = cline()


def check_expect(blast_output: str, min_e: float) -> list:
    """Checks a blast output file for alignments with e-value < min_e, returning a list
    of sequence IDs for the sequences that pass the condition

    Args:
        blast_output (str): path to the blast xml output file
        min_e (float): the minimum e-value thresshold

    Returns:
        list: a list of sequence IDs for sequences with e-value < min_e
    """
    ids_list = []

    with open(blast_output, "r") as blast_output:
        records = NCBIXML.parse(blast_output)

        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < min_e:
                        ids_list.append(record.query.split(" ")[0])

    return ids_list


def check_valid_seqs(input_path: str, ids_list: list, output_path: str) -> None:
    """Checks an input fasta file containing a collection of sequences for those that
    match the IDs in ids_list, writing to the output file only those which don't

    Args:
        input_path (str): path to the input fasta file
        ids_list (list): list of sequence IDs
        output_path (str): path to the output fasta file
    """
    with open(output_path, "a") as outfile:
        for seq in SeqIO.parse(input_path, "fasta"):
            if not (seq.id in ids_list):
                outfile.write(seq.format("fasta"))


def sanitize_args() -> None:
    """Checks if the user provided an argument before running"""
    if len(sys.argv) < 3:
        sys.exit("Run this script with the following args: [query_path] [subject_path]")


def main():
    # Change the min e-value here
    e_thresshold = 10 ** -10
    current_dir = os.getcwd()
    infile = os.path.join(current_dir, sys.argv[1])
    subject_path = os.path.join(current_dir, sys.argv[2])
    # Edit this string if you want to change the name of the Blast output
    blast_outname = f'{infile.split("/")[-1].split(".")[0]}_blast.xml'
    blast_outpath = os.path.join(current_dir, blast_outname)

    perform_blast(infile, e_thresshold, subject_path, blast_outname)
    ids_list = check_expect(blast_outpath, e_thresshold)

    # Edit this string if you want to change the name of the fasta output
    fasta_outname = f'{infile.split("/")[-1].split(".")[0]}_blast_filtered.fasta'
    fasta_outpath = os.path.join(current_dir, fasta_outname)
    check_valid_seqs(infile, ids_list, fasta_outpath)


if __name__ == "__main__":
    sanitize_args()
    main()