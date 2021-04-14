"""
Extracts a set of genes from a set of input genbank annotation files (GBFF),
generating a fasta file for each gene for all species or a file for each species
with all genes.
"""
import argparse
from collections import defaultdict
import os
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_text_file(path: str) -> list[str]:
    """Loads a line-separated file containing one or more strings to a list.

    Arguments:
        path (str): the path to the input file.
    Returns:
        list[str]: a list of the strings contained in the input file.
    """
    assert os.path.exists(path), "Please provide a valid path to an existing file"
    assert os.path.isfile(path), "Please, provide a path to a file, not a directory"

    with open(path, "r") as text_file:
        strings = [each.strip() for each in text_file.readlines()]

    return strings


def find_genes_in_annotation() -> tuple(str, list[SeqRecord]):
    pass


def find_genes_wrapper(
    gbff_folder: str, genes: list[str], exclude: list = []
) -> dict[str, list[SeqRecord]]:
    """Wrapper to load gbff files and find the wanted genes in them, while being aware
    of the strings to exclude.

    Arguments:
        gbff_folder (str): path to the directory containing one or more gbff files.
        genes (list[str]): a list of strings to search in the annotation.
        exclude (list): a list of strings to exclude from the search. Defaults to "[]".
    Returns:
        dict[str, list[SeqRecord]]: a dictionary where keys are the name of the species
        to whom the annotation belongs, and values are a list of SeqRecord objects that
        match any string from the list of genes and that do not contain any strings
        from the list of exclusions.
    """

    results = defaultdict(list)

    for gbff_filename in os.listdir(gbff_folder):
        if gbff_filename.endswith(".gbff"):
            gbff_fullpath = os.path.join(gbff_folder, gbff_filename)

            records = SeqIO.parse(gbff_fullpath, format="genbank")

            species, found_genes = find_genes_in_annotation(
                records=records, genes=genes, exclude=exclude
            )

            results[species].extend(found_genes)

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Retrieves a set of genes from a set of Genbank annotation"
            " files (.gbff files)"
        )
    )

    parser.add_argument(
        "gbff_path",
        metavar="G",
        type=str,
        help="path to the input genbank annotation files",
    )
    parser.add_argument(
        "genes_list",
        metavar="L",
        type=str,
        help=(
            "path to a text file with a list of strings to search the annotations"
            ", line-separated"
        ),
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        nargs="?",
        default=sys.stdout,
        help="name for the output file (by default outputs to STDOUT)",
    )
    parser.add_argument(
        "-e",
        "--exclude",
        type=str,
        help=(
            "path to a line-separated text file containing strings to exclude"
            " from the search"
        ),
    )

    args = parser.parse_args()

    try:
        genes = load_text_file(args.genes_list)

        if args.exclude:
            exclusions = load_text_file(args.exclude)
        else:
            exclusions = []

        found_genes = find_genes_wrapper(
            gbff_folder=args.gbff_path, genes=genes, exclusions=exclusions
        )
    except Exception as err:
        raise err


if __name__ == "__main__":
    main()
