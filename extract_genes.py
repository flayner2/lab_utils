"""
Extracts a set of genes from a set of input genbank annotation files (GBFF),
generating a fasta file for each gene for all species or a file for each species
with all genes.
"""
import argparse
from collections import defaultdict
import os
import sys
from typing import Generator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_text_file(path: str) -> list[str]:
    """Loads a line-separated file containing one or more strings to a list.

    Arguments:
        path (str): the path to the input file.
    Returns:
        list[str]: a list of the strings contained in the input file.
    """
    assert os.path.exists(path), "Please provide a valid path to an existing file."
    assert os.path.isfile(path), "Please, provide a path to a file, not a directory."

    with open(path, "r") as text_file:
        strings = [each.strip() for each in text_file.readlines()]

    return strings


def find_genes_in_annotation(
    records: Generator, genes: list[str], exclude: list[str] = []
) -> tuple[str, list[SeqRecord]]:
    """Finds any term from a list of terms in the features of a SeqRecord object,
    excluding any term from a list of exclusions, and returns a tuple with the name of
    the originating species and a list of SeqRecord objects for all features that pass
    those checks.

    Arguments:
        records (Generator): a Generator object which yields SeqRecord objects. Usually
        this is the return of the `SeqIO.parse()` method from BioPython's Bio module.
        genes (list[str]): a list of strings to find in the SeqRecord's description.
        exclude (list[str]): an optional list of strings used to rule out any SeqRecord
        that contains any of them. Defaults to "[]".
    Returns:
        tuple[str, list[SeqRecord]]: a 2-tuple with the name of the organism and a list
        of SeqRecords that contain any term from `genes` and do not contain any term
        from `exclude`.
    """
    assert len(genes) > 0, "Received an empty list of genes. Check your input."

    first_record = next(records)
    species_name = first_record.annotations["organism"]
    first_record_features = [
        SeqRecord(
            seq=feature.extract(first_record.seq),
            id=feature.qualifiers["gene"][0],
            name=feature.qualifiers["gene"][0],
            description=feature.qualifiers["product"][0],
        )
        for feature in first_record.features
        if (
            (feature.type.lower() == "CDS".lower() or feature.type == "rRNA".lower())
            and any(
                gene.lower() in feature.qualifiers["product"][0].lower()
                for gene in genes
            )
            and not any(
                exclusion.lower() in feature.qualifiers["product"][0].lower()
                for exclusion in exclude
            )
        )
    ]

    result = [
        SeqRecord(
            seq=feature.extract(record.seq),
            id=feature.qualifiers["gene"][0],
            name=feature.qualifiers["gene"][0],
            description=feature.qualifiers["product"][0],
        )
        for record in records
        for feature in record.features
        if (
            (
                feature.type.lower() == "CDS".lower()
                or feature.type.lower() == "rRNA".lower()
            )
            and any(
                gene.lower() in feature.qualifiers["product"][0].lower()
                for gene in genes
            )
            and not any(
                exclusion.lower() in feature.qualifiers["product"][0].lower()
                for exclusion in exclude
            )
        )
    ]

    result.extend(first_record_features)

    return species_name, result


def find_genes_wrapper(
    gbff_folder: str, genes: list[str], exclude: list[str] = []
) -> dict[str, list[SeqRecord]]:
    """Wrapper to load gbff files and find the wanted genes in them, while being aware
    of the strings to exclude.

    Arguments:
        gbff_folder (str): path to the directory containing one or more gbff files.
        genes (list[str]): a list of strings to search in the annotation.
        exclude (list[str]): a list of strings to exclude from the search.
        Defaults to "[]".
    Returns:
        dict[str, list[SeqRecord]]: a dictionary where keys are the name of the species
        to whom the annotation belongs, and values are a list of SeqRecord objects that
        match any string from the list of genes and that do not contain any strings
        from the list of exclusions.
    """
    assert len(genes) > 0, "Received an empty list of genes. Check your input."
    assert os.path.exists(
        gbff_folder
    ), "Please provide a valid path to an existing file."
    assert os.path.isdir(
        gbff_folder
    ), "Please, provide a path to a directory, not a file."

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


def write_dict_out(in_dict: dict[str, list[SeqRecord]], out: bool) -> None:
    """Uses SeqIO to write out a dictionary where the keys are strings with a species
    name and the values are a list of SeqRecord objects. Writes in FASTA format.
    """
    if out:
        with open(out, "w+") as outfile:
            pass


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
        action="store_true",
        help=(
            "whether an output file should be created for each species (uses species"
            " name automatically) (by default outputs to STDOUT)"
        ),
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
            gbff_folder=args.gbff_path, genes=genes, exclude=exclusions
        )

        write_dict_out(in_dict=found_genes, out=args.out)
    except Exception as err:
        raise err


if __name__ == "__main__":
    main()
