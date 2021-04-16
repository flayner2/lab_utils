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


def _parse_gene_block(
    gene_line: str, rest: list[str], path: str, gene_number: int
) -> tuple[str, list[str]]:
    """Parses a valid gene block from a config file.

    Arguments:
        gene_line (str): the line where the gene block starts.
        rest (list[str]): the following lines up to the end of the file.
        path (str): path to the config file. Used for debugging.
        gene_number (int): line number for the start of the gene block. Used for
        debugging only.

    Returns:
        tuple[str, list[str]]: a 2-tuple with a string representing the gene block
        index and a list of all terms for that gene block.
    """
    genes = []

    split_line = gene_line.split()

    assert len(split_line) > 1, (
        f"Error in {path}\n\nline {gene_number}: {gene_line}\n\n You probably forgot"
        " to add a space between the gene block and the index name."
    )

    index = split_line[-1]

    for line in rest:
        if line.startswith("!") or (line == "\n"):
            pass
        elif line.startswith("#"):
            break
        else:
            genes.append(line.strip())

    return index, genes


def parse_config_file(path: str) -> dict[str, list[str]]:
    """Loads a config file with a specified format. Check the example files.

    Arguments:
        path (str): the path to the input file.
    Returns:
        dict[str, list[str]]: a dictionary with a string representing an index as keys
        and a list of strings representing the terms (genes) as values.
    """
    assert os.path.exists(path), "Please provide a valid path to an existing file."
    assert os.path.isfile(path), "Please, provide a path to a file, not a directory."

    config = defaultdict(list)

    with open(path, "r") as config_file:
        lines = config_file.readlines()

        for number, line in enumerate(lines):
            line = line.strip().lower()

            if line.startswith("!") or (line == "\n"):
                continue
            elif line.startswith("#"):
                if "gene" in line:
                    try:
                        index, genes = _parse_gene_block(
                            gene_line=line,
                            rest=lines[number + 1 :],
                            path=path,
                            gene_number=number,
                        )

                        config[index].extend(genes)
                    except IndexError:
                        print("Found empty gene block at the end of file. Ignoring.")
                else:
                    sys.exit(
                        (
                            f"Error in {path}:\n\nline {number}: {line}\n\nFound a '#'"
                            " but no gene block follows. Remeber to name and index"
                            " your gene blocks."
                        )
                    )
            else:
                continue

    return config


def find_genes_in_annotation(
    records: Generator, genes: dict[str, list[str]], exclude: dict[str, list[str]] = {}
) -> tuple[str, list[tuple[str, SeqRecord]]]:
    """Finds any term from a list of terms in the features of a SeqRecord object,
    excluding any term from a list of exclusions, and returns a tuple with the name of
    the originating species and a list of SeqRecord objects for all features that pass
    those checks.

    Arguments:
        records (Generator): a Generator object which yields SeqRecord objects. Usually
        this is the return of the `SeqIO.parse()` method from BioPython's Bio module.
        genes (dict[str, list[str]]): a dict of strings to find in the SeqRecord's
        description.
        exclude (dict[str, list[str]]): an optional dict of strings used to rule out
        any SeqRecord that contains any of them. Defaults to "{}".
    Returns:
        tuple[str, list[tuple[str, SeqRecord]]]: a 2-tuple with the name of the species
        and a list of 2-tuples whith the gene index and the SeqRecords that contain any
        term from `genes` and do not contain any term from `exclude`.
    """
    assert len(genes) > 0, "Received an empty list of genes. Check your input."

    first_record = next(records)
    species_name = first_record.annotations["organism"]

    first_record_features = [
        (
            index,
            SeqRecord(
                seq=feature.extract(first_record.seq),
                id=feature.qualifiers["gene"][0],
                name=feature.qualifiers["gene"][0],
                description=f"{feature.qualifiers['product'][0]} for {species_name}",
            ),
        )
        for feature in first_record.features
        for index, block in genes.items()
        if (
            (
                feature.type.lower() == "CDS".lower()
                or feature.type.lower() == "rRNA".lower()
            )
            and any(
                (gene in feature.qualifiers["product"][0].lower()) for gene in block
            )
            and not any(
                exclusion in feature.qualifiers["product"][0].lower()
                for exclusion in exclude.get(index, [])
            )
        )
    ]

    result = [
        (
            index,
            SeqRecord(
                seq=feature.extract(record.seq),
                id=feature.qualifiers["gene"][0],
                name=feature.qualifiers["gene"][0],
                description=f"{feature.qualifiers['product'][0]} for {species_name}",
            ),
        )
        for record in records
        for feature in record.features
        for index, block in genes.items()
        if (
            (
                feature.type.lower() == "CDS".lower()
                or feature.type.lower() == "rRNA".lower()
            )
            and any(gene in feature.qualifiers["product"][0].lower() for gene in block)
            and not any(
                exclusion in feature.qualifiers["product"][0].lower()
                for exclusion in exclude.get(index, [])
            )
        )
    ]

    result.extend(first_record_features)

    return species_name, result


def find_genes_wrapper(
    gbff_folder: str, genes: dict[str, list[str]], exclude: dict[str, list[str]] = {}
) -> dict[str, list[tuple[str, SeqRecord]]]:
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


def write_dict_out(
    in_dict: dict[str, list[tuple[str, SeqRecord]]],
    out: str = "",
    for_species: bool = False,
) -> None:
    """Uses SeqIO to write out a dictionary where the keys are strings with a species
    name and the values are a list of SeqRecord objects. Writes in FASTA format.

    Arguments:
        in_dict (dict[str, list[tuple[str, SeqRecord]]]): a dictionary with the species
        name as keys and a list of tuples of gene identifier and SeqRecords as values.
        out (str): the path where the program should create the output files for each
        gene, with the file names generated from the gene identifier. If this is an
        empty string, prints to STDOUT. Defaults to `""`.
        for_species (bool): whether the output files should be created for each species
        instead of each gene. Only works if an output path is also provided. Defaults
        to `False`.
    """
    if out:
        if for_species:
            only_genes = []
            for species, index_and_genes in in_dict.items():
                for _, gene in index_and_genes:
                    only_genes.append(gene)

                filename = f"{species.replace(' ', '_')}_genelist.fasta"
                path = os.path.join(out, filename)

                with open(path, "w+") as outfile:
                    SeqIO.write(only_genes, outfile, "fasta")
        else:
            for species, index_and_genes in in_dict.items():
                for index, gene in index_and_genes:
                    filename = f"{index}_list.fasta"
                    path = os.path.join(out, filename)

                    with open(path, "a+") as outfile:
                        outfile.write(f"{gene.format('fasta')}\n")
    else:
        for _, index_and_genes in in_dict.items():
            for _, gene in index_and_genes:
                print(f"{gene.format('fasta')}")


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
        help="path to the directory with genbank annotation files",
    )
    parser.add_argument(
        "genes_list",
        metavar="L",
        type=str,
        help=(
            "path to a config file with the genes to search the annotation. See the"
            " example"
        ),
    )
    parser.add_argument(
        "-o",
        "--out",
        nargs="?",
        metavar="O",
        type=str,
        default=os.getcwd(),
        help=(
            "where the output files should be created for each gene (uses gene"
            " identifiers automatically). Defaults to the current directory. If not"
            " used, the program will output to STDOUT."
        ),
    )
    parser.add_argument(
        "-e",
        "--exclude",
        type=str,
        help=(
            "path to a config file with the terms to exclude from the search of some"
            " of the genes. See the example"
        ),
    )
    parser.add_argument(
        "-s",
        "--species",
        action="store_true",
        help=(
            "if used in conjunction with `-o`, the ouput files will be create for"
            " each species instead of each gene, using the species name as the file"
            " name by default."
        ),
    )

    args = parser.parse_args()

    try:
        genes = parse_config_file(args.genes_list)

        if args.exclude:
            exclusions = parse_config_file(args.exclude)
        else:
            exclusions = {}

        found_genes = find_genes_wrapper(
            gbff_folder=args.gbff_path, genes=genes, exclude=exclusions
        )

        write_dict_out(in_dict=found_genes, out=args.out)
    except Exception as err:
        raise err


if __name__ == "__main__":
    main()
