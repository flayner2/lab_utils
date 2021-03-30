from Bio import SeqIO

import sys
import csv
from collections import defaultdict


def read_bed(path: str) -> dict[str, list]:
    """
    A BED file is nothing more than a tsv file so this simple function just reads a
    tsv file with a specific conformation, that is, the output of bedtools intersect
    with the -wo flag.

    Arguments:
        path (str): a path to a BED file.

    Returns:
        dict[str, list]: a dictionary containing information about each contig
        on a BED file.
    """
    result = defaultdict(list)

    with open(path, "r") as bed_file:
        csv_reader = csv.reader(bed_file, delimiter="\t")

        locus = ""
        contig = ""
        op_len = 0

        for row in csv_reader:
            if locus == "" and contig == "":
                locus = row[7]
                contig = row[3]
                scaffold = row[0]
                op_len += int(row[-1])
            elif contig == row[3]:
                if locus == row[7]:
                    op_len += int(row[-1])
                elif result.get(contig) is not None:
                    for each_hsp in result[contig]:
                        if each_hsp["locus"] == row[7]:
                            each_hsp["op_len"] += int(row[-1])
                else:
                    each_row = {}
                    each_row["scaffold"] = row[0]
                    each_row["op_len"] = int(row[-1])
                    each_row["locus"] = row[7]

                    result[contig].append(each_row)
            else:
                each_row = {}
                each_row["scaffold"] = scaffold
                each_row["op_len"] = op_len
                each_row["locus"] = locus

                result[contig].append(each_row)

                locus = row[7]
                contig = row[3]
                op_len = int(row[-1])

    return result


def read_fasta(path: str) -> list[dict]:
    """
    Uses Biopython to read a FASTA file containing contigs that correspond to those
    from a BED file.

    Arguments:
        path: path to a FASTA file.

    Returns:
        list[dict]: a list of dictionaries containing the id of the contig and its
        length.
    """
    results = []

    for record in SeqIO.parse(path, "fasta"):
        each_record = {}
        each_record["contig"] = record.id
        each_record["len"] = len(record.seq)

        results.append(each_record)

    return results


def build_output(
    bed: list[dict], fasta: list[dict], outname: str = "output.tsv"
) -> None:
    """
    Write the information about contig coverages to a file.

    Arguments:
        bed (list[dict]): a list of dictionaries containing information about a set of
        BED records.
        fasta (list[dict]): a list of dictionaries containing information about a set of
        FASTA records.
        outname (str): a name for the output file. Defaults to "output.tsv".
    """
    with open(outname, "w") as outfile:
        outfile.write("scaffold\tcontig\tlocus\top_len\tcontig_len\tcoverage\n")

        for contig in fasta:
            for b_record in bed[contig["contig"]]:

                coverage = b_record["op_len"] / contig["len"]
                outfile.write(
                    (
                        f"{b_record['scaffold']}\t{contig['contig']}\t"
                        f"{b_record['locus']}\t{b_record['op_len']}\t{contig['len']}\t"
                        f"{coverage}\n"
                    )
                )


def main() -> None:
    try:
        bed_records = read_bed(sys.argv[1])
        fasta_records = read_fasta(sys.argv[2])
        if len(sys.argv) > 3:
            build_output(bed_records, fasta_records, sys.argv[3])
        else:
            build_output(bed_records, fasta_records)
    except IndexError:
        sys.exit(
            (
                "Didn't provide valid paths to at least one of the expected files or"
                " the ouptut file name."
            )
        )
    except Exception as e:
        sys.exit(f"Error occurred while trying to run: {e}")


if __name__ == "__main__":
    main()
