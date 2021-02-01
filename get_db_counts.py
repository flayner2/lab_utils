from Bio import Entrez
from Bio.Entrez.Parser import StringElement


# Change to your valid email
USR_EMAIL = "flayner5@gmail.com"


def query(base: str, query_term: str, email: str = USR_EMAIL) -> StringElement:
    """Queries the specified NCBI databases for sequence counts using
    the specified query terms

    Args:
        base (str): the target NCBI database to be queried against
        query_term (str): the Entrez-style query
        email (str, optional): a valid e-mail address. Defaults to USR_EMAIL.

    Returns:
        StringElement: a string-like object containing the sequence counts
    """
    Entrez.email = email

    handle = Entrez.esearch(db=base, term=query_term)
    result = Entrez.read(handle)

    return result["Count"]


def build_counts_dict(tax_id: str) -> dict[str, dict[str, StringElement]]:
    """Builds the specified queries, perform the queries and returns a
    dictionary of dictionaries where each key is a string containing the taxon
    id for a taxon and each value is a dictionary containing an identifier for
    each count as a string and the actual counts as a StringElement

    Args:
        tax_id (str): the taxon id for a specified taxon

    Returns:
        dict[str, dict[str, StringElement]]: a dictionary containing the taxon
        id as key and a dictionary of [str, StringElement] containing database
        identifiers and the sequence counts themselves
    """
    # These are the Entrez-style queries
    # Change these, remove some or add new ones
    mrna_query = f"{tax_id}[orgn] AND biomol_mrna[PROP]"
    refseq_query = f"{tax_id}[orgn] AND refseq[filter]"
    est_query = f"{tax_id}[orgn] AND is_est[filter]"
    sra_illumina_query = (
        f'{tax_id}[orgn] AND ("strategy rip seq"[Properties] OR'
        ' ("strategy other"[Properties] AND "library selection cDNA"[Properties]))'
    )

    counts_dict = dict()

    # Change each dict key and each query correspondingly
    counts_dict["mrna_counts"] = query(base="nucleotide", query_term=mrna_query)
    counts_dict["refseq_counts"] = query(base="nucleotide", query_term=refseq_query)
    counts_dict["est_counts"] = query(base="nucleotide", query_term=est_query)
    counts_dict["sra_illumina_rnaseq_counts"] = query(
        base="sra", query_term=sra_illumina_query
    )

    return counts_dict


def write_out(file: str, collection: dict[str, dict[str, StringElement]]) -> None:
    """Writes the counts to an output .csv file

    Args:
        file (str): the path to the output file
        collection (dict[str, dict[str, StringElement]]): the dictionary
        containing the database counts
    """
    with open(file, "w") as outfile:
        # Change these to be the desired column names
        colnames = ",".join(
            ["tax_id", "num_ests", "num_mrna", "num_refseq", "num_sra_illumina_rnaseq"]
        )
        outfile.write(f"{colnames}\n")

        for k, v in collection.items():
            # Change the accessed keys accordingly
            outfile.write(
                (
                    f'{k},{v["est_counts"]},{v["mrna_counts"]},'
                    f'{v["refseq_counts"]},'
                    f'{v["sra_illumina_rnaseq_counts"]}\n'
                )
            )


def main() -> None:
    # Change this to be the input path to your taxon ids file
    filepath = "./taxids.txt"

    tax_collection = dict()

    with open(filepath, "r") as file:

        tax_ids = file.readlines()

        for tax_id in tax_ids:
            tax_id = f"txid{tax_id.strip()}"
            tax_collection[tax_id] = build_counts_dict(tax_id)

    # Change this to be the path to your output file
    result_file = "./query_result.csv"

    write_out(file=result_file, collection=tax_collection)


if __name__ == "__main__":
    main()
