from Bio import Entrez
from Bio.Entrez.Parser import ListElement

# Change this to your valid email
USR_EMAIL = "flayner5@gmail.com"


class Taxon:
    """Defines a Taxon object containing a name, a taxon id, a count of the
    available EST sequences and a list of all EST sequence ids. Attributes
    should be accessed via the defined getters and setters;
    """

    def __init__(
        self, name: str, tax_id: str, est_count: int, est_list: list = []
    ) -> None:
        self.name = name
        self.tax_id = tax_id
        self.est_count = est_count
        self.est_list = est_list

    def get_name(self) -> str:
        return self.name

    def get_tax_id(self) -> str:
        return self.tax_id

    def get_est_count(self) -> int:
        return self.est_count

    def get_est_list(self) -> list:
        return self.est_list

    def set_name(self, name: str) -> None:
        self.name = name

    def set_tax_id(self, tax_id: str) -> None:
        self.tax_id = tax_id

    def set_est_count(self, est_count: int) -> None:
        self.est_count = est_count

    def set_est_list(self, est_list: list) -> None:
        self.est_list = est_list
        self.est_list.sort()


def fetch_est_seq(id_list: list, email: str = USR_EMAIL) -> str:
    """Retrieves the sequence for a particular EST id in fasta format

    Args:
        id_list (list): the EST sequence ids list
        email (str, optional): a valid e-mail address. Defaults to USR_EMAIL.

    Returns:
        str: the EST sequence in fasta format, with its identifier
    """
    Entrez.email = email

    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta")
    result = handle.read()

    return result


def search_id(query_term: str, retmax: int, email: str = USR_EMAIL) -> ListElement:
    """Searches NCBI Nucleotide for all available EST sequence ids for a taxon

    Args:
        query_term (str): the Entrez-style query for each taxon
        retmax (int): the total count of EST sequences for that taxon
        email (str, optional): a valid e-mail address. Defaults to USR_EMAIL.

    Returns:
        ListElement: a list-like object with all EST sequence ids for a taxon
    """
    Entrez.email = email

    handle = Entrez.esearch(db="nucleotide", term=query_term, retmax=retmax)
    result = Entrez.read(handle)

    return result["IdList"]


def retrieve_est_ids(taxon: Taxon) -> None:
    """Wrapper to retrieve the EST ids and construct the ids list for a Taxon

    Args:
        taxon (Taxon): a valid Taxon object
    """
    est_query = f"txid{taxon.get_tax_id()}[orgn] AND is_est[filter]"
    est_ids = search_id(query_term=est_query, retmax=taxon.get_est_count())

    taxon.set_est_list(est_ids)


def build_seq_and_write(taxon: Taxon) -> None:
    """Concatenates all sequences from a particular Taxon on a single file
    and writes them to it

    Args:
        taxon (Taxon): a valid Taxon object with a non-empty EST ids list
    """
    print(f"Fetching for taxon: {taxon.get_name()}")
    file = f"{taxon.get_name()}_ests.fasta"

    with open(file, "w") as outfile:
        est_list = taxon.get_est_list()

        if len(est_list) >= 10000:
            while len(est_list) >= 10000:
                seqs = fetch_est_seq(est_list[:10000])
                outfile.write(seqs)
                del est_list[:10000]

            seqs = fetch_est_seq(est_list)
            outfile.write(seqs)
        else:
            seqs = fetch_est_seq(est_list)
            outfile.write(seqs)


def main() -> None:
    # Change this to be a list of your taxon objects.
    # The taxon name is gonna be used to generate the output file names.
    # You can get the EST counts by manually querying NCBI Nucleotide or
    # by using `get_db_counts.py`
    taxa = [
        Taxon("Apis_mellifera", "7460", 169511),
        Taxon("Solenopsis_invicta", "13686", 22883),
        Taxon("Polistes_canadensis", "91411", 43),
    ]

    for taxon in taxa:
        retrieve_est_ids(taxon)
        build_seq_and_write(taxon)


if __name__ == "__main__":
    main()
