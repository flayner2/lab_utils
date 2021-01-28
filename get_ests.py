from Bio import Entrez
from Bio.Entrez.Parser import ListElement


USR_EMAIL = 'flayner5@gmail.com'


class Taxon():
    def __init__(self, name: str, tax_id: str, est_count: int, est_list: list = []) -> None:
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


def fetch_est_seq(id: str, email: str = USR_EMAIL) -> str:
    Entrez.email = email

    handle = Entrez.efetch(db='nucleotide', id=id, rettype='fasta')
    result = handle.read()

    return result


def search_id(query_term: str, retmax: int, email: str = USR_EMAIL) -> ListElement:
    Entrez.email = email

    handle = Entrez.esearch(db='nucleotide', term=query_term, retmax=retmax)
    result = Entrez.read(handle)

    return result['IdList']


def retrieve_est_ids(taxon: Taxon) -> None:
    est_query = f'txid{taxon.get_tax_id()} [orgn] AND is_est[filter]'
    est_ids = search_id(query_term=est_query, retmax=taxon.get_est_count())

    taxon.set_est_list(est_ids)


def build_seq_and_write(taxon: Taxon) -> None:
    file = f'{taxon.get_name()}_ests.fasta'

    with open(file, 'w') as outfile:
        for id in taxon.get_est_list():
            each_seq = fetch_est_seq(id)
            outfile.write(each_seq)


def main() -> None:
    taxa = [Taxon('Apis_mellifera', '7460', 169511), Taxon(
        'Solenopsis_invicta', '13686', 22883), Taxon('Polistes_canadensis', '91411', 43)]

    for taxon in taxa:
        retrieve_est_ids(taxon)
        build_seq_and_write(taxon)


if __name__ == '__main__':
    main()
