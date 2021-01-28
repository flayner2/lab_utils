from Bio import Entrez
from Bio.Entrez.Parser import StringElement


USR_EMAIL = 'flayner5@gmail.com'


def query(base: str, query_term: str, email: str = USR_EMAIL) -> StringElement:
    Entrez.email = email

    handle = Entrez.esearch(db=base, term=query_term)
    result = Entrez.read(handle)

    return result['Count']


def build_counts_dict(tax_id: str) -> dict:
    mrna_query = tax_id + '[orgn] AND biomol_mrna[PROP]'
    refseq_query = tax_id + '[orgn] AND refseq[filter]'
    est_query = tax_id + '[orgn] AND is_est[filter]'
    sra_illumina_query = tax_id + ('[orgn] AND ("biomol rna"[Properties] AND'
                                   ' "platform illumina" [Properties]')

    counts_dict = dict()

    counts_dict['mrna_counts'] = query(
        base='nucleotide', query_term=mrna_query)
    counts_dict['refseq_counts'] = query(
        base='nucleotide', query_term=refseq_query)
    counts_dict['est_counts'] = query(base='nucleotide', query_term=est_query)
    counts_dict['sra_illumina_rnaseq_counts'] = query(
        base='sra', query_term=sra_illumina_query)

    return counts_dict


def write_out(file: str, collection: dict) -> None:
    with open(file, 'w') as outfile:
        outfile.write(
            'tax_id,num_ests,num_mrna,num_refseq,num_sra_illumina_rnaseq\n')

        for k, v in collection.items():
            outfile.write((f'{k},{v["est_counts"]},{v["mrna_counts"]},'
                           f'{v["refseq_counts"]},'
                           f'{v["sra_illumina_rnaseq_counts"]}\n'))


def main() -> None:
    filepath = './taxids.txt'

    tax_collection = dict()

    with open(filepath, 'r') as file:

        tax_ids = file.readlines()

        for tax_id in tax_ids:
            tax_id = 'txid' + tax_id.strip()
            tax_collection[tax_id] = build_counts_dict(tax_id)

    result_file = './query_result.csv'

    write_out(file=result_file, collection=tax_collection)


if __name__ == '__main__':
    main()
