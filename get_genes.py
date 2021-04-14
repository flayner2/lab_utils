import json

from Bio import Entrez
import xmltodict

Entrez.email = "flayner5@gmail.com"
qry = Entrez.esearch(db="gene", term="Apis mellifera[orgn] AND 18s")
record = Entrez.read(qry)
gene_id = record["IdList"][0]

n_qry = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
print(
    json.dumps(
        xmltodict.parse(n_qry)["Entrezgene-Set"]["Entrezgene"]["Entrezgene_comments"],
        indent=2,
    )
)
n_qry.close()
