import os
from Bio import SeqIO

for filename in os.listdir("../data/genomes/"):
    if filename.endswith(".gbff"):
        in_name = "../data/genomes/" + filename
        outname = "../data/fna/" + filename.replace(".gbff", ".fna")
        print("Converting the file " + str(filename))
        SeqIO.convert(in_name, "genbank", outname, "fasta")

