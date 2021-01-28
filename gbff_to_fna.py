import os
from Bio import SeqIO

"""
Converts a set of .gbff files to .fna files using Biopython.
Example run inside LGM: `python3 gbff_to_fna.py`
"""

# Change this to the input path containing your files
input_path = '../data/genomes/'
# Chage this to the desired output path. It doesn't need to already be created
output_path = '../data/fna/'

for filename in os.listdir(input_path):
    if filename.endswith('.gbff'):
        in_name = f'{input_path}{filename}'
        outname = f'{output_path}{filename.replace(".gbff", ".fna")}'
        print(f'Converting the file {str(filename)}')
        SeqIO.convert(in_name, 'genbank', outname, 'fasta')
