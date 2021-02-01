# 1. Introduction

A simple repository to store all the utility scripts I use for my daily work on the LAB. Free to use, modify and distribute.

# 2. Description for each script

Each script is documented on itself. If you have any doubts on how to run each script, check the internal documentation. Here I provide a brief description of each script, in alphabetical order:

- `create_dict.sh`: create a dictionary file for KOMODO2 for a specific database from a set of Interproscan output files;
- `create_symlinks.sh`: creates a set of shorter-named symbolic links from a set of KOMODO2 input files (the parsed id/descriptor files). This is useful for trimming the multiple file extensions and shortening the species names for better visual output/tree branch naming;
- `extract_ORF.pl`: exctract all filtered and quality-controlled Open Reading Frames from a set of genome files (e.g. `.gbff` annotation files);
- `gbff_to_fna.py`: converts a set of Genbank format `.gbff` files to `.fna` files;
- `get_db_counts.py`: queries a set of NCBI databases to fetch the sequence counts for the specified types of sequences. Queries are hard-coded inside the script and correspond to Entrez-style queries. All queries are based off of taxon ids, which are informed through an input `taxids.txt`;
- `get_ests.py`: retrieves all EST sequences for a certain taxon specified by a taxon id from NCBI Nucleotide and outputs a file with all EST sequences for that particular taxon. Each taxon and the respective name, taxon id and EST counts are hard-coded. You can get the EST counts using `get_db_counts.py`;
- `get_longest_seq.pl`: gets the longest coding sequence for each locus from a fasta format ORFs file. It runs on a single file and outputs to the STDOUT, so its better used with a `for` loop and a redirection (`>`) operator;
- `get_prot_ids.sh`: extracts the protein ids from a fasta file where this information is describred (in the sequence identifier). Path to input files is the only argument. You could also copy the command and run from the command line directly;
- `get_seq_efetch.pl`: uses Entrez's efetch to retrieve all protein sequences from a set of protein identifiers. Identifiers are passed by a line-separated file. Runs on a single file each time so it's better used with a `for` loop. No need to redirect the output since output files are generated automatically;
- `hymenoptera_data.py`: connects to the Google Docs API to load data from a particular Google sheet into memmory, as a Pandas DataFrame;

# 3. External resources

- [TaIGa - Taxonomy Information Gatherer](https://github.com/flayner2/taiga): utility program coded in Python. Used to retrieve taxonomical information for a set of taxa. Accepts multiple forms of inputs, and produces a .csv file with all taxonomical information available from NCBI Taxonomy for the input taxa;
- [taiga-bio](https://github.com/flayner2/taiga-bio): a Python package version of TaIGa. Works roughly the same as the above, but allows you to manipulate TaIGa's results from a custom Python script, as well as use TaIGa's functions as standalone functions. Also available on [PyPi](https://pypi.org/project/taiga-bio/).

# 4. Final regards

Live long, stay safe and respect people. Thanks!
