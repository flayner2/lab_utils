# 1. Introduction

A simple repository to store all the utility scripts I use for my daily work on the LAB. Free to use, modify and distribute.

# 2. Description for each script

Each script is documented on itself. If you have any doubts on how to run each script, check the internal documentation. Here I provide a brief description of each script, in alphabetical order:

- `blast_2_bed`: a program written in Rust that converts a BLAST `.tsv` output generated with `-outfmt` 6 or 7 to a BED file. Multiple alignments for the same query sequence are numbered as something like `Query1_1`, `Query1_2` and so forth. This is an executable so you can just run it directly. To run this, do `./blast_2_bed < [input_file]` or `cat [input_file] | blast_2_bed`, since `blast_2_bed` expects an input from the `STDIN`. The source code for `blast_2_bed` can be found [here](https://github.com/flayner2/blast_2_bed);
- `create_df_from_sheets.py`: connects to the Google Docs API to load data from a particular Google sheet into memmory, as a Pandas DataFrame. Check the [reference](https://developers.google.com/sheets/api/quickstart/python) for additional information on dependencies, keys, permissions and how to modify the code for your needs. You will also need [python-dotenv](https://pypi.org/project/python-dotenv/), [pandas](https://pypi.org/project/pandas/) and a `.env` file like the `.env.example` one provided within this repository. This is not meant to be ran directly. Instead, import this from another Python script and execute it like the example:

  ```python
  from create_df_from_sheets import load_data

  data = load_data()
  ```
- `create_dict.sh`: create a dictionary file for KOMODO2 for a specific database from a set of Interproscan output files;
- `create_symlinks.sh`: creates a set of shorter-named symbolic links from a set of KOMODO2 input files (the parsed id/descriptor files). This is useful for trimming the multiple file extensions and shortening the species names for better visual output/tree branch naming;
- `extract_genes.py`: retrieves a list of strings from the features of a collection of Genbank annotation files (`.gbff` files). Check the help with `python extract_genes.py -h` to understand how to run. A couple of example files are distributed alongside this repo: `example_genes.txt` and `example_exclude.txt`. Check those files to understand the config format;
- `extract_ORF.pl`: exctract all filtered and quality-controlled Open Reading Frames from a set of genome files (e.g. `.gbff` annotation files);
- `filter_seqs_by_blast_expect.py`: filters a set of sequences encoded in a FASTA file based on similarity to a set of target sequences using local BLAST. Produces a set of intermediate `xml` files with information about the BLAST alignments. Final output is a FASTA file with all sequences that didn't match any of the subject sequences with an e-value lower than 10e-10. Run with `python filter_seqs_by_blast_expect.py [path_to_input] [path_to_subject]`. Change output file names and the minimum e-value trhough the code, the comments will guide you;
- `filter_seqs_by_size.py`: filters a set of sequences encoded in a FASTA file based on a length thresshold provided by the user. Outputs the valid sequences to `STDOUT`. Run with `python filter_seqs_by_size.py [path_to_fasta] [min_seq_len: integer]`;
- `gb_2_bed.py`: converts a set of Genbank (`.gbff`) files into a set of BED files with coordinates for the genes identified in the genome files. Run this with `python gb_2_bed.py [genomes_directory_path] [output_directory_path]`, where `[genomes_directory_path]` is a path to a diretory containing one or more `.gbff` files and `[output_directory_path]` is a path to an existing folder where the output files will be saved. The output files are named automatically based on the name of the taxon the genome belongs to;
- `gbff_to_fna.py`: converts a set of Genbank format `.gbff` files to `.fna` files;
- `get_db_counts.py`: queries a set of NCBI databases to fetch the sequence counts for the specified types of sequences. Queries are hard-coded inside the script and correspond to Entrez-style queries. All queries are based off of taxon ids, which are informed through an input `taxids.txt`;
- `get_ests.py`: retrieves all EST sequences for a certain taxon specified by a taxon id from NCBI Nucleotide and outputs a file with all EST sequences for that particular taxon. Each taxon and the respective name, taxon id and EST counts are hard-coded. You can get the EST counts using `get_db_counts.py`;
- `get_longest_seq.pl`: gets the longest coding sequence for each locus from a fasta format ORFs file. It runs on a single file and outputs to `STDOUT`, so its better used with a `for` loop and a redirection (`>`) operator;
- `get_prot_ids.sh`: extracts the protein ids from a fasta file where this information is describred (in the sequence identifier). Path to input files is the only argument. You could also copy the command and run from the command line directly;
- `get_seq_efetch.pl`: uses Entrez's efetch to retrieve all protein sequences from a set of protein identifiers. Identifiers are passed by a line-separated file. Runs on a single file each time so it's better used with a `for` loop. No need to redirect the output since output files are generated automatically.

# 3. External resources

- [TaIGa - Taxonomy Information Gatherer](https://github.com/flayner2/taiga): utility program coded in Python. Used to retrieve taxonomical information for a set of taxa. Accepts multiple forms of inputs, and produces a .csv file with all taxonomical information available from NCBI Taxonomy for the input taxa;
- [taiga-bio](https://github.com/flayner2/taiga-bio): a Python package version of TaIGa. Works roughly the same as the above, but allows you to manipulate TaIGa's results from a custom Python script, as well as use TaIGa's functions as standalone functions. Also available on [PyPi](https://pypi.org/project/taiga-bio/).

# 4. Final regards

Live long, stay safe and respect people. Thanks!
