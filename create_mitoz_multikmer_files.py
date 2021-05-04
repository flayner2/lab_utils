import sys
from collections import defaultdict

infile = sys.argv[1]
outpath = sys.argv[2]

with open(infile, "r") as input_file:
    lines = [line.strip() for line in input_file.readlines()]
    genes = defaultdict(list)

    for line in lines:
        if "Gene_name" in line:
            lines = lines[lines.index(line) + 1 :]
            break

    for line in lines:
        if line.startswith("#") or line.startswith("---"):
            continue
        elif not line:
            break
        else:
            each_col = line.split()
            seq_id = each_col[0]
            gene_name = each_col[6]
            genes[seq_id].append(gene_name)

    with open(f"{outpath}/found_genes.txt", "w+") as output_file:
        for key, vals in genes.items():
            to_write = f"{key} {' '.join(vals)}"
            output_file.write(f"{to_write}\n")

    for line in lines:
        if "Potential missing genes" in line:
            lines = lines[lines.index(line) + 1 :]
            break

    with open(f"{outpath}/missing_genes.txt", "w+") as output_file:
        to_write = ""

        for line in lines:
            if line.startswith("#") or line.startswith("---"):
                continue
            if not line:
                break
            else:
                each_col = line.split()
                to_write += f"{each_col[0]} "

        output_file.write(f"{to_write}\n")
