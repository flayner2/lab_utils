import pprint
import sys

from Bio import Entrez, SeqIO
from Bio.Entrez.Parser import StringElement


def _parse_file(handle: str) -> list[str]:
    with open(handle, "r") as infile:
        return [line.strip() for line in infile.readlines()]


def get_ids(species: list[str], email: str) -> dict[str, StringElement]:
    results = {}
    pp = pprint.PrettyPrinter()

    Entrez.email = email

    for each_sp in species:
        handle = Entrez.esearch(
            db="nucleotide",
            term=f'{each_sp}[orgn] AND ("complete sequence"[All Fields] OR "complete genome"[All Fields]) AND mitochondrion[filter]',
            sort="Date Modified",
        )
        record = Entrez.read(handle)
        handle.close()

        try:
            results[each_sp] = record["IdList"][0]
        except KeyError:
            sys.exit(f"No id found for {each_sp}")
        except IndexError:
            print(f"Problem indexing ids for {each_sp}. Take a look at the dict:")
            pp.pprint(record)
            sys.exit()

    return results


def main() -> None:
    try:
        infile = sys.argv[1]
        email = sys.argv[2]

        species = _parse_file(infile)

        ids = get_ids(species, email)

        for each_name, each_id in ids.items():
            n_handle = Entrez.efetch(
                db="nucleotide", id=each_id, retmode="text", rettype="gb"
            )
            test = SeqIO.read(n_handle, "genbank")
            print(f"{each_name}: {len(test.seq)}")
            # SeqIO.write(test, "test.gbff", "genbank")
            n_handle.close()
    except IndexError:
        sys.exit(
            (
                "Please provide a path to an input file and an email: python3 "
                "download_mitogenomes.py {input_file} {email}"
            )
        )


if __name__ == "__main__":
    main()
