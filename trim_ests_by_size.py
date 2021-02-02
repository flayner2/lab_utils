from Bio import SeqIO
from Bio.SeqIO import SeqRecord
import sys


def filter_by_size(seq: SeqRecord, n: int) -> bool:
    """Checks if the length of a SeqRecord is greater than a minimun thresshold 'n'

    Args:
        seq (SeqRecord): a Biopython SeqRecord
        n (int): the min thresshold value

    Returns:
        bool: True if len(seq) > n, else False
    """
    return len(seq) > n


def main() -> None:
    infile = sys.argv[1]
    f_format = "fasta"
    min_size = 50

    for seq in SeqIO.parse(infile, f_format):
        if filter_by_size(seq, min_size):
            print(seq.format(f_format))


def test_filter_by_size() -> None:
    """Tests the function filter_by_size. Run with pytest"""
    valid_seq = SeqRecord(seq="ACTGCTG", id="valid")
    invalid_seq = SeqRecord(seq="ACTG", id="invalid")
    min_size = 5

    assert filter_by_size(valid_seq, min_size) is True
    assert filter_by_size(invalid_seq, min_size) is False


if __name__ == "__main__":
    main()
