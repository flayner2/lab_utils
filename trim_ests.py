from Bio import SeqIO
from Bio.SeqIO import SeqRecord


def filter_by_size(seq: SeqRecord, n: int) -> bool:
    """Checks if the length of a SeqRecord is greater than a minimun thresshold 'n'

    Args:
        seq (SeqRecord): a Biopython SeqRecord
        n (int): the min thresshold value

    Returns:
        bool: True if len(seq) > n, else False
    """
    return len(seq) > n


def main():
    for seq in SeqIO.parse("./Polistes_canadensis_ests.fasta", "fasta"):
        pass


def test_filter_by_size():
    valid_seq = SeqRecord(seq="ACTGCTG", id="valid")
    invalid_seq = SeqRecord(seq="ACTG", id="invalid")
    min_size = 5

    assert filter_by_size(valid_seq, min_size) == True
    assert filter_by_size(invalid_seq, min_size) == False


if __name__ == "__main__":
    main()