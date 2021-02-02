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


for seq in SeqIO.parse("./Polistes_canadensis_ests.fasta", "fasta"):
