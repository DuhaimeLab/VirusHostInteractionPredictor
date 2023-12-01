from Bio import SeqIO


def read_sequence(path):
    """Return the sequence as a string.
    Take path where file reside as input.
    """
    for record in SeqIO.parse(path, "fasta"):
        return str(record.seq)


def read_headers(path):
    """ """
    result = []
    for record in SeqIO.parse(path, "fasta"):
        result.append(record.id)
    return result
