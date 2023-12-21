'''Pytest for reading sequences.'''

from vhip.mlmodel.read_sequence import read_headers


def test_read_headers():
    '''Test to read header of a fasta file.'''
    headers = ["NZ_CP065712.1"]
    filename = "tests/datatests/test_sequence.fasta"
    res = read_headers(filename)
    assert res == headers
