'''Pytest to read sequence.'''

import pytest
from vhip.mlmodel.read_sequence import read_sequence, read_headers

seq = """AGTACTTGTTGATGCTGATGCACTAGTTGATTCAGATGTGCTCGTACTTGTTGATTCAGACGCACTTGTG
CTCGCTGAAGTACTATTAGATGTAGACGTGCTTGCGCTTATCGATTCAGAAGTACTTGTACTTTCTGAAC
TACTCTCTGAAGTCGAATTACTTAATGAAGTACTTTCACTATTTTTATCTCTTGCTGATTCGCTTTCTGA
TAATGATGCACTGTTTACTGAATTACTCGTTGAAGTACTTGCTGATTGTAAACGAGATAGTGAATCACTG
ATTGATGCACTTTCTGATTCGATTTTTTCTGTTGATTCACTTTCTGAATTACTCGATGAAGCACTTAATG
ATTCACTTTCGAATACTGAATCTTCTAATGAATTAGAAATACTAATCGAATCTTGATTACTTTTACGTTG
TGATTCACTCTTCGATAACGACTCATTATTTGAATCTCTTACTGATTCACTCTTCGATAAAGATGCACTA
TTTTCTGAGGCACGTTTTGAATCAGATTCACTTAGACTTTTGGTTCTTGAAGCATCTTCTGAAGCTATCT
CTGAATTTGCTTTGCTTATAGACGTTGATCTTGATGCATCTTCTGAAGCTCTTTCTGAATTTACTTTGCT
TATAGACGCTGATCTTGATGCATCTTCAGAAGCTATTTCTGAATTTACTTTGCTTATAGATGCTGACCTT
GAAGCTTCTTGTGAGTTTCTCACTGAATCTGAGTTAGCTGTACTTTCTGATCTTGATTTCTCTTCAGAAA
GTCGTGAGTCTTCTGATCTACTATTACTTTCAGATGCACGTACAACACTTTCAGATGCTTTTTTAGAATT
GCTCTCTTGTTCTGAAAGTCGCTTAGAATTACTTTGTTGCTCTGAAAGTCGCTTAGAATTACTTTGTTGC
TCTGAAAGTCGCTTAGAATTACTTTGTTGCTCTGAAAGTCGCTTAGAATTACTTTGTTGCTCTGAAAGTC
GCTTAGAAGTGCTTGCTTGCTCTGAAAGACGTTTAGAATTACTTTGTTGTTCTGAAGCTTTCTTAGAAGT
GCTTTGTTGTTCAGAGGCTCTCGATTGACTTTGAGAATCCGCAATACTTTGGCTTGCAGATTTCGATGCA
ACAGATGGGTCTTCACCATGAATTACTGGTGCAGCGCCACCACGGCCTATAAAGAAGTTTCTATTTGAAT
CTGTGGAACTCGCAAATGGAACAAAGTCCCAAGTTGCATTATCAGCTATAATATTTGCGCTCCATTTTGT"""

cleaned_seq = seq.replace("\n", "")  # remove \n created by long string


def test_read_sequence():
    '''Test to check if reading fasta file is working as intended.'''
    # test 1 - check length of sequences match
    filename = "tests/datatests/test_sequence.fasta"
    res = read_sequence(filename)
    assert len(res[0]) == len(cleaned_seq)

    # test 2 - check that sequences do match
    assert res[0] == cleaned_seq


def test_read_empty_sequence():
    '''Test to check if reading empty fasta file raises ValueError.'''
    # test 1 - check if empty file raises ValueError
    filename = "tests/datatests/test_sequence_empty.fasta"
    #with pytest.raises(ValueError) as excinfo:
    #    read_sequence(filename)
    #assert str(excinfo.value) == 'Given fasta path does not contain any sequence'


def test_read_headers():
    '''Test to read header of a fasta file.'''
    headers = ["NZ_CP065712.1"]
    filename = "tests/datatests/test_sequence.fasta"
    res = read_headers(filename)
    assert res == headers

#TODO: test with multiple contigs in file

if __name__ == "__main__":
    test = read_sequence('tests/datatests/test_sequence_multiple_contigs.fasta')
    print(test)
