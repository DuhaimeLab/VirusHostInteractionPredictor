"""Pytest to read sequence."""

from fileinput import filename
from vhip.mlmodel.read_sequence import read_headers, read_sequence, read_gene_products

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
    """Test to check if reading fasta file is working as intended."""
    # test 1 - check length of sequences match
    filename = "tests/datatests/test_sequence.fasta"
    res = read_sequence(filename)
    assert len(res[0]) == len(cleaned_seq)

    # test 2 - check that sequences do match
    assert res[0] == cleaned_seq


def test_read_headers():
    """Test to read header of a fasta file."""
    # Test genome with one contig
    headers_1 = ["NZ_CP065712.1"]
    filename_1 = "tests/datatests/test_sequence.fasta"
    res = read_headers(filename_1)
    assert res == headers_1

    # Test annotated gene file
    headers_2 = ["ABDEAL_00005","ABDEAL_00010","ABDEAL_00015"]
    filename_2 = "tests/datatests/test_annotated_genes.ffn"
    res = read_headers(filename_2)
    assert res == headers_2


def test_read_sequence_many_contigs():
    """Test to read a fasta file with many contigs."""
    filename = "tests/datatests/test_sequence_multiple_contigs.fasta"
    res = read_sequence(filename)
    assert len(res) == 2


def test_read_gene_products():
    """Test to read gene products from an annotated fasta gene file."""
    products = ["Chromosomal replication initiator protein DnaA","Beta sliding clamp","S4 domain-containing protein YaaA"]
    filename = "tests/datatests/test_annotated_genes.ffn"
    res = read_gene_products(filename)
    assert res == products
