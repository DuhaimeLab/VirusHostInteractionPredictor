from vip.mlmodel.features.util.read_sequence import read_sequence

seq = '''AGTACTTGTTGATGCTGATGCACTAGTTGATTCAGATGTGCTCGTACTTGTTGATTCAGACGCACTTGTG
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
CTGTGGAACTCGCAAATGGAACAAAGTCCCAAGTTGCATTATCAGCTATAATATTTGCGCTCCATTTTGT'''

cleaned_seq = seq.replace('\n', '') #remove \n created by long string

def test_read_sequence():
    # test 1 - check length of sequences match
    filename = 'vip/tests/datatests/test_sequence.fasta'
    res = read_sequence(filename)
    assert len(res) == len(cleaned_seq)

    # test 2 - check that sequences do match
    assert res == cleaned_seq


