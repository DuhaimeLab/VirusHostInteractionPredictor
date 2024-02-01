'''Pytest for genomes_features module.'''
import pytest
import numpy as np
from vhip.mlmodel.genomes_features import KmerProfile


def test_KmerProfile_generate_kmer_words():
    '''Test code to generate k-mer words from given length.'''
    # Kmer words generation
    seq = "ATCG"
    profile = KmerProfile(seq, 3)
    assert (len(profile.generate_kmer_words(1)) == 4)
    assert (len(profile.generate_kmer_words(2)) == 16)
    assert (len(profile.generate_kmer_words(3)) == 64)
    assert (len(profile.generate_kmer_words(6)) == 4096)


def test_KmerProfile_generate_profile():
    '''Test code to generate k-mer profiles from DNA sequences.'''
    # GC content calculation test #1
    seq = "ATCG"
    expected_GCcontent = 50
    profile = KmerProfile(seq, 1)
    profile.generate_profile()
    assert profile.GCcontent == expected_GCcontent

    # GC content calculation test #2
    seq = "AAAAAAA"
    expected_GCcontent = 0
    profile = KmerProfile(seq, 1)
    profile.generate_profile()
    assert profile.GCcontent == expected_GCcontent

    # GC content calculation test #3
    seq = 'GGGGGGGGGGGGG'
    expected_GCcontent = 100
    profile = KmerProfile(seq, 1)
    profile.generate_profile()
    assert profile.GCcontent == expected_GCcontent

    # Kmer profile generation - checking presence of a k3-mer in profile
    seq = "ATCG"
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    k3mer = 'TTC'
    assert k3mer in profile.kmer_words

    # Kmer profile generation - checking presence of a k6-mer in profile
    seq = 'ATCG'
    profile = KmerProfile(seq, 6)
    profile.generate_profile()
    k6mer = 'ATCGGG'
    assert k6mer in profile.kmer_words

    # Kmer profile generation - checking kmer3-word count is correct
    seq = "ATCG"
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    assert sum(profile.profile_counts) == 2
    assert len(profile.profile_counts) == 64
    assert isinstance(profile.profile_counts, np.ndarray)

    # Kmer profile generation - checking kmer6-word count is correct
    seq = 'TTGTCTGCTGTATC'
    profile = KmerProfile(seq, 6)
    profile.generate_profile()
    assert sum(profile.profile_counts) == 9
    assert len(profile.profile_counts) == 4096
    assert isinstance(profile.profile_counts, np.ndarray)

    # Profile is correctly generated even if non-nucleotide characters are present
    seq = "AATCCGGNNKKZZNGG"
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    assert sum(profile.profile_counts) == 5


def test_KmerProfile_empty_string():
    '''Test for empty string as input to KmerProfile.'''
    # Empty string given as input
    seq = ""
    with pytest.raises(ValueError) as excinfo:
        profile = KmerProfile(seq, 3)
        profile.generate_profile()
    assert str(excinfo.value) == "seq cannot be an empty string"



