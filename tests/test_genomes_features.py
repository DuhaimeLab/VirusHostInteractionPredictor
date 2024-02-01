'''Pytest for genomes_features module.'''

import numpy as np
from vhip.mlmodel.genomes_features import KmerProfile, d2Distance


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

    seq = "ATCG"
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    assert sum(profile.profile_counts) == 2

    # test 4 (profile generated even if N are present in sequence)
    seq = "AATCCGGNNNGG"
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    assert sum(profile.profile_counts) == 5

    # test 5 (empty string)
    seq = ""
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
