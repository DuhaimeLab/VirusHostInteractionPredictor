"""Pytest for d2distances measurements of k-mer profiles."""

import numpy as np
from vhip.mlmodel.genomes_features import KmerProfile, d2Distance


def test_d2Distance_distance():
    """Test distance measurement of k-mer profiles."""
    # test 1 (distance value is meaningful)
    seq1 = ["ATCCTGAGTA", "ATCCTGGGGCACGGTGCG"]
    seq2 = ["CCAGGCCTGA"]

    seq1_profile = KmerProfile(seq1, 6)
    seq1_profile.generate_profile()
    seq2_profile = KmerProfile(seq2, 6)
    seq2_profile.generate_profile()

    test = d2Distance(seq1_profile, seq2_profile)
    test.distance()
    dist = 0.5004 # distance for the two small sequences above
    assert np.round(test.dist, 4) == dist #pyright: ignore[reportGeneralTypeIssues]

    # test 2 (sequences are the same so distance is 0)
    seq1 = ["ATTCCTGGAGTGACCGTGATGA"]
    seq2 = ["ATTCCTGGAGTGACCGTGATGA"]

    seq1_profile = KmerProfile(seq1, 3)
    seq1_profile.generate_profile()
    seq2_profile = KmerProfile(seq2, 3)
    seq2_profile.generate_profile()

    test = d2Distance(seq1_profile, seq2_profile)
    test.distance()
    dist = 0  # distance is 0 since sequences are completely similar
    assert test.dist == dist

    # test 3 (edge case where k used for the sequences do not match)
    seq1 = ["ATCCTGAGTA"]
    seq2 = ["CCAGGCCTGA"]

    seq1_profile = KmerProfile(seq1, 3)
    seq1_profile.generate_profile()
    seq2_profile = KmerProfile(seq2, 6)
    seq2_profile.generate_profile()

    test = d2Distance(seq1_profile, seq2_profile)
    test.distance()
    assert test.dist is None
