from vip.mlmodel.features.genomes_features import d2Distance, KmerProfile
import numpy as np

def test_d2Distance_init():
    # test 1 (distance value is meaningful)
    seq1 = 'ATCCTGAGTA'
    seq2 = 'CCAGGCCTGA'

    seq1_profile = KmerProfile(seq1, 6)
    seq1_profile.generate_profile()
    seq2_profile = KmerProfile(seq2, 6)
    seq2_profile.generate_profile()

    test = d2Distance(seq1_profile, seq2_profile)
    test.distance()
    dist = 0.4997 # distance for the two small sequences above
    assert np.round(test.dist, 4) == dist

    # test 2 (sequences are the same so distance is 0)
    seq1 = 'ATTCCTGGAGTGACCGTGATGA'
    seq2 = 'ATTCCTGGAGTGACCGTGATGA'

    seq1_profile = KmerProfile(seq1, 3)
    seq1_profile.generate_profile()
    seq2_profile = KmerProfile(seq2, 3)
    seq2_profile.generate_profile()

    test = d2Distance(seq1_profile, seq2_profile)
    test.distance()
    dist = 0 # distance is 0 since sequences are completely similar
    assert test.dist == dist

    # test 3 (edge case where k used for the sequences do not match)
    seq1 = 'ATCCTGAGTA'
    seq2 = 'CCAGGCCTGA'

    seq1_profile = KmerProfile(seq1, 3)
    seq1_profile.generate_profile()
    seq2_profile = KmerProfile(seq2, 6)
    seq2_profile.generate_profile()

    test = d2Distance(seq1_profile, seq2_profile)
    test.distance()
    assert test.dist == None



