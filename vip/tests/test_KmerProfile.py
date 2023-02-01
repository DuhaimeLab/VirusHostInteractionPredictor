from vip.mlmodel.features.genomes_features import KmerProfile

def test_KmerProfile_GCcontent():
    # test 1
    seq1 = 'ATCG'
    k1 = 1
    expected_GCcontent = 50
    profile = KmerProfile(seq1, k1)
    profile.generate_profile()
    assert profile.GCcontent == expected_GCcontent

    # test 2
    seq2 = 'AAAAAAA'
    k2 = 1
    expected_GCcontent = 0
    profile = KmerProfile(seq2, k2)
    profile.generate_profile()
    assert profile.GCcontent == expected_GCcontent
    

def test_KmerProfile_kmers():
    # test 1
    seq1 = 'ATCG'
    k1 = 3

    profile = KmerProfile(seq1, k1)
    profile.generate_profile()

    assert sum(profile.profile_counts) == 2


