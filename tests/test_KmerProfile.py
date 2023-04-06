from vip.mlmodel.genomes_features import KmerProfile, d2Distance

def test_KmerProfile_generate_profile():
    # test 1 (GC content calculation)
    seq = 'ATCG'
    expected_GCcontent = 50
    profile = KmerProfile(seq, 1)
    profile.generate_profile()
    assert profile.GCcontent == expected_GCcontent

    # test 2 (GC content calculation)
    seq = 'AAAAAAA'
    expected_GCcontent = 0
    profile = KmerProfile(seq, 1)
    profile.generate_profile()
    assert profile.GCcontent == expected_GCcontent
    
    # test 3 (profile generated)
    seq = 'ATCG'
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    assert sum(profile.profile_counts) == 2

    # test 4 (profile generated even if N are present in sequence)
    seq = 'AATCCGGNNNGG'
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    assert sum(profile.profile_counts) == 5

    # test 5 (empty string)
    seq = ''
    profile = KmerProfile(seq, 3)
    profile.generate_profile()
    






