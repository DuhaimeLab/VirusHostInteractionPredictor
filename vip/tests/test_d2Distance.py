from vip.mlmodel.features.genomes_features import d2Distance

# test that d2Distance returns an error because k length used
# between the two different sequences are different
def d2Distance_init():
    # sequence
    seq1 = 'ATCCTGAGTA'
    seq2 = 'CCAGGCCTGA'

    # generate k-mer profile
    seq1_profile = KmerProfile(seq1, 3)
    seq1_profile.generate_profile()
    seq2_profile = KmerProfile(seq2, 6)
    seq2_profile.generate_profile()

