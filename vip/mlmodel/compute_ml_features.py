'''
Description goes here
'''


from features.genomes_features import KmerProfile, d2Distance

testseq = 'ATCCAGTGGTAGACCAGTGGGGGGGGGGGGGGGGGGGGGG'


test = KmerProfile(testseq, 1)
test.generate_profile()

