'''
Description goes here
'''

from features.genomes_features import KmerProfile, d2Distance
from pairs.pairs import Pairs, determine_pairs
import os


print(os.getcwd())


virus_directory_path = './data/sequences/viruses/'
host_directory_path = './data/sequences/hosts/'

test = determine_pairs(virus_directory_path, host_directory_path)

for virus, host in test:
    print(virus, host)
