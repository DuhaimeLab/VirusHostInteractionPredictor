'''
Description goes here
'''

from features.genomes_features import KmerProfile, d2Distance
from pairs.pairs import Pairs, determine_pairs, list_files

from Bio import SeqIO
import os

def read_sequence(dir, filename):
    path = dir + filename
    for record in SeqIO.parse(path, 'fasta'):
        return str(record.seq)


virus_directory_path = './data/sequences/viruses/'
host_directory_path = './data/sequences/hosts/'

virus_files = list_files(virus_directory_path) 
host_files = list_files(host_directory_path)
files = virus_files + host_files

GCcontent = dict.fromkeys(files)
k3profiles = dict.fromkeys(files)
k6profiles = dict.fromkeys(files)

print('.....generating k-mer profiles.....')

# Generate k-mer profiles for viruses and their potential hosts
for virus in virus_files:
    seq = read_sequence(virus_directory_path, virus)

    seq_profile = KmerProfile(seq, k=1)
    seq_profile.generate_profile()
    GCcontent[virus] = seq_profile.GCcontent

    seq_profile = KmerProfile(seq, k=3)
    seq_profile.generate_profile()
    k3profiles[virus] = seq_profile

    seq_profile = KmerProfile(seq, k=6)
    seq_profile.generate_profile()
    k6profiles[virus] = seq_profile

for host in host_files: 
    seq = read_sequence(host_directory_path, host)

    seq_profile = KmerProfile(seq, k=1)
    seq_profile.generate_profile()
    GCcontent[host] = seq_profile.GCcontent

    seq_profile = KmerProfile(seq, k=3)
    seq_profile.generate_profile()
    k3profiles[host] = seq_profile

    seq_profile = KmerProfile(seq, k=6)
    seq_profile.generate_profile()
    k6profiles[host] = seq_profile


print('.....computing GC difference and distances between k-mer profiles.....')


# Calculate features for each pair of interest
pairs = []
count = 0
pairs_to_test = determine_pairs(virus_directory_path, host_directory_path)
for virus, host in pairs_to_test:
    
    current_pair = Pairs(virus, host)
    
    k3distance = d2Distance(k3profiles[virus], k3profiles[host])
    k3distance.distance()

    k6distance = d2Distance(k6profiles[virus], k6profiles[host])
    k6distance.distance()

    current_pair.GCdifference = GCcontent[virus] - GCcontent[host]
    current_pair.k3dist = k3distance.dist
    current_pair.k6dist = k6distance.dist

    # print message
    count += 1
    if (count % 100) == 0:
        print(f'Progress -- {round(count / len(pairs_to_test) * 100, 1)}%')


    



# ATTIC 

'''
def prepend(list, str):
    str += '{0}'
    list = [str.format(i) for i in list]
    return list
'''