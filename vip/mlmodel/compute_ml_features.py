'''
This script calls methods and functions to compute signals of coevolutions. 
Those signals are necessary for virus0-host predictions. 
'''

from Bio import SeqIO
import os

from features.genomes_features import KmerProfile, d2Distance, HomologyMatch
from pairs.pairs import Pairs, determine_pairs, list_files
from util.read_sequence import read_sequence, read_headers



print('.....setup.....')

# define location of viruses and hosts of interest
virus_directory_path = './data/sequences/viruses/'
host_directory_path = './data/sequences/hosts/'
virus_files = list_files(virus_directory_path) 
host_files = list_files(host_directory_path)
files = virus_files + host_files

# determine all possible virus-host pairs
pairs_to_test = determine_pairs(virus_directory_path, host_directory_path)

# initialize dictionaries to store k-mer profiles and GC content
GCcontent = dict.fromkeys(files)
k3profiles = dict.fromkeys(files)
k6profiles = dict.fromkeys(files)

# blastn and spacers information
blastn_path = './vip/tests/datatests/StaphStudy_virusvhosts_blastn.tsv'
spacer_path = None
#blastn_cnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
#                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']



print('.....processing blastn and spacers files.....')

# get headers for each fasta files
header_filename = {}

for virus in virus_files:
    path = virus_directory_path + virus
    headers = read_headers(path)
    for header in headers:
        header_filename[header] = virus

for host in host_files:
    path = host_directory_path + host
    headers = read_headers(path)
    for header in headers:
        header_filename[header] = host

# build dictionary 
blastn = {}

with open(blastn_path, 'r') as f:
    lines = [line.rstrip() for line in f]
    for line in lines:
        split = line.split('\t')
        # get accessions (blastn naming system)
        virus_accession = split[0]
        host_accession = split[1]
        # switch to my naming system (based on filenames)
        virus = header_filename[virus_accession]
        host = header_filename[host_accession]
        # add virus host relation to dictionary
        if virus not in blastn:
            blastn[virus] = [host]
        elif host not in blastn[virus]:
            blastn[virus].append(host)


print(blastn)

print('.....determine if homology match exist between virus and host.....')

homology_match = HomologyMatch(blastn, None)

for virus, host in pairs_to_test:
    current_pair = Pairs(virus, host)
    current_pair.blastn_hit = homology_match.check_blastn(virus, host)
    


print('.....generating k-mer profiles.....')

# generate k-mer profiles for viruses and their potential hosts
for virus in virus_files:
    path = virus_directory_path + virus
    seq = read_sequence(path)

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
    path = host_directory_path + host
    seq = read_sequence(path)

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

# calculate features for each pair of interest
count = 0
for virus, host in pairs_to_test:
    
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


print(current_pair.k3dist)
print(current_pair.k6dist)
print(current_pair.GCdifference)
print(current_pair.blastn_hit)

# ATTIC 

'''
def prepend(list, str):
    str += '{0}'
    list = [str.format(i) for i in list]
    return list
'''