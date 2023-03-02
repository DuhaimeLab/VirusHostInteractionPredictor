'''
This script calls methods and functions to compute signals of coevolutions. 
Those signals are necessary for virus0-host predictions. 
'''

from Bio import SeqIO
import os

from features.genomes_features import KmerProfile, d2Distance
from pairs.pairs import Pairs, determine_pairs, list_files
# TODO: find a fix for absolute import rather than doing relative import
from vip.util.read_sequence import read_sequence


virus_directory_path = './data/sequences/viruses/'
host_directory_path = './data/sequences/hosts/'

virus_files = list_files(virus_directory_path) 
host_files = list_files(host_directory_path)
files = virus_files + host_files

GCcontent = dict.fromkeys(files)
k3profiles = dict.fromkeys(files)
k6profiles = dict.fromkeys(files)

<<<<<<< HEAD
# blastn and spacers information
blastn_path = './vip/tests/datatests/blastnNahantCollection_phagevhost.tsv'
spacer_path = './vip/tests/datatests/blastnNahantCollection_phagevspacers.tsv'
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
    current_pair.blastn_hit = homology_match.match(virus, host)




=======
>>>>>>> parent of acbbaba (refactored code)
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
count = []
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


<<<<<<< HEAD
=======
    



>>>>>>> parent of acbbaba (refactored code)
# ATTIC 

'''
def prepend(list, str):
    str += '{0}'
    list = [str.format(i) for i in list]
    return list
'''