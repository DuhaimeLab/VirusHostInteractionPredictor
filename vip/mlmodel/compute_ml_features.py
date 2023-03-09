'''
This script calls methods and functions to compute signals of coevolutions. 
Those signals are necessary for virus0-host predictions. 
'''

from dataclasses import dataclass
from Bio import SeqIO
import itertools
import os

from features.genomes_features import KmerProfile, d2Distance, HomologyMatch
from util.read_sequence import read_sequence, read_headers


@dataclass
class Pairs:

    virus: str
    host: str
    interaction: int = None

    # genome level features
    GCdifference: float = None
    k3dist: float = None
    k6dist: float = None

    homology_hit: bool = False

    # TODO: gene level features



class ComputeFeatures:
    '''
    '''

    def __init__(self, virus_directory, host_directory, ext='fasta', pairs_dict=None) -> None:
        self.virus_directory = virus_directory
        self.host_directory = host_directory
        self.ext = ext

        if pairs_dict:
            self.pairs_dict = pairs_dict
        else:
            self.pairs_dict = None

    
    def add_blastn_files(self, blastn_path, spacer_path):
        '''
        '''
        self.blastn_path = blastn_path
        self.spacer_path = spacer_path

    
    def setup(self, custom_pairs = False):
        '''
        '''

        print('SETUP - ...indexing fasta filenames for viruses and hosts...')
        self.list_files()

        print('SETUP - ...initialize all pairs...')
        if custom_pairs:
            pass
        else:
            self.determine_pairs()

        print('SETUP - ...getting fasta headers...')
        self.get_headers()

        print('SETUP - ...process blastn and spacers output...')
        self.process_blastn()
        self.process_spacers()
        self.homology_matches = HomologyMatch(self.blastn, self.spacers)

        print('SETUP - ...calculate GC content and k-mer profiles...')
        self.generate_kmer_profiles()


    def list_files(self):
        '''
        '''
        self.virus_filenames = [f for f in os.listdir(self.virus_directory) if f.endswith('.' + self.ext)]
        self.host_filenames = [f for f in os.listdir(self.host_directory) if f.endswith('.' + self.ext)]
        self.all_files = self.virus_filenames + self.host_filenames


    def determine_pairs(self):
        '''
        '''
        total_interactions = len(self.virus_filenames) * len(self.host_filenames)
        print(f'-------> There are {len(self.virus_filenames)} viral sequences')
        print(f'-------> There are {len(self.host_filenames)} host sequences')
        print(f'-------> Total number of interactions: {total_interactions}')

        # determine all virus-host pair possible (every host is going to be considered for every virus of interest)
        virus_inter = list(itertools.chain.from_iterable(itertools.repeat(x, len(self.host_filenames)) for x in self.virus_filenames))
        host_inter = self.host_filenames * len(self.virus_filenames)
        pairs = list(zip(virus_inter, host_inter))

        # create list of Pairs
        self.pairs = []
        for pair in list(zip(virus_inter, host_inter)):
            virus = pair[0]
            host = pair[1]
            self.pairs.append(Pairs(virus, host))


    def custom_pairs(self):
        '''
        '''
        pass


    def get_headers(self):
        '''
        '''
        header_filename = {}

        for virus in self.virus_filenames:
            path = self.virus_directory + virus
            headers = read_headers(path)
            for header in headers:
                header_filename[header] = virus
        
        for host in self.host_filenames:
            path = self.host_directory + host
            headers = read_headers(path)
            for header in headers:
                header_filename[header] = host

        self.headers = header_filename
            
    
    def process_blastn(self):
        '''
        '''
        self.blastn = {}

        with open(self.blastn_path, 'r') as f:
            lines = [line.rstrip() for line in f]
            for line in lines:
                split = line.split('\t')
                # get accessions (blastn naming system)
                virus_acc = split[0]
                host_acc = split[1]
                # switch to naming by filenames
                virus = self.headers[virus_acc]
                host = self.headers[host_acc]
                # add virus host relation to dictionary
                if virus not in self.blastn:
                    self.blastn[virus] = [host]
                elif host not in self.blastn[virus]:
                    self.blastn[virus].append(host)
    

    def process_spacers(self):
        '''
        '''
        self.spacers = {}

        with open(self.spacer_path, 'r') as f:
            lines = [line.rstrip() for line in f]
            for line in lines:
                split = line.split('\t')
                # get filename for virus
                virus_acc = split[0]
                virus = self.headers[virus_acc]
                # get filename for host
                tmp = split[1].split('_')
                host_partial = tmp[0] + '_' + tmp[1]
                host = list(filter(lambda x: x.startswith(host_partial), self.host_filenames))
                if virus not in self.spacers:
                    self.spacers[virus] = host
                elif host not in self.spacers[virus]:
                    self.spacers[virus].append(host[0])
    

    def generate_kmer_profiles(self):
        '''
        '''

        self.GCcontent = dict.fromkeys(self.all_files)
        self.k3profiles = dict.fromkeys(self.all_files)
        self.k6profiles = dict.fromkeys(self.all_files)

        for virus in self.virus_filenames:
            path = self.virus_directory + virus
            seq = read_sequence(path)

            seq_profile = KmerProfile(seq, k=1)
            seq_profile.generate_profile()
            self.GCcontent[virus] = seq_profile.GCcontent

            seq_profile = KmerProfile(seq, k=3)
            seq_profile.generate_profile()
            self.k3profiles[virus] = seq_profile

            seq_profile = KmerProfile(seq, k=6)
            seq_profile.generate_profile()
            self.k6profiles[virus] = seq_profile
        
        for host in self.host_filenames:
            path = self.host_directory + host
            seq = read_sequence(path)

            seq_profile = KmerProfile(seq, k=1)
            seq_profile.generate_profile()
            self.GCcontent[host] = seq_profile.GCcontent

            seq_profile = KmerProfile(seq, k=3)
            seq_profile.generate_profile()
            self.k3profiles[host] = seq_profile

            seq_profile = KmerProfile(seq, k=6)
            seq_profile.generate_profile()
            self.k6profiles[host] = seq_profile


    def run(self):
        '''
        '''

        print('RUN - ...checking for homology matches...')
        for pair in self.pairs:
            pair.homology_hit = self.homology_matches.match(pair.virus, pair.host)
        
        print('RUN - ...compute k-mer distances and GC differences...')

        count = 0

        for pair in self.pairs:
            count += 1
            print(f'-------> current pair: {pair.virus} | {pair.host} || Percent done: {round(count / len(self.pairs) * 100, 1)}')
            k3distance = d2Distance(self.k3profiles[pair.virus], self.k3profiles[pair.host])
            k3distance.distance()
            pair.k3dist = k3distance.dist

            k6distance = d2Distance(self.k6profiles[pair.virus], self.k6profiles[pair.host])
            k6distance.distance()
            pair.k6dist = k6distance.dist

            pair.GCdifference = self.GCcontent[pair.virus] - self.GCcontent[pair.host]


    def save_features(self):
        '''
        '''
        pass


#TODO: To transfer code below to tests
'''
virus_directory_path = './test_set/virus_sequences/'
host_directory_path = './test_set/host_sequences/'

blastn_path = './test_set/StaphStudy_virusvhosts.tsv'
spacer_path = './test_set/StaphStudy_virusvspacers_blastn.tsv'


test = ComputeFeatures(virus_directory_path, host_directory_path)
test.add_blastn_files(blastn_path, spacer_path)
test.setup()
test.run()

'''