'''
Description goes here
'''

from dataclasses import dataclass
import itertools
import os


@dataclass
class Pairs:

    virus: str
    host: str
    interaction: int = None

    # genome level features
    GCdifference: float = None
    k3dist: float = None
    k6dist: float = None

    spacers_hit: bool = False
    blastn_hit: bool = False

    # TODO: gene level features



def determine_pairs(virus_directory, host_directory):
    virus_filenames = list_files(virus_directory)
    host_filenames = list_files(host_directory)

    total_interactions = len(virus_filenames) * len(host_filenames)
    print(f'There are {len(virus_filenames)} viral sequences')
    print(f'There are {len(host_filenames)} host sequences')
    print(f'Total number of interactions: {total_interactions}')

    # determine all virus-host pair possible (every host is going to be considered for every virus of interest)
    virus_inter = list(itertools.chain.from_iterable(itertools.repeat(x, len(host_filenames)) for x in virus_filenames))
    host_inter = host_filenames * len(virus_filenames)
    return list(zip(virus_inter, host_inter))



def list_files(path, ext='fasta'):
    return [f for f in os.listdir(path) if f.endswith('.' + ext)]


