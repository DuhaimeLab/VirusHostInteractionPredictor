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
    interaction: int

    # genome level features
    GCdifference: float
    k3dist: float
    k6dist: float

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

    virus_inter = list(itertools.chain.from_iterable(itertools.repeat(x, len(host_filenames)) for x in virus_filenames))
    host_inter = host_filenames * len(virus_filenames)
    return list(zip(virus_inter, host_inter))



def list_files(path, ext='fasta'):
    return [f for f in os.listdir(path) if f.endswith('.' + ext)]


