'''
Description goes here
'''

from dataclasses import dataclass

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

    

