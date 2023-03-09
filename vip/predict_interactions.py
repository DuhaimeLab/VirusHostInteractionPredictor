import pandas as pd
import joblib
import importlib

from mlmodel.compute_ml_features import *


class PredictInteractions(ComputeFeatures):
    '''
    '''

    def __init__(self, virus_directory, host_directory, ext='fasta', pairs_dict=None) -> None:
        super().__init__(virus_directory, host_directory, ext, pairs_dict)
        self.model = None
    
    
    def load_model(self, path):
        '''
        '''

        self.model = None






'''
virus_directory_path = './test_set/virus_sequences/'
host_directory_path = './test_set/host_sequences/'

blastn_path = './test_set/StaphStudy_virusvhosts.tsv'
spacer_path = './test_set/StaphStudy_virusvspacers_blastn.tsv'


test = PredictInteractions(virus_directory_path, host_directory_path)
test.add_blastn_files(blastn_path, spacer_path)
test.setup()
test.run()
'''