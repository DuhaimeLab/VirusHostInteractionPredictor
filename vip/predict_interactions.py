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

        self.model = joblib.load(path)

    
    def predict(self):
        '''
        '''

        # call method to transfer dataclass pairs into a dataframe
        # this is a requirement for scikit-learn
        self.convert_to_dataframe()

        # TODO: check model was loaded
        
        # run predictions
        print('MODEL - ...making predictions...')
        predictions = self.model.predict(self.features_df)
        probabilities = self.model.predict_proba(self.features_df)






    
    def convert_to_dataframe(self):
        '''
        '''

        pairs = []
        
        k3dist = []
        k6dist = []
        GCdiff = []
        Homology = []

        for pair in self.computed_pairs:
            virus_host = str(pair.virus + ':' + pair.host)
            pairs.append(virus_host)

            k3dist.append(pair.k3dist)
            k6dist.append(pair.k6dist)
            GCdiff.append(pair.GCdifference)
            Homology.append(int(pair.homology_hit))
        
        self.features_df = pd.DataFrame(list(zip(pairs, GCdiff, k3dist, k6dist, Homology)),
                                        columns = ['pairs', 'GCdiff', 'k3dist', 'k6dist', 'Homology'])
        self.features_df = self.features_df.set_index('pairs')





if __name__ == '__main__':
    virus_directory_path = './test_set/virus_sequences/'
    host_directory_path = './test_set/host_sequences/'

    blastn_path = './test_set/StaphStudy_virusvhosts.tsv'
    spacer_path = './test_set/StaphStudy_virusvspacers_blastn.tsv'

    model_path = './vip/gbrt.pkl'


    test = PredictInteractions(virus_directory_path, host_directory_path)
    test.add_blastn_files(blastn_path, spacer_path)
    test.load_model(model_path)
    test.do_setup()
    test.run_parallel(6)
    test.predict()

