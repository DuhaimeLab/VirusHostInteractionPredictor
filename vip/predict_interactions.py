import pandas as pd
import joblib
import importlib

from .mlmodel.compute_ml_features import *


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
        self.predictions = self.model.predict(self.features_df)
        self.scores = self.model.predict_proba(self.features_df)[:, 1]

        print('MODEL - ...predictions are done!...')
    

    def save_predictions(self, filename):
        '''
        '''

        self.features_df['Predictions'] = self.predictions
        self.features_df['Scores'] = self.scores

        self.features_df.to_csv(filename, sep='\t')




