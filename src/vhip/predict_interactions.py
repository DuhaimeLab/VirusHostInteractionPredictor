import skops.io as sio

from .mlmodel.compute_ml_features import *


class PredictInteractions(ComputeFeatures):
    """Add methods to ComputeFeatures related to making predictions.

    :param virus_directory: Path to the directory of viruses filenames. 1 fasta file = 1 unique virus
    :type virus_directory: str
    :param host_directory: Path to the directory of host filenames. 1 fasta file = 1 unique host
    :type host_directory: str
    :param ext: Extension of fasta filenames. It assumes it is .fasta.
    :type ext: str
    :param model: Pathway of model to be loaded
    :type model: str
    """

    def __init__(
        self, virus_directory, host_directory, ext="fasta", pairs_of_interest=None
    ) -> None:
        super().__init__(virus_directory, host_directory, ext, pairs_of_interest)
        self.model = None

    def load_model(self, path) -> None:
        """Load machine learning model."""
        self.model = sio.load(path, trusted=True)

    def load_model_user(self, ml_model):
        """Take user model."""
        self.model = ml_model

    def predict(self) -> None:
        """Uses machine learning and uses signals of co-evolution to predict if virus
        infect host. Also calculate score for each prediction.
        """
        # call method to transfer dataclass pairs into a dataframe
        # this is a requirement for scikit-learn
        self.convert_to_dataframe()

        # TODO: check model was loaded

        # run predictions
        print("MODEL - ...making predictions...")
        self.predictions = self.model.predict(self.features_df)
        self.scores = self.model.predict_proba(self.features_df)[:, 1]

        print("MODEL - ...predictions are done!...")

    def save_predictions(self, filename):
        """Save predictions as a tsv file.

        :param filename: Name of output filename.
        :type filename: str
        """
        self.features_df["Predictions"] = self.predictions
        self.features_df["Scores"] = self.scores

        self.features_df.to_csv(filename, sep="\t")
