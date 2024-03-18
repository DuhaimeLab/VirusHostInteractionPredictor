"""Predict class."""

from sklearn.ensemble import (  # pyright: ignore[reportMissingTypeStubs]
    GradientBoostingClassifier,
)

from .mlmodel.compute_ml_features import ComputeFeatures


class PredictInteractions(ComputeFeatures):
    """Predict virus-host ecological interactions.

    Args:
        virus_directory (str): Pathway to the virus directory
        host_directory (str): Pathway to the host directory
        ext (str): Extension used for the fasta file. Default to fasta
        model (str): Pathway to model to be loaded
    """

    def __init__(
        self, virus_directory: str, host_directory: str, ext: str = "fasta"
    ) -> None:
        """Initialize class variables."""
        super().__init__(virus_directory, host_directory, ext)
        self.model: GradientBoostingClassifier

    def predict(self) -> None:
        """Make interaction prediction for each virus-host pair."""
        # call method to transfer dataclass pairs into a dataframe
        # this is a requirement for scikit-learn
        self.convert_to_dataframe()

        # TODO: check model was loaded

        # run predictions
        print("MODEL - ...making predictions...")
        self.predictions = self.model.predict(self.features_df)  # pyright: ignore
        self.scores = self.model.predict_proba(self.features_df)[:, 1]  # pyright: ignore

        print("MODEL - ...predictions are done!...")

    def save_predictions(self, filename: str):
        """Save predictions as a tsv file.

        Args:
            filename (str): Filename to be used when saving predictions
        """
        self.features_df["Predictions"] = self.predictions  # pyright: ignore
        self.features_df["Scores"] = self.scores  # pyright: ignore

        self.features_df.to_csv(filename, sep="\t")  # pyright: ignore
