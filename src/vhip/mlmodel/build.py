"""Build machine learning model de novo.

In case the skops file is unusable, this file contains all the information needed to retrain VHIP.
"""

import pandas as pd  # pyright: ignore[reportMissingTypeStubs]
from sklearn.ensemble import (  # pyright: ignore[reportMissingTypeStubs]
    GradientBoostingClassifier,
)
from sklearn.model_selection import (  # pyright: ignore[reportMissingTypeStubs]
    train_test_split,  # pyright: ignore[reportUnknownVariableType]
)
from sklearn.utils import (  # pyright: ignore[reportMissingTypeStubs]
    resample,  # pyright: ignore[reportUnknownVariableType]
)


class BuildModel:
    """Class to re-build the model from host range data."""

    def __init__(self, training_data_path: str) -> None:
        """Initialize class variables.

        Args:
            training_data_path (str): Pathway to the training/testing host range dataset.
        """
        # load training data
        data = pd.read_csv(training_data_path)  # pyright: ignore[reportUnknownMemberType]
        data = data.set_index("pairs")  # pyright: ignore[reportUnknownMemberType]

        # downsample non-infection data
        noninf_messages = data[data["infection"] == "NoInf"]  # pyright: ignore
        inf_messages = data[data["infection"] == "Inf"]  # pyright: ignore
        ntarget_samples = int(len(inf_messages) + 0.50 * len(inf_messages))  # pyright: ignore
        noninf_downsample = resample(  # pyright: ignore[reportUnknownVariableType]
            noninf_messages, replace=True, n_samples=ntarget_samples
        )
        tmp = pd.concat([noninf_downsample, inf_messages])  # pyright: ignore

        # select relevant rows for machine learning model
        self.ml_input = tmp[["GCdiff", "k3dist", "k6dist", "Homology"]]  # pyright: ignore
        print(
            "The dataframe is made of {} rows and {} columns!".format(
                self.ml_input.shape[0],  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
                self.ml_input.shape[1],  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
            )
        )

        # retrieve target (what we want to predict)
        ml_target = tmp["infection"]  # pyright: ignore
        self.ml_target = ml_target.map({"NoInf": 0, "Inf": 1})  # pyright: ignore

        # parameters resulting in best performing ml model.
        # they were determined by running a grid search
        self.default_max_depth = 15
        self.default_learning_rate = 0.75
        self.default_loss = "exponential"

    def build(self):
        """Build machine learning model de novo."""
        X_train, X_test, y_train, y_test = train_test_split(  # pyright: ignore[reportUnknownVariableType]
            self.ml_input,  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
            self.ml_target,  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
            random_state=5,
            test_size=0.3,
            train_size=0.7,  # pyright: ignore
        )

        self.gbrt = GradientBoostingClassifier(
            random_state=42,
            max_depth=self.default_max_depth,
            learning_rate=self.default_learning_rate,
            loss=self.default_loss,  # pyright: ignore
        )

        self.gbrt.fit(X_train, y_train)  # pyright: ignore

        print(
            "Accuracy on training set: {:.3f}".format(self.gbrt.score(X_train, y_train))  # pyright: ignore
        )
        print("Accuracy on test set: {:.3f}".format(self.gbrt.score(X_test, y_test)))  # pyright: ignore

        return self.gbrt
