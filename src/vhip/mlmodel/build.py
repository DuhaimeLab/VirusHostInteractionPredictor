"""

"""

from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.utils import resample

import pandas as pd


class BuildModel:
    def __init__(self, training_data_path) -> None:
        # load training data
        data = pd.read_csv(training_data_path)
        data = data.set_index("pairs")

        # downsample non-infection data
        noninf_messages = data[data["infection"] == "NoInf"]
        inf_messages = data[data["infection"] == "Inf"]
        ntarget_samples = int(len(inf_messages) + 0.50 * len(inf_messages))
        noninf_downsample = resample(
            noninf_messages, replace=True, n_samples=ntarget_samples
        )
        tmp = pd.concat([noninf_downsample, inf_messages])

        # select relevant rows for machine learning model
        self.ml_input = tmp[["GCdiff", "k3dist", "k6dist", "Homology"]]
        print(
            "The dataframe is made of {} rows and {} columns!".format(
                self.ml_input.shape[0], self.ml_input.shape[1]
            )
        )

        # retrieve target (what we want to predict)
        ml_target = tmp["infection"]
        self.ml_target = ml_target.map({"NoInf": 0, "Inf": 1})

        # parameters resulting in best performing ml model.
        # they were determined by running a grid search
        self.default_max_depth = 15
        self.default_learning_rate = 0.75
        self.default_loss = "exponential"

    def build(self):
        X_train, X_test, y_train, y_test = train_test_split(
            self.ml_input, self.ml_target, random_state=5, test_size=0.3, train_size=0.7
        )

        self.gbrt = GradientBoostingClassifier(
            random_state=42,
            max_depth=self.default_max_depth,
            learning_rate=self.default_learning_rate,
            loss=self.default_loss,
        )

        self.gbrt.fit(X_train, y_train)

        print(
            "Accuracy on training set: {:.3f}".format(self.gbrt.score(X_train, y_train))
        )
        print("Accuracy on test set: {:.3f}".format(self.gbrt.score(X_test, y_test)))

        return self.gbrt
