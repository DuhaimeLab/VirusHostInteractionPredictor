import pandas as pd

import os

print(os.getcwd())

### load and format input for ml model ###
data = pd.read_csv("data/ml_input.csv")
data = data.set_index("pairs")

ml_input = data[["GCdiff", "k3dist", "k6dist", "Homology"]]
ml_target = data["infection"]
ml_target = ml_target.map({"NoInf": 0, "Inf": 1})

### Determine best hyperparameters
jter = 100
iter = 100

train_frac = 0.7
test_frac = 1 - train_frac

# setting different parameters
learning_rates = [0.1, 0.2, 0.3, 0.5, 0.75]
loss = ["log_loss", "exponential"]
max_depth = [3, 5, 6, 7, 8]


# split data into train+validation and test set
X_trainval, X_test, y_trainval, y_test = train_test_split(
    ml_data, ml_target, random_state=0, test_size=testing_size, train_size=training_size
)
