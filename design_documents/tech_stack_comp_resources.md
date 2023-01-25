# Technology used in the project and computing requirements

The need of using high computer resources depends on the scale the user would like to apply the model. 
On small dataset, VIP can be ran locally. 

For large datasets (dozens of potential hosts, thousands of viruses), signals of coevolution can be time consuming to run. For that reason, I strongly
suggest using HPC (i.e., Great Lakes at the University of Michigan). 

The training and testing of the machine learning model was also done on Great Lakes due to the need of training thousands of models to determine
the best hyperparameters. 
Approximation: it takes about 2 minute to create 1 model considering the size of the data and the number of features considered. To explore the feature
space, typically 1000s models are generated, so 2000 minutes are needed. That is ~33.3 hours. 
