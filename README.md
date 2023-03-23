# Virus Infection Predictor (VIP)


## Introduction 

**Goal**: Predict hosts for virus of interest from sequence data  
**Input**: Sequences for viruses and hosts of interest and the blast results between viruses and hosts, and between viruses and spacers. 
**Output**: Predictions and score for each virus-host pair. This can be used to plot virus-host infection network. 

## To run the tool 

The script file to make the prediction is here: `/vip/predict_interactions.py`. User will need to change the inputs to fit for their projects. 

In addition, I strongly recommend creating a conda environment (I personally named vip). The list of required packages can be found in the `requirements.txt` file. 

## Example

To showcase how this tool work, I made a file named `example.ipynb` to showcase the setup needed to run the tool. 

