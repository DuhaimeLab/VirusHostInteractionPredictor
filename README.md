# Virus Infection Predictor (VIP)


## Introduction 

**Goal**: Predict hosts for virus of interest from sequence data (fasta files) 
**Input**: Sequences for viruses and hosts of interest and the blast results between viruses and hosts, and between viruses and spacers. 
**Output**: Predictions and score for each virus-host pair. This can be used to plot virus-host infection network. 


## Requirements

Before the tool can be used, you will need to:
1. Create a virus and a host directory. Each fasta file should only represent 1 host/virus.
2. Run blastn between viruses and hosts. Also run CRISPRCasFinder on hosts sequences, then blastn between viruses and the spacers (output of CRISPRCasFinder). 
3. Create a conda environment named `vip` with the modules/packages that are in the `requirements.txt` file. 

An example test set is provided. For the example, the directories are already created (`/test_set`) and the blastn have already been computed. 


## To run the tool 

To showcase how this tool work, I made a file named `example.ipynb` to showcase the setup needed to run the tool. 

To use it for your own purposes, you will need to modify the user input code chunk. 

