# Virus Host Infection Predictor (VHIP)

[![Test](https://github.com/DuhaimeLab/VirusHostInteractionPredictor/actions/workflows/test.yml/badge.svg?branch=main&event=push)](https://github.com/DuhaimeLab/VirusHostInteractionPredictor/actions/workflows/test.yml)
[![Coverage Status](https://coveralls.io/repos/github/DuhaimeLab/VirusHostInteractionPredictor/badge.svg?branch=main)](https://coveralls.io/github/DuhaimeLab/VirusHostInteractionPredictor?branch=main)

## Introduction

VHIP is a machine-learning model that predict virus-microbe interactions (i.e., infection or non-infection) from genomic sequences of viruses and microbes of interest. It leverages virus-microbe coevolution signal that are extracted from genomic sequences. VHIP was trained on lab-verified virus-microbe pairs collected from literature and the NCBI virus database. 

You can find more information about the philosophy and performance of VHIP here: https://www.biorxiv.org/content/10.1101/2023.11.03.565433v1 

## Installation 

This module will be made available as a conda environment. To be updated once it is up on conda. 

Also will be available with `pip install virushostinteractionpredictor`. 


## Inputs

The inputs for VHIP are sequences for the viruses and hosts of interest. Each file should represent an unique virus/host species. In addition, a blastn between the virus and host sequences, and between viruses and host spacers are needed. If there are no results from the blastn, then the files can be empty. 

## Example

An example is included in the `example` folder. It includes virus and microbes sequences, how they should be organized (separate folders), and example blastn files. Make sure you are able to run VHIP on the example before applying to a new dataset. 


## Citation

https://www.biorxiv.org/content/10.1101/2023.11.03.565433v1