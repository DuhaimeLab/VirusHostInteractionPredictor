# Virus Host Infection Predictor (VHIP)

[![Test](https://github.com/DuhaimeLab/VirusHostInteractionPredictor/actions/workflows/test.yml/badge.svg?branch=main&event=push)](https://github.com/DuhaimeLab/VirusHostInteractionPredictor/actions/workflows/test.yml)

## Introduction

**Goal**: Predict hosts for virus of interest from sequence data (fasta files)  
**Input**: Sequences for viruses and hosts of interest and the blast results between viruses and hosts, and between viruses and spacers.  
**Output**: Predictions and score for each virus-host pair. This can be used to plot virus-host infection network.

## Requirements

Before the tool can be used, you will need to:

1. Create a virus and a host directory. The virus and host sequences need to be in separate folders. Each fasta file should only represent 1 host or virus.
2. Run blastn between viruses and hosts. Also run CRISPRCasFinder on hosts sequences, then blastn between viruses and the spacers (output of CRISPRCasFinder).
3. Create a conda environment named `vhip` with the modules/packages that are in the `requirements.txt` file.

## Tutorial

A tutorial is included in the `tutorials/` folder, with a test dataset.
The test dataset include sequences for hosts, viruses, and the blastn files (viruses against spacers, and viruses against host sequences).
