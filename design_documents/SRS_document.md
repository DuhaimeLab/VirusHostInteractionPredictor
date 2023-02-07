# Software requirement specification

Tool: Virus Infection Predictor (VIP)
Author: Eric Bastien

# I. Introduction

This tool predict ecological interactions between a set of viruses (either infection bacteria and/or archaea) and a set of hosts solely based on genomic sequences. This is possible by leveraging signals of coevolution/adaptation between the viruses and their hosts that is embedded in DNA sequences. 

# II. Overall description

Input - sequences of viruses and hosts of interest, and the output of blastn between the viruses and hosts, and the viruses and CRISPR-Cas spacers. 

Output - a dataframe containing the signals of coevolution for each possible virus-host pair, and the predicted ecological interaction along with a score representing the strength of that prediction. 

# III. Requirements

To apply the model on a small dataset (less than 50,000 pairs to be tested), this tool can be ran locally. On any dataset larger than 50k pairs, we recommend using a HPC cluster. 

To re-train and re-test the machine learning model, a HPC cluster is needed. 


 
