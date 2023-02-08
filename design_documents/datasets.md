# Datasets

The training and testing of the machine learning model uses the following datasets:
- Host range data (lab-tested knowledge of which virus infect which host). This data came from 3 major sources (NCBI viral database, the Nahant Collection, 
and a large scale host range experiment against Staphylococcus)
- Sequences of viruses and hosts corresponding to the viruses and hosts in the host range dataset

There are a total of 8,849 pairs in the host range dataset. 

## Training and testing set to build machine learning model

The computation of coevolution signals between viruses and hosts are computationally expensive. To reduce time to develop the code, only the host range 
from the Staphylococcus study and their respective sequences for the viruses and hosts were used. 

Once the code is validated, the machine learning model will be tested and trained on the entire dataset. 

The smaller set data is stored in the `/data/` folder. 
