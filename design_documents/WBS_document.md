Activity 1 - Set up folder structure
- [x] Task 1.1 - Make list of functionalities
- [x] Task 1.2 - Broadly design where functions should reside 

Activity 2 - Rewrite basic functionality of VIP 
- [x] Task 2.1 - Keep code from previous VIP iteration that are needed
- [x] Task 2.2 - Implement OOP for predictions
- [ ] Task 2.2 - Write scripts to parse output of CRISPRCasFinder
- [ ] Task 2.3 - Write script to barn blastn output
- [ ] Task 2.4 - Write main python script to be ran by user (this script will calculate all features needed by the model and output virus-host predictions) 

Activity 3 - Implement novel gene features
- [ ] Task 3.1 - Determine pipeline for de novo genome annotation so it is consistent across the board
- [ ] Task 3.2 - Design potential new features at the gene level that might encode signals of coevolution
	- [ ] Task 3.2.1 - Check if # of shared tRNAs between virus and host encode a signal (viruses that encode tRNAs have to go it from somewhere) 
	- [ ] Task 3.2.2 - Check if tRNAs encoded by virus that compensate for the codon bias encode a signal for prediction
	- [ ] Task 3.2.3 - Check codon bias differences between virus and host pairs encode a signal
	- [ ] Task 3.2.4 - Determine if certain viral genes (i.e., structural genes) are under heavier pressure in minimizing codon bias difference (because they are translated more often)
- [ ] Task 3.3 - Implement features and test their usefulness for predicting ecological virus-host interactions
- [ ] Task 3.4 - Decide which new features to keep and include in the pipeline 

Activity 4 - Computing network properties
- [x] Task 4.1 - Implement algorithm to calculate nestedness of network
- [ ] Task 4.2 - Implement algorithm to calculate modularity of network 

Activity 5 - Plots/visualization of data
- [x] Task 5.1 - Code to transform output of machine learning model into input for Gephi for network visualization purposes
- [ ] Task 5.2 - Re-write code to plot/visualize features of coevolution
