R codes and documentation for Shi, Wei and Huang (SIM2020)
In this paper, we proposed landmark analysis of a survival dataset with longitudinal biomarker using FPCA. In order to demonstrate the advantages of our model (LMFPCA), 
we compared it with the traditional landmark model (LMonly) using real datasets and simulated data. Please refer to our paper for details.
Here we use the simulated data to show how to apply both approaches in data analysis. In this folder, the following files are included:   

1. R codes for generating a survival dataset with longitudinal biomarker: dataset_sim.R 
2. Simulated data folder: dataset_gen. It consists of four files: phi.RData, tc_mu.RData, tc_FPC.RData are eigen functions, means and FPC scores estimated from the Framingham Heart Study data for FPCA-based simulation.
   dat_sim.csv is an example of simulated datasets, which is used for the codes: LMonly.R and LMFPCA.R.
3. R codes for running the traditional landmark analysis model on simulated data: LMonly.R.
4. R codes for running the model of landmark analysis coupled with FPCA on simulated data: LMFPCA.R.
5. R codes for simulation study to compare the performance of both models, LMonly and LMFPCA, in term of Brier scores and AUC values: simulation-LMonly-LMFPCA.R  

