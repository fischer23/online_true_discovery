# online_true_discovery

This repository contains the R code for reproducing the results in the paper "Online closed testing with e-values". In the aforementioned paper, we introduced new online algorithms that provide simultaneous lower bounds for the number of discoveries in data-adaptive rejection sets.

Files:

boosted_e_vals.R:            Contains the functions for boosting e-values. 

SeqE-Guard:                  Contains the functions for our general algorithms for online true discovery guarantee.

generate_data_tau:           Generates the simulated data for Section 3.2 and saves it in the file "TD_data_taus.rda" in the results folder.

generate_data:               Applies "boosted_e_vals.R" and "SeqE-Guard" to generate the simulated data for Section 4 and saves it in the file "TD_data.rda" in the results folder.

generate_plots:              Uses the data in "results/TD_data.rda" and "results/TD_data_taus.rda" to create plots 1-4 of the paper.
