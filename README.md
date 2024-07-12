# online_true_discovery

This repository contains the R code for reproducing the results in the paper "Online true discovery control with e-values". In the aforementioned paper, we introduced new online algorithms that provide simultaneous lower bounds for the number of discoveries in data-adaptive rejection sets.

Files:

boosted_e_vals.R:            Contains the functions for boosting e-values. 

SeqE-Guard:                  Contains the functions for our general algorithms for online true discovery guarantee.

generate_data:               Applies "boosted_e_vals.R" and "SeqE-Guard" to generate the simulated data and saves it in the file "true_disc_prop.rda" in the results folder.

generate_plots:              Uses the data in "results/true_disc_prop.rda" to create plots 1-3 of the paper.
