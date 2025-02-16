# online_true_discovery

This repository contains the R code for reproducing the results in the paper "Admissible online closed testing must employ e-values". In the aforementioned paper, we proved via the (online) closure principle that online multiple testing procedures with true discovery guarantee must employ anytime-valid tests for the intersection hypotheses. Motivated by the recently established importance of test martingales for anytime-valid testing, we introduced new a online true discovery algorithm based on sequential e-values.

Files:

boosted_e_vals.R:            Contains the functions for boosting e-values. 

SeqE-Guard:                  Contains the functions for our general algorithms for online true discovery guarantee.

generate_data_tau:           Generates the simulated data for Section 3.3 and saves it in the file "TD_data_taus.rda" in the results folder.

generate_data:               Applies "boosted_e_vals.R" and "SeqE-Guard" to generate the simulated data for Section 4 and the supplementary material and saves it in the file "TD_data.rda" in the results folder.

generate_plots:              Uses the data in "results/TD_data.rda" and "results/TD_data_taus.rda" to create plots 4-7 and S.1-S.2 of the paper.
