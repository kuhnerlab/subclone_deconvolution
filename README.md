# Subclone Deconvolution

In our analysis of data from multiple samples taken from patients with Barrett's Esophagus, we discovered patterns of shared mutations that could be explained by positing the existence of distinguishable subclones within sequenced samples.  Unfortunately, the existence of extensive copy number variation and genome doubling prevented us from using other programs to construct phylogenetic trees from our data.  We instead created subclone phylogenies by hand, but used various programs and scripts to prepare and collate that data for our hand analysis, and then more programs and scripts to collect, quantify, and validate the trees we created.

This repository is a collection of those programs and scripts.  The data itself is not directly downloadable, as it contains genetic information from living individuals who have not consented to having their genetic information disclosed in this way, but will be made available to researchers of Barretts Esophagus as part of the database of Genotypes and Phenotypes (dbGaP) in July of 2020.

Required data files/directories (not included):
* P01CA91955-WGS80-Full-Pilot-Samples.csv List of all samples and information about the patient and location they came from.
* calling_evidence_odds.tsv	File that calculates the odds of a given sample being diploid or tetraploid.  Last column is the final call.
* lucian_from_kanika.csv	List of every called mutation in every sample.
* noninteger_processed_CNs/	Directory containing CNV calls for every sample.
* signature_probabilities/all_muts_with_signature_probabilities.txt	File with mutation signature assignments per mutation.


To run and create everything, run the python scripts in the following order:

* branch_lengths_and_signature_calculator.py
* assign_branch_signatures.py
* make_branch_signature_bar_plots.py
* make_trees_with_sigpngs.py
* calculate_signature_trajectories.py

