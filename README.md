# Subclone Deconvolution

In our analysis of data from multiple samples taken from patients with Barrett's Esophagus, we discovered patterns of shared mutations that could be explained by positing the existence of distinguishable subclones within sequenced samples.  Unfortunately, the existence of extensive copy number variation and genome doubling prevented us from using other programs to construct phylogenetic trees from our data.  We instead created subclone phylogenies by hand, but used various programs and scripts to prepare and collate that data for our hand analysis, and then more programs and scripts to collect, quantify, and validate the trees we created.

This repository is a collection of those programs and scripts.  The data itself is not directly downloadable, as it contains genetic information from living individuals who have not consented to having their genetic information disclosed in this way, but will be made available to researchers of Barretts Esophagus as part of the database of Genotypes and Phenotypes (dbGaP) in July of 2020.


branch_lengths_and_signature_calculator.py
assign_branch_signatures.py
make_branch_signature_bar_plots.py
make_trees_with_sigpngs.py
calculate_signature_trajectories.py

