#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:58:00 2020

@author: lpsmith

This script takes as input one file that contains the assigned signature
probabilities (from a separate by-patient analysis of all unique mutations
found in that patient), as well as the by-branch lists of mutations that
were output from 'branch_lengths_and_signature_calculator.py', and combines
them to calculate the mutation signatures on every branch in the trees.
These are then written to an output file.
"""

from __future__ import division
from os import walk
import numpy as np

onlysomepatients = False
somepatients = ["483"]
#somepatients = ["387", "483", "609", "611", "626", "631", "635", "652", "852"] #double split patients

#input files and directories:
signature_file = "BE_signatures\BE_patient_based_corrected\Mutation_with_signature_probability_annotated_patient_based.txt"
sigdir = "branch_signatures/"

outfilename = "all_signatures.tsv"

def readSigProbs():
    """
    Read in the signature probabilities from signature_file, and store in
    a dictionary.
    """
    sigs = []
    sigprobs = {}
    sigfile = open(signature_file, "r")
    for line in sigfile:
        lvec = line.rstrip().split("\t")
        if "Patient" in line:
            sigs = lvec[4:]
            continue
        (patient, chrom, pos, refalt) = lvec[0:4]
        patient = patient.split("_")[1]
        (ref, alt) = refalt.split(">")
        if patient not in sigprobs:
            sigprobs[patient] = {}
        sigprobs[patient][(chrom, pos)] = lvec[4:]
        for n, prob in enumerate(sigprobs[patient][(chrom, pos)]):
            sigprobs[patient][(chrom, pos)][n] = float(prob)
    sigfile.close()
    return (sigprobs, sigs)


def hasSix(patient):
    """
    Ten patients have six samples instead of four.  These patients are 
    sometimes analyzed two ways: once with 'all' and once without.  This
    returns 'True' if they're one of the ten.
    """
    if patient in ("55", "59", "126", "184", "381", "478", "609", "635", "865", "909"):
        return True
    return False


sigprobs, sigs = readSigProbs()

outfile = open(sigdir + outfilename, "w")
outfile.write("Patient")
outfile.write("\tBranch")
for sig in sigs:
    outfile.write("\t" + sig)
outfile.write("\tOther\n")

sigfiles = []
for __, _, files in walk(sigdir):
    sigfiles += files

for file in sigfiles:
    if "mutations" not in file:
        continue
    filesigs = np.zeros(len(sigs))
    fvec = file.split(".")
    label = fvec[0]
    patient = label.split("_")[0]
    if "all" not in label and hasSix(label):
        print("Skipping", label, "so we only analyze the full version.")
        continue
    if onlysomepatients and patient not in somepatients:
        continue
    tips = fvec[1:-2]
    sigfile = open(sigdir + file, "r")
    nmuts = 0
    for line in sigfile:
        if "Chrom" in line:
            continue
        (chrom, pos, alt) = line.rstrip().split("\t")
        nsamples = 0
        sigset = np.zeros(len(sigs))
        for tip in tips:
            onetip = tip.split("-")[0].split("_")[0]
            if (chrom, pos) in sigprobs[patient]:
                probs = sigprobs[patient][(chrom, pos)]
                sigset = [sigset[i]+probs[i] for i in range(len(sigset))]
                nsamples += 1
        assert(nsamples > 0)
        for n in range(len(sigset)):
            sigset[n] = sigset[n]/nsamples
        filesigs = [filesigs[i]+sigset[i] for i in range(len(filesigs))]
        nmuts += 1
    sigfile.close()
    remainder = nmuts
    outfile.write(label)
    outfile.write("\t")
    tipstr = ""
    for tip in tips:
        tipstr += tip + "."
    tipstr = tipstr[:-1]
    outfile.write(tipstr)
    for sigprob in filesigs:
        outfile.write("\t" + str(sigprob/nmuts))
        remainder = remainder - sigprob
    outfile.write("\t" + str(remainder/nmuts))
    outfile.write("\n")

outfile.close()
