#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:58:00 2020

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
from ete3 import Tree, TreeStyle, TextFace, AttrFace
import numpy as np

import csv

import lucianSNPLibrary as lps
#Alternative for importing the 'lucianSNPLibrary':
#import imp
#lps = imp.load_source("lps","/home/mkkuhner/Papers/phylo/lucianSNPLibrary.py")

onlysomepatients = False
somepatients = ["483"]
#somepatients = ["387", "483", "609", "611", "626", "631", "635", "652", "852"] #double split patients

#input files and directories:
signature_file = "signature_probabilities/all_muts_with_signature_probabilities.txt"
sigdir = "branch_signatures/"

outfilename = "all_signatures.tsv"

def readSigProbs():
    sigs = []
    sigprobs = {}
    sigfile = open(signature_file, "r")
    for line in sigfile:
        lvec = line.rstrip().split("\t")
        if "Sample" in line:
            sigs = lvec[4:]
            continue
        (sample, chrom, pos, refalt) = lvec[0:4]
        (ref, alt) = refalt.split(">")
        if sample not in sigprobs:
            sigprobs[sample] = {}
        sigprobs[sample][(chrom, pos)] = lvec[4:]
        for n, prob in enumerate(sigprobs[sample][(chrom, pos)]):
            sigprobs[sample][(chrom, pos)][n] = float(prob)
    sigfile.close()
    return (sigprobs, sigs)





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
    if "all" in label:
        print("Skipping", label, "because T3 mutations not analyzed.")
        continue
    if onlysomepatients and label not in somepatients:
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
            if (chrom, pos) in sigprobs[onetip]:
                probs = sigprobs[onetip][(chrom, pos)]
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
