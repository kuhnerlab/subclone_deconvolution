#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:53:12 2020

@author: lpsmith

Because the 'optim' files must be created by hand, there is a high probability
of typos and other mistake creeping in.  However, while it cannot be
generated automatically, it can be *validated* automatically, which this
program does.  It checks:
    * That the VAF ratio is a valid ratio for the given call
    * That the VAF values and 'nmut' values are present in the 'histpeak'
      files that were used to create the optim files and trees.
"""
#


from __future__ import division
from os import walk


VAFdir = "histpeaks"
optimdir = "final_optim/"

onlysomepatients = False
somepatients = ["483"]

def createSixToFourMap():
    """
    For patients with six samples, we did our analysis in two ways:  once with
     all six samples, and once with just the four T1/T2 samples.  In the latter
     case, the 'histpeak' inputs had different letter codes for the samples
     than the six-sample 'histpeak' inputs.  However, to simplify things, the
     optim files used the naming scheme from the six-sample inputs.  This
     small routine simply creates a conversion from one to the other, so that
     the four-sample optim files can be correctly validated against the 
     'histpeak' inputs.
    """
    ret = {}
    ret["126"] = {}
    ret["126"]["b"] = "a"
    ret["126"]["c"] = "b"
    ret["126"]["d"] = "c"
    ret["126"]["e"] = "d"
    ret["55"] = {}
    ret["55"]["a"] = "a"
    ret["55"]["b"] = "b"
    ret["55"]["c"] = "c"
    ret["55"]["f"] = "d"
    ret["609"] = {}
    ret["609"]["a"] = "a"
    ret["609"]["b"] = "b"
    ret["609"]["e"] = "c"
    ret["609"]["f"] = "d"
    ret["635"] = {}
    ret["635"]["a"] = "a"
    ret["635"]["b"] = "b"
    ret["635"]["c"] = "c"
    ret["635"]["f"] = "d"
    ret["865"] = {}
    ret["865"]["a"] = "a"
    ret["865"]["b"] = "b"
    ret["865"]["e"] = "c"
    ret["865"]["f"] = "d"
    ret["909"] = {}
    ret["909"]["a"] = "a"
    ret["909"]["b"] = "b"
    ret["909"]["c"] = "c"
    ret["909"]["f"] = "d"
    return ret


sixtofour = createSixToFourMap()


optimfiles = []
for __, _, files in walk(optimdir):
    optimfiles += files

for file in optimfiles:
    if "input" not in file:
        continue
    patient = file.split("_")[0]
    if onlysomepatients and patient not in somepatients:
        continue

    allsets = {}
    oneset = {}
    setkey = set()
    for line in open(optimdir + "/" + file, "r"):
        if line[0] == "#":
            continue
        if line=="\n":
            setkey = list(setkey)
            setkey.sort()
            groupid = ""
            for letter in setkey:
                groupid += letter
            allsets[groupid] = oneset
            oneset = {}
            setkey = set()
            continue
        lvec = line.rstrip().split()
        (letter, vaf, nmuts, tipgroup) = lvec[0:4]
        if "only4" in file and patient in sixtofour:
            letter = sixtofour[patient][letter]
        call = "1,1"
        if len(lvec)>4:
            call = lvec[4]
        calls = call.split("_")
        if len(calls) > 1:
            call = calls[2]
        call = "(" + call + ")"
        call = eval(call)
        tipvec = tipgroup.split("_")
        if "x" not in tipvec[-1]:
            if tipvec[-1][-1] != str(sum(call)):
                print("Possible discrepancy in file", file, ":  tip group alleles don't match call total.")
                print(line)
        if len(calls) > 1:
            call = "(" + calls[0] + ")"
            call = eval(call)
        setkey.add(letter)
        lettercall = letter + "-" + str(call)
        if lettercall not in oneset:
            oneset[lettercall] = []
        truth = False
        if nmuts=="0":
            truth = True #These won't be found in the spreadsheet
        oneset[lettercall].append([vaf, nmuts, line, truth])
    setkey = list(setkey)
    setkey.sort()
    groupid = ""
    for letter in setkey:
        groupid += letter
    allsets[groupid] = oneset
    peakfilename = VAFdir
    if "only4" in file:
        peakfilename += "_4"
    peakfilename += "/" + patient + "_report.tsv"
    peakfile = open(peakfilename, "r")
    headers = []
    for line in peakfile:
        lvec = line.rstrip().split("\t")
        if "Patient" in line:
            headers = lvec
            continue
        lettercall = lvec[0]
        for n, header in enumerate(headers):
            if header not in allsets:
                continue
            if lettercall not in allsets[header]:
                continue
            vafs = lvec[n].split(":")
            nmuts = lvec[n+1].split(":")
            for n, vaf in enumerate(vafs):
                vaf = eval(vaf)
                vaf = vaf*100
                vaf = str(round(vaf))
                mut = nmuts[n]
                for n, (o_vaf, o_nmuts, line, __) in enumerate(allsets[header][lettercall]):
                    if vaf == o_vaf and o_nmuts == mut:
                        allsets[header][lettercall][n][3] = True
                    elif vaf == o_vaf or o_nmuts == mut:
                        print("Potential mismatch. Histfile:", vaf, mut, "Optim:", o_vaf, o_nmuts)
    for header in allsets:
        for lettercall in allsets[header]:
            for __, __, line, truth in allsets[header][lettercall]:
                if truth == False:
                    print("VAF and nmuts not found in", file, ":")
                    print(line)
                
            
        
        
