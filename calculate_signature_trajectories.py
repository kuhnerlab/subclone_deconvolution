#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:54:07 2020

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
from ete3 import Tree, TreeStyle, TextFace, AttrFace

import csv
import numpy as np
import random

import lucianSNPLibrary as lps
#Alternative for importing the 'lucianSNPLibrary':
#import imp
#lps = imp.load_source("lps","/home/mkkuhner/Papers/phylo/lucianSNPLibrary.py")

onlysomepatients = False
somepatients = ["1005"]
#somepatients = ["387", "483", "609", "611", "626", "631", "635", "652", "852"] #double split patients

#input files and directories:
treedir = "branch_lengths/"
sigfile = "branch_signatures/all_signatures.tsv"



class Treeset:
    def __init__(self, treelist):
        self.treelist = treelist
        self.branchlist_nosplits = []
        self.branchlist_splits = []
        nosplit_index = 0
        split_index = 0
        indexmap = {}
        for tree in treelist:
            for branch in tree.traverse():
                if hasattr(branch, "sigs"):
                    if "-" in branch.name:
                        branch.nosplit_index = -1
                        if "-1" in branch.name:
                            self.branchlist_splits.append(branch);
                            branch.split_index = split_index
                            sampleid = branch.name.split("-")[0]
                            indexmap[sampleid] = split_index
                            split_index += 1
                    else:
                        self.branchlist_nosplits.append(branch)
                        self.branchlist_splits.append(branch)
                        branch.split_index = split_index
                        branch.nosplit_index = nosplit_index
                        split_index += 1
                        nosplit_index += 1
                else:
                    branch.nosplit_index = -1
                    branch.split_index = -1
        for tree in treelist:
            for tip in tree:
                if "-" in tip.name and "-1" not in tip.name:
                    sampleid = tip.name.split("-")[0]
                    tip.split_index = indexmap[sampleid]
        #Just to check:
        for n, branch in enumerate(self.branchlist_splits):
            assert(branch.split_index == n)
        for n, branch in enumerate(self.branchlist_nosplits):
            assert(branch.nosplit_index == n)
        for tree in self.treelist:
            for branch in tree.traverse():
                assert(hasattr(branch, "split_index"))
                assert(hasattr(branch, "nosplit_index"))

    def getUpWithSigs(self, branch):
        up = branch.up
        if up and not(hasattr(up, "sigs")):
            return self.getUpWithSigs(up)
        return up
    
    def calculateAvgAll(self):
        allsigs = []
        for branch in self.branchlist_splits:
            allsigs.append(branch.sigs)
        return np.average(allsigs, axis=0)
    
    def calculateAvgNosplits(self):
        allsigs = []
        for branch in self.branchlist_nosplits:
            allsigs.append(branch.sigs)
        return np.average(allsigs, axis=0)
    
    def calculateOneUpTrajectory(self, splitvec=None):
        if not splitvec:
            splitvec = list(range(len(self.branchlist_splits)))
        alldiffs = []
        for tree in alltrees[label]:
            for branch in tree.traverse():
                if not (hasattr(branch, "sigs")):
                    continue
                upbranch = self.getUpWithSigs(branch)
                if upbranch:
                    branchsigs = self.branchlist_splits[splitvec[branch.split_index]].sigs
                    upsigs = self.branchlist_splits[splitvec[upbranch.split_index]].sigs
                    alldiffs.append(np.subtract(branchsigs, upsigs))
        return np.average(alldiffs, axis=0)
    
    def calculateUpTrajectories(self):
        ntries = 100000
        canonical = self.calculateOneUpTrajectory()
        nlarger = np.zeros(len(canonical))
        nsmaller = np.zeros(len(canonical))
        nsame = np.zeros(len(canonical))
        significance = np.zeros(len(canonical))
        for n in range(ntries):
            splitvec = list(range(len(self.branchlist_splits)))
            random.shuffle(splitvec)
            diffrand = self.calculateOneUpTrajectory(splitvec=splitvec)
            for n in range(len(canonical)):
                if diffrand[n] < canonical[n]:
                    nsmaller[n] += 1
                elif diffrand[n] > canonical[n]:
                    nlarger[n] += 1
                else:
                    nsame[n] += 1
        for n in range(len(significance)):
            significance[n] = (min(nlarger[n], nsmaller[n]) + nsame[n])/ntries
        return canonical, significance
    
    def calculateAverageAndTrajectories(alltrees):
        sigavgs = {}
        diffavgs = {}
        for label in alltrees:
            allsigs = []
            alldiffs = []
            for tree in alltrees[label]:
                for branch in tree.traverse():
                    if hasattr(branch, "sigs"):
                        allsigs.append(branch.sigs)
                    if hasattr(branch, "sigdiff"):
                        alldiffs.append(branch.sigdiff)
            sigavgs[label] = np.average(allsigs, axis=0)
            diffavgs[label] = np.average(alldiffs, axis=0)
        return sigavgs, diffavgs
    
        



def readTrees():
    """
    Read in the hand-generated Newick trees for all patients.  In many cases,
    multiple trees are present for a single patient, because fewer than 100
    mutations arose on the common ancestor of all the samples' subclones.
    """
    alltrees = {}
    treefiles = []
    for __, _, files in walk(treedir):
        treefiles += files
    for treefile in treefiles:
        if "newick" not in treefile:
            continue
        label = treefile.split(".")[0]
        if onlysomepatients and label not in somepatients:
            continue
        newicks = []
        for line in open(treedir + treefile, "r"):
            if "#" in line:
                continue
            newicks.append(line.rstrip())
    #    halfindexes = range(round(len(newicks)/2))
    #    newicks = numpy.delete(newicks, halfindexes)
        trees = []
        for newick in newicks:
            tree = Tree(newick)
            trees.append(tree)
        alltrees[label] = trees
    return alltrees

def readAndStoreSignatures(alltrees):
    for line in open(sigfile):
        lvec = line.rstrip().split()
        if "Patient" in line:
            sigids = lvec[2:-1]
            continue
        sigvals = lvec[2:-1]
        for n, valstr in enumerate(sigvals):
            sigvals[n] = float(valstr)
        label = lvec[0]
        if onlysomepatients and label not in somepatients:
            continue
        trees = alltrees[label]
        tipids = lvec[1].split(".")
        found = False
        if len(tipids)==1:
            for tree in trees:
                for tip in tree:
                    if tipids[0] in tip.name:
                        tip.sigs = sigvals
                        found = True
            pass
        else:
            for tree in trees:
                try:
                    branch = tree.get_common_ancestor(tipids)
                    if hasattr(branch, "sigs"):
                        assert(False)
                    branch.sigs = sigvals
                    found = True
                except:
                    continue
        assert(found)
    #Now put labels on the branches
#    for label in alltrees:
        
    return sigids


def makeTreesets(alltrees):
    treesets = {}
    for label in alltrees:
        if "all" in label:
            #LS DEBUG:  REMOVE THIS WHEN WE GET T3 DATA!!
            continue
        treesets[label] = Treeset(alltrees[label])
    return treesets




    


alltrees = readTrees()
sigids = readAndStoreSignatures(alltrees)
treesets = makeTreesets(alltrees)

avgfile = open("average_sigs.tsv", "w")
avgfile.write("Patient")
for sigid in sigids:
    avgfile.write("\t" + sigid)
avgfile.write("\n")

difffile = open("sig_diffs.tsv", "w")
difffile.write("Patient")
for sigid in sigids:
    difffile.write("\t" + sigid)
    difffile.write("\t" + sigid + " significance")
difffile.write("\n")


for label in treesets:
    #calculate the average signature per patient:
    avgfile.write(label)
    treeset = treesets[label]
    sigavgs = treeset.calculateAvgAll()
    for sigavg in sigavgs:
        avgfile.write("\t" + str(sigavg))
    avgfile.write("\n")
    print(str(sigavgs))
    
    #Calculate the trajectory of every signature in the patient, i.e.
    # test whether it goes up or down as you go from root to tip.
    #
    #The 'significance' is the percentage of randomized trials
    # that produced the same or a more extreme value for that signature.
    # A value of '1' means all the values were the same.
    difffile.write(label)
    diffavgs, diff_significance = treeset.calculateUpTrajectories()
    print(str(diffavgs))
    print(str(diff_significance))
    for n in range(len(diffavgs)):
        difffile.write("\t" + str(diffavgs[n]))
        difffile.write("\t" + str(diff_significance[n]))
    difffile.write("\n")
        
    
avgfile.close()
difffile.close()