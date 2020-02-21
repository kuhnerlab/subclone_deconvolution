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
        branches = lvec[1].split(".")
        found = False
        if len(branches)==1:
            for tree in trees:
                for tip in tree:
                    if branches[0] in tip.name:
                        tip.sigs = sigvals
                        found = True
            pass
        else:
            for tree in trees:
                try:
                    branch = tree.get_common_ancestor(branches)
                    if hasattr(branch, "sigs"):
                        assert(False)
                    branch.sigs = sigvals
                    found = True
                except:
                    continue
        assert(found)
    return sigids

def getUpWithSigs(branch):
    up = branch.up
    if up and not(hasattr(up, "sigs")):
        return getUpWithSigs(up)
    return up

def calculateUpTrajectories(alltrees):
    for label in alltrees:
        for tree in alltrees[label]:
            for branch in tree.traverse():
                if not (hasattr(branch, "sigs")):
                    continue
                upbranch = getUpWithSigs(branch)
                if upbranch:
                    branch.sigdiff = np.subtract(branch.sigs, upbranch.sigs)

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



alltrees = readTrees()
sigids = readAndStoreSignatures(alltrees)
calculateUpTrajectories(alltrees)
sigavgs, diffavgs = calculateAverageAndTrajectories(alltrees)