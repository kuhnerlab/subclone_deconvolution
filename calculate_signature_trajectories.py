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

splitstoo = False
ntries = 10000

#output files:
avgout = "average_sigs.tsv"
sigdiffs = "sig_diffs.tsv"
pairdists = "pair_distances.tsv"

if splitstoo:
    avgout = "average_sigs_splitstoo.tsv"
    sigdiffs = "sig_diffs_splitstoo.tsv"
    pairdists = "pair_distances_splitstoo.tsv"


class Treeset:
    def __init__(self, treelist, splits=True):
        self.treelist = treelist
        self.branchlist = []
        self.allbranches = []
        self.ntries = ntries
        index = 0
        indexmap = {}
        for tree in treelist:
            for branch in tree.traverse():
                self.allbranches.append(branch)
                if hasattr(branch, "sigs"):
                    if "-" in branch.name:
                        if splits:
                            if "-1" in branch.name:
                                self.branchlist.append(branch);
                                branch.index = index
                                sampleid = branch.name.split("-")[0]
                                indexmap[sampleid] = index
                                index += 1
                        else:
                            branch.index = -1
                    else:
                        self.branchlist.append(branch)
                        branch.index = index
                        index += 1
                else:
                    branch.index = -1
        if splits:
            for tree in treelist:
                for tip in tree:
                    if "-" in tip.name and "-1" not in tip.name:
                        sampleid = tip.name.split("-")[0]
                        tip.index = indexmap[sampleid]
        #Just to check:
        for n, branch in enumerate(self.branchlist):
            assert(branch.index == n)
        for tree in self.treelist:
            for branch in tree.traverse():
                assert(hasattr(branch, "index"))

    def getUpWithSigs(self, branch):
        up = branch.up
        if up and not(hasattr(up, "sigs")):
            return self.getUpWithSigs(up)
        return up
    
    def calculateAvg(self):
        allsigs = []
        for branch in self.branchlist:
            allsigs.append(branch.sigs)
        return np.average(allsigs, axis=0)
    
    def calculateOneUpTrajectory(self, splitvec=None):
        if not splitvec:
            splitvec = list(range(len(self.branchlist)))
        alldiffs = []
        for tree in alltrees[label]:
            for branch in tree.traverse():
                if not (hasattr(branch, "sigs")):
                    continue
                upbranch = self.getUpWithSigs(branch)
                if upbranch:
                    branchsigs = self.branchlist[splitvec[branch.index]].sigs
                    upsigs = self.branchlist[splitvec[upbranch.index]].sigs
                    alldiffs.append(np.subtract(branchsigs, upsigs))
        return np.sum(alldiffs, axis=0), len(alldiffs)

    def calculateUpTrajectories(self):
        canonical, ndiffs = self.calculateOneUpTrajectory()
        nlarger = np.zeros(len(canonical))
        nsmaller = np.zeros(len(canonical))
        nsame = np.zeros(len(canonical))
        significance = np.zeros(len(canonical))
        shuffled = []
        for n in range(self.ntries):
            splitvec = list(range(len(self.branchlist)))
            random.shuffle(splitvec)
            diffrand, __ = self.calculateOneUpTrajectory(splitvec=splitvec)
            shuffled.append(diffrand)
            for n in range(len(canonical)):
                if diffrand[n] < canonical[n]:
                    nsmaller[n] += 1
                elif diffrand[n] > canonical[n]:
                    nlarger[n] += 1
                else:
                    nsame[n] += 1
        for n in range(len(significance)):
            significance[n] = significanceTwoTailed(nlarger[n], nsmaller[n], nsame[n], self.ntries)
        return canonical, significance, shuffled, ndiffs

    def getAbsDistance(self, b1sigs, b2sigs):
        tot = 0
        for n in range(len(b1sigs)):
            tot += abs(b1sigs[n] - b2sigs[n])
        return tot
    
    def getSqDistance(self, b1sigs, b2sigs):
        tot = 0
        for n in range(len(b1sigs)):
            tot += pow(b1sigs[n] - b2sigs[n], 2)
        return tot
        

    def sortOnePairDistance(self, splitvec=None):
        if not splitvec:
            splitvec = list(range(len(self.branchlist)))
        absdistances = {}
        sqdistances = {}
        for d in [absdistances, sqdistances]:
            d["ancestral"] = []
            d["skew"] = []
        for b1 in range(len(self.allbranches)):
            branch1 = self.allbranches[b1]
            if (branch1.index==-1):
                continue
            for b2 in range(b1+1, len(self.allbranches)):
                branch2 = self.allbranches[b2]
                if (branch2.index==-1):
                    continue
                tree = branch1.get_tree_root()
                disttype = "skew"
                try:
                    common = tree.get_common_ancestor(branch1, branch2)
                    if common==branch1 or common==branch2:
                        disttype = "ancestral"
                    else:
                        disttype = "skew"
                except:
                    pass
                b1sigs = self.branchlist[splitvec[branch1.index]].sigs
                b2sigs = self.branchlist[splitvec[branch2.index]].sigs
                absdistances[disttype].append(self.getAbsDistance(b1sigs, b2sigs))
                sqdistances[disttype].append(self.getSqDistance(b1sigs, b2sigs))
        distances = {}
        npairs = {}
        for stype in absdistances:
            distances[stype + "_abs"] = np.sum(absdistances[stype])
            npairs[stype + "_abs"] = len(absdistances[stype])
        for stype in sqdistances:
            distances[stype + "_sq"] = np.sum(sqdistances[stype])
            npairs[stype + "_sq"] = len(absdistances[stype])
        return distances, npairs

    def sortPairDistances(self):
        canonical, npairs = self.sortOnePairDistance()
        shuffled = {}
        nlarger = {}
        nsmaller = {}
        nsame = {}
        significance = {}
        for stype in canonical:
            nlarger[stype] = 0
            nsmaller[stype] = 0
            nsame[stype] = 0
            significance[stype] = 0
            shuffled[stype] = []
        for n in range(self.ntries):
            splitvec = list(range(len(self.branchlist)))
            random.shuffle(splitvec)
            distrand, __ = self.sortOnePairDistance(splitvec=splitvec)
            for stype in canonical:
                if distrand[stype] < canonical[stype]:
                    nsmaller[stype] += 1
                elif distrand[stype] > canonical[stype]:
                    nlarger[stype] += 1
                else:
                    nsame[stype] += 1
                shuffled[stype].append(distrand[stype])
        for stype in significance:
            significance[stype] = significanceTwoTailed(nlarger[stype], nsmaller[stype], nsame[stype], self.ntries)
        return canonical, significance, shuffled, npairs

    def getNBranches(self):
        return len(self.branchlist)

def significanceOneTailed(nlarger, nsmaller, nsame, ntries):
    same_rnd = random.randint(0,nsame)
    pos = min(nlarger, nsmaller) + same_rnd
    if pos*2>ntries:
        pos = ntries-pos
    return pos/ntries

def significanceTwoTailed(nlarger, nsmaller, nsame, ntries):
    return 2*significanceOneTailed(nlarger, nsmaller, nsame, ntries)




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


def makeTreesets(alltrees, splits=True):
    treesets = {}
    for label in alltrees:
        if "all" in label:
            #LS DEBUG:  REMOVE THIS WHEN WE GET T3 DATA!!
            continue
        treesets[label] = Treeset(alltrees[label], splits=splits)
    return treesets




    


alltrees = readTrees()
sigids = readAndStoreSignatures(alltrees)
treesets = makeTreesets(alltrees, splits=splitstoo)

avgfile = open(avgout, "w")
avgfile.write("Patient")
for sigid in sigids:
    avgfile.write("\t" + sigid)
avgfile.write("\n")

difffile = open(sigdiffs, "w")
difffile.write("Patient")
for sigid in sigids:
    difffile.write("\t" + sigid)
    difffile.write("\t" + sigid + " significance")
difffile.write("\n")

stypes = ("ancestral_abs", "skew_abs", "ancestral_sq", "skew_sq")

pairfile = open(pairdists, "w")
pairfile.write("Patient")
for stype in stypes:
    pairfile.write("\t" + stype)
    pairfile.write("\t" + stype + " significance")
pairfile.write("\n")

pairs_canonical_overall = {}
shuffled_pairs_overall = {}
for stype in stypes:
    pairs_canonical_overall[stype] = {}
    pairs_canonical_overall[stype]["average"] = 0
    pairs_canonical_overall[stype]["sq_av"] = 0
    pairs_canonical_overall[stype]["sum"] = 0
    
    shuffled_pairs_overall[stype] = {}
    shuffled_pairs_overall[stype]["average"] = np.zeros(ntries)
    shuffled_pairs_overall[stype]["sq_av"] = np.zeros(ntries)
    shuffled_pairs_overall[stype]["sum"] = np.zeros(ntries)


canon_trajects = {}
canon_trajects["average"] = np.zeros(10)
canon_trajects["sum"] = np.zeros(10)

shuffled_trajects = {}
shuffled_trajects["average"] = np.zeros((ntries,10))
shuffled_trajects["sum"] = np.zeros((ntries,10))

for label in treesets:
    #calculate the average signature per patient:
    print("Processing", label)
    avgfile.write(label)
    treeset = treesets[label]
    sigavgs = treeset.calculateAvg()
    for sigavg in sigavgs:
        avgfile.write("\t" + str(sigavg))
    avgfile.write("\n")
    
    #Calculate the trajectory of every signature in the patient, i.e.
    # test whether it goes up or down as you go from root to tip.
    #
    #The 'significance' is the percentage of randomized trials
    # that produced the same or a more extreme value for that signature.
    # A value of '1' means all the values were the same.
    difffile.write(label)
    canonical_traj, signif_traj, shuffled_traj, ndiffs = treeset.calculateUpTrajectories()
    for nsig in range(len(canonical_traj)):
        difffile.write("\t" + str(canonical_traj[nsig]))
        difffile.write("\t" + str(signif_traj[nsig]))
        canon_trajects["average"][nsig] += canonical_traj[nsig] / ndiffs
        canon_trajects["sum"][nsig] += canonical_traj[nsig]
        for repeat in range(treeset.ntries):
            shuffled_trajects["average"][repeat][nsig] += shuffled_traj[repeat][nsig]/ndiffs
            shuffled_trajects["sum"][repeat][nsig] += shuffled_traj[repeat][nsig]
    difffile.write("\n")
    
    #Calculate the average distance between signatures on the three types of branch pairs
    pairfile.write(label)
    pairsums, pair_significance, shuffled, npairs = treeset.sortPairDistances()
    for stype in pairsums:
        pairfile.write("\t" + str(pairsums[stype]))
        pairfile.write("\t" + str(pair_significance[stype]))
        pairs_canonical_overall[stype]["average"] += pairsums[stype] / npairs[stype]
        pairs_canonical_overall[stype]["sq_av"] += pairsums[stype] / np.sqrt(npairs[stype])
        pairs_canonical_overall[stype]["sum"] += pairsums[stype]
        for repeat in range(treeset.ntries):
            shuffled_pairs_overall[stype]["average"][repeat] += shuffled[stype][repeat] / npairs[stype]
            shuffled_pairs_overall[stype]["sq_av"][repeat] += shuffled[stype][repeat] / np.sqrt(npairs[stype])
            shuffled_pairs_overall[stype]["sum"][repeat] += shuffled[stype][repeat]
    pairfile.write("\n")
        
    
avgfile.close()
difffile.close()
pairfile.close()


print("Overall significance of pairwise comparisons:")
for stype in stypes:
    print("Averages:")
    nsmaller  = 0
    nlarger = 0
    nsame = 0
    for repeat in range(ntries):
        if shuffled_pairs_overall[stype]["average"][repeat] < pairs_canonical_overall[stype]["average"]:
            nsmaller += 1
        elif shuffled_pairs_overall[stype]["average"][repeat] > pairs_canonical_overall[stype]["average"]:
            nlarger += 1
        else:
            nsame += 1
    significance = significanceOneTailed(nlarger, nsmaller, nsame, ntries)
    print(stype, ": canonical distance = ", str(pairs_canonical_overall[stype]["average"]), "avg shuffled distance =",  str(np.average(shuffled_pairs_overall[stype]["average"])), "Significance:", str(significance))

    print("\nSquare root 'average':")
    nsmaller  = 0
    nlarger = 0
    nsame = 0
    for repeat in range(ntries):
        if shuffled_pairs_overall[stype]["sq_av"][repeat] < pairs_canonical_overall[stype]["sq_av"]:
            nsmaller += 1
        elif shuffled_pairs_overall[stype]["sq_av"][repeat] > pairs_canonical_overall[stype]["sq_av"]:
            nlarger += 1
        else:
            nsame += 1
    significance = significanceOneTailed(nlarger, nsmaller, nsame, ntries)
    print(stype, ": canonical distance = ", str(pairs_canonical_overall[stype]["sq_av"]), "avg shuffled distance =",  str(np.average(shuffled_pairs_overall[stype]["sq_av"])), "Significance:", str(significance))

    print("\nSums:")
    nsmaller  = 0
    nlarger = 0
    nsame = 0
    for repeat in range(ntries):
        if shuffled_pairs_overall[stype]["sum"][repeat] < pairs_canonical_overall[stype]["sum"]:
            nsmaller += 1
        elif shuffled_pairs_overall[stype]["sum"][repeat] > pairs_canonical_overall[stype]["sum"]:
            nlarger += 1
        else:
            nsame += 1
    significance = significanceOneTailed(nlarger, nsmaller, nsame, ntries)
    print(stype, ": canonical distance = ", str(pairs_canonical_overall[stype]["sum"]), "avg shuffled distance =",  str(np.average(shuffled_pairs_overall[stype]["sum"])), "Significance:", str(significance))


print("\nOverall significance of trajectory comparisons:")
print("Sums:")
for nsig in range(10):
    nsmaller  = 0
    nlarger = 0
    nsame = 0
    print("\nResults for signature", sigids[nsig], ":")
    for repeat in range(ntries):
        if shuffled_trajects["sum"][repeat][nsig] < canon_trajects["sum"][nsig]:
            nsmaller += 1
        elif shuffled_trajects["sum"][repeat][nsig] > canon_trajects["sum"][nsig]:
            nlarger += 1
        else:
            nsame += 1
    significance = significanceOneTailed(nlarger, nsmaller, nsame, ntries)
    print("Total canonical distance = ", str(canon_trajects["sum"][nsig]), "shuffled total distance =",  str(np.average(shuffled_trajects["sum"], axis=0)[nsig]), "Significance:", str(significance))

print("Average across patients:")
for nsig in range(10):
    nsmaller  = 0
    nlarger = 0
    nsame = 0
    print("\nResults for signature", sigids[nsig], ":")
    for repeat in range(ntries):
        if shuffled_trajects["average"][repeat][nsig] < canon_trajects["average"][nsig]:
            nsmaller += 1
        elif shuffled_trajects["average"][repeat][nsig] > canon_trajects["average"][nsig]:
            nlarger += 1
        else:
            nsame += 1
    significance = significanceOneTailed(nlarger, nsmaller, nsame, ntries)
    print("Average canonical distance = ", str(canon_trajects["average"][nsig]), "shuffled avg distance =",  str(np.average(shuffled_trajects["average"], axis=0)[nsig]), "Significance:", str(significance))

