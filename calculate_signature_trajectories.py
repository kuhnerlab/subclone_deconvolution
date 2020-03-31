#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:54:07 2020

@author: lpsmith
"""

from __future__ import division
from os import walk
from ete3 import Tree

import numpy as np
import random

onlysomepatients = False
somepatients = ["1005"]
#somepatients = ["387", "483", "609", "611", "626", "631", "635", "652", "852"] #double split patients

#input files and directories:
treedir = "branch_lengths/"
sigfile = "branch_signatures/all_signatures.tsv"

ntries = 100000

#output files:
avgout = "average_sigs.tsv"
sigdiffs = "sig_diffs.tsv"
pairdists = "pair_distances.tsv"
tvifilename = "tip_vs_internal_dist.tsv"


class Treeset:
    """
    Hey, an actual python class!  w00t!
    So: this class takes and stores a list of trees, and provides various
    ways of analyzing the mutation signatures on those trees' branches.
    When created, you declare whether you want to include the 'split tips' 
    in your analysis, i.e. mutation signatures from mutations that were
    private to a particular sample, but which showed up as multiple tips 
    in the trees.  While developing these algorithms, it wasn't entirely
    clear what we wanted to do with them; in the end, we almost always
    just drop them.  Dropping them is what happens by default.
    """
    def __init__(self, treelist, splits=False):
        """
        Set up a Treeset with a list of trees, either dropping the 'split
        tips' (the default; splits=False) or including them.
        """
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
        """
        Finds the nearest ancestor of the passed-in branch that actually
        has mutation signatures associated with it.  (In at least one case,
        there's an internal branch of length 0 with no mutations associated
        with it.)
        """
        up = branch.up
        if up and not(hasattr(up, "sigs")):
            return self.getUpWithSigs(up)
        return up
    
    def calculateAvg(self):
        """
        Returns the average signature profile of all branches in these
        trees.
        """
        allsigs = []
        for branch in self.branchlist:
            allsigs.append(branch.sigs)
        return np.average(allsigs, axis=0)
    
    def excludeFrom(self, vec, excludes):
        """
        Once we know that one signature's trajectory is significant, we then
        want to remove it from subsequent analyses so that we can see if
        other signatures also change their trajectory, without being 
        enhanced or depressed by the main signal.
        """
        if len(excludes)==0:
            return vec
        newvec = []
        for n in range(len(vec)):
            if n not in excludes:
                newvec.append(vec[n])
        return np.divide(newvec, np.sum(newvec))
    
    def calculateOneUpTrajectory(self, excludes, splitvec=None):
        """
        Calculate the difference between a child mutation signature and its
        parent's mutation signature.  If 'splitvec' is not set, it's created
        as a straight 0-based vector. The splitvec is used to shuffle the
        branches when performing the analysis: you grab the signatures from
        branchlist[splitvec[branch.index]] instead of 'branch'.  Because a 
        branch's index is its position in branchlist, if 'splitvec' is
        [0, 1, 2, ...], you get 'branch' back.  In all other cases, you
        get a random branch with a signature.
        """
        if not splitvec:
            splitvec = list(range(len(self.branchlist)))
        alldiffs = []
        for tree in self.treelist:
            for branch in tree.traverse():
                if not (hasattr(branch, "sigs")):
                    continue
                upbranch = self.getUpWithSigs(branch)
                if upbranch:
                    branchsigs = self.branchlist[splitvec[branch.index]].sigs
                    upsigs = self.branchlist[splitvec[upbranch.index]].sigs
                    branchsigs = self.excludeFrom(branchsigs, excludes)
                    upsigs = self.excludeFrom(upsigs, excludes)
                    alldiffs.append(np.subtract(branchsigs, upsigs))
        return np.sum(alldiffs, axis=0), len(alldiffs)

    def calculateUpTrajectories(self, excludes):
        """
        Calculate the canonical 'up trajectories', and then compare that
        to a randomized tree, where every branch gets a random mutation
        signature instead of its own.
        """
        canonical, ndiffs = self.calculateOneUpTrajectory(excludes)
        nlarger = np.zeros(len(canonical))
        nsmaller = np.zeros(len(canonical))
        nsame = np.zeros(len(canonical))
        significance = np.zeros(len(canonical))
        shuffled = []
        for n in range(self.ntries):
            splitvec = list(range(len(self.branchlist)))
            random.shuffle(splitvec)
            diffrand, __ = self.calculateOneUpTrajectory(excludes, splitvec=splitvec)
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
        """
        Calculate the sum of the absolute differences between the two vectors.
        """
        tot = 0
        for n in range(len(b1sigs)):
            tot += abs(b1sigs[n] - b2sigs[n])
        return tot
    
    def getSqDistance(self, b1sigs, b2sigs):
        """
        Calculate the sum of the squares of the differences between the
        two vectors.
        """
        tot = 0
        for n in range(len(b1sigs)):
            tot += pow(b1sigs[n] - b2sigs[n], 2)
        return tot
        

    def sortOnePairDistance(self, splitvec=None):
        """
        Calculate the pairwise distances between every pair of mutation
        signatures on the trees, and store them either as 'ancestral' (when
        one is the direct ancestor of the other) or 'skew' (all other cases,
        including when one branch is on one tree, and the other branch is
        on a completely different tree.)
        As in calculateOneUpTrajectory(), splitvec is used so that the 
        same routine can be used to calculate the true set of distances
        and then can be re-used when shuffling the mutation signatures on the
        tree and re-calculating everything. It relies on the fact that 
        branch.index is equal that branch's position in self.branchlist.
        """
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
        """
        Calculate and sort the pairwise distances between every pair of 
        mutation signatures in the trees, and then compare the canonical
        calculations to the re-sorted distances you get when every signature
        on the tree is shuffled.
        """
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
        """
        Returns the number of brances in the tree with mutation signatures.
        """
        return len(self.branchlist)

    def calcOneTipVsInternalSigs(self, splitvec=None):
        """
        Calculates the overall difference between the tip branch mutation
        signatures and the internal branch mutation signatures.  As in other
        similar functions 'splitvec' allows us to use the same function 
        whether we're calculating the true distance, or whether we're
        calculating what the distance would be if we shuffled all the
        signatures in the trees.
        """
        if not splitvec:
            splitvec = list(range(len(self.branchlist)))
        tipSigTotal = np.zeros(10)
        branchSigTotal = np.zeros(10)
        ntips = 0
        nbranches = 0
        for branch in self.branchlist:
            if branch.is_leaf():
                tipSigTotal = np.add(tipSigTotal, self.branchlist[splitvec[branch.index]].sigs)
                ntips += 1
            else:
                branchSigTotal = np.add(branchSigTotal, self.branchlist[splitvec[branch.index]].sigs)
                nbranches += 1
        if ntips==0 or nbranches==0:
            return np.zeros(10)
        tipSigs = np.divide(tipSigTotal, ntips)
        branchSigs = np.divide(branchSigTotal, nbranches)
        return np.subtract(tipSigs, branchSigs)
    
    def calcTipVsInternalSigs(self):
        """
        Calculates the difference between the tips and the internal branches,
        and then re-calculates this difference with the branch mutation
        signatures shuffled, to determine significance.
        """
        canonical = self.calcOneTipVsInternalSigs()
        tvi_shuffled = []
        significances = []
        for repeat in range(self.ntries):
            splitvec = list(range(len(self.branchlist)))
            random.shuffle(splitvec)
            shuff_dist = self.calcOneTipVsInternalSigs(splitvec)
            tvi_shuffled.append(shuff_dist)
        for nsig in range(len(canonical)):
            nlarger = 0
            nsmaller = 0
            nsame = 0
            for repeat in range(self.ntries):
                if tvi_shuffled[repeat][nsig] < canonical[nsig]:
                    nsmaller += 1
                elif tvi_shuffled[repeat][nsig] > canonical[nsig]:
                    nlarger += 1
                else:
                    nsame += 1
            significance = significanceTwoTailed(nlarger, nsmaller, nsame, self.ntries)
            if significance < 0:
                assert(False)
            significances.append(significance)
        #print(str(significances))
        return canonical, significances, tvi_shuffled
            

def significanceOneTailed(nlarger, nsmaller, nsame, ntries):
    """
    Given a number of observations that are smaller, larger, or the same
    as a given value, calculates the significance of that observation,
    assuming we were expecting the observed value to be in the direction
    it actually was.  This would be super straightforward except for the
    number of ties:  it turns out that both the conservative option (where
    you count the ties against the signficance) and the aggressive option
    (where you count the ties in your favor) is wrong: you must include a
    random number of ties. This is more obvious when you run this algorithm
    with nothing but ties:  the result should be significant exactly 5%
    of the time, which you only get if you choose randomly.
    """
    same_rnd = random.randint(0,nsame)
    pos = min(nlarger, nsmaller) + same_rnd
    if pos*2>ntries:
        pos = ntries-pos
    return pos/ntries

def significanceTwoTailed(nlarger, nsmaller, nsame, ntries):
    """
    When you didn't have an expectation of what direction the true
    value should be in relation to the average, you double the 
    one-tailed test significance.
    """
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
    """
    Reads in the mutation signature file, and assigns the signatures there
    to the branches of the trees.
    """
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
        
    return sigids


def makeTreesets(alltrees, splits=True):
    """
    Takes the trees (which have already been annotated with the mutation
    signatures) and creates Treesets from them.
    """
    treesets = {}
    for label in alltrees:
#        if "all" in label:
#            #LS DEBUG:  REMOVE THIS WHEN WE GET T3 DATA!!
#            continue
        if label + "_all" in alltrees:
            #Skips 4-sample trees if there are 6-sample trees available.
            print("Skipping", label, "since there's a 6-tip version available.")
            continue
        treesets[label] = Treeset(alltrees[label], splits=splits)
    return treesets


def calculateAverages(treesets, sigids):
    """
    Calculates the average signature per patient.
    """
    print("\nCalculating averages...")
    avgfile = open(avgout, "w")
    avgfile.write("Patient")
    for sigid in sigids:
        avgfile.write("\t" + sigid)
    avgfile.write("\n")
    for label in treesets:
        #calculate the average signature per patient:
        avgfile.write(label)
        treeset = treesets[label]
        sigavgs = treeset.calculateAvg()
        for sigavg in sigavgs:
            avgfile.write("\t" + str(sigavg))
        avgfile.write("\n")
    avgfile.close()


def calculateTrajectories(treesets, oldsigids, excludes=()):
    """
    Calculates how much mutations signatures change as you move from the root
    of the tree to the tips.
    """
    sigids = []
    for n in range(len(oldsigids)):
        if n not in excludes:
            sigids.append(oldsigids[n])
        
    print("\n\nCalculating trajectories...")
    difffile = open(sigdiffs, "w")
    difffile.write("Patient")
    for sigid in sigids:
        difffile.write("\t" + sigid)
        difffile.write("\t" + sigid + " significance")
    difffile.write("\n")
    
    canon_trajects = {}
    canon_trajects["average"] = np.zeros(10)
    canon_trajects["sum"] = np.zeros(10)
    
    shuffled_trajects = {}
    shuffled_trajects["average"] = np.zeros((ntries,10))
    shuffled_trajects["sum"] = np.zeros((ntries,10))

    for label in treesets:
#        print("Processing", label, "for trajectories.")
        treeset = treesets[label]
        #Calculate the trajectory of every signature in the patient, i.e.
        # test whether it goes up or down as you go from root to tip.
        #
        #The 'significance' is the percentage of randomized trials
        # that produced the same or a more extreme value for that signature.
        # A value of '1' means all the values were the same.
        difffile.write(label)
        canonical_traj, signif_traj, shuffled_traj, ndiffs = treeset.calculateUpTrajectories(excludes)
        for nsig in range(len(canonical_traj)):
            difffile.write("\t" + str(canonical_traj[nsig]))
            difffile.write("\t" + str(signif_traj[nsig]))
            canon_trajects["average"][nsig] += canonical_traj[nsig] / ndiffs
            canon_trajects["sum"][nsig] += canonical_traj[nsig]
            for repeat in range(treeset.ntries):
                shuffled_trajects["average"][repeat][nsig] += shuffled_traj[repeat][nsig]/ndiffs
                shuffled_trajects["sum"][repeat][nsig] += shuffled_traj[repeat][nsig]
        difffile.write("\n")
    difffile.close()
    print("\nOverall significance of trajectory comparisons:")
    print("Sums:")
    for nsig in range(len(sigids)):
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
        print("Total canonical distance = ", str(canon_trajects["sum"][nsig]), "shuffled total distance =",  str(np.average(shuffled_trajects["sum"], axis=0)[nsig]), "\nSignificance:", str(significance))
    
    print("Average across patients:")
    for nsig in range(len(sigids)):
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
        print("Average canonical distance = ", str(canon_trajects["average"][nsig]), "shuffled avg distance =",  str(np.average(shuffled_trajects["average"], axis=0)[nsig]), "\nSignificance:", str(significance))

def calculateAncestralVsSkew(treesets):
    """
    Calculate whether the signatures are significantly different in ancestral
    lineages vs. skew lineages.
    """
    print("\n\nCalculating ancestral vs. skew...")
    #stypes = ("ancestral_abs", "skew_abs", "ancestral_sq", "skew_sq")
    
#    pairfile = open(pairdists, "w")
#    pairfile.write("Patient")
#    for stype in stypes:
#        pairfile.write("\t" + stype)
#        pairfile.write("\t" + stype + " significance")
#    pairfile.write("\n")
    
    pairs_canonical_overall_diff = 0
    pairs_canonical_overall_anc = 0
    pairs_canonical_overall_skew = 0
    shuffled_pairs_overall_diff = np.zeros(ntries)
    shuffled_pairs_overall_anc = np.zeros(ntries)
    shuffled_pairs_overall_skew = np.zeros(ntries)
    for label in treesets:
        #Calculate the average distance between signatures on the three types of branch pairs
        #pairfile.write(label)
        print("Processing", label, "for ancestral vs. skew analysis.")
        treeset = treesets[label]
        pairsums, pair_significance, shuffled, npairs = treeset.sortPairDistances()
#            pairfile.write("\t" + str(pairsums[stype]))
#            pairfile.write("\t" + str(pair_significance[stype]))
        pairs_canonical_overall_diff += (pairsums["ancestral_abs"]/npairs["ancestral_abs"]) - (pairsums["skew_abs"] / npairs["skew_abs"])
        pairs_canonical_overall_anc += (pairsums["ancestral_abs"]/npairs["ancestral_abs"])
        pairs_canonical_overall_skew += (pairsums["skew_abs"] / npairs["skew_abs"])
        for repeat in range(treeset.ntries):
            shuffled_pairs_overall_diff[repeat] += (shuffled["ancestral_abs"][repeat]/npairs["ancestral_abs"]) - (shuffled["skew_abs"][repeat] / npairs["skew_abs"])
            shuffled_pairs_overall_anc[repeat] += (shuffled["ancestral_abs"][repeat]/npairs["ancestral_abs"]) 
            shuffled_pairs_overall_skew[repeat] += (shuffled["skew_abs"][repeat] / npairs["skew_abs"])
 #       pairfile.write("\n")
 #   pairfile.close()
    print("Overall significance of pairwise comparisons of ancestral vs. skew (diff, ancestral, then skew):")
    for (canon, shuffled) in ((pairs_canonical_overall_diff, shuffled_pairs_overall_diff), (pairs_canonical_overall_anc, shuffled_pairs_overall_anc), (pairs_canonical_overall_skew, shuffled_pairs_overall_skew)):
        nsmaller  = 0
        nlarger = 0
        nsame = 0
        for repeat in range(ntries):
            if shuffled[repeat] < canon:
                nsmaller += 1
            elif shuffled[repeat] > canon:
                nlarger += 1
            else:
                nsame += 1
        significance = significanceOneTailed(nlarger, nsmaller, nsame, ntries)
        print("Canonical = ", str(canon), "Avg shuffled =",  str(np.average(shuffled)), "\nSignificance:", str(significance))


def calculateTipVsInternal(treesets):
    """
    Calculate whether the signatures are significantly different in the tips 
    vs. in the internal branches.
    """
    print("\n\nCalculating differences between tip and internal signatures.")
    canon_tvis = []
    shuffled_tvis = np.zeros((ntries,10))
    
    tvifile = open(tvifilename, "w")
    tvifile.write("Patient")
    for sigid in sigids:
        tvifile.write("\t" + sigid + " diff")
        tvifile.write("\t" + sigid + " significance")
    tvifile.write("\n")
    for label in treesets:
#        print("Processing", label, "for tip vs. internal analysis.")
        treeset = treesets[label]
        ###Calculate the differences in tip signatures vs. branch signatures:
        tvi, tvi_sig, tvi_shuff = treeset.calcTipVsInternalSigs()
        canon_tvis.append(tvi)
        shuffled_tvis = np.add(shuffled_tvis, tvi_shuff)
        tvifile.write(label)
        for nsig in range(len(sigids)):
            tvifile.write("\t" + str(tvi[nsig]))
            tvifile.write("\t" + str(tvi_sig[nsig]))
        tvifile.write("\n")
    tvifile.close()
    print("Overall significance of tip vs. internal comparisons:")
    print("(Note:  average shuffled values are always zero.)")
    canon_tvi_avg = np.average(canon_tvis, axis=0)
    for nsig in range(10):
        nlarger = 0
        nsmaller = 0
        nsame = 0
        for repeat in range(ntries):
            if shuffled_tvis[repeat][nsig] > canon_tvi_avg[nsig]:
                nlarger += 1
            elif shuffled_tvis[repeat][nsig] < canon_tvi_avg[nsig]:
                nsmaller += 1
            else:
                nsame += 1
        significance = significanceTwoTailed(nlarger, nsmaller, nsame, ntries)
        print(sigids[nsig], 
              " average = ", 
              str(canon_tvi_avg[nsig]), 
              "stdev shuffled =", 
              str(np.std(shuffled_tvis, axis=1)[nsig]), 
              "significance = ", 
              str(significance))
        print("Larger", str(nlarger), "Smaller", str(nsmaller), "Same", str(nsame))
    
    


#First, do all analyses that require no split tips in them.
alltrees = readTrees()
sigids = readAndStoreSignatures(alltrees)
treesets_nosplits = makeTreesets(alltrees, splits=False)

calculateAverages(treesets_nosplits, sigids)
#calculateTrajectories(treesets_nosplits, sigids)
#calculateTrajectories(treesets_nosplits, sigids, excludes=[7])
calculateAncestralVsSkew(treesets_nosplits)


#Next, do all analyses that can include the split tips.  Have to re-read
# the trees themselves, since the signatures were already stored on them
# in the first 'readAndStoreSignatures'.
alltrees = readTrees()
sigids = readAndStoreSignatures(alltrees)
treesets_splits = makeTreesets(alltrees, splits=True)
calculateTipVsInternal(treesets_splits)
