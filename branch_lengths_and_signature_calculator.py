# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:43:41 2018

@author: Lucian

This program takes a list of mutants, the hand-created Newick trees, and the
hand-created 'optim' files, and assigns shared mutations to the internal
branches of the tree, for two purposes:  to determine branch lengths, and
to group mutations by branch so that their mutation signatures can be
analyzed.  In some cases, mutations can be used to determine branch length,
but not mutation signatures.
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
from ete3 import Tree, TreeStyle, TextFace, AttrFace

import csv

import lucianSNPLibrary as lsl
#Alternative for importing the 'lucianSNPLibrary':
#import imp
#lsl = imp.load_source("lsl","/path/to/lucianSNPLibrary.py")

#Flags for changing what the program outputs:
showTrees = True
needDoubleSplitInfo = False

onlysomepatients = False
somepatients = ["184"]
#somepatients = ["387", "483", "609", "611", "626", "631", "635", "652", "852"] #double split patients

#input files and directories:
mutation_file = "all_SNVs.csv"
treedir = "final_trees/"
optimdir = "final_optim/"

#output directories:
lengthdir = "branch_lengths/"
sigdir = "branch_signatures/"


if not path.isdir(lengthdir):
    mkdir(lengthdir)

if not path.isdir(sigdir):
    mkdir(sigdir)

def getCNVCall(patient, sample, chrom, pos, CNVs):
    """
    Takes a given patient, sample, chromosome, and position, and 
    returns the CNV call present in the 'CNVs' object.
    """
    if patient not in CNVs:
        assert(False)
        return (-1, -1)
    if sample not in CNVs[patient]:
        assert(False)
        return (-1, -1)
    if chrom not in CNVs[patient][sample]:
        print("No chromosome", str(chrom), "found.")
        assert(False)
        return (-1, -1)
    for (start, end, call) in CNVs[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return call
    return (-1, -1)

def isDeleted(patient, sample, chrom, pos, deletions):
    """
    Takes a given patient, sample, chromosome, and position, and 
    returns whether there is a deletion present there.
    """
    if patient not in deletions:
        return False
    if sample not in deletions[patient]:
        return False
    if chrom not in deletions[patient][sample]:
        return False
    for (start, end) in deletions[patient][sample][chrom]:
        if start <= pos and end >= pos:
            return True
    return False

def getSampleStatuses():
    """
    Takes the file that indicates various information about the samples and
    packages it up as a python dictionary.
    """
    statuses = {}
    with open("P01CA91955-WGS80-Full-Pilot-Samples.csv", "r") as csvfile:
        for lvec in csv.reader(csvfile):
            if "RandomID" in lvec:
                continue
            sample = lvec[6]
            statuses[sample] = {}
            statuses[sample]["patient"] = lvec[0]
            statuses[sample]["prog"] = lvec[1]
            statuses[sample]["gender"] = lvec[2]
            time = lvec[8]
            if time=="LongFollowUp":
                time = "T3"
            if time=="Index":
                time = "T2"
            if sample == "23521":
                time = "T2"
            statuses[sample]["time"] = time
    return statuses

def readMutations():
    """
    Reads in the mutations from mutation_file and saves the VAFs as a
    dictionary of the form mutation[patient][chr][pos][alt][sample].
    """
    mutations = {}
    with open(mutation_file, 'r') as csvfile:
        for lvec in csv.reader(csvfile):
            if "DNANum" in lvec[0]:
                continue
            (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
            if (is_snv=="f"):
                continue
            if (is_2p=="f"):
                continue
            if ("N" in sample):
                continue
            if sample not in samplePatientMap:
                continue
            refcnt = int(lvec[-2])
            bafcnt = int(lvec[-1])
            VAF = bafcnt/(refcnt+bafcnt)
            patient = samplePatientMap[sample][0]
    #        if patient not in patientSampleMap:
    #            patientSampleMap[patient] = set()
    #        patientSampleMap[patient].add(sample)
            if onlysomepatients and patient not in somepatients:
                continue
            if patient not in mutations:
                mutations[patient] = {}
                print("Reading mutations from patient", patient)
            if chr not in mutations[patient]:
                mutations[patient][chr] = {}
            pos = int(pos)
            if pos not in mutations[patient][chr]:
                mutations[patient][chr][pos] = {}
            if alt not in mutations[patient][chr][pos]:
                mutations[patient][chr][pos][alt] = {}
            mutations[patient][chr][pos][alt][sample] = VAF
    return mutations

def getCodeToSampleMap():
    """
    Since the 'optim' input files refer to samples by a one-letter code, but
    everything else uses the full sample ID (a five-digit number), this
    creates a dictionary that maps the code to the sample IDs, by patient.
    The one-letter code is the letter corresponding to the sorted IDs
    of the samples.
    """
    lists = {}
    codes = {}
    alpha = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    with open("P01CA91955-WGS80-Full-Pilot-Samples.csv", "r") as csvfile:
        for lvec in csv.reader(csvfile):
            if "RandomID" in lvec:
                continue
            sample = lvec[6]
            patient = lvec[0]
            time = lvec[8]
            if time=="LongFollowUp":
                time = "T3"
            if time=="Index":
                time = "T2"
            if sample == "23521":
                time = "T2"
            if time in ["T1", "T2", "T3"]:
                if patient not in lists:
                    lists[patient] = []
                lists[patient].append(sample)
    for patient in lists:
        lists[patient].sort()
        codes[patient] = {}
        for n in range(len(lists[patient])):
            codes[patient][alpha[n]] = lists[patient][n]
    return codes

def readOptimInputs(codeToSampleMap):
    """
    Converts the 'optim' input files to a dictionary that can be used in
    the rest of the program.
    
    The 'optim' files were originally created to calculate cell fractions from
    the mutation data, but turned out to be useful for this purpose as well.
    """
    allsets = {}
    labelSamples = {}
    optimfiles = []
    for __, _, files in walk(optimdir):
        optimfiles += files
    for file in optimfiles:
        if "input" not in file:
            continue
        fvec = file.split("_")
        patient = fvec[0]
        if onlysomepatients and patient not in somepatients:
            continue
        tag = ""
        if patient + "_input_only4" in optimfiles and len(fvec)<3:
            tag = "_all"
        label = patient + tag
        labelSamples[label] = set()
        allsets[label] = {}
        #print(file, label)

        oneset = set()
        inputs = {}
        for line in open(optimdir + file, "r"):
            lvec = line.rstrip().split(" ")
            if len(lvec)==1:
                if len(oneset) > 0:
                    oneset = list(oneset)
                    oneset.sort()
                    oneset = tuple(oneset)
                    allsets[label][oneset] = inputs
                oneset = set()
                inputs = {}
                continue
            if "#" in line:
                continue
            if lvec[0] not in codeToSampleMap[patient]:
                assert(lvec[1] == "32" and lvec[2] == "799")
                lvec[0] = "c"
                lvec[3] = "x4"
            sample = codeToSampleMap[patient][lvec[0]]
            oneset.add(sample)
            labelSamples[label].add(sample)
            vaf = int(lvec[1])
            nmuts = int(lvec[2])
            tips = lvec[3]
            call = "1,1"
            if len(lvec)>4:
                call = lvec[-1]
            assert call != ''
            tipvec = tips.split("_")
            tipset = set()
            for entry in tipvec:
                if "x" in entry:
                    tipset.add(entry)
            if sample not in inputs:
                inputs[sample] = {}
            tipset = list(tipset)
            tipset.sort()
            tipset = tuple(tipset)
            if tipset not in inputs[sample]:
                inputs[sample][tipset] = {}
                
            if call not in inputs[sample][tipset]:
                inputs[sample][tipset][call] = {}
                inputs[sample][tipset][call]["vaf"] = vaf
                inputs[sample][tipset][call]["nmuts"] = nmuts
            else:
                inputs[sample][tipset][call]["vaf"] = [inputs[sample][tipset][call]["vaf"], vaf]
                inputs[sample][tipset][call]["nmuts"] = [inputs[sample][tipset][call]["nmuts"], nmuts]
        if len(oneset) > 0:
            oneset = list(oneset)
            oneset.sort()
            oneset = tuple(oneset)
            allsets[label][oneset] = inputs
        if len(allsets[label]) == 1:
            print(label, "might not have spacing in it.")
    return allsets, labelSamples

def getOptimInputSplits(allsets):
    """
    Takes the output of 'readOptimInputs' and analyzes them for places where
    a single partition exists where the VAFs of two samples are split.  It is
    necessary to find these partitions because it is much more complicated to
    assign these mutations to branches:  the VAFs of the mutation in multiple
    samples much be considered.
    """
    inputSplits = {}
    for label in allsets:
        for group in allsets[label]:
            for sample in allsets[label][group]:
                if len(allsets[label][group][sample]) > 1:
                    #assert len(allsets[label][group][sample]) == 2 ##Actually get 3 for 360!
                    if (label, group) not in inputSplits:
                        inputSplits[(label, group)] = {}
                    else:
                        print("Double splits for", label)
                    inputSplits[(label, group)][sample] = {}
                    for tips in allsets[label][group][sample]:
                        for call in allsets[label][group][sample][tips]:
                            vaf = allsets[label][group][sample][tips][call]["vaf"]
                            nmuts = allsets[label][group][sample][tips][call]["nmuts"]
                            if call not in inputSplits[(label, group)][sample]:
                                inputSplits[(label, group)][sample][call] = {}
                            inputSplits[(label, group)][sample][call][tips] = (vaf, nmuts)
        
    return inputSplits

def readSampleTrees(label):
    """
    Read in the hand-generated Newick trees for all patients.  In many cases,
    multiple trees are present for a single patient, because fewer than 100
    mutations arose on the common ancestor of all the samples' subclones.
    """
    treefile = treedir + label + ".newick"
    newicks = []
    for line in open(treefile, "r"):
        if "#" in line:
            continue
        newicks.append(line.rstrip())
#    halfindexes = range(round(len(newicks)/2))
#    newicks = numpy.delete(newicks, halfindexes)
    trees = []
    samples = set()
    for newick in newicks:
        tree = Tree(newick)
        for branch in tree.traverse():
            branch.dist = 0
        for branch in tree:
            samples.add(branch.name.split('-')[0].split('_')[0])
        trees.append(tree)
        if onlysomepatients:
            print("Tree(s) for", label)
            print(tree)
    return trees, samples

def isSplit(group, inputSplits, label):
    """
    Check if a particular partition (group) needs to be analyzed differently,
    since its inputs are split.
    """
    if (label, group) in inputSplits:
        return True
    return False

def getBranch(trees, tipnames):
    """
    Given a set of trees from a single patient, and a set of tip names,
    get the common ancestor of those tips.  Note that this ancestor may
    be the ancestor of more tips than are present in 'tipnames': this
    is actually common due to 'dropout' mutations, where shared mutations
    were observed in some, but not all, of the samples that inherit from
    the originating branch.
    
    Note that if there was only one tip from a particular sample, the
    code that generates 'tipnames' only knows that the tip it's looking
    for is the first one, and calls it 'sample-1'.  However, in the
    newick tree itself, the sample will just be called 'sample'.
    """
    tipbranches = []
    for tree in trees:
        for tip in tree:
            for tipname in tipnames:
                if tipname in tip.name:
                    tipbranches.append(tip)
                elif "-1" in tipname:
                    if tip.name.split("_")[0] + "-1" == tipname:
                        tipbranches.append(tip)
        if len(tipbranches) == len(tipnames):
            if len(tipnames) == 1:
                return tipbranches[0]
            else:
                return tree.get_common_ancestor(tipbranches)
        if (len(tipbranches) != 0):
            print("Unable to find the cluster:", str(tipnames))
            print("In the trees:")
            for treeB in trees:
                print(treeB)
        assert(len(tipbranches) == 0)
        tipbranches.clear()
    print("Unable to find the cluster:", str(tipnames))
    print("In the trees:")
    for tree in trees:
        print(tree)
    assert(False)
    return None

def getBranchFrom(trees, xtipnames, samplename):
    """
    In the optim file, the particular tips associated with each sample are
    labeled as 'x1', 'x2' etc.  This converts this nomenclature to the
    form used for the actual newick trees, i.e. 'sample_1', etc.
    """
    tipnames = []
    for xtip in xtipnames:
        tipnames.append(combineSampleAndTipLabel(samplename, xtip))
    return getBranch(trees, tipnames)

def getClosestBranch(trees, tipset, basebranch, samplename):
    """
    This function finds the common ancestor of the tipset, but then moves
     rootward in the tree to get as close to the 'basebranch' as possible
     without actually becoming the basebranch.  This is necessary in some
     situations when a split is found where other descendants are also in
     the picture.
    """
    branch = getBranchFrom(trees, tipset, samplename)
    isAncestor = False
    for testbranch in basebranch.traverse():
        if testbranch == branch:
            isAncestor = True
            break
    assert(isAncestor)
    while(True):
        ancestor = branch.up
        if ancestor==basebranch:
            return branch
        else:
            branch = ancestor


def getSigFileFrom(sigfilename, files):
    """
    We need to write to various files all at once.  This function creates
    a file if we need it, and stores it in 'files' so that it can be
    written to as needed later.
    """
    if sigfilename not in files:
        sigfile = open(sigdir + sigfilename, "w")
        files[sigfilename] = sigfile
        sigfile.write("Chrom\tpos\talt\n")
    return files[sigfilename]

def getSigFile(branch, label, files):
    """
    Here, we figure out what the name of the file is that we need to write
    to, then calls getSigFileFrom() to get the actual file handle.
    """
    tipnames = []
    for tip in branch:
        tipnames.append(tip.name)
    tipnames.sort()
    sigfilename = label
    for tipname in tipnames:
        sigfilename += "." + tipname
    sigfilename += ".mutations.tsv"
    return getSigFileFrom(sigfilename, files)

def closeFiles(files):
    """
    Closes all the file handles in the dictionary.
    """
    for filename in files:
        files[filename].close()

def isDeletedInNonGroupSample(group, treesamples, label, chrom, pos):
    """
    Here, we know that we have a position that was found to be mutant in 
    several samples, but we don't know if it was deleted in other samples.
    This checks to see if this is the case.
    """
    patient = label.split('_')[0]
    for sample in treesamples:
        if sample not in group:
            #Check to see if it was deleted in that sample:
            if isDeleted(patient, sample, chrom, pos, deletions):
                return True
    return False

def saveTips(group, label, trees, files, chrom, pos, alt):
    """
    This writes out mutations that are known to be in a tip, but not
    necessarily *which* tip, to a file for mutation signature analysis
    purposes.
    """
    assert(len(group)==1)
    sigfilename = label
    for tipname in group:
        sigfilename += "." + tipname
    sigfilename += ".mutations.tsv"
    sigfile = getSigFileFrom(sigfilename, files)
    sigfile.write(chrom)
    sigfile.write("\t" + str(pos))
    sigfile.write("\t" + alt)
    sigfile.write("\n")
    tips = []
    for tree in trees:
        for tip in tree:
            if tipname in tip.name:
                tips.append(tip)
    ntips = len(tips)
    for tip in tips:
        tip.dist += 1/ntips
    

def saveNonSplit(label, trees, treesamples, group, allsets, files, chrom, pos, alt):
    """
    In this scenario, we've found a mutation and we know exactly what branch it
    should be on, since there was no 'split' of mutations in this cluster.
    So, we can just save it to the appropriate signature file, and add it to
    the branch length.
    """
    if group not in allsets[label]:
        #At some point, we might save these somewhere, but for now they're just considered errors.
        return False
    if isDeletedInNonGroupSample(group, treesamples, label, chrom, pos):
        return False
    tipnames = set()
    for sample in allsets[label][group]:
        for tips in allsets[label][group][sample]:
            for tip in tips:
                tipnames.add(combineSampleAndTipLabel(sample, tip))
    #print(label)
    assert(len(tipnames)>0)
    branch = getBranch(trees, tipnames)
    branch.dist += 1
    sigfile = getSigFile(branch, label, files)
    sigfile.write(chrom)
    sigfile.write("\t" + str(pos))
    sigfile.write("\t" + alt)
    sigfile.write("\n")
    return True

def combineSampleAndTipLabel(sample, tiplabel):
    """
    The 'tiplabel' is (from the optim file) in the form 'x1', 'x2', etc.
    This converts this to a newick-format tip label, which is the sample
    name, a dash, and the number.  So, "x1" is converted to something
    like "23466-1".
    """
    nox = tiplabel.replace("x", "-")
    return sample + nox

def saveSingleSplit(label, group, branch, chrom, pos, alt, treesamples, files):
    """
    Once we've determined which particular branch a mutation belongs in
    (after figuring out which branch of the split was more appropriate),
    this function writes it to the appropriate signature file, and increments
    the branch length.
    """
    if isDeletedInNonGroupSample(group, treesamples, label, chrom, pos):
        return False
    #print(label)
    branch.dist += 1
    sigfile = getSigFile(branch, label, files)
    sigfile.write(chrom)
    sigfile.write("\t" + str(pos))
    sigfile.write("\t" + alt)
    sigfile.write("\n")
    return True
    

def saveForSplits(splits, doublesplits, group, sampleVAFs, chrom, pos, alt, CNVs, patient, splitPartitions):
    """
    When going through all the mutations the first time, this function is
     called when a 'split' is found:  when one or more of the samples in the
     partition has a split VAF with two peaks.  Here, we just save that
     information into 'splits' and 'doublesplits' so we can deal with it 
     later.
    """
    for keySample in splitPartitions:
        call = getCNVCall(patient, keySample, chrom, pos, CNVs)
        strcall = str(call[0]) + "," + str(call[1])
        if len(splitPartitions) > 1:
            #Save double splits separately.
            if group not in doublesplits:
                doublesplits[group] = {}
            if keySample not in doublesplits[group]:
                doublesplits[group][keySample] = {}
            if strcall not in doublesplits[group][keySample]:
                doublesplits[group][keySample][strcall] = []
            doublesplits[group][keySample][strcall].append((sampleVAFs[keySample], chrom, pos, alt))
        else:
            if group not in splits:
                splits[group] = {}
            if strcall not in splits[group]:
                splits[group][strcall] = []
            splits[group][strcall].append((sampleVAFs[keySample], chrom, pos, alt))

def getSplitBranchesFrom(trees, tipsets, keySample):
    """
    Find the ancestral branch of each of the tipsets, each of which will be
     assigned some of the mutations assigned to that split.
    """
    assert(len(tipsets) > 1)
    sortedtips = sorted(tipsets)
    shortesttipset = sortedtips[0]
    for tipset in sortedtips:
        if len(tipset) < len(shortesttipset):
            shortesttipset = tipset

    #We need the *next-shortest* tipset, not the longest.  We have three sets, at one point.
    longertipset = shortesttipset
    for tipset in sortedtips:
        if len(tipset) > len(shortesttipset):
            if len(longertipset) == len(shortesttipset):
                longertipset = tipset
            elif len(tipset) < len(longertipset):
                longertipset = tipset

    branches = {}
    for tipset in sortedtips:
        if tipset == shortesttipset:
            continue
        branches[tipset] = getBranchFrom(trees, tipset, keySample)
    branches[shortesttipset] = getClosestBranch(trees, shortesttipset, branches[longertipset], keySample)
#    for tipset in branches:
#        print(keySample, str(tipset), ":")
#        print(branches[tipset])
        
    return branches

def divideLengthFromSplit(lengthOnly, branches):
    """
    When we can't assign a mutation to a particular branch, but we know they
     belong to a subset of branches, we add those mutations to the given
     branches proportionally to their existing branch length (which is the
     the sum of the mutations we *can* assign them.)
    """
    totlen = 0
    for tipset in branches:
        totlen += branches[tipset].dist

    for tipset in branches:
        branch = branches[tipset]
        branch.dist += lengthOnly * (branch.dist/totlen)

def saveSplits(splits, trees, treesamples, files, inputSplits, label):
    """
    After we've assigned all the mutations to the 'splits' object, this routine
     goes through them and parcels out each mutation to the branch it's most
     likely to belong to.  If it's impossible to assign a mutation to a
     particular branch, it can't add it to any branch's mutation signature,
     but its length is still added proportionally to the branches where it
     might be.
    """
    for group in splits:
        lengthOnly = 0
        
        #First, find the key branches
        tipsets = set()
        for call in splits[group]:
            splitinfo = inputSplits[(label, group)]
            keySample = list(splitinfo.keys())[0]
            if call not in splitinfo[keySample]:
                continue
            splitinfo = splitinfo[keySample][call]
            for tips in splitinfo:
                tipsets.add(tips)
        assert(len(tipsets) > 1)
        branches = getSplitBranchesFrom(trees, tipsets, keySample)
        for call in splits[group]:
            vafvec = splits[group][call]
            splitinfo = inputSplits[(label, group)]
            keySample = list(splitinfo.keys())[0]
            if call not in splitinfo[keySample]:
                #We didn't say how to split this call, so the only thing we can do is 
                # add it to the length, but not the signature file
                
                #However, we still need to check if the site has been deleted in a different sample
                for (vaf, chrom, pos, alt) in vafvec:
                    if not isDeletedInNonGroupSample(group, treesamples, label, chrom, pos):
                        lengthOnly += 1
                continue
            splitinfo = splitinfo[keySample][call]
            vnts = []
            muttot = 0
            for tips in splitinfo:
                tipsets.add(tips)
                (vaf, nmut) = splitinfo[tips]
                if isinstance(vaf, list):
                    #The call is unbalanced, leading to two different VAFs for the same set of tips.
                    for i in range(len(vaf)):
                        vnts.append((vaf[i], nmut[i], tips))
                        muttot += nmut[i]
                else:
                    vnts.append((vaf, nmut, tips))
                    muttot += nmut
            vafvec.sort()
            vnts.sort()
            currindex = 0
            currend = 0
            for (vaf, nmut, tips) in vnts:
                currend += round((nmut/muttot)*len(vafvec))
                while currindex < currend and currindex < len(vafvec):
                    (__, chrom, pos, alt) = vafvec[currindex]
                    saveSingleSplit(label, group, branches[tips], chrom, pos, alt, treesamples, files)
                    currindex += 1
            if currend < len(vafvec):
                assert len(vafvec)-currend == 1
                #just in case we don't sum to 100%
                (__, chrom, pos, alt) = vafvec[-1]
                saveSingleSplit(label, group, branches[tips], chrom, pos, alt, treesamples, files)
        divideLengthFromSplit(lengthOnly, branches)
#        if len(group)==4:
#            assert(False)

def divideLengthAmongBranches(divlen, branches, families):
    """
    When a double split mutation can't be assigned to a particular branch,
     its length is instead assigned proportionally to all of the potential
     branches where it might have been.
    """
    totlen = 0
    for tiplists in families:
        totlen += len(families[tiplists])
    for tiplists in families:
        branch = branches[tiplists]
        branch.dist += len(families[tiplists])
        branch.dist += (divlen * len(families[tiplists])/totlen)
    

def ensureEntireGroupIsChildOfBranch(listbasebranch, group, basebranch):
    """
    In exactly one patient, there's a double split where the common ancestor
     of the smallest subgroup is actually *not* the branch you want, because
     there's a non-split tip that needs to be in there as well.  This routine
     ensures that every sample in the group is represented in the returned
     common ancestor.
    """
    if basebranch == listbasebranch:
        return basebranch
    foundset = set()
    for tip in listbasebranch:
        for sample in group:
            if sample in tip.name:
                foundset.add(sample)
    for sample in group:
        if sample not in foundset:
            print("Branch found that doesn't include entire group")
            print(str(foundset))
            print(str(group))
            print(listbasebranch)
            return ensureEntireGroupIsChildOfBranch(listbasebranch.up, group, basebranch)
            
    return listbasebranch

def saveDoubleSplitFamilies(families, trees, treesamples, files, label, lengthOnlyPositions, group):
    """
    Parcel out mutations among all the appropriate branches, to add to the
     signature lists, or, if that is impossible, to the branch lengths.
    """
    alltips = set()
    tipsets = {}
    divlen = len(lengthOnlyPositions)
    for tiplists in families:
        onetipset = set()
        for sample, tips in tiplists:
            for tip in tips:
                alltips.add(combineSampleAndTipLabel(sample, tip))
                onetipset.add(combineSampleAndTipLabel(sample, tip))
        tipsets[tiplists] = onetipset
    basebranch = getBranch(trees, alltips)
    badtiplists = []
    branches = {}
    for tiplists in tipsets:
        listbasebranch = getBranch(trees, tipsets[tiplists])
        if listbasebranch == basebranch and alltips != tipsets[tiplists]:
            divlen += len(families[tiplists])
            badtiplists.append(tiplists)
        else:
            listbasebranch = ensureEntireGroupIsChildOfBranch(listbasebranch, group, basebranch)
            branches[tiplists] = listbasebranch
    for tiplists in badtiplists:
        del tipsets[tiplists]
        del families[tiplists]
    for tiplists in branches:
        branch = branches[tiplists]
        sigfile = getSigFile(branch, label, files)
        for (chrom, pos, alt) in families[tiplists]:
            sigfile.write(chrom)
            sigfile.write("\t" + str(pos))
            sigfile.write("\t" + alt)
            sigfile.write("\n")
    divideLengthAmongBranches(divlen, branches, families)


def removeDeletedPositions(group, treesamples, label, positions):
    """
    Remove all the mutations from positions that were deleted in samples that
     did not have the mutation.
    """
    badpos = []
    for (chrom, pos, alt) in positions:
        if isDeletedInNonGroupSample(group, treesamples, label, chrom, pos):
            badpos.append((chrom, pos, alt))
            
    for location in badpos:
        if isinstance(positions, set):
            positions.remove(location)
        else:
            del positions[location]


def sortDoubleSplits(doublesplits, trees, treesamples, files, inputSplits, label):
    """
    Take the 'doublesplits' result of the initial sorting algorithm, and parse
    them out to the appropriate branches, or, if the particular branch
    cannot be found, at least divvy the lengths appropriately.
    """
    #For the double splits, they're not stored as vectors, so we have to store them.
    for group in doublesplits:
        positions = {}
        lengthOnlyPositions = set()
        for keySample in doublesplits[group]:
            for call in doublesplits[group][keySample]:
                vafvec = doublesplits[group][keySample][call]
                #print(keySample, str(call), str(len(vafvec)))
                splitinfo = inputSplits[(label, group)][keySample]
                if call not in splitinfo:
                    #We didn't say how to split this call, so the only thing we can do is 
                    # add it to the length, but not the signature file
                    for (vaf, chrom, pos, alt) in vafvec:
                        location = (chrom, pos, alt)
                        if location in positions:
                            del positions[location]
                        lengthOnlyPositions.add(location)
                    continue
                splitinfo = splitinfo[call]
                vnts = []
                muttot = 0
                for tips in splitinfo:
                    (vaf, nmut) = splitinfo[tips]
                    if isinstance(vaf, list):
                        #The call is unbalanced, leading to two different VAFs for the same set of tips.
                        for i in range(len(vaf)):
                            vnts.append((vaf[i], nmut[i], tips))
                            muttot += nmut[i]
                    else:
                        vnts.append((vaf, nmut, tips))
                        muttot += nmut
                vafvec.sort()
                vnts.sort()
                currindex = 0
                currend = 0
                for (vaf, nmut, tips) in vnts:
                    currend += round((nmut/muttot)*len(vafvec))
                    while currindex < currend and currindex < len(vafvec):
                        (__, chrom, pos, alt) = vafvec[currindex]
                        currindex += 1
                        location = (chrom, pos, alt)
                        if location in lengthOnlyPositions:
                            continue
                        if location not in positions:
                            positions[location] = []
                        positions[location].append((keySample, tips))
                if currend < len(vafvec):
                    assert len(vafvec)-currend == 1
                    #just in case we don't sum to 100%
                    (__, chrom, pos, alt) = vafvec[-1]
                    if location not in positions:
                        positions[location] = []
                    positions[location].append((keySample, tips))
        families = {}
        
        removeDeletedPositions(group, treesamples, label, positions)
        removeDeletedPositions(group, treesamples, label, lengthOnlyPositions)
            
        for location in positions:
            positions[location].sort()
            alltips = tuple(positions[location])
            if alltips not in families:
                families[alltips] = []
            families[alltips].append(location)
        saveDoubleSplitFamilies(families, trees, treesamples, files, label, lengthOnlyPositions, group)
        if needDoubleSplitInfo:
            print("Double splits for", label, "group", str(group))
            tiplists = list(families.keys())
            tiplists.sort()
            for alltips in tiplists:
                print(str(alltips), "\t", str(len(families[alltips])))
            print("Unusable:\t", len(lengthOnlyPositions))
#    foo()

def writeTrees(trees, label):
    """
    Now that we have branch lengths, write out new Newick trees with the
    branch lengths included.  Note that ete3 does not seem to write out
    the length of the root automatically, so we add it in by hand.
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.branch_vertical_margin = 20
    ts.scale = 0.05
    outtrees = open(lengthdir + label + ".newick", "w")
    for index, tree in enumerate(trees):
        #print(tree)
        #If needed:  round the branch lengths
#        for branch in tree.traverse():
#            branch.dist = round(branch.dist)
        line = tree.write(format=1)
        rootlen = round(tree.get_tree_root().dist)
        line = line[:-1] + ":" + str(rootlen) + ";"
        outtrees.write(line + "\n")
        if (showTrees):
            for branch in tree.traverse():
                if branch.name != "":
                    name_face = AttrFace("name", fsize=30)
                    branch.add_face(name_face, column=0, position="branch-right")
                if branch.dist != 0:
                    if branch.name == "":
                        text_face = TextFace(str(round(branch.dist, 2)), fsize=30)
                        branch.add_face(text_face, column=1, position="branch-top")
                else:
                    print("Branch found of length zero but not a tip in tree", label)
                    print(branch)
            filename = label
            if len(trees)>1:
                filename += "_" + str(index)
            filename += ".png"
            tree.render(lengthdir + filename, tree_style = ts)
            #foo()
    outtrees.close()
    #foo()
    

def sortMutations(mutations, allsets, inputSplits, deletions, CNVs, labelSamples):
    """
    Loop through all the mutations for all the samples in all the patients in
    the study, and parcel them out to three groups:
        * The simple ones that can be assigned to a single branch.
        * The 'split' ones that can be assigned to one of two branches.
        * The 'double split' ones that can be assigned to one of three branches.
    While doing this, we write out three things:
        * The new tree with branch lengths.
        * 'mutation signature' files, which are lists (one per branch) of
          mutations for subsequent mutation signature analysis
        * A summary file that lists how many mutations were stored and how
          many had to be skipped.
    """
    savesort = open(sigdir + "saved_or_skipped.tsv", "w")
    savesort.write("Tree\tsaved\tskipped\tsplit\n")
    for label in allsets:
        print("Working on", label)
        trees, treesamples = readSampleTrees(label)
        patient = label.split("_")[0]
        samples = labelSamples[label]
        muts = mutations[patient]
        splits = {}
        doublesplits = {}
        files = {}
        used = 0
        skipped = 0
        split = 0
        for chrom in muts:
            for pos in muts[chrom]:
                for alt in muts[chrom][pos]:
                    sampleVAFs = muts[chrom][pos][alt]
                    group = list(sampleVAFs.keys())
                    toRemove = []
                    for sample in group:
                        if sample not in samples:
                            toRemove.append(sample)
                    for sample in toRemove:
                        group.remove(sample)
                    if len(group) < 2:
                        if len(group)==1:
                            saveTips(group, label, trees, files, chrom, pos, alt)
                        continue
                    group.sort()
                    group = tuple(group)
                    if isSplit(group, inputSplits, label):
                        saveForSplits(splits, doublesplits, group, sampleVAFs, chrom, pos, alt, CNVs, patient, inputSplits[(label, group)])
                        split += 1
                    else:
                        if saveNonSplit(label, trees, treesamples, group, allsets, files, chrom, pos, alt):
                            used += 1
                        else:
                            skipped += 1
        sortDoubleSplits(doublesplits, trees, treesamples, files, inputSplits, label)
        saveSplits(splits, trees, treesamples, files, inputSplits, label)
        closeFiles(files)
#        print("For", label, ",", str(used), "non-split mutations, and", str(skipped), "skipped mutations.")
        savesort.write(label)
        savesort.write("\t" + str(used))
        savesort.write("\t" + str(skipped))
        savesort.write("\t" + str(split))
        savesort.write("\n")
        writeTrees(trees, label)
    savesort.close()

codeToSampleMap = getCodeToSampleMap()

#Read in the 'optim' input files, and parse them to figure out where the splits are.
allsets, labelSamples = readOptimInputs(codeToSampleMap)
inputSplits = getOptimInputSplits(allsets)

#Read in the deletions/CNVs
(patientSampleMap, samplePatientMap) = lsl.getPatientSampleMap()
deletions, CNVs = lsl.loadDeletionsAndCNVs(samplePatientMap)

#Read in the mutations
mutations = readMutations()

sortMutations(mutations, allsets, inputSplits, deletions, CNVs, labelSamples)
