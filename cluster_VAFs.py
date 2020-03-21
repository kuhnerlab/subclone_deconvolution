# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 12:43:41 2018

@author: Lucian

This script is the start of the subclonal deconvolution analysis.  It takes
 as input a file with all the SNVs ('all_SNVs.csv'), information about where
 deletions were found (from the 'noninteger_processed_CNs/' directory), and
 outputs several files in a new 'VAFclusters/' directory; one per patient.

During the analysis, all mutations are read in from the input file, and
 classified according to what samples that mutation was found in.

Deletions complicate matters, because if the site of a mutation has been
 deleted in a sample where it was not observed, we don't know whether
 it was not observed because that sample's lineage never had it, or whether
 that sample's linage had it at one point, but then lost it due to the
 deletion.  Accordingly, we flag that sample as having an 'unknown' status
 for that mutation, and cluster it in a completely different cluster:
 Normally, if a mutation is only seen in samples A and B, that mutation's
 partition is ("A", "B").  But if that site is deleted in C, it instead
 is classified into the ("A", "B", "-C") partition.

These partitions ended up never being used in subsequent analyses, and
 were dropped entirely.  But the classification scheme is still here
 in case a future analysis includes them afer all.
"""


from __future__ import division
from os import path
from os import mkdir

import csv

import lucianSNPLibrary as lsl

#For quicker runs/testing:
onlysomepatients = False
somepatients = ["160"]

mutation_file = "all_SNVs.csv"
outdir = "VAFclusters/"

if not path.isdir(outdir):
    mkdir(outdir)

def writeAllSampleVAFs(mutations, patientSampleMap, deletions):
    """
    Takes the sorted mutations and a list of the deletions, and outputs
     cluster information to files; one per patient.
     
    'mutations' is a dictionary of the form:
     mutations[patient][sample][chr][pos][alt][sample] = VAF
    """
    for patient in mutations:
        #Collect a set of clusters
        print("Writing data for patient", patient)
        clustercount = {}
        #As we collect clusters, ensure that the deleted sites
        # are accounted for properly.
        for chr in mutations[patient]:
            for pos in mutations[patient][chr]:
                for alt in mutations[patient][chr][pos]:
                    cluster = list(mutations[patient][chr][pos][alt].keys())
                    cluster.sort()
                    for sample in patientSampleMap[patient]:
                        if sample in cluster:
                            continue
                        if isDeleted(patient, sample, chr, pos, deletions):
                            cluster.append("-" + sample)
                            mutations[patient][chr][pos][alt]["-" + sample] = ""
                    cluster = tuple(cluster)
                    if cluster not in clustercount:
                        clustercount[cluster] = 0
                    clustercount[cluster] += 1
        clusterlist = []
        for cluster in clustercount:
#            if clustercount[cluster] > 40:
                clusterlist.append(cluster)
        for sample in patientSampleMap[patient]:
            patientVAFs = open(outdir + patient + "_" + sample + "_VAFs.tsv", "w")
            patientVAFs.write("Patient")
            patientVAFs.write("\tSample")
            patientVAFs.write("\tchr")
            patientVAFs.write("\tpos")
            patientVAFs.write("\talt")
            for cluster in clusterlist:
                if sample in cluster:
                    patientVAFs.write("\t" + str(cluster))
            patientVAFs.write("\n")
            
            for chr in mutations[patient]:
                posvec = list(mutations[patient][chr].keys())
                posvec.sort()
                for pos in posvec:
                    for alt in mutations[patient][chr][pos]:
                        outstr = patient
                        outstr += ("\t" + sample)
                        outstr += ("\t" + chr)
                        outstr += ("\t" + str(pos))
                        outstr += ("\t" + alt)
                        writeVAF = False
                        for cluster in clusterlist:
                            if sample in cluster and sample in mutations[patient][chr][pos][alt]:
                                #We have to write something: VAF or a tab.
                                outstr += "\t"
                                if set(cluster) == set(mutations[patient][chr][pos][alt].keys()):
                                    outstr += str(mutations[patient][chr][pos][alt][sample])
                                    writeVAF = True
                        if writeVAF:
                            patientVAFs.write(outstr + "\n")
        patientVAFs.close()

def isDeleted(patient, sample, chrom, pos, deletions):
    """
    Returns True if the position is deleted in the given sample, 
     False if not.
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


mutations = {}
(__, samplePatientMap) = lsl.getPatientSampleMap()
patientSampleMap = {}

with open(mutation_file, 'r') as csvfile:
    for lvec in csv.reader(csvfile):
        if "DNANum" in lvec[0]:
            continue
        (sample, __, __, chr, pos, ref, alt, is_snv, is_2p) = lvec[0:9]
        if (is_snv=="f"):
            continue
        if (is_2p=="f"):
            continue
    #    if ("N" in sample):
    #        continue
        refcnt = int(lvec[-2])
        bafcnt = int(lvec[-1])
        VAF = bafcnt/(refcnt+bafcnt)
        patient, ploidy = samplePatientMap[sample]
        if patient not in patientSampleMap:
            patientSampleMap[patient] = set()
        patientSampleMap[patient].add(sample)
        if onlysomepatients and patient not in somepatients:
            continue
        if patient not in mutations:
            mutations[patient] = {}
            print("Reached patient", patient)
        if chr not in mutations[patient]:
            mutations[patient][chr] = {}
        pos = int(pos)
        if pos not in mutations[patient][chr]:
            mutations[patient][chr][pos] = {}
        if alt not in mutations[patient][chr][pos]:
            mutations[patient][chr][pos][alt] = {}
        mutations[patient][chr][pos][alt][sample] = VAF

deletions = lsl.loadDeletions(samplePatientMap)
writeAllSampleVAFs(mutations, patientSampleMap, deletions)

