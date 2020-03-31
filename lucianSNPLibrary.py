# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:34:39 2016

@author: lpsmith
"""

#Useful functions for manipulating log2r data

from __future__ import division
import numpy
import matplotlib.pyplot as plt
from os.path import isfile

patientSampleMapFile = "calling_evidence_odds.tsv"
CNdir = "noninteger_processed_CNs/"

def BiweightKernel(t):
    """
    A bounded kernel suitable for smoothing a histogram.
    """
    if (abs(t) > 1.0):
        return 0.0
    return (15.0/16.0) * pow((1.0-pow(t,2)),2)

def addKernelToHistogram(val, weight, histogram, kernelwidth, numpoints, binwidth):
    """
    The heart of the histogram smoother:  takes a given value (val), a weight
    for that value (weight), and adds a kernel centered at that value, scaled
    by the weight, to the histogram.
    """
    digits = -int(numpy.floor(numpy.log10(abs(binwidth))))
    start = round(val-kernelwidth, digits)
    end = round(val+kernelwidth, digits)
    num = start
    while (num<=end):
        newval = BiweightKernel((val-num)/kernelwidth)
        if (histogram.get(num) == None):
            histogram[num] = weight*newval/(numpoints*kernelwidth)
        else:
            histogram[num] += weight*newval/(numpoints*kernelwidth)

        num = round(num+binwidth, 3)


def getKernelWidth(data, binwidth):
    """
    Estimates a reasonable kernel width, based on the quartiles of the data
    """
    if (len(data) <= 1):
        return 2.78 * 20*binwidth
    highquart, lowquart = numpy.percentile(data, [75, 25])
    minsig = min(numpy.std(data), (highquart-lowquart)/1.34)
    minsig = max(minsig, 20*binwidth)
    return 2*minsig*pow(len(data), -0.2)

def combineHistograms(newhist, fullhist, num, weight=1):
    """
    Takes two histograms and merges them (assumes a total of 'num' histograms
    to merge, with a local weight of 'weight').
    """
    for val in newhist:
        if (fullhist.get(val) == None):
            fullhist[val] = weight*newhist[val]/num
        else:
            fullhist[val] += weight*newhist[val]/num

def createPrintAndSaveHistogram(data, filename, binwidth, xdata="log2r", axis=(), show=True, savefig=False):
    """
    Given a set of input data (data), treated as a set of data points, 
    and a given binwidth, create a smoothed histogram.  The binwidth shouldn't
    be too small, or this function will run forever.  Can be set to display
    the resulting histogram, or not, and to write out the data for the smoothed
    histogram, if desired.
    """
    if len(data) == 0:
        print("Cannot create a histogram with no data for file ", filename)
        return
    kw = getKernelWidth(data, binwidth)
    hist = {}
    if (show):
        print("Histogram for", filename, "with", len(data), "datapoints:")
    for val in data:
        addKernelToHistogram(val, 1, hist, kw, len(data), binwidth)
    if (show or savefig):
        plt.plot(list(hist.keys()), list(hist.values()), "ro")
        paxis = list(plt.axis())
        for a in range(0,len(axis)):
            paxis[a] = axis[a]
        plt.axis(paxis)
        if savefig:
            if "png" not in filename:
                plt.savefig(filename + ".png")
            else:
                plt.savefig(filename)
        if show:
            plt.show()
        plt.close()
    if filename != "" and "png" not in filename:
        outfile = open(filename + ".tsv", "w")
        outfile.write(xdata + "\tprobability\tnumpoints:\t" + str(len(data)) + "\tkernel width:\t" + str(kw) + "\n")
        for data in hist.keys():
            outfile.write(str(data) + "\t" + str(hist[data]) + "\n")
    return hist

def getPatientSampleMap():
    """
    Reads the (globally-set) patientSampleMapFile, and organizes its data for
    use by other functions.  Creates a dictionary of patients to samples,
    and another dictionary of samples to patient and ploidy.  (Ploidy here
    refers to the fact that pASCAT will return estimated CNVs based on an
    assumption of near-diploid copy number or near-tetraploid copy number,
    and the decision about which to use must be arrived at separately.)
    The file we typically use contains all the indirect data we use to 
    guess whether a given sample's diploid or tetraploid results are better,
    culminating in the next-to-last column the estimated odds of the diploid 
    solution being better than the tetraploid solution, and the actually-last
    column being the final by-hand call of what the best solution should be.
    """
    s2p = {}
    p2s = {}
        
    callfile = open(patientSampleMapFile, "r")
    for line in callfile:
        if "Patient" in line:
            continue
        lvec = line.rstrip().split()
        (patient, sample) = lvec[0:2]
        ploidy = lvec[-1]
        if ploidy=="Unknown":
            odds = float(lvec[-2])
            if odds > 0.5:
                ploidy = "Diploid"
            else:
                ploidy = "Tetraploid"
        if ploidy=="Diploid":
            ploidy = "diploid"
        if ploidy=="Tetraploid":
            ploidy="tetraploid"
        s2p[sample] = [patient, ploidy]
        if patient not in p2s:
            p2s[patient] = []
        p2s[patient].append(sample)
    return p2s, s2p

def loadDeletions(samplePatientMap, CNdir=CNdir):
    """
    Searches the CNVs in CNdir (default defined above), and returns a 
    dictionary by patient, sample, and chromosome of the deletions (both
    single and double deletions) for all patients and samples (as listed
    in the samplePatientMap).
    """
    deletions = {}
    for sample in samplePatientMap:
        patient, ploidy = samplePatientMap[sample]
        filename = CNdir + patient + "_" + sample + "_g500_" + ploidy + "_nonint_CNs.txt"
        if not isfile(filename):
            filename = CNdir + patient + "_" + sample + "_g550_" + ploidy + "_nonint_CNs.txt"
        for line in open(filename, "r"):
            if "patient" in line:
                continue
            lvec = line.rstrip().split()
            if lvec[7] == "0" or lvec[8] == "0":
                (chrom, start, end) = lvec[2:5]
                start = int(start)
                end = int(end)
                if patient not in deletions:
                    deletions[patient] = {}
                if sample not in deletions[patient]:
                    deletions[patient][sample] = {}
                if chrom not in deletions[patient][sample]:
                    deletions[patient][sample][chrom] = []
                deletions[patient][sample][chrom].append((start, end))
    return deletions


def loadDeletionsAndCNVs(samplePatientMap, CNdir=CNdir):
    """
    Searches the CNVs in CNdir (default defined above), and returns a 
    dictionary by patient, sample, and chromosome of the deletions (both
    single and double deletions) for all patients and samples (as listed
    in the samplePatientMap), and a separate dictionary of all calls for
    all segments in the sample.  This includes wild-type "1,1"-style
    calls, so if a site is not found in any of the segments in the returned
    dictionary, this means that the number of copies of that site is unknown.
    """
    deletions = {}
    CNVs = {}
    for sample in samplePatientMap:
        (patient, ploidy) = samplePatientMap[sample]
        filename = CNdir + patient + "_" + sample + "_g500_" + ploidy + "_nonint_CNs.txt"
        if not isfile(filename):
            filename = CNdir + patient + "_" + sample + "_g550_" + ploidy + "_nonint_CNs.txt"
        for line in open(filename, "r"):
            if "patient" in line:
                continue
            lvec = line.rstrip().split()
            if lvec[7] == "0" or lvec[8] == "0":
                (chrom, start, end) = lvec[2:5]
                start = int(start)
                end = int(end)
                if patient not in deletions:
                    deletions[patient] = {}
                if sample not in deletions[patient]:
                    deletions[patient][sample] = {}
                if chrom not in deletions[patient][sample]:
                    deletions[patient][sample][chrom] = []
                deletions[patient][sample][chrom].append((start, end))
            intA = lvec[7]
            intB = lvec[8]
            if intA=="NA" or intB=="NA":
                continue
            intA = int(intA)
            intB = int(intB)
            if intB<intA:
                temp = intA
                intA = intB
                intB = temp
            call = (intA, intB)
            (chrom, start, end) = lvec[2:5]
            start = int(start)
            end = int(end)
            if patient not in CNVs:
                CNVs[patient] = {}
            if sample not in CNVs[patient]:
                CNVs[patient][sample] = {}
            if chrom not in CNVs[patient][sample]:
                CNVs[patient][sample][chrom] = []
            CNVs[patient][sample][chrom].append((start, end, call))
    return deletions, CNVs






