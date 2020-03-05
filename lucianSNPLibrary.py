# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:34:39 2016

@author: lpsmith
"""

#Useful functions for manipulating log2r data

from __future__ import division
import numpy
import matplotlib.pyplot as plt
import random
from os import walk
from os.path import isfile

def BiweightKernel(t):
    if (abs(t) > 1.0):
        return 0.0
    return (15.0/16.0) * pow((1.0-pow(t,2)),2)

def addKernelToHistogram(val, weight, histogram, kernelwidth, numpoints, binwidth):
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
    if (len(data) <= 1):
        return 2.78 * 20*binwidth
    highquart, lowquart = numpy.percentile(data, [75, 25])
    minsig = min(numpy.std(data), (highquart-lowquart)/1.34)
    minsig = max(minsig, 20*binwidth)
    return 2*minsig*pow(len(data), -0.2)

def combineHistograms(newhist, fullhist, num, weight=1):
    for val in newhist:
        if (fullhist.get(val) == None):
            fullhist[val] = weight*newhist[val]/num
        else:
            fullhist[val] += weight*newhist[val]/num

def createPrintAndSaveHistogram(data, filename, binwidth, xdata="log2r", axis=(), show=True, savefig=False):
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

def createPrintAndSaveMultipleHistograms(data, filename, binwidth, xdata="log2r", ydata = [], axis=(), show=True):
    if len(data) == 0:
        print("Cannot create a histogram with no data for file ", filename)
        return
    hists = []
    for onehist in data:
        if len(onehist) == 0:
            print("Cannot create a histogram with no data for a subset of the data in", filename)
            continue
        kw = getKernelWidth(onehist, binwidth)
        hist = {}
        if (show):
            print("One histogram in", filename, "with", len(onehist), "datapoints:")
        for val in onehist:
            addKernelToHistogram(val, 1, hist, kw, len(onehist), binwidth)
        if (show):
            plt.plot(hist.keys(), hist.values(), "ro")
            paxis = list(plt.axis())
            for a in range(0,len(axis)):
                paxis[a] = axis[a]
            plt.axis(paxis)
            plt.show()
            plt.close()
        hists.append(hist)
    outfile = open(filename, "w")
    outfile.write(xdata)
    for n in range(len(data)):
        if (n < len(ydata)):
            outfile.write("\t" + ydata[n])
        else:
            outfile.write("\tset" + str(n))
        outfile.write(" (" + str(len(data[n])) + " pts)")
    outfile.write("\n")
    allkeys = set()
    for hist in hists:
        allkeys = allkeys | set(hist.keys())

    for key in allkeys:
        outfile.write(str(key))
        for hist in hists:
            if key in hist:
                outfile.write("\t" + str(hist[key]))
            else:
                outfile.write("\t")
        outfile.write("\n")

def getPeaks(hist, binwidth):
    low = min(hist.keys())
    high = max(hist.keys())
    npoints = (high-low)/binwidth
    window = min(20, npoints)

def saveScatterPlot(data, filename, labels):
    if len(data) == 0:
        print("Cannot create a scatterplot with no data for file ", filename)
        return
    outfile = open(filename, "w")
    outfile.write(labels + "\n")
    for dataline in data:
        for n in range(len(dataline)):
            if (n>0):
                outfile.write("\t")
            outfile.write(str(dataline[n]))
        outfile.write("\n")
    outfile.close()

def getChromosomeMax(chr):
    if (chr==1):
        return 249300000
    elif (chr==2):
        return 243000000
    elif (chr==3):
        return 197950000
    elif (chr==4):
        return 191100000
    elif (chr==5):
        return 180900000
    elif (chr==6):
        return 170950000
    elif (chr==7):
        return 159200000
    elif (chr==8):
        return 146300000
    elif (chr==9):
        return 141050000
    elif (chr==10):
        return 135400000
    elif (chr==11):
        return 135100000
    elif (chr==12):
        return 133780000
    elif (chr==13):
        return 115150000
    elif (chr==14):
        return 107300000
    elif (chr==15):
        return 102500000
    elif (chr==16):
        return 90100000
    elif (chr==17):
        return 81100000
    elif (chr==18):
        return 78100000
    elif (chr==19):
        return 59100000
    elif (chr==20):
        return 62980000
    elif (chr==21):
        return 48100000
    elif (chr==22):
        return 51670000
    elif (chr==23):
        return 155300000
    elif (chr==24):
        return 28800000
    else:
        print("Unknown chromosome #", chr)

def getChromosomeCentromere(chr):
    if (chr==1):
        return (135380000, 142530000)
    elif (chr==2):
        return (92300000, 95330000)
    elif (chr==3):
        return (90500000, 93500000)
    elif (chr==4):
        return (49640000, 52330000)
    elif (chr==5):
        return (46400000, 49420000)
    elif (chr==6):
        return (58620000, 61900000)
    elif (chr==7):
        return (57850000, 61060000)
    elif (chr==8):
        return (43800000, 46800000)
    elif (chr==9):
        return (47200000, 65500000)
    elif (chr==10):
        return (39020000, 41650000)
    elif (chr==11):
        return (51580000, 54400000)
    elif (chr==12):
        return (34720000, 37850000)
    elif (chr==13):
        return (0, 19020000)
    elif (chr==14):
        return (0, 18200000)
    elif (chr==15):
        return (0, 20010000)
    elif (chr==16):
        return (35250000, 46400000)
    elif (chr==17):
        return (22300000, 25300000)
    elif (chr==18):
        return (15410000, 18500000)
    elif (chr==19):
        return (24500000, 27700000)
    elif (chr==20):
        return (26400000, 29400000)
    elif (chr==21):
        return (11300000, 14300000)
    elif (chr==22):
        return (0, 14670000)
    elif (chr==23):
        return (58420000, 61730000)
    elif (chr==24):
        return (10400000, 13200000)
    else:
        print("Unknown chromosome #", chr)

def getLengthFrom(chr, min, max):
    if (max=="inf"):
        max = getChromosomeMax(chr)
    else:
        max = int(max)
    min = int(min)
    return max-min

def getGlobalPosition(chr, min, max):
    prevlength = 0;
    for c in range(0,chr-1):
        prevlength += getChromosomeMax(c+1)
    if (max=="inf"):
        max = getChromosomeMax(chr)
    else:
        max = int(max)
    min = int(min)
    return prevlength + (max+min)/2

def WhichArm(chr, start, end):
    (cstart, cend) = getChromosomeCentromere(chr)
    if (end <= cstart):
        return -1
    if (start >= cend):
        return 1
    return 0

def getSNPLabels1M(zeroes):
    # read the probeset file, which correlates name to position.
    infilename = "HumanOmni1-Quad_v1-0_H_SNPlist.txt"
    infile = open(infilename,"r")
    #duplicates = open("duplicates.txt", "w")

    labels = {}
    rev_labels = {}
    for line in infile:
        if line.find("Chr") != -1:
            continue
        (id, __, chr, pos) = line.rstrip().split()
        if chr=="X":
            chr = "23"
        if chr=="XY":
            chr = "23"
        if chr=="Y":
            chr = "24"
        try:
            int(chr)
        except ValueError:
            #print("problematic chr: " + chr)
            continue
        if (not zeroes):
            if (pos == "0"):
                continue
            if (chr == "0"):
                continue
        if (chr, pos) not in rev_labels:
            rev_labels[(chr, pos)] = id
        elif "cnvi" in rev_labels[(chr, pos)]:
            rev_labels[(chr, pos)] = id
        #else:
            #print("Two SNPs with the same position:", id, "and", rev_labels[chr, pos], "both map to", chr, ",", pos)
            #duplicates.write(id + "\t" + rev_labels[chr, pos] + "\t" + chr + "\t" + pos + "\n")
            #continue
            #UPDATE:  it appears that Partek uses all SNPs, even if they map to the same location.
        labels[id] = (chr, pos)
    infile.close()
    return labels, rev_labels

def getSNPLabels2_5M(zeroes):
    #infilename = "CN_raw_data_25M/REI_12051_B01_SOM_WGS_443samples_12Dec2016_Partek_Partek.annotation.txt"
    infilename = "probe_sets/InfiniumOmni2-5-8v1-3_A1.annotated.txt"
    infile = open(infilename,"r")
    #duplicates = open("duplicates.txt", "w")

    labels = {}
    rev_labels = {}
    for line in infile:
        if line.find("Name") != -1:
            continue
        line = line.rstrip().split()

        if (len(line) >= 4):
            #(id, __, chr, pos) = line[0:4]
            (id, chr, pos) = line[0:3]
            if chr=="X":
                chr = "23"
            if chr=="XY":
                chr = "23"
            if chr=="Y":
                chr = "24"
            try:
                int(chr)
            except ValueError:
                #print("problematic chr: " + chr)
                continue
            if (not zeroes):
                if (pos == "0"):
                    continue
                if (chr == "0"):
                    continue
            if (chr, pos) not in rev_labels:
                rev_labels[(chr, pos)] = id
            elif "cnvi" in rev_labels[(chr, pos)]:
                rev_labels[(chr, pos)] = id
            #else:
                #print("Two SNPs with the same position:", id, "and", rev_labels[chr, pos], "both map to", chr, ",", pos)
                #duplicates.write(id + "\t" + rev_labels[chr, pos] + "\t" + chr + "\t" + pos + "\n")
                #continue
                #UPDATE:  it appears that Partek uses all SNPs, even if they map to the same location.
            labels[id] = (chr, pos)
    infile.close()
    return labels, rev_labels

def getSNPLabels2_5Mv12(zeroes):
    infilename = "HumanOmni2-5-8-v1-2-A-Gene-Annotation-File.txt"
    infile = open(infilename,"r")
    #duplicates = open("duplicates.txt", "w")

    labels = {}
    rev_labels = {}
    for line in infile:
        if line.find("Name") != -1:
            continue
        line = line.rstrip().split()

        (id, chr, pos) = line[0:3]
        if chr=="X":
            chr = "23"
        if chr=="XY":
            chr = "23"
        if chr=="Y":
            chr = "24"
        try:
            int(chr)
        except ValueError:
            #print("problematic chr: " + chr)
            continue
        if (not zeroes):
            if (pos == "0"):
                continue
            if (chr == "0"):
                continue
        if (chr, pos) not in rev_labels:
            rev_labels[(chr, pos)] = id
        elif "cnvi" in rev_labels[(chr, pos)]:
            rev_labels[(chr, pos)] = id
        #else:
            #print("Two SNPs with the same position:", id, "and", rev_labels[chr, pos], "both map to", chr, ",", pos)
            #duplicates.write(id + "\t" + rev_labels[chr, pos] + "\t" + chr + "\t" + pos + "\n")
            #continue
            #UPDATE:  it appears that Partek uses all SNPs, even if they map to the same location.
        labels[id] = (chr, pos)
    infile.close()
    return labels, rev_labels

def getSNPLabelsAll(zeroes):
    # read the probeset file, which correlates name to position.
    infilename = "probe_set_1_25_all.txt"
    infile = open(infilename,"r")
    #duplicates = open("duplicates.txt", "w")

    labels = {}
    rev_labels = {}
    discard_duplicates = ["kgp3675773", "rs1443", "kgp11488181", "kgp1946361", "rs33945777", "rs3752426", "kgp468181", "rs1518004", "kgp6538154", "kgp4057395", "rs9258315", "kgp11800323", "kgp4698400", "rs2523832", "kgp6575555", "kgp5212761", "401753", "rs45514397", "kgp5430667", "rs12315", "kgp3529433", "kgp4408472", "rs2070232", "rs1786171", "VG15S12103", "rs9420882", "rs7148", "rs45589133", "rs9332385", "rs45530145", "rs467225", "rs10010325", "rs13306723", "kgp3668772", "kgp17138187"]
    for line in infile:
        line = line.rstrip().split()
        if (len(line) == 3):
            (id, chr, pos) = line[0:3]
            if chr=="X":
                chr = "23"
            if chr=="Y":
                chr = "24"
            try:
                int(chr)
            except ValueError:
                #print("problematic chr: " + chr)
                continue
            if id in discard_duplicates:
                chr = "0"
                pos = "0"
            if (not zeroes):
                if (pos == "0"):
                    continue
                if (chr == "0"):
                    continue
            if (chr, pos) not in rev_labels:
                rev_labels[(chr, pos)] = id
            elif "cnvi" in rev_labels[(chr, pos)]:
                rev_labels[(chr, pos)] = id
            #else:
                #print("Two SNPs with the same position:", id, "and", rev_labels[chr, pos], "both map to", chr, ",", pos)
                #duplicates.write(id + "\t" + rev_labels[chr, pos] + "\t" + chr + "\t" + pos + "\n")
                #continue
                #UPDATE:  it appears that Partek uses all SNPs, even if they map to the same location.
            labels[id] = (chr, pos)
    infile.close()
    return labels, rev_labels

def getSNPLabelsAveraged(zeroes):
    # read the probeset file, which correlates name to position.
    infilename = "probe_set_1_25_averaged_overlap.txt"
    infile = open(infilename,"r")
    #duplicates = open("duplicates.txt", "w")

    labels = {}
    rev_labels = {}
    for line in infile:
        line = line.rstrip().split()
        if (len(line) == 3):
            (id, chr, pos) = line[0:3]
            if chr=="X":
                chr = "23"
            if chr=="Y":
                chr = "24"
            try:
                int(chr)
            except ValueError:
                #print("problematic chr: " + chr)
                continue
            if (not zeroes):
                if (pos == "0"):
                    continue
                if (chr == "0"):
                    continue
            if (chr, pos) not in rev_labels:
                rev_labels[(chr, pos)] = id
            elif "cnvi" in rev_labels[(chr, pos)]:
                rev_labels[(chr, pos)] = id
            #else:
                #print("Two SNPs with the same position:", id, "and", rev_labels[chr, pos], "both map to", chr, ",", pos)
                #duplicates.write(id + "\t" + rev_labels[chr, pos] + "\t" + chr + "\t" + pos + "\n")
                #continue
                #UPDATE:  it appears that Partek uses all SNPs, even if they map to the same location.
            labels[id] = (chr, pos)
    infile.close()
    return labels, rev_labels

def getSNPLabels(is1M, zeroes):
    if (is1M):
        return getSNPLabels1M(zeroes)
    return getSNPLabels2_5M(zeroes)

def addPopsToList(samepops, point, sample, pops):
    if not(sample in samepops):
        samepops[sample] = {}
    pops_t = tuple(pops)
    if not (pops_t in samepops[sample]):
        samepops[sample][pops_t] = []
    samepops[sample][pops_t].append(point)


def getPoints(data):
    allsegments = set()
    for segments in data:
        for segment in segments[1]:
            start = (segment[0], segment[1])
            end   = (segment[0], segment[2])
            allsegments.add(start)
            allsegments.add(end)
    return sorted(list(allsegments))

def getSamePops(data, points):
    samepops = {}
    for segments in data:
        sample = segments[0]
        seg = segments[1]
        pnum = 0
        snum = 0
        while (pnum < len(points) and snum < len(seg)):
            point = points[pnum]
            chr   = seg[snum][0]
            start = seg[snum][1]
            end   = seg[snum][2]
            pops  = seg[snum][3]
            if chr < point[0]:
                snum = snum+1
            elif chr > point[0]:
                pnum = pnum+1
            elif point[1] == 0 and start == 0:
                #include
                addPopsToList(samepops, point, sample, pops)
                pnum = pnum+1
            elif point[1] < start:
                pnum = pnum+1
            elif point[1] >= end:
                snum = snum+1
            else:
                addPopsToList(samepops, point, sample, pops)
                pnum = pnum+1
    return samepops

def calculateOverlap(data, points):
    samepops = getSamePops(data,points)
    pairs = {}
    for sample in samepops:
        for pops in samepops[sample]:
            points = samepops[sample][pops]
            for p in range(0, len(points)):
                for q in range(p+1, len(points)):
                    pair = (points[p][0], points[p][1], points[q][0], points[q][1])
                    if not(pair in pairs):
                        pairs[pair] = []
                    pairs[pair].append(sample)
    return pairs


def randomizePopAssignment(data):
    for segpairs in data:
        for n in range(0,len(segpairs[1])):
            s = random.randint(n,len(segpairs[1])-1)
            if segpairs[1][n][3] == segpairs[1][s][3]:
                continue
            temppop = segpairs[1][n][3]
            segpairs[1][n][3] = segpairs[1][s][3]
            segpairs[1][s][3] = temppop
    return


def tenvec(nrandom):
    return [nrandom, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

def collectRandomPairs(randompairs, pairs, nrandom):
    for pair in pairs:
        newval = len(pairs[pair])
        if not(pair in randompairs):
            randompairs[pair] = tenvec(nrandom)
        randompairs[pair][0] += -1
        randompairs[pair][newval] += 1

def combinePopsInList(combinedpops, point, sample, pops):
    if not(point in combinedpops):
        combinedpops[point] = {}
    combinedpops[point][sample] = pops

def allZeroes(pops):
    for pop in pops:
        if pop != 0:
            return False
    return True

def sufficientOverlap(seg1, seg2, mindiff):
    if (abs(seg1[1] - seg2[1]) < mindiff and \
        abs(seg1[2] - seg2[2]) < mindiff and \
        seg1[1] < seg2[2] and seg1[2] > seg2[1]):
        return True
    if allZeroes(seg2[3]) and seg1[1] >= seg2[1] and seg1[2] <= seg2[2]:
        return True
    return False

def getBestOverlap(armlist, numsamples, mindiff):
    allMutantSegments = []
    for segpair in armlist:
        (sample, segment) = segpair
        (chr, start, end, pops) = segment
        if not(allZeroes(pops)):
            allMutantSegments.append([segment, []])
    #print("All samples: ", allsamples, ", all mutant segments", allMutantSegments)
    for proposedSegment in allMutantSegments:
        #print("proposed segment", proposedSegment)
        for segpair in armlist:
            (sample, segment) = segpair
            (chr, start, end, pops) = segment
            if sufficientOverlap(proposedSegment[0], segment, mindiff):
                valid = True
                for oldsegcheck in proposedSegment[1]:
                    if sample == oldsegcheck[0]:
                        print("Two segments from sample", sample, "both seem to overlap.")
                        print("First = ", oldsegcheck)
                        print("New   = ", segpair)
                        print("Match = ", proposedSegment[0])
                        valid = False
                if (valid):
                    proposedSegment[1].append(segpair)
                #print("added-to proposed segment", proposedSegment)
    for proposedSegment in allMutantSegments:
        if len(proposedSegment[1]) == numsamples:
            return proposedSegment
        if len(proposedSegment[1]) > numsamples:
            print("What???", proposedSegment, numsamples)
    return []

def AddSegmentToPerArm(sample, segment, perarm):
    (chr, start, end, pops) = segment
    end = WhichArm(chr, start, end)
    lowpair= (chr, "low")
    highpair = (chr, "high")
    if end == 0:
        if allZeroes(pops):
            if not(lowpair in perarm):
                perarm[lowpair] = []
            perarm[lowpair].append((sample, segment))
            if not (highpair in perarm):
                perarm[highpair] = []
            perarm[highpair].append((sample, segment))
        #else:
            #print(pops)
    elif end == 1:
        if not (highpair in perarm):
            perarm[highpair] = []
        perarm[highpair].append((sample, segment))
    elif end == -1:
        if not(lowpair in perarm):
            perarm[lowpair] = []
        perarm[lowpair].append((sample, segment))
    return

def OneSegmentPerArm(segmentlist, numsamples, mindiff):
    perarm = {}
    for segpair in segmentlist:
        (sample, segments) = segpair
        for segment in segments:
            AddSegmentToPerArm(sample, segment, perarm)
    bestsegs = []
    for armlist in perarm:
        #print(armlist)
        bestOverlap = getBestOverlap(perarm[armlist], numsamples, mindiff)
        if (len(bestOverlap) > 0):
            bestsegs.append(bestOverlap)
    #print(len(bestsegs), "best segments: ", bestsegs)
    return bestsegs

def CalculateAverageDiseq(allpops, probones, numsegs, numpops):
    diseqs = []
    for n in range(0, numsegs-1):
        for s in range(n+1, numsegs):
            numoneones = 0
            for i in range(0,numpops):
                if (allpops[n][i] == 1 and allpops[s][i] == 1):
                    numoneones += 1
            diseqs.append(abs((numoneones/numpops) - (probones[n]*probones[s])))
    return numpy.average(diseqs)

def SwapFour(shuffled, goodpops, goodcols):
    if len(goodpops)==0 or len(goodcols)==0:
#        numpy.nothing()
        return
    unfound = True
    while (unfound):
        p1 = random.choice(goodpops)
        s1 = random.choice(goodcols)
        c1 = shuffled[p1][s1]
        pset = []
        pairs = []
        for x in goodpops:
            if x==p1:
                continue
            if shuffled[x][s1] != c1:
                pset.append(x)
        for x in goodcols:
            if x==s1:
                continue
            if shuffled[p1][x] != c1:
                for y in pset:
                    pairs.append((y, x))
        while (unfound and len(pairs)>0):
            pairnum = random.randint(0,len(pairs)-1)
            (p2, s2) = pairs[pairnum]
            if shuffled[p2][s2] != c1:
                del pairs[pairnum]
                continue
            unfound = False
            c2 = shuffled[p1][s2]
            shuffled[p1][s1] = c2
            shuffled[p2][s2] = c2
            shuffled[p1][s2] = c1
            shuffled[p2][s1] = c1
        unfound = False

def GetVariedRowsAndColumns(shuffled):
    goodpops = []
    for pop in range(0,len(shuffled)):
        psum = numpy.sum(shuffled[pop])
        if psum == 0 or psum == len(shuffled[0]):
#            print("original pass: ", pop, "no good in ")
#            print(shuffled)
            continue
        goodpops.append(pop)
    goodcols = []
    for col in range(0,len(shuffled[0])):
        csum = numpy.sum(shuffled[:,col])
        if csum == 0 or csum == len(shuffled):
#            print("original pass: ", col, "no good in ")
#            print(shuffled)
            continue
        goodcols.append(col)
    (goodpops, goodcols) = GetInternallyVariedRowsAndColumns(shuffled, goodpops, goodcols)
    return (goodpops, goodcols)

def GetInternallyVariedRowsAndColumns(shuffled, goodpops, goodcols):
    badpops = []
    for pop in goodpops:
        psum = 0
        for col in goodcols:
            psum += shuffled[pop][col]
        if psum == 0 or psum == len(goodpops):
            badpops.append(pop)
            continue
    badcols = []
    for col in goodcols:
        csum = 0
        for pop in goodpops:
            csum += shuffled[pop][col]
        if csum == 0 or csum == len(goodcols):
            badcols.append(col)
    if (len(badpops) > 0 or len(badcols) > 0):
#        print("Found some extra bad rows or columns", badpops, badcols)
#        numpy.nothing()
        for bp in badpops:
            goodpops.remove(bp)
        for bc in badcols:
            goodcols.remove(bc)
        goodpops, goodcols = GetInternallyVariedRowsAndColumns(shuffled, goodpops, goodcols)
    return (goodpops, goodcols)

def ShufflePops(allpops, goodpops, goodcols):
    shuffled = allpops
    for n in range(0,100):
        SwapFour(shuffled, goodpops, goodcols)
    return shuffled

def FindFourGametes(allpops):
    foundfour = 0
    notfour = 0
    for n in range(0, len(allpops)-1):
        gametes1 = allpops[n]
        for s in range(n+1, len(allpops)):
            gametes2 = allpops[s]
            if (len(gametes1) != len(gametes2)):
                print("Error!")
                return (0,0)
            g4test = set()
            g4test.add("zerozero") #Because we know that the ancestral allele was WT.
            for g in range(0,len(gametes1)):
                g1 = gametes1[g]
                g2 = gametes2[g]
                if (g1==1):
                    if (g2==1):
                        g4test.add("oneone")
                    else:
                        g4test.add("onezero")
                else:
                    if (g2==1):
                        g4test.add("zeroone")
                    else:
                        g4test.add("zerozero")
            if len(g4test) > 3:
                foundfour += 1
            else:
                notfour += 1
    return (foundfour, notfour)

def GetJustPops(segments):
    justpops = []
    for segment in segments:
        justpops.append(segment[3])
    return numpy.array(justpops)

def getNumSegmentsPerSample():
    directory = "expands_full_input/"

    flist = []
    for (_, _, f) in walk(directory):
        flist += f

    seglengths = {}
    segvec = []
    for f in flist:
        if f.find("BAF") == -1:
            continue
        (patient, sample, __) = f.split("_")
        numsegs = 0;
        baffile = open(directory + f, "r")
        for line in baffile:
            numsegs += 1
        seglengths[(patient, sample)] = numsegs
        segvec.append(numsegs)

    avgseg = numpy.average(segvec)
    return (seglengths, avgseg)

def fixNameForR(name):
    ret = name.replace("#", "").replace("_", "")
    return ret

def resegmentLOHes(segvec):
    retvec = segvec[:]
    retvec.sort()
    for n in range(0, len(retvec)-1):
        line1 = retvec[n]
        line2 = retvec[n+1]
        if line1[0] == line2[0]:
            if line1[2] > line2[1]:
                l1start = line1[1]
                l1end = line1[2]
                l2start = line2[1]
                l2end = line2[2]
                if l1start == l2start:
                    if l1end == l2end:
                        #Just keep the LOH segment; remove the other one entirely.
                        if line1[4] == "LOH":
                            del retvec[n+1]
                            return resegmentLOHes(retvec)
                        del retvec[n]
                        return resegmentLOHes(retvec)
                    if l1end > l2end:
                        if line1[4] == "LOH":
                            #The LOH segment is bigger: remove the other entirely.
                            del retvec[n+1]
                            return resegmentLOHes(retvec)
                        #The LOH segment is smaller: remove it from the other and continue.
                        line1[1] = line2[2]
                        return resegmentLOHes(retvec)
                    if l1end < l2end:
                        if line2[4] == "LOH":
                            #The LOH segment is bigger: remove the other entirely.
                            del retvec[n]
                            return resegmentLOHes(retvec)
                        #The LOH segment is smaller: remove it from the other and continue.
                        line2[1] = line1[2]
                        return resegmentLOHes(retvec)
                else: #l1start < l2start
                    if l1end == l2end:
                        if line1[4] == "LOH":
                            #The LOH segment is bigger: remove the other entirely.
                            del retvec[n+1]
                            return resegmentLOHes(retvec)
                        #The LOH segment is smaller: remove it from the other and continue.
                        line1[2] = line2[1]
                        return resegmentLOHes(retvec)
                    if l1end > l2end:
                        if line1[4] == "LOH":
                            #The LOH segment is bigger: remove the other entirely.
                            del retvec[n+1]
                            return resegmentLOHes(retvec)
                        #The LOH segment is smaller and embedded inside the larger segment: make the exisitng segment the left side, and add a brand new segment for the right.
                        newline = line1[:]
                        newline[1] = line2[2]
                        retvec.append(newline)
                        line1[2] = line2[1]
                        return resegmentLOHes(retvec)
                    if l1end < l2end:
                        #They overlap: find which one is LOH and shorten the other one.
                        if line1[4] == "LOH":
                            line2[1] = l1end
                            return resegmentLOHes(retvec)
                        line1[2] = l2start
                        return resegmentLOHes(retvec)
    return retvec


def validateSegments(BAFs_by_sample, validation_output, id, failfile, summaryfile):
    #assumes that any double deletion is already removed.
    output = {}
    samples = list(BAFs_by_sample.keys())
    sampairs = set()
    summatches = 0
    sumantimatches = 0
    sumfails = 0
    samplefails = {}
    for sample in samples:
        samplefails[sample] = [0, 0, 0]
    for n in range(0,len(samples)-1):
        sample1 = samples[n]
        for m in range(n+1, len(samples)):
            sample2 = samples[m]
            matches = []
            for segname in BAFs_by_sample[sample1]:
                if not(segname in BAFs_by_sample[sample2]):
                    continue
                match = 0
                antimatch = 0
                for pos in BAFs_by_sample[sample1][segname]:
                    if not(pos in BAFs_by_sample[sample2][segname]):
                        continue
                    baf1 = BAFs_by_sample[sample1][segname][pos]
                    baf2 = BAFs_by_sample[sample2][segname][pos]
#                    if 0.4 < baf1 and baf1 < 0.6:
#                        continue
#                    if 0.4 < baf2 and baf2 < 0.6:
#                        continue
                    if baf1 < 0.5 and baf2 < 0.5:
                        match += 1
                    elif baf1 > 0.5 and baf2 > 0.5:
                        match += 1
                    elif baf1 < 0.5 and baf2 > 0.5:
                        antimatch += 1
                    elif baf1 > 0.5 and baf2 < 0.5:
                        antimatch += 1
                    #print(baf1, baf2, match, antimatch)
                matches.append((match, antimatch))
                nmax = match + antimatch
                if nmax >= 2:
                    if match/nmax > 0.9:
                        samplefails[sample1][0] += 1
                        samplefails[sample2][0] += 1
                    elif antimatch/nmax > 0.9:
                        samplefails[sample1][1] += 1
                        samplefails[sample2][1] += 1
                    else:
                        samplefails[sample1][2] += 1
                        samplefails[sample2][2] += 1
                if not(segname in output):
                    output[segname] = {}
                output[segname][sample1+"_"+sample2] = (match, antimatch)
                sampairs.add(sample1 + "_" + sample2)
    for sample in samples:
        summaryfile.write(id + "\t" + sample + "\t" + str(samplefails[sample][0]) + "\t" + str(samplefails[sample][1]) + "\t" +  str(samplefails[sample][2])  + "\n")
    outfile = open(validation_output + id + "_segvalidate.txt", "w")
    outfile.write("patient\tchr\tstart\tend\tmax_nBAFs\tmatches\tantimatches\tfails")
    sampairs = sorted(sampairs)
    for sampair in sampairs:
        outfile.write("\t" + sampair + " match\t" + sampair + " anti-match")
    outfile.write("\n")
    for segname in output:
        outfile.write(id + "\t" + segname[0] +"\t" + str(segname[1]) +"\t" + str(segname[2]))
        outline = ""
        nummatches = 0
        numantimatches = 0
        numfails = 0
        maxBAFs = 0
        for sampair in sampairs:
            if sampair in output[segname]:
                match = output[segname][sampair][0]
                antimatch = output[segname][sampair][1]
                outline += "\t" + str(match) + "\t" + str(antimatch)
                numBAFs = match+antimatch
                if maxBAFs<numBAFs:
                    maxBAFs = numBAFsod
                if numBAFs<=1:
                    continue
                elif (match/numBAFs) > .9:
                    #print("Match", match, antimatch, numBAFs, match/numBAFs)
                    nummatches += 1
                elif (antimatch/numBAFs) > .9:
                    numantimatches += 1
                    #print("Antimatch", match, antimatch, numBAFs, match/numBAFs)
                else:
                    numfails += 1
                    (sample1, sample2) = sampair.split("_")
                    failfile.write(id + "\t" + segname[0] + "\t" + str(segname[1]) +"\t" + str(segname[2]) + "\t" + sample1)
                    for pos in BAFs_by_sample[sample1][segname]:
                        failfile.write("\t" + str(BAFs_by_sample[sample1][segname][pos]))
                    failfile.write("\n")
                    failfile.write(id + "\t" + segname[0] + "\t" + str(segname[1]) +"\t" + str(segname[2]) + "\t" + sample2)
                    for pos in BAFs_by_sample[sample2][segname]:
                        failfile.write("\t" + str(BAFs_by_sample[sample2][segname][pos]))
                    failfile.write("\n")
                    #print("Fail", match, antimatch, numBAFs, match/numBAFs)
            else:
                outline += "\t--\t--"
        outfile.write("\t" + str(maxBAFs) + "\t" + str(nummatches) + "\t" + str(numantimatches) + "\t" + str(numfails) + outline + "\n")
        summatches += nummatches
        sumantimatches += numantimatches
        sumfails += numfails
    outfile.close()
    summaryfile.write(id + "\tall\t" + str(summatches) + "\t" + str(sumantimatches) + "\t" + str(sumfails) + "\n")



#From Mary:
def readblock(data, keywd):
  # find chromosomes
  collection = []
  found = False
  done = False
  for line in data:
    if done: break     # we're past desired section

    line = line.rstrip()
    line = line.rstrip(",")
    line = line.replace(" ","")
    line = line.replace('"','')

    # the line consists of quoted comma-delimited entries with spaces
    # between them; a whole section is enclosed in parentheses.

    if keywd in line:  # start of section
      line = line.split("(")
      line = line[-1]
      found = True
      line = line.rstrip().split(",")
      for entry in line:
        collection.append(entry)
      found = True
      continue

    if found and ")" not in line:   # middle of section
      line = line.split(",")
      for entry in line:
        collection.append(entry)
      continue

    if found:
      line = line.split(")")
      line = line[0]
      if line != "":     # test for the lonely parenthesis
        line = line.split(',')
        for entry in line:
          collection.append(entry)
      done = True

  return collection

def collatepASCATOutput(infiledir, infilename, outfile,markerlocations):

  # we assume that "outfile" exists and is open for writing,
  # and that "markerlocations" has been created by a previous
  # call to "makeprobedict()"

  # read ASCAT segmentation data from dump of the ascat.output object
  # THIS FILE WON'T EXIST if P-ASCAT failed due to missing data, so we
  # return an error code to skip it
  (pid, sid) = infilename.split("_")[0:2]
  try:
    segdata = open(infiledir + infilename,"r").readlines()
  except:
    print("Unable to open file",infiledir + infilename)
    assert False
    return False

  startsnps = readblock(segdata,"track_start")
  endsnps = readblock(segdata,"track_end")
  Araw = readblock(segdata,"track_nAraw")
  Braw = readblock(segdata,"track_nBraw")
  Aint = readblock(segdata,"track_nAint")
  Bint = readblock(segdata,"track_nBint")
#  print(len(startsnps))
#  print(len(endsnps))

  # find the positions corresponding to start and end probes
  chr = []
  startpos = []
  endpos = []
  prevchr = 0
  prevend = 0
  for num in range(len(startsnps)):
    st = startsnps[num]
    en = endsnps[num]
#  for st, en in zip(startsnps, endsnps):
    startloc = markerlocations[st]
    ch = startloc[0]
    stp = startloc[1]
    endloc = markerlocations[en]
    assert ch == endloc[0]
    enp = endloc[1]
    chr.append(ch)
    startpos.append(stp)
    endpos.append(enp)
    if (ch == prevchr and float(prevend)>float(stp)):
        print("ERROR: overlapping segments, possibly due to pASCAT mislabeling segments.")
        assert False
        return False
    prevchr = ch
    prevend = enp

  for (ch, st, en, ar, br, ai, bi) in zip(chr,startpos,endpos,Araw,Braw,Aint,Bint):
    outline = pid + "\t" + sid + "\t"               # patient and sample
    outline += ch + "\t" + str(st) + "\t" + str(en) + "\t"    # segment position
    if ar=="NA" or br=="NA":
        outline += ar + "\t" + br + "\t"                # non-integer CN
    else:
        outline += str(float(ar)) + "\t" + str(float(br)) + "\t"                # non-integer CN
    outline += ai + "\t" + bi + "\n"                # integer CN
    outfile.write(outline)

  return True


def getCanonicalAscatCallsFor(patient):
    bestfilename = "best_analyses/" + patient + "_best.tsv"
    if not isfile(bestfilename):
        print("No such file", bestfilename)
        return []
    else:
        print("Reading best analysis file", bestfilename)
        best_file = open(bestfilename, "r")
        ret = set()
        for line in best_file:
            if "Patient" in line:
                continue
            lvec = line.split()
            (patient, sample, constraint, compare, accuracy, gamma) = lvec[0:6]
            if compare=="by_length" and constraint != "Xiaohong" and constraint != "overall":
                if constraint == "eight":
                    #Check to make sure that the tetraploid solution wasn't identical:
                    eightploidyfile = "gamma_test_output/pASCAT_input_g" + gamma + "/" + constraint + "/" + patient + "_fcn_ascat_ploidy.txt"
                    if not isfile(eightploidyfile):
                        continue
                    eightploidy = "0"
                    tetraploidy = "0"
                    for line in open(eightploidyfile):
                        if patient + "_" + sample in line:
                            eightploidy = line.split()[1]
                    tetraploidyfile = "gamma_test_output/pASCAT_input_g" + gamma + "/tetraploid/" + patient + "_fcn_ascat_ploidy.txt"
                    if isfile(tetraploidyfile):
                        for line in open(tetraploidyfile):
                            if patient + "_" + sample in line:
                                tetraploidy = line.split()[1]
                    if eightploidy == "0" or eightploidy == tetraploidy:
                        #found the same ploidy or no solution
                        continue
                ret.add((sample, gamma, constraint))
        return ret
        
    print("Unknown patient")
    return []

def getBestPloidyFor(patient, sample, challenge=False):
    oddsfile = "calling_evidence_odds.tsv"
    if challenge:
        oddsfile = "calling_evidence_challenge_inc_odds.tsv"
    for line in open(oddsfile, "r"):
        lvec = line.rstrip().split("\t")
        if lvec[0] != patient:
            continue
        if lvec[1] != sample:
            continue
        #first check if the hutch made a decision about this sample:
        if lvec[-1] == "Diploid":
            return "diploid"
        if lvec[-1] == "Tetraploid":
            return "tetraploid"
        #if there was only one solution, we go with that solution.
        if "Diploid only" in line:
            return "diploid"
        if "Tetraploid only" in line:
            return "tetraploid"
        #If not, we use the odds calculation:
        odds = float(lvec[-2])
        if odds >= 0.5:
            return "diploid"
        return "tetraploid"
    if challenge:
        print("Unknown patient/sample", patient, sample)
        assert(False)
        return "diploid"
    else:
        return "unknown"

def getSingleGammaCallsFor(patient, gamma):
    if gamma == "500":
        gamma = getGammaFor(patient)
    gammadir = "gamma_test_output/pASCAT_input_g" + gamma + "/"
    ret = set()
    for ploidy in ["diploid", "tetraploid", "eight"]:
        ploidyfile = gammadir + ploidy + "/" + patient + "_fcn_ascat_ploidy.txt"
        if isfile(ploidyfile):
            for line in open(ploidyfile):
                if patient not in line:
                    continue
                (pat_sam, ploidy_val) = line.split();
                sample = pat_sam.split("_")[1].replace('"','')
                if isfile(gammadir + ploidy + "/" + patient + "_" + sample + "_raw_segments.txt"):
                    ret.add((sample, gamma, ploidy))
    return ret

def readBalancedCalls(balanced_dir, patient, sample):
    balfile = open(balanced_dir + patient + "_" + sample + "_balanced_calls.tsv", "r")
    balcalls = {}
    for line in balfile:
        if "Chr" in line:
            continue
        (fpatient, fsample, chr, start, end, call) = line.rstrip().split()
        if chr not in balcalls:
            balcalls[chr] = {}
        balcalls[chr][(int(start), int(end))] = call
    return balcalls

def getPatientInfo():
    patientfile = open("Patient_status.txt", "r")
    progressions = {}
    sexes = {}
    for line in patientfile:
        if "Patient" in line:
            continue
        (patient, progression, sex) = line.split()
        progressions[patient] = progression
        sexes[patient] = sex
    return (progressions, sexes)

def getGammaFor(patient):
    #From WGS
    if patient=="772":
        return "550"
    #Can't get a tetraploid solution for this for any gamma:
#    if patient=="194":
#        return "1000"

    #Used to be on this list until we got updated flow data
#    if patient=="37":
#        return "450"
#    if patient=="891":
#        return "650"
    return "500"

def getPatientSampleMap(challenge=False, dipvtet_file=""):
    s2p = {}
    p2s = {}
    if dipvtet_file=="":
        dipvtet_file = "calling_evidence_challenge_inc_odds.tsv"
        if not challenge:
            dipvtet_file = "calling_evidence_odds.tsv"
        
    callfile = open(dipvtet_file, "r")
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

def loadDeletions(samplePatientMap, CNdir="noninteger_processed_CNs/"):
    deletions = {}
    for sample in samplePatientMap:
        patient, ploidy = samplePatientMap[sample]
        #ploidy = getBestPloidyFor(patient, sample)
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


def loadDeletionsAndCNVs(samplePatientMap, CNdir="noninteger_processed_CNs/"):
    deletions = {}
    CNVs = {}
    for sample in samplePatientMap:
        (patient, ploidy) = samplePatientMap[sample]
        ploidy = getBestPloidyFor(patient, sample)
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






