#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:45:47 2020

@author: lpsmith
"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
from ete3 import Tree, TreeStyle, TextFace, AttrFace
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import csv

import lucianSNPLibrary as lps
#Alternative for importing the 'lucianSNPLibrary':
#import imp
#lps = imp.load_source("lps","/home/mkkuhner/Papers/phylo/lucianSNPLibrary.py")

onlysomepatients = False
somepatients = ["483"]
#somepatients = ["387", "483", "609", "611", "626", "631", "635", "652", "852"] #double split patients


sigfilename = "branch_signatures/all_signatures.tsv"

pngdir = "signature_pngs/"


if not path.isdir(pngdir):
    mkdir(pngdir)

colormap = {
    'SBS1': '#046764',
    'SBS2': '#70ffc4',
    'SBS3': '#72869f',
    'SBS5': '#fdf415',
    'SBS13': '#47a003',
    'SBS17a': '#fedaff',
    'SBS17b': '#a306a7',
    'SBS18': '#c90c43',
    'SBS34': '#fd5180',
    'SBS40': '#ffa963',
}


sigbars = pd.read_csv(sigfilename, sep="\t")
sigbars.set_index('Branch', inplace=True)
sigbars.head()

col = [colormap[x] for x in sigbars.loc[:, 'SBS1':'SBS40'].columns]

pids = list(set(sigbars.Patient))
pids.sort()
for pid in pids:
    branches = sigbars[sigbars.Patient==pid].index
    for branch in branches:
        f = plt.figure(figsize=(18,3))
        ax = f.add_subplot(111)
        ax.margins(0)
        sigs = sigbars.loc[[branch], :].iloc[::-1]
        sigs.loc[:, 'SBS1':'Other'].plot(kind='barh', stacked=True, ax=ax, color=col, width=1.2, linewidth=0)
        #ax.legend(bbox_to_anchor=(1.1,0.5))
        plt.legend(bbox_to_anchor=(1.1, 1., 1., .102), loc=2, ncol=1, borderaxespad=2.)
        plt.legend().remove()
        plt.axis('off')
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
#        plt.savefig(pngdir + "legend.png", bbox_inches='tight', pad_inches=0)
        plt.savefig(pngdir + str(pid) + "." + str(branch) + ".signature.png", bbox_inches = 'tight', pad_inches = 0)
#        plt.show()
#        plt.clear()

