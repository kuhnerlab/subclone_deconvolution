#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:26:42 2020

@author: lpsmith

This script takes the trees from branch_lengths_and_signature_calculator.py,
plus the bar plots from make_branch_signature_bar_plots.py, and creates ete3
tree pictures with the signature bars printed over every branch.

"""

from __future__ import division
from os import walk
from os import path
from os import mkdir
from ete3 import Tree, TreeStyle, TextFace, AttrFace, ImgFace

onlysomepatients = False
somepatients = ["42"]

treedir = "branch_lengths/"
sigpngdir = "signature_pngs/"
treesigdir = "trees_with_signatures/"

if not path.isdir(treesigdir):
    mkdir(treesigdir)


def loadBranchSignaturePngs():
    """
    Load all the bar plot pictures, and use the filenames to store what
    patient and set of tips they belong to.
    """
    sigpngs = {}
    pngfiles = []
    for __, _, files in walk(sigpngdir):
        pngfiles += files
    
    for file in pngfiles:
        fvec = file.split(".")
        patient = fvec[0]
        if patient not in sigpngs:
            sigpngs[patient] = {}
        tipset = fvec[1:-2]
        tipset.sort()
        tipset = tuple(tipset)
        sigpngs[patient][tipset] = file
    return sigpngs


sigpngs = loadBranchSignaturePngs()


treefiles = []
for __, _, files in walk(treedir):
    treefiles += files


ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = False
ts.branch_vertical_margin = 20
ts.scale = 0.05


for file in treefiles:
    if "newick" not in file:
        continue
    if "all" in file:
        continue
    label = file.split(".")[0]
    if onlysomepatients and label not in somepatients:
        continue
    trees = []
    for line in open(treedir + file, "r"):
        if "#" in line:
            continue
        tree = Tree(line.rstrip())
        trees.append(tree)
    
    for index, tree in enumerate(trees):
        print(tree)
        for branch in tree.traverse():
            branch.dist = round(branch.dist)
        line = tree.write(format=1)
        rootlen = round(tree.get_tree_root().dist)
        line = line[:-1] + ":" + str(rootlen) + ";"
        for branch in tree.traverse():
            if branch.name != "":
                name_face = AttrFace("name", fsize=30)
                branch.add_face(name_face, column=3, position="branch-right")
                pngfile = sigpngs[label][tuple(branch.name.split("-")[0].split("_")[0:1])]
                if "-" in branch.name:
                    brackets = TextFace("** ", fsize=30)
                    branch.add_face(brackets, column=2)
                image_face = ImgFace(sigpngdir + pngfile, width=400)
                branch.add_face(image_face, column=1, position="branch-right")
            elif branch.dist != 0:
                tipset = []
                for tip in branch:
                    tipset.append(tip.name)
                tipset.sort()
                tipset = tuple(tipset)
                pngfile = sigpngs[label][tipset]
                image_face = ImgFace(sigpngdir + pngfile, width=400)
                branch.add_face(image_face, column=0, position="branch-top")
                
        filename = label + "_signatures"
        if len(trees)>1:
            filename += "_" + str(index)
        filename += ".png"
        tree.render(treesigdir + filename, tree_style = ts)
