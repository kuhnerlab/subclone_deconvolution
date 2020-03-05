# Count origins and classify them

import sys
if len(sys.argv} != 2:
  print "USAGE:  python origins.py treedir"
  exit()

inputdir = sys.argv[1]

legal_timepoints = ["T1","T2","T3"]


# read in P/NP data
pdict = {}
for line in open("TS_statuses.tsv","r"):
  line = line.rstrip().split()
  # header
  if line[] == "Patient":
    pid_index = line.index("Patient")
    sid_index = line.index("Sample")
    prog_index = line.index("Prog")
    tp_index = line.index("Time")
    continue
  # non-header
  pid = line[pid_index]
  sid = line[sid_index]
  if sid.endswith("N"):  continue   # normal samples
  prog = line[prog_index]
  tp = line[tp_index]
  if tp not in legal_timepoints:  continue  
  if pid not in pdict:
    pdict[pid] = prog

# read in the trees and count origins
import os
rawfiles = []
for root,dirs,files in os.walk(inputdir):
  for file in files:
    if file.endswith(".newick"):
      rawfiles.append(file)
   
trees4 = []
trees_all = []
for pid in pdict.keys():
  if pid+".newick" not in rawfiles:
    print "File not found:",pid+".newick"
    exit()
  trees4.append(pid+".newick")
  if pid+"_all.newick" in rawfiles:
    trees_all.append(pid+"_all.newick")
  else:
    trees_all.append(pid+".newick")
