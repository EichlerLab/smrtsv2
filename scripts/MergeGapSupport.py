#!/usr/bin/env python


import argparse

import numpy as np
import math
import sys


import Tools

ap = argparse.ArgumentParser(description="Print gap support from output of PrintGaps.py.")
ap.add_argument("--table", help="Input tabular file, stdin reads from cl.", default="/dev/stdin")
ap.add_argument("--out", help="Output file, stdout implies stdout", default="/dev/stdout")
ap.add_argument("--overlap", help="Required overlap consistency", type=float, default=0.60)
ap.add_argument("--minSupport", help="Min overlapping clusters", type=int, default=2)
args = ap.parse_args()
outFile = open(args.out, 'w')
inFile = open(args.table)


prevLine = ""
prevVals = None
prevCoords = None
lineIndex = 0
for line in inFile.readlines():
    vals = line.split()
    (chrom, start, end ) = (vals[0], int(vals[1]), int(vals[2]))
#    vals[3] = 1
    overlap = 0
    if (prevCoords is not None):
        overlap = (2.0*Tools.Overlap((start,end), (prevCoords[1], prevCoords[2])))/ (end - start + prevCoords[2] - prevCoords[1])
    try:
        if (prevCoords is not None and chrom == prevCoords[0] and overlap > args.overlap):
            vals = [chrom, #0
                    min(start, prevCoords[1]), #1
                    max(end, prevCoords[2]), #2
                    (float(vals[3])*int(vals[4]) + float(prevVals[3])*int(prevVals[4])) / (int(vals[4]) + int(prevVals[4])), #3
                    int(vals[4])+ int(prevVals[4]), #4
                    vals[5], #5
                    ";".join([prevVals[6], vals[6]]), #6
                    ";".join([prevVals[7], vals[7]]), #7
                    ";".join([prevVals[8], vals[8]])] #8

        else:
            if (prevVals is not None):
                #            print str(prevVals)
                outFile.write("\t".join([str(i) for i in prevVals]) + "\n")
    except:
        print "ERROR at line " + str(lineIndex)
        sys.exit(0)
    prevVals = vals
    prevCoords = (vals[0], int(vals[1]), int(vals[2]))
    lineIndex += 1

if (outFile is not sys.stdout):
    outFile.close()

