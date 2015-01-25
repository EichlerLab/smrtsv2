#!/usr/bin/env python

import sys
import Tools
import pdb
import argparse

ap = argparse.ArgumentParser(description="Remove overlapping entries, with options on how to determine overlap.")
ap.add_argument("bed", help="Input bed")
ap.add_argument("bedout", help="Output bed")
ap.add_argument("--leftjustify", help="Only look to see if left (lower) coordinate overlaps, using window.", default=False, action='store_true')
ap.add_argument("--rightjustify", help="Only look to see if right (upper) coordinate overlaps, using window.", default=False, action='store_true')
ap.add_argument("--exact", help="Just remove lines with this fraction of overlap.", default=None, type=float)
ap.add_argument("--window", help="Add this to coordinates to see if they overlap.", type=int, default=20)
ap.add_argument("--addLength", help="Add length in this column to the start", type=int, default=None)
ap.add_argument("-v", help="Opposite (like grep -v), print what would have been removed from ovelrap.",action='store_true', default=False)
ap.add_argument("--sourceIndex", help="Make sure that overlapping elements are not from the same source contig.  This specifies the index to check.", type=int, default=None)

args = ap.parse_args()

inFile = open(args.bed)
outFile = open(args.bedout, 'w')

pc = ""
ps = 0
pe = 0
ovpps = 0
ovppe = 0

maxLine = ""

def Print(l,o,t):
    
        o.write(l)

prevSource = None
curSource  = None
ovp=False
line = None
for line in inFile:
    vals = line.split()
    cc = vals[0]
    cs = int(vals[1])
    ce = int(vals[2])
    if (args.addLength is not None):
        ce += int(vals[args.addLength])
    if (args.sourceIndex is not None):
        curSource = vals[args.sourceIndex]
    ovpcs = cs
    ovpce = ce
    if (args.leftjustify):
        ovpce = cs
    if (args.rightjustify):
        ovpcs = ce
        
    ovp = int(Tools.Overlap( (ovpcs - args.window, ovpce+args.window ), (ovpps - args.window, ovppe+args.window)))
#    print str(ovp) + "\t" + str(ovpcs - args.window) + " " + str(ovpce+args.window ) + " " + str(ovpps - args.window) + " " + str(ovppe+args.window) 
    doPrint = True
    if (ovp == 0):
        ovp = False
    else:
        ovp = True

    #
    # Check to see if overlapping calls must be from separate source alignments.
    #
    if (ovp == True and args.sourceIndex is not None):
        if (prevSource == curSource):
            ovp = False

    prevSource = curSource
        
    if (args.exact is not None):
        if (ovp == False):
            outFile.write(line)
        elif (cs != ps and ce != pe):
            outFile.write(line)
            
    else:        
        if (ovp == False):
            if (args.v == False):
                outFile.write(maxLine)
            pc = cc
            ps = cs
            pe = ce

            maxLine = line
        else:
            if (ce - cs > pe - ps):
                # current line is longer than the previous, keep that one,
                if (args.v == True):
                    outFile.write(maxLine)
                maxLine = line
                # record this entry as the one that should be checked for overlaps.
                pc = cc
                ps = cs
                pe = ce
            else:
                if (args.v == True):
                    outFile.write(line)
                
    ovppe = ovpce
    ovpps = ovpcs

if (line is not None and ovp == False and args.v == False):
    outFile.write(line)
