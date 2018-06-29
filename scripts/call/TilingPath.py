#!/usr/bin/env python3

import argparse
import intervaltree
import sys
import os
import inspect

# Append smrtsvlib to path
smrtsv_base = os.path.abspath(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe()))
    ))),
)

sys.path.append(smrtsv_base)

from smrtsvlib import tools

ap = argparse.ArgumentParser(description="Given a SAM file of locally aligned contigs, determine which contigs best overlap the genome segmented by contig alignment positions. ")
ap.add_argument("sam", help="Alignments of the tiled assembly.")
ap.add_argument("--minContigLength", help="Minimum contig length.", default=0,type=int)
args = ap.parse_args()

samFile = open(args.sam)

chromIntervals = {}
chromLengths = {}
lineNumber = 0
chroms = []

for line in samFile:
    vals = line.split()

    if line[0] == '@':
        if vals[0] == "@SQ":
            sn = vals[1].split(":")[1]
            chromLen = int(vals[2].split(":")[1])

            if chromLen < args.minContigLength:
                continue

            chromLengths[sn] = chromLen
            chromIntervals[sn] = intervaltree.IntervalTree()
            chroms.append(sn)

        # Header: Stop parsing the line
        continue

    samEntry = tools.SAMEntry(line)

    if samEntry.title is not None and samEntry.tLen > 0:
        chromIntervals[samEntry.tName].addi(samEntry.tStart, samEntry.tEnd, samEntry.title)

for chrom in chroms:
    intvs = chromIntervals[chrom]
    chromPos = {0: 1, chromLengths[chrom]: 1}

    for intv in intvs.items():
        chromPos[intv[0]] = 1
        chromPos[intv[1]] = 1

    chromPos = sorted(chromPos.keys())

    i = 0
    j = i
    nCondense = 0

    # Ignore chromosomes with no intervals (i.e., no contigs mapped to this
    # chromosome).
    if len(chromPos) == 2:
        continue

    for i in range(0, len(chromPos)-1):

        midPoint = (chromPos[i] + chromPos[i+1])/2

        # Ignore intervals that are too small.
        if chromPos[i+1] - chromPos[i] < 2:
            continue

        ovpMidPoint = intvs.search(midPoint)
        opt = None
        optDist = 0

        if len(ovpMidPoint) == 0:
            continue

        for intv in ovpMidPoint:
            dist = abs(((intv[1]+intv[0])/2) - midPoint)
            if opt is None or dist < optDist:
                optDist = dist
                opt = intv[2]

        print(chrom + "\t" + str(chromPos[i]) + "\t" + str(chromPos[i+1]) + "\t" + opt + "\t" + str(optDist))
