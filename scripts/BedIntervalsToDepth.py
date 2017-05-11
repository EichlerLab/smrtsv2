#!/usr/bin/env python

import argparse
import Tools
import sys
import numpy as np

ap = argparse.ArgumentParser(description="Use a bed file to compute depth.")
ap.add_argument("bed", help="Input bed file.")
ap.add_argument("genome", help="Genome file with .fai")
ap.add_argument("--bin", help="Bin size", type=int, default=100)
ap.add_argument("--out", help="Output file.", default="/dev/null")

args = ap.parse_args()

bed = open(args.bed)

fai = Tools.read_fai_file(args.genome + ".fai")

outFile = open(args.out, 'w')

bins = {}
for k in fai.keys():
    l = fai[k][0]
    bins[k] = np.zeros(l / args.bin + 1)

lineNumber = 0

for line in bed:
    vals = line.split()
    chrom = vals[0]
    start = int(vals[1])
    end = int(vals[2])

    if chrom in bins:
        b = bins[chrom]
        startBin = start / args.bin
        endBin = end / args.bin
        b[startBin:endBin+1] += 1

    lineNumber += 1

    if lineNumber % 100000 == 0:
        sys.stderr.write("processed {} \n".format(lineNumber))

for k in bins.keys():

    b = bins[k]

    for i in range(0,len(b)):
        if b[i] > 0:

            outFile.write("{}\t{}\t{}\t{}\n".format(k, i*args.bin, (i+1)*args.bin, b[i]))

outFile.close()
