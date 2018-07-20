#!/usr/bin/env python
import argparse
import sys
ap = argparse.ArgumentParser(description="Modify a bed file")
ap.add_argument("bed", help="Input bed file.")
ap.add_argument("out", help="Output file")
ap.add_argument("--leftjustify", help="Print left side of bed coordinates plus delta", type=int, default=None)
ap.add_argument("--rightjustify", help="Print right side of bed coordinates only.", type=int, default=None)
ap.add_argument("--leftslop", help="Subtract this from the left side of every coordinate", type=int, default=0)
ap.add_argument("--rightslop", help="Add this to the right side of every coordinate", type=int, default=0)
args = ap.parse_args()

outFile = open(args.out, 'w')
inFile  = open(args.bed, 'r')
for line in inFile:
    vals = line.split()
    vals[1] = int(vals[1])
    vals[2] = int(vals[2])
    if (args.leftjustify is not None):
        vals[2] = int(vals[1]) + args.leftjustify
    elif (args.rightjustify is not None):
        vals[1] = int(vals[2]) - args.rightjustify

    vals[1] = max(0, vals[1] - args.leftslop)
    vals[2] += args.rightslop
        
    outFile.write(vals[0] + "\t" + str(vals[1]) + "\t" + str(vals[2]) + "\t" + "\t".join(vals[3:]) + "\n")


