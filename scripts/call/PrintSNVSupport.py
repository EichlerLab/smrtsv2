#!/usr/bin/env python

import sys
import argparse
ap = argparse.ArgumentParser(description="Print support for snv or indel.")
ap.add_argument("bed", help="Input bed file.")
ap.add_argument("out", help="Output bed file.")

args = ap.parse_args()

bed = open(args.bed)
out = open(args.out, 'w')


merge = False
support = []
prevVals = None
first = True
for line in bed:
    vals = line.split()

    merge = False
    if (first == False and vals[0:5] == prevVals[0:5]):
        support.append(vals[5])
        merge = True

    
    if (first != True and merge == False):
        v = prevVals[0:5] + [str(len(support))] + [";".join(support)]
        
        out.write("\t".join(v) + "\n")
        support = []

    prevVals = vals;
    if (merge == False):
        support.append(vals[5])
    first = False

if (merge == False):
    v = prevVals[0:5] + [str(len(support))] + [";".join(support)]
    out.write("\t".join(v) + "\n")
    
