#!/usr/bin/env python

import argparse
import intervaltree
import sys
ap = argparse.ArgumentParser(description="Print gaps from a gaps.short.bed file, but cross reference by which contig is in the minimal tiling path of contigs")
ap.add_argument("gaps", help="Bed file of gaps.short.bed, the 8th column needs to be the contig name.")
ap.add_argument("tiling", help="Tiling path.  The 4th column is the contig neame.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")

args = ap.parse_args()


gaps = open(args.gaps)
tiling = open(args.tiling)
out = open(args.out, 'w')


chromTiling = {}
for line in tiling:
    vals = line.split()
    if (vals[0] not in chromTiling):
        chromTiling[vals[0]] = intervaltree.IntervalTree()
    chromTiling[vals[0]].addi(int(vals[1]), int(vals[2]), vals[3])


for line in gaps:
    vals = line.split()
    if (vals[0] not in chromTiling):
        continue
    intvs = chromTiling[vals[0]]
    ovps = intvs.search(int(vals[1]))
    if (len(ovps)== 0):
        sys.stderr.write("Missing for " + line)
        continue
    if (len(ovps) > 1):
        sys.stderr.write( "Multiple overlaps for " + line)
        continue
    for o in ovps:
        if (o[2] == vals[7]):
            print line.strip()
    



