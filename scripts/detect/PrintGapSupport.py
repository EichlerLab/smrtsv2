#!/bin/env python

import argparse
import numpy as np
import sys

ap = argparse.ArgumentParser(description="Print gap support from output of PrintGaps.py.")

ap.add_argument("table",
                help="Input tabular file")
ap.add_argument("out",
                help="Output file, stdout implies stdout")
ap.add_argument("--overlap", type=float, default=0.60,
                help="Required overlap consistency")
ap.add_argument("--minSupport", type=int, default=2,
                help="Min overlapping clusters")
ap.add_argument('--header', default=False, action='store_true',
                help='Print BED header line in output.')

args = ap.parse_args()

if args.out == "stdout":
    outFile = sys.stdout
else:
    outFile = open(args.out, 'w')

inFile = open(args.table, 'r')
prevChrom = ""
intv = ()
prevRead = ""

def Overlap(a, b, pct):

    if a[1] < b[0] or a[0] > b[1]:
        return False

    span = max(a[1], b[1]) - min(a[0], b[0])
    overlap = min(a[1], b[1]) - max(a[0], b[0])

    if float(overlap) / span >= pct:
        return True

intv = None
prevChrom = ""
prevRead  = ""
prevOp    = ""
strings   = []
names     = []
tsds      = []

# Write header
if args.header:
    outFile.write('#CHROM\tPOS\tEND\tMEAN_LEN\tSUPPORT\tSVTYPE\tSEQ_L\tREAD_L\tTSD_L\n')

for line in inFile:

    line = line.strip()

    if not line or line.startswith('#'):
        continue

    vals = line.split()

    if len(vals) == 7:
        vals.append(vals[6])
        vals[7] = "notsd"

    curChrom = vals[0]

    try:
        read = '/'.join(vals[7].split('/')[0:2])

    except:
        sys.stderr.write("Error joining %s\n" % vals)
        raise

    op = vals[3]

    intvStart = int(vals[1])
    intvEnd = int(vals[2])

    curIntv = (int(intvStart), int(intvEnd))

    if intv is not None and \
            Overlap(intv, curIntv, args.overlap) and \
            curChrom == prevChrom and \
            op == prevOp:

        if read != prevRead:
            intv = (min(intv[0], curIntv[0]), max(intv[1], curIntv[1]))
            support += 1
            strings.append(vals[5])
            tsds.append(vals[6])
            names.append(vals[7])

    else:
        uniqueSupport = len(np.unique(names))

        if intv is not None and uniqueSupport >= args.minSupport:
            meanLength = np.mean(np.asarray([len(seq) for seq in strings]))

            outFile.write(
                intvChrom + "\t" +
                str(intv[0]) + "\t" +
                str(intv[1]) + "\t" +
                str(meanLength) + "\t" +
                str(uniqueSupport) + "\t" +
                intvOp + "\t" +
                ';'.join(strings) + '\t' +
                ';'.join(names) + '\t' +
                ';'.join(tsds) + "\n"
            )

        support = 1
        intv = curIntv
        intvOp = op
        strings = [vals[5]]
        tsds = [vals[6]]

        names = [vals[7]]
        intvChrom = curChrom

    prevChrom = curChrom
    prevRead = read
    prevOp = op

if outFile != sys.stdout:
    outFile.close()
