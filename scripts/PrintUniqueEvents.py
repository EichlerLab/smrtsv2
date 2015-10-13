#!/usr/bin/env python

#chr1	9230436	9230711	insertion	275	AAAAAATTAGCCAGGTGTGGTGGTGGTGCCTGTAGTCCAGCTACTAGGAGGCTATGGAGGAGATGGAGAACAGGAGCAGAGGTGAGTGAGCGAGATCACGCACTGCACCCAGCTGGGTGACAAAAAAAACCATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGTAATGGCTGGTGCGGTGCTCACACTGTAATCCAGCACTTGAGCCAAGAGGTGATCATAGGTCAGCATCCAGACAGCTGGCAAATAGTAAACCTGTCTCTACTAAAAATA	AluYh9:INC,AluYb:INC	AAAAATTA
import sys
import argparse

ap = argparse.ArgumentParser(description="Print insertions based on how many annotations there are.")
ap.add_argument("infile", help="input bed file, 7th column is the annotation.")
ap.add_argument("--complex", help="Print complex insertions.", default=False, action='store_true')
ap.add_argument("--exact", help="Print items matching exactly this number", type=int, default=None)
ap.add_argument("--max", help="Print items matching at most this number", type=int, default=None)
ap.add_argument("--min", help="Print items matching at least this number", type=int, default=None)
ap.add_argument("--prefix", help="Count events labeled by this prefix, e.g. L1P", default=None)
ap.add_argument("--minPrefix", help="Count events with at least this many matching", type=int, default=None)
ap.add_argument("--maxPrefix", help="Count events with at most this many matching", type=int, default=None)
ap.add_argument("--maxNotPrefix", help="Maximum number of events not matching this prefix.", type=int, default=None)
ap.add_argument("--minNotPrefix", help="Minimum number of events not matching this prefix.", type=int, default=None)
ap.add_argument("--remainder", help="Print anything that does not match to this file.", default=None)
ap.add_argument("--maxSTR", help="Count STRs separately, maximum allowed", default=None, type=int)
ap.add_argument("--minSTR", help="Minimum number of STRs to keep.", default=None, type=int)
ap.add_argument("--secondary", help="Count secondary elements", default=None)
ap.add_argument("--minSecondary", help="Require at least this many secondary, for example, AluY is primary, and Alu is secondary.", default=None, type=int)
ap.add_argument("--maxSecondary", help="Require no more than this many secondary.", default=None, type=int)
ap.add_argument("-v", help="Complement of any previous logic.", action='store_true', default=False)

args = ap.parse_args()

inFile =  open(args.infile, 'r')

if (args.remainder is not None):
    remainder = open(args.remainder, 'w')
else:
    remainder = None
import pdb

for line in inFile:
    vals = line.split()
    ann = vals[6]
    nann = len(ann.split(';'))
    annotations = ann.split(';')


    doPrint = False
    
    if (True):
        nPrefix = 0
        nNotPrefix = 0
        nSTR = 0
        nSecondary = 0
        lp = 0
        if (args.prefix is not None):
            lp = len(args.prefix)
        ls = 0
        if (args.secondary is not None):
            ls = len(args.secondary)

        for a in annotations:
            if (args.prefix is not None and a[0:lp] == args.prefix):
                nPrefix +=1
            elif (args.secondary is not None and a[0:ls] == args.secondary):
                nSecondary += 1
            elif (a.find(")n") != -1 or a.find("_rich") != -1 or a.find("-rich") != -1):
                nSTR += 1
                if (args.minNotPrefix is None and args.maxNotPrefix is None):
                    nNotPrefix += 1
            else:
                nNotPrefix += 1

        if ((args.max == None or len(annotations) < args.max) and
            (args.min == None or len(annotations) >= args.min) and
            (args.minPrefix == None or nPrefix >= args.minPrefix) and
            (args.maxNotPrefix == None or nNotPrefix <= args.maxNotPrefix) and
            (args.maxPrefix == None or nPrefix <= args.maxPrefix) and
            (args.minNotPrefix == None or nNotPrefix >= args.minNotPrefix) and
            (args.minSTR == None or nSTR >= args.minSTR) and
            (args.maxSTR == None or nSTR <= args.maxSTR) and
            (args.minSecondary == None or nSTR >= args.minSecondary) and
            (args.maxSecondary == None or nSTR <= args.maxSecondary) ):
            doPrint = True

    if ((args.complex == False and ((args.exact is not None and nann == args.exact) or (args.max is not None and nann <= args.max))) or (nann > 1 and args.complex)):
        doPrint = True


    if (args.v == True):
        if (doPrint):
            doPrint = False
        else:
            doPrint = True
    if (doPrint):
        print line.strip()
    elif (remainder is not None):
        remainder.write(line)
