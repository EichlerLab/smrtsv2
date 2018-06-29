#!/usr/bin/env python

import Bio.SeqIO
import Bio.SeqRecord
import Bio.Seq
import sys
import argparse

#format
# 0     1       2       3   4   5
#chr10	603788	603835	ins	47	gggggagggggcaggggcaggcaggcaggcagggggagggggcaggg

ap = argparse.ArgumentParser(description="Print gap sequences to fasta files.")
ap.add_argument("bed", help="Input bed file. Position of gap sequence defaults to 5")
ap.add_argument("insertion", help="Output insertion file.")
ap.add_argument("--deletion", help="Output deletion file.", default=None)
ap.add_argument("--seqidx", help="Index of sequence. (5)", default=5, type=int)
ap.add_argument("--indelidx", help="Index of insertion/deletion annotaiton. (3)", type=int, default=3)
ap.add_argument("--unmask", help="Print all upper case", default=False, action='store_true')

args = ap.parse_args()
bedFile = open(args.bed, 'r')
insertionsFile = open(args.insertion, 'w')
if (args.deletion is not None):
    deletionsFile=open(args.deletion, 'w')

for line in bedFile:
    vals = line.split()
    seqstr = vals[args.seqidx].split(';')[0]

    namestr = '/'.join(vals[0:3])
    seq = Bio.Seq.Seq(seqstr)
    if (args.unmask):
        seq = seq.upper()
        
    rec = Bio.SeqRecord.SeqRecord(seq, id=namestr,name="",description="")
    if (vals[args.indelidx] == "insertion"):
        Bio.SeqIO.write(rec, insertionsFile, "fasta")
    if (vals[args.indelidx] == "deletion"):
        if (args.deletion is not None):
            Bio.SeqIO.write(rec, deletionsFile, "fasta")
        else:
            Bio.SeqIO.write(rec, insertionsFile, "fasta")
            
insertionsFile.close()
if (args.deletion is not None):
    deletionsFile.close()
