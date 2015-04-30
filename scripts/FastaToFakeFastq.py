#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

import sys
if (len(sys.argv) < 2):
    print "usage: FastaToFakeFastq.py input.fasta output.fastq"
    sys.exit(1)

inFile = open(sys.argv[1],'r')
outFile = open(sys.argv[2], 'w')
index = 0
for rec in SeqIO.parse(inFile, "fasta"):
    rec.letter_annotations['phred_quality'] =  [ 40 ]*len(rec.seq)
    qual = 'I'*len(rec.seq)
    outFile.write("@" + str(index) + "\n")
    outFile.write(str(rec.seq) + "\n")
    outFile.write("+" + str(index) + "\n")
    outFile.write(qual + "\n")
    index+=1

outFile.close()
