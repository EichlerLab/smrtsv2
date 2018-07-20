#!/usr/bin/env python

import argparse
from Bio import SeqIO

ap = argparse.ArgumentParser(description='Print gap sequences to FASTA files.')

ap.add_argument('bedin', help='Input bed file.')
ap.add_argument('bedout', help='Output bed file.')
ap.add_argument('dotout', help='RepeatMasker file.out annotation file.')
ap.add_argument('maskedout', help='Masked output file.', default=None)
ap.add_argument('--seqidx', help='Index of gap sequence.', default=6, type=int)

args = ap.parse_args()

bedFileIn = open(args.bedin, 'r')
bedFileOut = open(args.bedout, 'w')

annotations = {}

dotoutFile = open(args.dotout, 'r')
maskedDict = {}

if args.maskedout is not None:
    maskedSequences = open(args.maskedout)
    maskedDict = SeqIO.to_dict(SeqIO.parse(maskedSequences, 'fasta'))


for i in range(3):
    dotoutFile.readline()

for line in dotoutFile:
    vals = line.split()
    name = vals[4]
    rep = vals[9]
    pre = vals[11]
    post = vals[13]
    pre = int(pre.replace('(', '').replace(')', ''))
    post = int(post.replace('(', '').replace(')', ''))

    if pre + post < 30:
        rep += ':FULL'
    else:
        rep += ':INC'

    if name not in annotations:
        annotations[name] = []

    annotations[name].append(rep)


for line in bedFileIn:
    vals = line.split()
    name = '/'.join(vals[0:3])

    if name in annotations:
        annotation = ';'.join(annotations[name])
    else:
        annotation = 'NONE'

    repeatContent = ''

    if 'name in maskedDict':
        vals[5] = str(maskedDict[name].seq)
        repeatContent = '\t{:2.2f}'.format(
            float(vals[5].count('a') + vals[5].count('c') + vals[5].count('g') + vals[5].count('t')) / len(vals[5])
        )

    bedFileOut.write(
        '\t'.join(vals[0:args.seqidx]) + '\t' +
        annotation + '\t' +
        '\t'.join(vals[args.seqidx:]) +
        repeatContent + '\n'
    )

bedFileOut.close()
