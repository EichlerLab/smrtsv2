#!/usr/bin/env python

import sys
import argparse
import Tools
import Align
from  Bio import SeqIO


ap = argparse.ArgumentParser(description="Print gaps in a SAM file.")
ap.add_argument("genome", help="Genome file with a .fai")
ap.add_argument("sam", help="Sam file of alignment.", nargs="+")
ap.add_argument("--onTarget", help="Assume the query encodes the position of the aligned sequence, and make sure at least the chromosomes match.", default=False, action='store_true')
ap.add_argument("--gapFree", help="Print sequences without gaps.", default=None)
ap.add_argument("--minContigLength", help="Only parse alignments from contigs this length or greater", default=0, type=int)
ap.add_argument("--minLength", help="Minimum gap length.", default=50, type=int)
ap.add_argument("--maxLength", help="Maximum gap length.", default=None, type=int)
ap.add_argument("--outFile", help="Print output here, default= stdout", default=None)
ap.add_argument("--context", help="Print surrounding context", default=0, type=int)
ap.add_argument("--condense", help="Pack indels if the matches separating them is less than this value.", default=10, type=int)
ap.add_argument("--tsd", help="Attempt to find Target Site Duplications at most this length", default=20, type=int)
ap.add_argument("--outsam", help="Write the modified condensed sam to a file.", default=None)
ap.add_argument("--minq", help="Minimal mapping quality to consider (10)",default=10,type=int)
ap.add_argument("--qpos", help="Write query position of gaps", default=False,action='store_true')
ap.add_argument("--snv", help="Print SNVs to this file.", default=None)
ap.add_argument("--contigBed", help="Print where contigs map.", default=None)
ap.add_argument("--status", help="Print how far along the alignments are.", default=False, action='store_true')
ap.add_argument("--blacklist", help="Exclude contigs on this list from callsets.", default=None)
ap.add_argument("--removeAdjacentIndels", help="Find instances of SNVs pushed into indels, in the format: NIXMND., and remove these operations.", default=False, action='store_true')
args = ap.parse_args()

genome = file(args.genome, 'r')
handle = open(args.genome, "r")



if (args.outFile is None):
    outFile = sys.stdout
else:
    outFile = open(args.outFile, 'w')

if (args.gapFree is not None):
    gapFree = open(args.gapFree, 'w')

if (args.contigBed is not None):
    contigBed = open(args.contigBed, 'w')

blacklist = {}
if (args.blacklist is not None):
    bl = open(args.blacklist)
    blacklist = { i.strip() : true for i in bl.readlines() }
    
fai = Tools.ReadFAIFile(args.genome + ".fai")

#genomeDict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

if (args.outsam is not None):
    outsam = open(args.outsam, 'w')

snvOut = None
if (args.snv is not None):
    snvOut = open(args.snv, 'w')

fai = Tools.ReadFAIFile(args.genome + ".fai")
genomeFile = open(args.genome, 'r')

M = 'M'
I = 'I'
D = 'D'
N = 'N'
S = 'S'
H = 'H'
P = 'P'
if (args.sam[0].find(".fofn") >= 0):
    fofnFile = open(args.sam[0])
    samFiles = [line.strip() for line in fofnFile.readlines()]
    args.sam = samFiles
lineNumber = 0
contextLength = 8
for samFileName in args.sam:
    samFile = open(samFileName)

    for line in samFile:
        lineNumber = lineNumber + 1
        sys.stderr.write(str(lineNumber) + "\n")
        if (line[0] == "@"):
            if (args.outsam is not None):
                outsam.write(line)
            continue
        if (len(line) <= 1):
            continue

        aln = Tools.SAMEntry(line)
        if (aln.title is None):
            continue

        if (args.onTarget == True):
            srcChrom = aln.title.split(":")[0]
            if (srcChrom != aln.tName):
                continue
        if (aln.flag & 256 != 0):
            continue
        if (aln.mapqv < args.minq):
            continue

        if (args.contigBed is not None):
            contigBed.write("{}\t{}\t{}\t{}\n".format(aln.tName, aln.tStart, aln.tStart + aln.tlen, aln.title))

        tPos = aln.tStart
        qPos = 0
        #
        # condense matches.
        #
        packedCigar = []
        i = 0
        i1 = 1
        niter = 0
        maxGap = 0
        maxGapType = 0
        packedOps = []
        packedLengths = []
        #print str(aln.ops)

        foundGap = False

        while (i < len(aln.lengths)):
            l = aln.lengths[i]
            op = aln.ops[i]
            j = i
            if (op == I or op == D):
                if (l > maxGap):
                    maxGap = l
                    maxGapType = op

            if (op == I or op == D and i < len(aln.ops) - 2 and aln.ops[i+2][0] == op):
                matchLen = 0
                gapLen   = 0
                while (j+2 < len(aln.ops) and aln.ops[j+2][0] == op and aln.ops[j+1][0] == M and aln.lengths[j+1] < args.condense):

                    matchLen += aln.lengths[j+1]
                    gapLen   += aln.lengths[j+2]
                    j+=2
                if (j > i):
                    newIndel = (op, l+gapLen)
                    newMatch = (M, matchLen)
                    packedOps.append(op)
                    packedLengths.append(l+gapLen)

                    packedOps.append(M)
                    packedLengths.append(matchLen)

                else:
                    packedLengths.append(l)
                    packedOps.append(op)

            else:
                packedLengths.append(l)
                packedOps.append(op)

            i = j + 1
            niter +=1
            if (niter > len(aln.ops)):
                sys.stderr.write("ERROR! too many interations.\n")

        for i in range(len(packedOps)):
            op = packedOps[i]
            l  = packedLengths[i]
            if (op == S):
                qPos += l

            if (op == N):
                tPos += l
                qPos += l
            if (op == M):
                # Inside match block (if op == M)
                if (args.snv is not None):
                    targetSeq = Tools.ExtractSeq((aln.tName, tPos-1,tPos+l-1), genomeFile, fai)
                    querySeq  = aln.seq[qPos:qPos+l]
                    nMis = 0

                    for mp in range(0,len(targetSeq)):
                        if (querySeq[mp].upper() != targetSeq[mp].upper() and targetSeq[mp].upper() != 'N' and querySeq[mp].upper() != 'N'):
                            nMis +=1
                            snvOut.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(aln.tName, tPos+mp, tPos+mp+1, targetSeq[mp], querySeq[mp], aln.title ))
                tPos += l
                qPos += l

            if (op == I):
                if (l > args.minLength and (args.maxLength is None or l < args.maxLength)):
                    #                    print "gap at " + str(i) + " " + str(op) + " " + str(l) + " " + str(qPos) + " qlen: " + str(len(aln.seq))
                    foundGap = True
                    chrName = aln.tName
                    #gapSeq = aln.seq[max(0,qPos-args.context):min(qPos+l+args.context, len(aln.seq))]
                    gapSeq = aln.seq[qPos:qPos+l]

                    context= aln.seq[qPos+l:min(qPos+l+args.context, len(aln.seq))]
                    if (context == "A"*len(context) or context == "T"*len(context)):
                        homopolymer="T"
                    else:
                        homopolymer="F"
                    tsd = "notsd"
                    if (len(gapSeq) == 0):
                        sys.stderr.write("ERROR, gap seq is of zero length\n")
                    if (args.tsd):
                        # try and find the target site duplications, this may be on either side of the alignemnt
                        tsdSuffix = gapSeq[-args.tsd:]
                        tsdSuffix = tsdSuffix.upper()
                        tsdPrefix = gapSeq[0:args.tsd]
                        tsdPrefix = tsdPrefix.upper()
                        targetPrefix = Tools.ExtractSeq((chrName, tPos-args.tsd,tPos), genomeFile, fai)
#                        targetPrefix = genomeDict[chrName].seq[tPos-args.tsd:tPos]

                        targetPrefix = targetPrefix.upper()
                        #targetSuffix = genomeDict[chrName].seq[tPos:tPos+args.tsd]
                        targetSuffix = Tools.ExtractSeq((chrName, tPos,tPos+args.tsd), genomeFile, fai)
                        targetSuffix = targetSuffix.upper()
                        (sp, ss, sScore) = Align.TSDAlign(tsdSuffix, targetPrefix, 'suffix')
                        (pp, ps, pScore) = Align.TSDAlign(tsdPrefix, targetSuffix, 'prefix')
                        if (sScore > pScore ):
                            tsd = ss
                        elif (pScore > sScore ):
                            tsd = ps

                    nucs = ['A', 'C', 'G', 'T']
                    fracs = [float(gapSeq.count(n))/(len(gapSeq)+1) for n in nucs]
                    doPrint = True
                    for frac in fracs:
                        if (frac > 0.85):
                            doPrint = False
                    if (doPrint):
                        outFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrName, tPos, tPos + l, "insertion", l, gapSeq, tsd, aln.title, qPos, qPos + l))
                        if (args.context > 0):
                            outFile.write("\t{}".format(homopolymer))
                        if (args.qpos):
                            outFile.write("\t{}\t{}".format(qPos, qPos + len(gapSeq)))
                        outFile.write("\n")


                qPos += l
            if (op == D):
                if (l > args.minLength and (args.maxLength is None or l < args.maxLength)):
                    foundGap = True
                    chrName = aln.tName
                    if (tPos > fai[chrName][0]):
                        sys.stderr.write("ERROR! tpos is past the genome end." + str(tPos) + " " + str(fai[chrName][0]) + "\n")
                    #delStart = max(tPos - args.context, 0)
                    #delEnd   = min(tPos + args.context + l, fai[chrName][0])
                    delStart = max(tPos - args.context, 0)
                    delEnd   = min(tPos + args.context + l, fai[chrName][0])
                    context= aln.seq[qPos+l:min(qPos+l+args.context, len(aln.seq))]
                    if (context == "A"*len(context) or context == "T"*len(context)):
                        homopolymer="T"
                    else:
                        homopolymer="F"

                    #delSeq = genomeDict[chrName].seq[delStart:delEnd].tostring()
                    delSeq = Tools.ExtractSeq([chrName, delStart, delEnd], genomeFile, fai)

                    outFile.write("{}\t{}\t{}\t{}\t{}\t{}\tno_tsd\t{}\t{}\t{}".format(chrName, tPos, tPos + l, "deletion", l, delSeq, aln.title, qPos, qPos + l))
                    if (args.context > 0):
                        outFile.write("\t{}".format(homopolymer))

                    if (args.qpos):
                        outFile.write("\t{}\t{}".format(qPos, qPos))
                    outFile.write("\n")

                tPos += l
            if (op == H):
                pass

        if (foundGap == False and args.gapFree is not None):
            gapFree.write(aln.tName + "\t" + str(aln.tStart) + "\t" + str(aln.tEnd) + "\t" + aln.title + "\n")

        if (args.outsam is not None):
            packedCigar= ''.join([str(v[0]) + str(v[1]) for v in zip(packedLengths, packedOps)])
            vals = line.split()
            packedLine = '\t'.join(vals[0:5]) + "\t" + packedCigar + '\t'.join(vals[6:]) + "\n"
            outsam.write(packedLine)

if (args.gapFree is not None):
    gapFree.close()

outFile.close()
if (args.outsam is not None):
    outsam.close()
