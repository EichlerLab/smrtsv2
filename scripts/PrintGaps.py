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
ap.add_argument("--condense", help="Pack indels if the matches separating them is less than this value.", default=0, type=int)
ap.add_argument("--tsd", help="Attempt to find Target Site Duplications at most this length", default=20, type=int)
ap.add_argument("--outsam", help="Write the modified condensed sam to a file.", default=None)
ap.add_argument("--minq", help="Minimal mapping quality to consider (10)",default=10,type=int)
ap.add_argument("--qpos", help="Write query position of gaps", default=False,action='store_true')
ap.add_argument("--snv", help="Print SNVs to this file.", default=None)
ap.add_argument("--nloc", help="Print locations of aligned N's here.", default=None)
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
    for line in bl:
        v = line.split()
        if (v[0] not in blacklist):
            blacklist[v[0]] = []
        if (len(v) > 1):
            blacklist[v[0]].append(int(v[1])+1)

fai = Tools.ReadFAIFile(args.genome + ".fai")

#genomeDict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

if (args.outsam is not None):
    outsam = open(args.outsam, 'w')

snvOut = None
if (args.snv is not None):
    snvOut = open(args.snv, 'w')

nLocOut = None
if (args.nloc is not None):
    nLocOut = open(args.nloc, 'w')

fai = Tools.ReadFAIFile(args.genome + ".fai")
genomeFile = open(args.genome, 'r')

M = 'M'
X = 'X'
E = '='
I = 'I'
D = 'D'
N = 'N'
S = 'S'
H = 'H'
P = 'P'
def IsMatch(c):
    return (c == M or c == X or c == E)

if (args.sam[0].find(".fofn") >= 0):
    fofnFile = open(args.sam[0])
    samFiles = [line.strip() for line in fofnFile.readlines()]
    args.sam = samFiles
lineNumber = 0
contextLength = 8
#import pdb
import re
coordRe = re.compile(".*(chr.*)\.(\d+)-(\d+).*")
for samFileName in args.sam:
    samFile = open(samFileName)

    for line in samFile:
        lineNumber = lineNumber + 1

        if (line[0] == "@"):
            if (args.outsam is not None):
                outsam.write(line)
            continue
        if (len(line) <= 1):
            continue

        aln = Tools.SAMEntry(line)
        if (aln.title is None):
            continue
	#
	# Use 0-based coordinate system
	#
        aln.tStart -=1
        aln.tEnd -=1
        if (args.onTarget == True):
            coordReMatch = coordRe.match(aln.title)
            if (coordReMatch is not None):
                coordMatchGroups = coordReMatch.groups()
                srcChrom = coordMatchGroups[0]
                srcStart = int(coordMatchGroups[1])
                srcEnd   = int(coordMatchGroups[2])
                if (srcChrom != aln.tName):
                    print "off target chromosome: " + srcChrom + " " + aln.tName
                    continue
                if (((srcStart >= aln.tStart and srcStart < aln.tEnd) or (srcEnd >= aln.tStart and srcEnd < aln.tEnd) or (srcStart < aln.tStart and srcEnd > aln.tEnd )) == False):
                    print "no overlap " + srcChrom + " " + str(srcStart) + " " + str(srcEnd) + " alignment: " + str(aln.tStart) + " "+ str(aln.tEnd)
                    continue

        if (aln.flag & 256 != 0):
            continue
        if (aln.mapqv < args.minq):
            continue

        if (args.contigBed is not None):
            contigBed.write("{}\t{}\t{}\t{}\n".format(aln.tName, aln.tStart, aln.tStart + aln.tlen, aln.title))

        if (args.minContigLength > len(aln.seq)):
            continue


        if (args.blacklist is not None):
            if (aln.title in blacklist):
                if (len(blacklist[aln.title]) == 0):
                    sys.stderr.write("Skipping " + aln.title + " in blacklist.\n")
                    continue
                else:
                    foundPos = False

                    for p in blacklist[aln.title]:
                        if int(aln.tPos) == p:
                            foundPos = True
                            break

                    if foundPos:
                        sys.stderr.write("Skipping " + aln.title + " in blacklist.\n")
                        continue

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
        #print str(aln.ops)

        foundGap = False

        if (args.removeAdjacentIndels):
            for i in range(1,len(aln.lengths)-1):
                if (aln.ops[i-1] != 'M' and
                    aln.ops[i+1] != 'M' and
                    aln.ops[i-1] != aln.ops[i+1] and
                    aln.ops[i] == 'M' and
                    aln.lengths[i-1] == aln.lengths[i+1] and
                    aln.lengths[i] < 4):
                    aln.lengths[i-1] = 0
                    aln.lengths[i+1] = 0

            newLengths = []
            newOps = []
            for i in range(0,len(aln.lengths)):
                if (aln.lengths[i] != 0):
                    newLengths.append(aln.lengths[i])
                    newOps.append(aln.ops[i])
            aln.lengths = newLengths
            aln.ops = newOps

        packedOps = []
        packedLengths = []
        i = 0
        if (args.condense > 0):
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
                    while (j+2 < len(aln.ops) and aln.ops[j+2][0] == op and IsMatch(aln.ops[j+1][0])  and aln.lengths[j+1] < args.condense):

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
                    print "ERROR! too many interations."
        else:
            packedOps = aln.ops
            packedLengths = aln.lengths

        for i in range(len(packedOps)):
            op = packedOps[i]
            oplen  = packedLengths[i]

            if (op == N or op == S):
                # Inside match block (if op == M)
                qPos += oplen
            if (IsMatch(op)):
                # Inside match block (if op == M)
                if (args.snv is not None):
                    targetSeq = Tools.ExtractSeq((aln.tName, tPos,tPos+oplen), genomeFile, fai)
                    querySeq  = aln.seq[qPos:qPos+oplen]
                    nMis = 0

                    for mp in range(0,len(targetSeq)):
                        if (querySeq[mp].upper() != targetSeq[mp].upper() and targetSeq[mp].upper() != 'N' and querySeq[mp].upper() != 'N'):
                            nMis +=1
                            snvOut.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(aln.tName, tPos+mp, tPos+mp+1, targetSeq[mp], querySeq[mp], aln.title, mp+qPos ))
                        if (args.nloc is not None and (targetSeq[mp].upper() == 'N' or querySeq[mp].upper() == 'N')):
                            nLocOut.write("{}\t{}\t{}\n".format(aln.tName, tPos+mp,tPos+mp+1));

                tPos += oplen
                qPos += oplen

            if (op == I):
                if (oplen >= args.minLength and (args.maxLength is None or oplen < args.maxLength)):

                    foundGap = True
                    chrName = aln.tName
                    #gapSeq = aln.seq[max(0,qPos-args.context):min(qPos+oplen+args.context, len(aln.seq))]
                    gapSeq = aln.seq[qPos:qPos+oplen]

                    context= aln.seq[qPos+oplen:min(qPos+oplen+args.context, len(aln.seq))]
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
                        if (tsd == ""):
                            print line
                        outFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrName, tPos, tPos + oplen, "insertion", oplen, gapSeq, tsd, aln.title, qPos, qPos + oplen))
                        if (args.context > 0):
                            outFile.write("\t{}".format(homopolymer))
                        if (args.qpos):
                            outFile.write("\t{}\t{}".format(qPos, qPos + len(gapSeq)))
                        outFile.write("\n")


                qPos += oplen
            if (op == D):
                if (oplen >= args.minLength and (args.maxLength is None or oplen < args.maxLength)):
                    foundGap = True
                    chrName = aln.tName
                    if (tPos > fai[chrName][0]):
                        sys.stderr.write("ERROR! tpos is past the genome end." + str(tPos) + " " + str(fai[chrName][0]) + "\n")
                    #delStart = max(tPos - args.context, 0)
                    #delEnd   = min(tPos + args.context + l, fai[chrName][0])
                    delStart = max(tPos - args.context, 0)
                    delEnd   = min(tPos + args.context + oplen, fai[chrName][0])
                    if (delEnd < delStart):
                        continue
                    context= aln.seq[qPos+oplen:min(qPos+oplen+args.context, len(aln.seq))]
                    if (context == "A"*len(context) or context == "T"*len(context)):
                        homopolymer="T"
                    else:
                        homopolymer="F"

                    #delSeq = genomeDict[chrName].seq[delStart:delEnd].tostring()
                    delSeq = Tools.ExtractSeq([chrName, delStart, delEnd], genomeFile, fai)

                    outFile.write("{}\t{}\t{}\t{}\t{}\t{}\tno_tsd\t{}\t{}\t{}".format(chrName, tPos, tPos + oplen, "deletion", oplen, delSeq, aln.title, qPos, qPos))
                    if (args.context > 0):
                        outFile.write("\t{}".format(homopolymer))

                    if (args.qpos):
                        outFile.write("\t{}\t{}".format(qPos, qPos))
                    outFile.write("\n")

                tPos += oplen
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
