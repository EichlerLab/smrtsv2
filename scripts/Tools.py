#!/usr/bin/env python
import os
import re
import sys


def GetKV(key, vals):
    lk = len(key)
    for v in vals:
        if (len(v) >= len(key) and v[0:lk] == key):
            return v[lk:]
    else:
        return None

def GetStrand(value):
    if (value & 16 != 0):
        return 1
    else:
        return 0


def IsPrimary(value):
    if (value & 256 != 0):
        return 0
    else:
        return 1

def CIGARToArrays(cigar):
    ops = []
    lengths = []
    i1,i2 = 0,0
    end = len(cigar)
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0,len(opVals))]
    ops = [opVals[i][1] for i in range(0,len(opVals))]
    return ops, lengths

def ParseReadTitle(title):
    values = title.split('/')
    if (len(values) >= 3):
        coords = values[-1].split('_')
        coords = (int(coords[0]), int(coords[1]))
    else:
        coords = None
    return (values[0], int(values[1]), coords)


class Reference:
    def __init__(self, name, length, index):
        self.name = name
        self.length = length
        self.index = index


class SAMHeader:
    def __init__(self):
        self.references = {}
        self.IDS = []
    def ReadHeader(self, samFile):
        line = samFile.readline()
        index = 0
        while (len(line)> 0 and line[0] == '@'):
            if (line[0:3] == "@SQ"):
                vals = line.split("\t")
                valkvs = [kv.split(":") for kv in vals[1:]]
                self.references[valkvs[0][1]] = Reference(valkvs[0][1], valkvs[1][1], index)
                index += 1
            if (line[0:3] == "@PG"):
                self.IDS.append(line.strip())
            line = samFile.readline()
        return line

def SAMToAccuracy(cigar, readseq, refseq):
    q = 0
    t = 0
    refseq = refseq.upper()
    nMatch = 0
    nMisMatch =0
    nIns = 0
    nDel = 0
    #M	BAM_CMATCH	0
    #I	BAM_CINS	1
    #D	BAM_CDEL	2
    #N	BAM_CREF_SKIP	3
    #S	BAM_CSOFT_CLIP	4
    #H	BAM_CHARD_CLIP	5
    #P	BAM_CPAD	6
    #=	BAM_CEQUAL	7
    #X	BAM_CDIFF	8
    alnStart = 0

    nq = 0
    nt = 0
    for i in range(0,len(cigar)):
        if (cigar[i][0] == 4):
            q += cigar[i][1]
            alnStart = cigar[i][1]
        if (cigar[i][0] == 0):
            for j in range(0,cigar[i][1]):
                if (t >= len(refseq) or q >= len(readseq)):
                    print "ERROR at " + "\t".join([str(j) for j in [t, len(refseq), q, len(readseq), i, len(cigar)]])
                    print cigar
                if (refseq[t] == readseq[q]):
                    nMatch+=1
                else:
                    nMisMatch+=1
                q+=1
                t+=1
            nq += cigar[i][1]
        if (cigar[i][0] == 1):
            q += cigar[i][1]
            nIns += cigar[i][1]
            nq+= cigar[i][1]
        if (cigar[i][0] == 2):
            t += cigar[i][1]
            nDel += cigar[i][1]
    return float(nMatch - (nMisMatch + nIns + nDel))/nMatch

class SAMEntry:
    def __init__(self, line):
        v = ParseSamLine(line)
        if (v is None):
            self.title = None
            return None
        else:
            (self.title, self.flag, self.tName, self.tPos, self.mapqv, self.qStart, self.qEnd, self.readlen, self.seq, self.tlen)  = (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])

        vals = line.split("\t")
        self.cigar = vals[5]
        self.ops, self.lengths = CIGARToArrays(self.cigar)
        self.strand = GetStrand(self.flag)
        self.tLen = 0

        for i in range(len(self.ops)):
            if (self.ops[i] in ('M', 'D', '=', 'X')):
                self.tLen += self.lengths[i]
        prefixSoftClip = 0
        suffixSoftClip = 0

        # find the first soft clip cigar
        l = len(self.ops)
        if (l > 1 and self.ops[0] == 'S'):
            prefixSoftClip = self.lengths[0]
        elif (l > 2 and self.ops[1] == 'S'):
            prefixSoftClip = self.lengths[1]


        # The SAM alignment is in the direction of the target, so
        # the soft clipped end is the beginning of the reverse
        # strand.
        l = len(self.ops)
        if (l > 2 and self.ops[-1] == 'S' and self.ops[-2] != 'S'):
            suffixSoftClip = self.lengths[-1]
        elif (l > 3 and self.ops[-2] == 'S' and self.ops[-3] != 'S'):
            suffixSoftClip = self.lengths[-2]
        if (self.strand == 1):
            tmp = prefixSoftClip
            prefixSoftClip = suffixSoftClip
            suffixSoftClip = tmp
        self.qStart += prefixSoftClip
        self.qEnd   -= suffixSoftClip
        self.tStart = self.tPos
        self.tEnd = self.tPos + self.tLen
        self.line = line
        self.fullReadLength = GetKV("XQ:i:",vals[11:])
        self.vals = vals
        if (self.fullReadLength is not None):
            self.fullReadLength = int(self.fullReadLength)

    def PrintIntervals(self, out):
        out.write("{},{}\t{},{}\n".format(self.tStart , self.tEnd, self.qStart, self.qEnd))

titlei = 0
flagi = 1
tnamei = 2
tposi = 3
mapqvi = 4
qstarti = 5
qendi = 6
readleni = 7
seqi = 8
tleni = 9


def ParseSamLine(line):
    try:
        vals = line.split("\t")
        #if (vals[6] == "*" and vals[7] == "0" and vals[8] == "0"):
        #    return None
        title = vals[0]
        flag = int(vals[1])
        tName = vals[2]
        tPos  = int(vals[3])
        mapqv = int(vals[4])
        seq   = vals[9]
        idx = 11
        readlen = len(vals[9])

        start = GetKV("XS:i:", vals[11:])
        if (start is not None):
            start = int(start)
        else:
            start = 0
        end   = GetKV("XE:i:", vals[11:])

        if (end is not None):
            end = int(end)
        else:
            end = len(seq)

        tLen = int(vals[8])
    except:
        print "Error parsing"
        print line
        return None


    #       0      1     2      3     4      5      6    7,     8      9
    return (title, flag, tName, tPos, mapqv, start, end, readlen, seq, tLen)

def BuildAlignOpStrings(ops, lengths, qPos, tPos, qSeq):
    qStr = ""
    tStr = ""

    for i in range(0,len(ops)):
        if (ops[i] == "S"):
            qPos += lengths[i]
        elif (ops[i] == "M"):

            for j in range(0,lengths[i]):
                if (qPos >= len(qSeq)):
                    print "error at " + str(i)  + " of " + str(len(ops))
                    continue

                qStr += qSeq[qPos]
                tStr += "M"
                qPos += 1
                tPos += 1
        elif (ops[i] == "I"):
            for j in range(0,lengths[i]):
                if (qPos >= len(qSeq)):
                    print "insertion error at " + str(i) + " of " + str(len(ops))
                qStr += qSeq[qPos]
                tStr += "-"
                qPos += 1
        elif (ops[i] == "D"):
            for j in range(0,lengths[i]):
                qStr += "-"
                tStr += "D"
                tPos += 1
    return (qStr, tStr)

def Overlap( a, b):
    v = (a,b)
    if (a[1] < b[1]):
        i,j = 0,1
    else:
        i,j = 1,0

    if (v[i][1] < v[j][0]):
        return 0.0

    if (v[i][0] < v[j][0]):
        overlap = v[i][1] - v[j][0]
    else:
        overlap = v[i][1] - v[j][0]

    return abs(float(overlap))


def GapBetweenIntervals(a,b):
    if a[1] < b[1]:
        return b[0] - a[1]
    else:
        return a[0] - b[1]



def ReadFAIFile(faiFileName):
    fai = {}
    faiFile = open(faiFileName, 'r')
    for line in faiFile:
        vals = line.split()
        fai[vals[0]] = [int(vals[i]) for i in range(1,len(vals))]
    return fai


def ParseRegionStr(regionStr):

    a = regionStr.split(':')
    if (len(a) != 2):
        return None
    b = a[1].split("-")
    if (len(b) != 2):
        return None
    b[0] = b[0].replace(',','')
    b[1] = b[1].replace(',','')
    return (a[0], int(b[0]), int(b[1]))

def FormatRegion(regions):
    if (type(regions) == list and len(regions) >= 3):
        chrom = regions[0]
        start = regions[1]
        end   = regions[2]
        start = start.replace(',','')
        end   = end.replace(',', '')
    elif (type(regions) == list and len(regions) == 1):
        val = ParseRegionStr(regions[0])
        if (val is None):
            return None
        else:
            (chrom, start,end) = (val[0], val[1], val[2])
    elif (type(regions) == str):
        if (regions.find(":") >= 0):
            val = ParseRegionStr(regions)
            if (val is None):
                return None
            else:
                (chrom, start,end) = (val[0], int(val[1]), int(val[2]))
        else:
            vals = regions.split()
            (chrom, start, end) = (vals[0], int(vals[1]), int(vals[2]))
    else:
        return None
    start = int(start)
    end   = int(end)
    return (chrom,start,end)

def BedToRegion(bedline):
    return str(bedline[0]) + ":" + str(bedline[1]) + "-" + str(bedline[2])


def AddPreSuf(region, fai, pre, suf):
    newRegion = (region[0], max(0, region[1] - pre), min(fai[region[0]][0], region[2] + suf))
    return newRegion

def AddSlop(region, fai, slop):
    newRegion = (region[0], max(0, region[1] - slop), min(fai[region[0]][0], region[2] + slop))
    return newRegion


def ExtractSeq(region, seqFile, fai):
    if (region[0] not in fai):
        sys.stderr.write(region[0] + " missing from index file.\n")
        sys.exit(0)
    chrStart   = fai[region[0]][1]
    seqLength  = fai[region[0]][2]
    lineLength = fai[region[0]][3]

    startLine, startLinePos = region[1] / seqLength, region[1] % seqLength
    endLine, endLinePos = region[2] / seqLength, region[2] % seqLength

    startFilePos = chrStart + startLine * lineLength + startLinePos
    endFilePos   = chrStart + endLine * lineLength  + endLinePos
    if (startFilePos < 0):
        sys.stderr.write("ERROR! seeking before 0\n")
        sys.exit(0)

    if endFilePos < startFilePos:
        sys.stderr.write("ERROR! End position is less than start position (%s:%s-%s)\n" % (region[0], region[1], region[2]))
        sys.exit(0)

    seqFile.seek(startFilePos)

    seqNewLine = seqFile.read(endFilePos - startFilePos)
    seq = seqNewLine.replace("\n","")
    return seq

def FindBasFile(barcode, zmw, dirname):
	suffix = ""
	if (zmw >= 0 and zmw <= 54493):
	    suffix = ".1.bax.h5"
	elif (zmw >= 54494 and zmw <= 108987):
	    suffix = ".2.bax.h5"
	else:
	    suffix = ".3.bax.h5"

	fileName = dirname + "/" + barcode + suffix
	if (os.path.exists(fileName)):
	    return fileName
	else:
	    fileName = dirname + "/" + barcode + ".bas.h5"
	    if (os.path.exists(fileName)):
	        return fileName
	    else:
	        return None

def GetSoftClip(ops, lengths):
    firstS = len(ops)
    lastS  = len(ops)
    for i in range(0,len(ops)):
        if (ops[i] == 'S'):
            firstS = i
            break
    for i in range(len(ops)-1, firstS, -1):
        if (ops[i] == 'S'):
            lastS = i
            break
    if (firstS < 2):
        preSoftClip = lengths[firstS]
    else:
        preSoftClip = 0

    if (lastS >= len(ops)-2 and lastS < len(ops)):
        sufSoftClip = lengths[lastS]
    else:
        sufSoftClip = 0

    return (preSoftClip, sufSoftClip)
