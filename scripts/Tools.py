#!/usr/bin/env python

import re


def GetKV(key, vals):

    lk = len(key)

    for v in vals:
        if len(v) >= len(key) and v[0:lk] == key:
            return v[lk:]

    else:
        return None


def GetStrand(value):
    if value & 16 != 0:
        return 1
    else:
        return 0


def CIGARToArrays(cigar):
    ops = []
    lengths = []
    i1, i2 = 0, 0
    end = len(cigar)
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    return ops, lengths


class Reference:
    def __init__(self, name, length, index):
        self.name = name
        self.length = length
        self.index = index


class SAMEntry:
    def __init__(self, line):
        v = ParseSamLine(line)
        if v is None:
            self.title = None
            return None
        else:
            (self.title, self.flag, self.tName, self.tPos, self.mapqv, self.qStart, self.qEnd, self.readlen, self.seq, self.tlen) = (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])

        vals = line.split("\t")
        self.cigar = vals[5]
        self.ops, self.lengths = CIGARToArrays(self.cigar)
        self.strand = GetStrand(self.flag)
        self.tLen = 0

        for i in range(len(self.ops)):
            if self.ops[i] in ('M', 'D', '=', 'X'):
                self.tLen += self.lengths[i]

        prefixSoftClip = 0
        suffixSoftClip = 0

        # find the first soft clip cigar
        l = len(self.ops)
        if l > 1 and self.ops[0] == 'S':
            prefixSoftClip = self.lengths[0]

        elif l > 2 and self.ops[1] == 'S':
            prefixSoftClip = self.lengths[1]


        # The SAM alignment is in the direction of the target, so
        # the soft clipped end is the beginning of the reverse
        # strand.
        l = len(self.ops)
        if l > 2 and self.ops[-1] == 'S' and self.ops[-2] != 'S':
            suffixSoftClip = self.lengths[-1]

        elif l > 3 and self.ops[-2] == 'S' and self.ops[-3] != 'S':
            suffixSoftClip = self.lengths[-2]

        if self.strand == 1:
            tmp = prefixSoftClip
            prefixSoftClip = suffixSoftClip
            suffixSoftClip = tmp

        self.qStart += prefixSoftClip
        self.qEnd   -= suffixSoftClip
        self.tStart = self.tPos
        self.tEnd = self.tPos + self.tLen
        self.line = line
        self.fullReadLength = GetKV("XQ:i:", vals[11:])
        self.vals = vals

        if self.fullReadLength is not None:
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

        title = vals[0]
        flag = int(vals[1])
        tName = vals[2]
        tPos = int(vals[3])
        mapqv = int(vals[4])
        seq = vals[9]
        idx = 11
        readlen = len(vals[9])

        start = GetKV("XS:i:", vals[11:])

        if start is not None:
            start = int(start)
        else:
            start = 0

        end = GetKV("XE:i:", vals[11:])

        if end is not None:
            end = int(end)
        else:
            end = len(seq)

        tLen = int(vals[8])

    except:
        print("Error parsing")
        print(line)
        return None


    #       0      1     2      3     4      5      6    7,     8      9
    return (title, flag, tName, tPos, mapqv, start, end, readlen, seq, tLen)


def Overlap(a, b):

    v = (a, b)

    if a[1] < b[1]:
        i, j = 0, 1

    else:
        i, j = 1, 0

    if v[i][1] < v[j][0]:
        return 0.0

    if v[i][0] < v[j][0]:
        overlap = v[i][1] - v[j][0]
    else:
        overlap = v[i][1] - v[j][0]

    return abs(float(overlap))


def read_fai_file(fai_file_name):
    """
    Get a dictionary of an FAI file (reference index).

    :param fai_file_name: Name of the file to read.

    :return: A dictionary keyed on the sequence names with one tuple of the record values per record (all integers).
    """

    fai_record_dict = {}

    with open(fai_file_name, 'r') as fai_file:

        for line in fai_file:

            line = line.strip()

            if not line:
                continue

            fai_record = line.split()

            fai_record_dict[fai_record[0]] = [int(x) for x in fai_record[1:]]

    return fai_record_dict


def parse_region_string(region):
    """
    Separate a region, such as "chr1:1000-2000", into a tuple of "name", "pos", and "end". The coordinates will
    be stripped of commas, if present, and cast as integers.

    :param region: Region string.

    :return: Tuple of "name", "pos", and "end" or `None` if the string could not be parsed.
    """

    a = region.split(':')

    if len(a) != 2:
        return None

    b = a[1].split("-")

    if len(b) != 2:
        return None

    b[0] = b[0].replace(',', '')
    b[1] = b[1].replace(',', '')

    return a[0], int(b[0]), int(b[1])


def FormatRegion(regions):

    if type(regions) == list and len(regions) >= 3:

        chrom = regions[0]
        start = regions[1]
        end = regions[2]

        start = start.replace(',', '')
        end = end.replace(',', '')

    elif type(regions) == list and len(regions) == 1:
        val = parse_region_string(regions[0])

        if (val is None):
            return None
        else:
            (chrom, start, end) = (val[0], val[1], val[2])

    elif type(regions) == str:
        if regions.find(":") >= 0:
            val = parse_region_string(regions)
            if val is None:
                return None

            else:
                (chrom, start, end) = (val[0], int(val[1]), int(val[2]))

        else:
            vals = regions.split()
            (chrom, start, end) = (vals[0], int(vals[1]), int(vals[2]))
    else:
        return None

    start = int(start)
    end = int(end)

    return chrom, start, end


def BedToRegion(bedline):
    return str(bedline[0]) + ":" + str(bedline[1]) + "-" + str(bedline[2])


def extract_sequence(region, sequence_file, fai):
    """
    Get a sequence region from an indexed sequence file.

    :param region: A tuple of the genomic coordinates with elements: name, pos, end.
    :param sequence_file: Open file to extract sequences from.
    :param fai: Dictionary of FAI records keyed by sequence name (see `read_fai_file()`).

    :return: Sequence in `seqFile` specified by `region`.
    """

    # Check FAI for region
    if region[0] not in fai:
        raise ValueError('Sequence {} is missing in the index'.format(region[0]))

    # Get coordinates and lengths
    chr_start = int(fai[region[0]][1])
    seq_len = int(fai[region[0]][2])
    line_len = int(fai[region[0]][3])

    region_pos = int(region[1])
    region_end = int(region[2])

    # Calculate psotions
    start_line = int(region_pos / seq_len)
    start_line_pos = region_pos % seq_len

    end_line = int(region_end / seq_len)
    end_line_pos = region_end % seq_len

    start_file_pos = chr_start + start_line * line_len + start_line_pos
    end_file_pos = chr_start + end_line * line_len + end_line_pos

    # Check file positions
    if start_file_pos < 0:
        raise ValueError('Region {0}:{1}-{2} attempts to seek before 0 in the sequence file'.format(*region))

    if end_file_pos < start_file_pos:
        raise ValueError(
            'Region {0}:{1}-{2} attempts to seek before the start position of its record in the '
            'sequence file'.format(*region)
        )

    # Read sequence
    sequence_file.seek(start_file_pos)

    return sequence_file.read(end_file_pos - start_file_pos).replace('\n', '')
