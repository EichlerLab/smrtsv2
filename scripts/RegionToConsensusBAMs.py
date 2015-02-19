#!/usr/bin/env python

import argparse
import sys
import Tools
import subprocess
import tempfile
import pdb
import os

def Call(command, logfile=None):
    print "calling " + command
    if (logfile is not None):
        if (isinstance(logfile, str)):
            outFile = open(logfile, 'w')
        elif (isinstance(logfile, file)):
            outFile = logfile
        else:
            print "bad input."
            sys.exit(1)

        retval = subprocess.call(command, shell=True, stdout=outFile)
        outFile.close()
    else:
        retval = subprocess.call(command.split())
    if (retval != 0):
        print "Error running " + command
        raise ValueError
    return 0

def DeleteFiles(tempfiles):
    for filename in tempfiles:
        command = "/bin/rm -f {}".format(filename)
        retval = subprocess.call(command.split())

ap = argparse.ArgumentParser(description="Create a consensus from high mapqv viewed from a bam file. To run this, "
                             "you must have samtools, quiver, cmph5tools.py in your path.")

ap.add_argument("bam_list", help="File with list of paths to BAMs generated from blasr SAM output that contains full quality values.")
ap.add_argument("--region", help="Region from the bam file.  If 3 options, assume bed format, otherwise, UCSC chr:start-end format.", nargs="+", default=None, required=True)

ap.add_argument("--reference", help="Use this reference file.", default=None)
ap.add_argument("--referenceWindow", help="Use this window from the reference, similar to region, can be 3 arguments or 1.", nargs="+", default=None)

ap.add_argument("--consensus", help="Write consensus to this file.  The special file name \"fromregion\" implies build the refrence name from the region", default="consensus.fa")

ap.add_argument("--minq", help="Minimum mapping quality (20)", default=20, type=int)
ap.add_argument("--tmpdir", help="Create temporary files here.", default=".")
ap.add_argument("--keeptemp", help="Do not delete temporary files.", default=False,action='store_true')
ap.add_argument("--delta", help="Increase/decrease region by delta for fetching reads.", default=0,type=int)
ap.add_argument("--p5c3", help="Force p5c5 chemistry flag.", action='store_true', default=False)

args = ap.parse_args()

regionSamFileName   = tempfile.mktemp(suffix=".sam", dir=args.tmpdir)
regionBasH5FileName = tempfile.mktemp(suffix=".bas.h5", dir=args.tmpdir)

if (args.region is None):
    print "Required argument region is missing."
    sys.exit(1)

bedRegion = Tools.FormatRegion(args.region)
print args.region
print bedRegion
expandedRegion = (bedRegion[0], max(0, bedRegion[1] - args.delta), int(bedRegion[2]) + args.delta)
region = Tools.BedToRegion(expandedRegion)

path="/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp"

if (os.path.exists(args.tmpdir) == False):
    os.makedirs(args.tmpdir)

tempfiles = []
try:
    with open(args.bam_list, "r") as fh:
        bam_list = [line.strip() for line in fh]

    tempfiles.append(regionSamFileName)
    Call("/net/eichler/vol4/home/jlhudd/src/bamtools-2.3.0/bin/bamtools filter -in {} -region {} -mapQuality \">={}\" | samtools view -".format(" -in ".join(bam_list), region.replace("-", ".."), args.minq), regionSamFileName)
    tempfiles.append(regionBasH5FileName)
    Call("{}/pbihdfutils/bin/samtobas {} {}".format(path, regionSamFileName, regionBasH5FileName))

    # build a reference if necessary
    if (args.referenceWindow is not None):
        tempReference = tempfile.mktemp(suffix=".reference.fasta", dir=args.tmpdir)
        refRegion = Tools.BedToRegion(Tools.FormatRegion(args.referenceWindow))
        tempfiles.append(tempReference)
        Call("samtools faidx {} {} ".format(args.reference, refRegion), tempReference)
        args.reference = tempReference

    alignmentSamFileName = tempfile.mktemp(suffix=".alignment.sam", dir=args.tmpdir)

    cmph5FileName = tempfile.mktemp(suffix=".cmp.h5", dir=args.tmpdir)

    tempfiles.append(alignmentSamFileName)

    Call("{}/alignment/bin/blasr {} {} -minAlignLength 1000 -sam -bestn 1 -nproc 6 -out {}".format(path, regionBasH5FileName, args.reference, alignmentSamFileName))
    tempfiles.append(cmph5FileName)
    Call("{}/pbihdfutils/bin/samtoh5 {} {} {}".format(path, alignmentSamFileName, args.reference, cmph5FileName))
    Call("cmph5tools.py sort --deep {}".format(cmph5FileName))
    Call("{}/pbihdfutils/bin/loadPulses {} {} -metrics InsertionQV,DeletionQV,SubstitutionQV,MergeQV,SubstitutionTag,DeletionTag".format(path, regionBasH5FileName, cmph5FileName))


    if (args.consensus == "fromregion"):
        tmpRegion = Tools.FormatRegion(args.region)
        args.consensus = "{}_{}-{}.fasta".format(tmpRegion[0], tmpRegion[1], tmpRegion[2])

    p5C3Flag = ""
    if (args.p5c3 == True):
        p5C3Flag = "-p P5-C3.AllQVsMergingByChannelModel"
    Call("samtools faidx {}".format(args.reference))
    Call("quiver -j 6 -r {} {} -o {} {}  ".format(args.reference, cmph5FileName, args.consensus, p5C3Flag))

except ValueError:
    if (args.keeptemp == False):
        DeleteFiles(tempfiles)
    sys.exit(1)

if (args.keeptemp == False):
    DeleteFiles(tempfiles)
