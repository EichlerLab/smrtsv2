#!/usr/bin/env python

import sys
import argparse
ap = argparse.ArgumentParser(description="Format fasta for falcon reading")
ap.add_argument("--fakename", help="Copy the name of the first movie into all for local assembly.", default=False, action='store_true')
args = ap.parse_args()

movieIndex = 1
firstMovieName = None
for line in sys.stdin.readlines():
    if (line[0] == '>'):
        title = line
    else:
        if (args.fakename):
            if (firstMovieName is None):
                vals = title.split("/")
                firstMovieName = vals[0][1:]
            print ">{}/{}/1_{}".format(firstMovieName, movieIndex, len(line))
        else:
            print title.strip()
        movieIndex +=1

        L=50
        i = 0
        while (i < len(line)):
            e = min(i+L,len(line))
            print line[i:e].strip()
            i+=L
