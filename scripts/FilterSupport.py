#!/usr/bin/env python

import argparse
import sys
ap = argparse.ArgumentParser(description="Filter coverage bed, 5th field is hard stop, 9th is average coverage.")
ap.add_argument("--me", help="Min event", default=5, type=int)
ap.add_argument("--Me", help="Max event", default=20, type=int)
ap.add_argument("--mc", help="Min coverage", default=5, type=int)
ap.add_argument("--Mc", help="Max coverage", default=25, type=int)
ap.add_argument("--mr", help="Min hard stop to coverage ratio", default=0.5, type=float)
ap.add_argument("--Mr", help="Max hard stop to coverage ratio", default=1.5, type=float)
ap.add_argument("--il", help="Length index.", default=3, type=int)
ap.add_argument("--ie", help="Event index.", default=4, type=int)
ap.add_argument("--ic", help="coverage index.", default=8, type=int)
ap.add_argument("--ml", help="Accept all at least this length, -1 = everything", default = 300, type=int)
args = ap.parse_args()


for line in sys.stdin.readlines():
    vals = line.split()
    length = float(vals[args.il])
    h = float(vals[args.ie])
    c = float(vals[args.ic])
    if (length >= args.ml and (h >= args.me and h <= args.Me and c >= args.mc and c <= args.Mc and h/c > args.mr and h/c < args.Mr)):
        print line.strip()
