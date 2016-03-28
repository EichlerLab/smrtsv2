#!/usr/bin/env python
"""
Recalculate total repeatmasked bases in the given variants based on the number
of lowercased ATCG bases.
"""
import sys
idx=11
if (len(sys.argv) == 2):
	idx = int(sys.argv[1])

while sys.stdin:
	line = sys.stdin.readline()
	if (line == ""):
		break
	v = line.split()
	nl = sum([v[5].count(i) for i in ['a','g','c','t']])
	t = len(v[5])
	v[11] = "{:2.2f}".format(float(nl)/t)
	sys.stdout.write('\t'.join(v) + '\n')
