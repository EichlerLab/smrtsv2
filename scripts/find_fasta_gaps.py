#!/bin/env python
"""
Report the positions of all gaps (N bases) in the given FASTA sequence(s).
"""
import argparse
from Bio import SeqIO


def _range(n):
    # range in python2 takes too much memory when iterating over a large set of integers. range() in python3 does what
    # xrange() did in python2, and xrange() is not in python3. To keep this script from eating too much memory and to
    # enable both major versions, _range() is defined.
    val = 0

    while val < n:
        yield val
        val += 1


def find_gaps(input_filename):
    # Load the original FASTA sequence.
    fasta = SeqIO.parse(input_filename, "fasta")

    for record in fasta:
        sequence = str(record.seq)
        gap_start = None  # Not in a gap

        for i in _range(len(record)):

            if sequence[i].upper() == "N":
                if gap_start is None:
                    gap_start = i
            else:
                if gap_start is not None:
                    print("\t".join(map(str, (record.id, gap_start, i))))
                    gap_start = None

        if gap_start is not None:
            print("\t".join(map(str, (record.id, gap_start, len(record)))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    args = parser.parse_args()

    find_gaps(args.input_file)
