#!/bin/env python
"""
Report the positions of all gaps (N bases) in the given FASTA sequence(s).
"""
import argparse
from Bio import SeqIO


def find_gaps(input_filename):
    # Load the original FASTA sequence.
    fasta = SeqIO.parse(input_filename, "fasta")

    for record in fasta:
        sequence = str(record.seq)

        for i in xrange(len(record)):
            if sequence[i].upper() == "N":
                print "\t".join(map(str, (record.id, i, i + 1)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    args = parser.parse_args()

    find_gaps(args.input_file)
