#!/usr/bin/env python3

import argparse
import pysam


if __name__ == "__main__":

    # Get arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('fastq')
    parser.add_argument('--max_mask_prop', type=float, default=0.5)

    args = parser.parse_args()

    with pysam.FastqFile(args.fastq) as fh:
        for record in fh:
            if record.sequence.count('N') < args.max_mask_prop * len(record.sequence):
                print('@%s' % record.name)
                print(record.sequence)
                print('+')
                print(record.quality)
