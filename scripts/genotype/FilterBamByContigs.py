#!/usr/bin/env python3

import argparse
import pysam

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(description='Extract contigs from a BAM file.')

    parser.add_argument('input_bam', help='BAM to filter')
    parser.add_argument('queries_to_keep', help='Read names to extract (one per line)')
    parser.add_argument('output_bam', help='filtered BAM')
    parser.add_argument('output_sizes', help='A table of contig names (col 1) and their sizes (col 2) with no header. '
                                             'Desigend to be input into bedtool slop.')
    args = parser.parse_args()

    # Get query names
    with open(args.queries_to_keep, 'r') as fh:
        queries = set([query.strip() for query in fh])

    # Get reads
    with pysam.AlignmentFile(args.input_bam, 'rb') as input_bam:
        with pysam.AlignmentFile(args.output_bam, 'wb', template=input_bam) as output_bam:
            with open(args.output_sizes, 'w') as output_sizes:
                for alignment in input_bam:
                    if alignment.query_name in queries:
                        output_bam.write(alignment)
                        output_sizes.write('{}\t{}\n'.format(alignment.query_name, alignment.query_length))
