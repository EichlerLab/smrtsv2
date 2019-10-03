#!/usr/bin/env python3

import argparse
import pysam


def filter_bam(input_file, queries, output_file, output_sizes_file, ref_filename):
    """
    Filter a BAM file.

    :param input_file: Input BAM file.
    :param queries: A set of query names to be extracted.
    :param output_file: Output SAM file.
    :param output_sizes_file: File query sizes are written to.
    :param ref_filename: Reference file name. Query sequences are aligned against this reference.
    """

    found_queries = set()

    with pysam.AlignmentFile(input_file, 'rb', reference_filename=ref_filename) as input_bam:
        with pysam.AlignmentFile(output_file, 'w', template=input_bam, reference_filename=ref_filename) as output_bam:
            with open(output_sizes_file, 'w') as output_sizes:
                for alignment in input_bam:
                    if alignment.query_name in queries:

                        # Check for multiple records with the same query
                        if alignment.query_name in found_queries:
                            raise ValueError('Query has multiple records in the input BAM: {}'.format(alignment.query_name))

                        found_queries.add(alignment.query_name)

                        output_bam.write(alignment)
                        output_sizes.write('{}\t{}\n'.format(alignment.query_name, alignment.query_length))

    # Check for missing queries
    missing_queries = queries - found_queries
    n_missing = len(missing_queries)

    if n_missing > 0:

        missing_str = ', '.join(sorted(list(missing_queries))[:3])

        if n_missing > 3:
            missing_str += '...'

        raise ValueError('Missing {} query sequences while filtering BAM: {}'.format(n_missing, missing_str))


def filter_sam(input_file, queries, output_file, output_sizes_file):
    """
    Filter a SAM file.

    :param input_file: Input SAM file.
    :param queries: A set of query names to be extracted.
    :param output_file: Output SAM file.
    :param output_sizes_file: File query sizes are written to.
    """

    found_queries = set()

    with open(input_file, 'r') as input_sam:
        with open(output_file, 'w') as output_sam:
            with open(output_sizes_file, 'w') as output_sizes:

                for line in input_sam:

                    if not line.strip():
                        continue

                    if line.startswith('@'):
                        output_sam.write(line)

                    tok = line.split('\t')

                    if tok[0] in queries:

                        # Check for multiple records with the same query
                        if tok[0] in found_queries:
                            raise ValueError('Query has multiple records in the input SAM: {}'.format(tok[0]))

                        found_queries.add(tok[0])

                        output_sam.write(line)
                        output_sizes.write('{}\t{}\n'.format(tok[0], len(tok[9])))

    # Check for missing queries
    missing_queries = queries - found_queries
    n_missing = len(missing_queries)

    if n_missing > 0:

        missing_str = ', '.join(sorted(list(missing_queries))[:3])

        if n_missing > 3:
            missing_str += '...'

        raise ValueError('Missing {} query sequences while filtering SAM: {}'.format(n_missing, missing_str))


if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(description='Extract contigs from an alignment file (SAM, BAM, or CRAM).')

    parser.add_argument('input_file', help='SAM or BAM to filter')
    parser.add_argument('queries_to_keep', help='Read names to extract (one per line)')
    parser.add_argument('output_file', help='filtered SAM')
    parser.add_argument('output_sizes', help='A table of contig names (col 1) and their sizes (col 2) with no header. '
                                             'Designed to be input into bedtools slop.')
    parser.add_argument('ref_filename', nargs='?', default=None, help='Reference for records are aligned against.')

    args = parser.parse_args()

    # Get query names
    with open(args.queries_to_keep, 'r') as fh:
        query_set = set([query.strip() for query in fh])

    if args.input_file.endswith('.bam') or args.input_file.endswith('.cram'):
        filter_bam(args.input_file, query_set, args.output_file, args.output_sizes, args.ref_filename)

    elif args.input_file.endswith('.sam'):
        filter_sam(args.input_file, query_set, args.output_file, args.output_sizes)

    else:
        raise RuntimeError('Unrecognized input file type (SAM, BAM, or CRAM): {}'.format(args.input_file))
