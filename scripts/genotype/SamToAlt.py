#!/usr/bin/env python3

"""
Translate a SAM file to a .alt file for BWA alt-aware mapping.
"""

import argparse
import re

CIGAR_MATCH_SET = {'M', '=', 'X'}


def cigar_seq_to_aln(in_file, out_file):
    """
    Translate sequence match operations (=, X) to alignment match operations (M). Neighboring match operations are
    coalesced.

    :param in_file: Open input file.
    :param out_file: Open output file.
    """

    line_count = 0

    # Process each record
    for line in in_file:

        line_count += 1

        line = line.strip()

        if not line:
            continue

        if line.startswith('@'):
            out_file.write(line)
            out_file.write('\n')
            continue

        # Initialize CIGAR transformation
        match_len = 0
        cigar_string = ''

        tok = line.split('\t')

        if len(tok) < 11:
            raise RuntimeError('Truncated SAM record on line {} ({}): Found {} fields (min = 11)'.format(
                line_count, tok[0], len(tok)
            ))

        # Process each CIGAR operation in the record
        for cigar_match in re.finditer('(\d+)([MIDNSHP=X])', tok[5]):

            cigar_count = cigar_match.group(1)
            cigar_op = cigar_match.group(2)

            if cigar_op in CIGAR_MATCH_SET:
                # Save match operations to be coalesced
                match_len += int(cigar_count)

            else:
                # Non-match operation, first output coalesced matches
                if match_len > 0:
                    cigar_string += '{}M'.format(match_len)
                    match_len = 0

                # Output current operation
                cigar_string += '{}{}'.format(cigar_count, cigar_op)

        # Output matches at the ends
        if match_len > 0:
            cigar_string += '{}M'.format(match_len)

        # Set new CIGAR and write
        tok[5] = cigar_string

        out_file.write('\t'.join(tok))
        out_file.write('\n')


# Main
if __name__ == '__main__':

    # Get arguments
    arg_parser = argparse.ArgumentParser(description='SAM to ALT (.alt) file for BWA alt-aware mapping.')

    arg_parser.add_argument('input',
                            help='SAM file of read alignments.')

    arg_parser.add_argument('output',
                            help='Output SAM file.')

    args = arg_parser.parse_args()

    # Open input and output files
    with open(args.input, 'r') as in_file:
        with open(args.output, 'w') as out_file:
            cigar_seq_to_aln(in_file, out_file)
