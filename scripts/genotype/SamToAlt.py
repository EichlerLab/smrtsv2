#!/usr/bin/env python3

"""
Translate a SAM file to a .alt file for BWA alt-aware mapping.
"""

import argparse
import pysam

# Declarations
CIGAR_DICT = {
    0: 'M',
    1: 'I',
    2: 'D',
    3: 'N',
    4: 'S',
    5: 'H',
    6: 'P',
    7: '=',
    8: 'X',
    9: 'B',
}


def cigar_seq_to_aln(in_file, out_file):
    """
    Translate sequence match operations (=, X) to alignment match operations (M). Neighboring match operations are
    coalesced.

    :param in_file: Open input pysam alignment file.
    :param out_file: Open output pysam alignment file.
    """

    # Process each record
    for record in in_file.fetch():

        match_len = 0
        cigar_string = ''

        # Process each CIGAR operation in the record
        for cigar_record in record.cigartuples:

            if cigar_record[0] in (0, 7, 8):
                # Save match operations to be coalesced
                match_len += cigar_record[1]

            else:
                # Non-match operation, first output coalesced matches
                if match_len > 0:
                    cigar_string += '{}M'.format(match_len)
                    match_len = 0

                # Output current operation
                cigar_string += '{}{}'.format(cigar_record[1], CIGAR_DICT[cigar_record[0]])

        # Output matches at the ends
        if match_len > 0:
            cigar_string += '{}M'.format(match_len)

        # Set new CIGAR and write
        record.cigarstring = cigar_string
        out_file.write(record)


# Main
if __name__ == '__main__':

    # Get arguments
    arg_parser = argparse.ArgumentParser(description='SAM to ALT (.alt) file for BWA alt-aware mapping.')

    arg_parser.add_argument('input',
                            help='SAM or BAM file of read alignments.')

    arg_parser.add_argument('output',
                            help='Output SAM or BAM file.')

    args = arg_parser.parse_args()

    # Get input and output file types
    if args.input.endswith('bam'):
        in_file_type = 'rb'
    else:
        in_file_type = 'r'

    # Get output file and type
    if args.output.endswith('bam'):
        out_file_type = 'wb'
    else:
        out_file_type = 'w'

    # Open input and output files
    with pysam.AlignmentFile(args.input, in_file_type) as in_file:
        with pysam.AlignmentFile(args.output, out_file_type, in_file) as out_file:
            cigar_seq_to_aln(in_file, out_file)
