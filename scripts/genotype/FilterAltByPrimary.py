#!/usr/bin/env python3

"""
Filter an ALT file (SAM format) for all contigs on a given primary contig
"""

import argparse


def filter_contigs(in_file, out_file, primary):
    """
    Filter contigs aligned to a primary contig.

    :param in_file: Open input file.
    :param out_file: Open output file.
    :param primary: Primary contig name.
    """

    for line in in_file:

        if not line.strip():
            continue

        tok = line.split('\t')

        # Process header line
        if tok[0].startswith('@'):

            # Only output sequence headers that match the primary contig
            if tok[0] == '@SQ':
                if tok[1].split(':')[1] == primary:
                    out_file.write(line)
            else:
                # Write all other headers
                out_file.write(line)

            continue

        # Process record
        if tok[2] == primary:
            out_file.write(line)

# Main
if __name__ == '__main__':

    # Get arguments
    arg_parser = argparse.ArgumentParser(
        description='Filter ALT (SAM format) for all contigs aligned to a given primary contig. All contigs and header'
                    'sequences that do not match the primary contig are removed.'
    )

    arg_parser.add_argument('input',
                            help='ALT (SAM) of all contigs.')

    arg_parser.add_argument('output',
                            help='Filtered ALT.')

    arg_parser.add_argument('primary',
                            help='Name of the primary contig.')

    args = arg_parser.parse_args()

    # Filter
    with open(args.input, 'r') as in_file:
        with open(args.output, 'w') as out_file:
            filter_contigs(in_file, out_file, args.primary)
