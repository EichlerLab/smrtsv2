#!/usr/bin/env python3

import argparse
import os
import pysam

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio import SeqIO


# Main
if __name__ == '__main__':

    # Get arguments
    parser = argparse.ArgumentParser(description="Format fasta for falcon reading")

    parser.add_argument('bam', help='BAM input file.')

    parser.add_argument('fasta', help='FASTA output file.')

    parser.add_argument('--fastq', help='Also write records as a FASTQ file if set.')

    parser.add_argument('--fakename',
                        default=False, action='store_true',
                        help="Copy the name of the first movie into all for local assembly."
                        )

    parser.add_argument('--qual', default=40, type=int,
                        help='Phred scaled quality score to assign to each base when writing a FASTQ file.')

    args = parser.parse_args()

    # Initialize record list
    record_list = list()

    # Read BAM file
    with pysam.AlignmentFile(args.bam, 'r') as sam_file:

        movie_name = None
        read_count = 0

        for read in sam_file.fetch():

            read_count += 1

            # Set fake name by taking the first read name and recyciling it with a ZMW equal to read_count
            # Read name format: moviename/zmw/start_stop
            if args.fakename:

                if movie_name is None:
                    movie_name = read.query_name.split('/')[0]

                read.query_name = '{}/{}/0_{}'.format(movie_name, read_count, len(read.query_sequence))

            record_list.append(
                SeqRecord(
                    Seq(read.query_sequence),
                    id=read.query_name,
                    description=''
                )
            )

    # Write FASTA
    os.makedirs(os.path.dirname(args.fasta), exist_ok=True)

    with open(args.fasta, 'w') as fasta_file:
        SeqIO.write(record_list, fasta_file, 'fasta')

    # Write FASTQ
    if args.fastq is not None:

        # Add fake quality scores
        for record in record_list:
            record.letter_annotations["phred_quality"] = [args.qual] * len(record)

        # Write
        with open(args.fastq, 'w') as fastq_file:
            SeqIO.write(record_list, fastq_file, 'fastq')
