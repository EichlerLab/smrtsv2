#!/usr/bin/env python3

"""
Get insert size differences over SV breakpoints on the reference.
"""

import argparse
import gc
import numpy as np
import os
import pandas as pd
import pysam

# Number of SVs to process before resetting pysam (close and re-open file). Avoids a memory leak in pysam.
PYSAM_RESET_INTERVAL = 1000


def get_insert_size_list(sv_record, bam_file, ref_flank):
    """
    Get a list of insert sizes over an SV.

    :param sv_record: SV record.
    :param bam_file: BAM file to extract reads from.
    :param ref_flank: Number of bases upstream and downstream of the SV that should be considered. This value should
        be larger than a reasonable upper-limit of the insert-size distribution.

    :return: A list of insert sizes for paired-end reads that mapped across the SV.
    """

    # Get window at each breakpoint
    bp_l_lower = max(0, sv_record['POS'] - ref_flank)
    bp_l_upper = sv_record['POS']

    bp_r_lower = sv_record['END']
    bp_r_upper = sv_record['END'] + ref_flank

    chr = sv_record['#CHROM']

    # Get insert sizes over the variant
    insert_size_list = list()

    for record in bam_file.fetch(str(chr), bp_l_lower, bp_l_upper):
        if record.next_reference_name == chr and bp_r_lower <= record.next_reference_start <= bp_r_upper:
            insert_size_list.append(abs(record.template_length))

    return insert_size_list


def get_insert_size_distribution(
        bam_file_name, sample_size=5e6, size_min=100, size_limit=10000, ref_filename=None
):
    """
    Get the expected insert size distribution.

    :param bam_file_name: Name of BAM file.
    :param sample_size: Number of records to retrieve.
    :param size_min: Lower limit on the insert size for outlier removal.
    :param size_limit: Upper limit on the insert size for outlier removal.
    :param ref_filename: Reference file name. Required for CRAM files.

    :return: A Series with elements "N", "MEAN", and "STDEV" describing the insert size sample.
    """

    n_records = 0
    size_list = list()

    with pysam.AlignmentFile(bam_file_name, 'r', reference_filename=ref_filename) as in_file:
        for record in in_file:

            if record.is_proper_pair and \
                    not record.is_supplementary and \
                    not record.is_secondary and \
                    size_min < abs(record.template_length) < size_limit:

                n_records += 1
                size_list.append(abs(record.template_length))

                if n_records >= sample_size:
                    break

    # Check for aligned reads.
    if n_records == 0:
        raise RuntimeError(
            'No primary aligned reads in {}: Cannot estimate the insert size distribution'.format(bam_file_name)
        )

    # Return
    return pd.Series((n_records, np.mean(size_list), np.std(size_list)), index=('N', 'MEAN', 'STDEV'))


# Main
if __name__ == '__main__':

    # Get arguments
    arg_parser = argparse.ArgumentParser(description='Get insert size deltas on the reference over the SV breakpoints.')

    arg_parser.add_argument('bam', help='BAM file of short read alignments.')

    arg_parser.add_argument('bed',
                            help='SV info BED file with columns "#CHROM", "POS", "END", "SVTYPE", "CONTIG", '
                                 '"CONTIG_START", and "CONTIG_END", including a header line.')

    arg_parser.add_argument('out',
                            help='Output file.')

    arg_parser.add_argument('--ref', nargs='?', default=None,
                            help='Reference for records are aligned against.')

    arg_parser.add_argument('--out_stats',
                            help='Output size distribution statistics used for z-scores.')

    arg_parser.add_argument('--mapq', type=int, default=20,
                            help='Minimum mapping quality of aligned reads.')

    arg_parser.add_argument('--ref_flank', type=int, default=5e3,
                            help='Search a window this far upstream of the breakpoints for a read with a paired end '
                                 'that maps this far downstream of the breakpoints.')

    arg_parser.add_argument('--force', '-f', action='store_true',
                            help='Force overwrite the output file if it exists.')

    arg_parser.add_argument('--sample_size', type=int, default=5e6,
                            help='Sample this many reads while estimating the insert size distribution.')

    arg_parser.add_argument('--size_limit', type=int, default=1500,
                            help='Filter insert sizes of this length or greater. Used for outlier control.')

    arg_parser.add_argument('--size_min', type=int, default=100,
                            help='Filter insert sizes less than this length. Used for outlier control.')

    arg_parser.add_argument('--z_threshold', type=float, default=1.5,
                            help='Z-score applied to insert sizes for an SV to determine if it may support an '
                                 'insertion or a deletion.')

    args = arg_parser.parse_args()

    # Check arguments
    if not os.path.isfile(args.bam):
        raise RuntimeError('Input BAM file does not exist or is not a regular file: {}'.format(args.bam))

    if args.mapq < 0:
        raise RuntimeError('Mapping quality is negative: {}'.format(args.mapq))

    if args.sample_size < 1:
        raise RuntimeError('Sample size is not positive: {}'.format(args.sample_size))

    if args.size_limit < 1:
        raise RuntimeError('Insert size limit is not positive: {}'.format(args.size_limit))

    if args.out is not None:
        args.out = args.out.strip()

        if not args.out:
            raise RuntimeError('Output file name is empty.')

        if not args.force and os.path.exists(args.out):
            raise RuntimeError('Output file exists and --force was not set: {}'.format(args.out))

    z_limit = abs(args.z_threshold)

    if z_limit < 0.000001:  # Check for 0 or effectively 0
        raise RuntimeError('Z-score limit is 0.')

    # Get variant info
    df_bed = pd.read_table(args.bed, header=0)

    # Get insert size distribution
    insert_stats = get_insert_size_distribution(
        args.bam, args.sample_size, args.size_min, args.size_limit, ref_filename=args.ref
    )

    insert_mean = np.float(insert_stats['MEAN'])
    insert_sd = np.float(insert_stats['STDEV'])

    if args.out_stats:
        insert_stats.to_csv(args.out_stats, sep='\t')

    # Open files and process
    bam_file_in = None

    with open(args.out, 'w') as out_file:
        out_file.write('INDEX\tN_INSERT\tINSERT_LOWER\tINSERT_UPPER\n')

        # Iterate over SV calls
        for index in range(df_bed.shape[0]):

            # Close and re-open pysam at set intervals (memory-leak work-around)
            if index % PYSAM_RESET_INTERVAL == 0:
                if bam_file_in is not None:
                    bam_file_in.close()

                gc.collect()

                bam_file_in = pysam.AlignmentFile(args.bam, 'r', reference_filename=args.ref)

            # Get record
            sv_rec = df_bed.iloc[index]

            # Get insert size stats
            insert_array = get_insert_size_list(
                sv_rec, bam_file_in, args.ref_flank
            )

            sv_rec['N_INSERT'] = len(insert_array)

            if sv_rec['N_INSERT'] > 0:
                insert_array = np.array(insert_array)  # To array
                insert_array = (insert_array - insert_mean) / insert_sd  # Z-score

                sv_rec['INSERT_LOWER'] = sum(insert_array < -z_limit) / sv_rec['N_INSERT']
                sv_rec['INSERT_UPPER'] = sum(insert_array > z_limit) / sv_rec['N_INSERT']

            else:
                sv_rec['INSERT_LOWER'] = 0
                sv_rec['INSERT_UPPER'] = 0

            # Write
            out_file.write(
                '{INDEX}\t{N_INSERT}\t{INSERT_LOWER}\t{INSERT_UPPER}\n'.format(**sv_rec)
            )
