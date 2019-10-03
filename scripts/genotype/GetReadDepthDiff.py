#!/usr/bin/env python3

import argparse
import gc
import numpy as np
import os
import pandas as pd
import pysam


# Number of SVs to process before resetting pysam (close and re-open file). Avoids a memory leak in pysam.
PYSAM_RESET_INTERVAL = 1000


def get_read_depth(df_subset, bam_file_name, mapq, ref_filename=None):
    """
    Get read depths over one or more breakpoints.

    :param df_subset: Subset dataframe with a column for contigs (first column) and one or more columns for the
        location of breakpoints to quantify.
    :param bam_file_name: Name of alignment file to query.
    :param mapq: Minimum mapping quality.

    :return: A Series with with one element for each row of `df_subset` containing the average of read depths over
        the breakpoints for each variant.
    """

    # Init pysam query count (for memory leak prevention)
    pysam_count = 0
    bam_file = pysam.AlignmentFile(bam_file_name, 'r', reference_filename=ref_filename)

    # Init dataframe
    df_subset = df_subset.copy()

    n_loc_cols = df_subset.shape[1] - 1  # Number of location columns; depth is averaged for each

    df_subset.columns = ['CONTIG'] + ['LOC_{}'.format(col) for col in range(n_loc_cols)]

    # Init count
    df_subset['N'] = np.zeros(df_subset.shape[0], np.float64)
    n_index = df_subset.shape[1] - 1

    # Count
    for subset_index in range(n_loc_cols):

        # Use numeric index, skip chromosome column
        subset_index += 1

        for row_index in range(df_subset.shape[0]):

            n_reads = 0

            # Get position
            contig = df_subset.iloc[row_index, 0]
            pos = df_subset.iloc[row_index, subset_index]

            # Reset pysam periodically (avoids memory leak)
            pysam_count += 1

            if pysam_count >= PYSAM_RESET_INTERVAL:
                if bam_file is not None:
                    bam_file.close()

                gc.collect()

                bam_file = pysam.AlignmentFile(bam_file_name, 'r', reference_filename=ref_filename)

                pysam_count = 0

            # Count
            for segment in bam_file.fetch(str(contig), pos, pos + 1):
                if segment.mapping_quality >= mapq and segment.is_proper_pair:
                    n_reads += 1

            df_subset.iloc[row_index, n_index] += n_reads

    # Return mean of depths (divide by the number of locations)
    return df_subset['N'] / n_loc_cols


def get_ref_contig_sizes(altref_file):
    """
    Get a Series of contigs lengths. Includes primary and alt contigs.

    :param altref_file: BED file of contig information where each record spans the whole contig. Must contain
        columns "#CHROM" and "END".

    :return: Series of contig lengths indexed by the contig name.
    """

    # Get reference chromosome sizes
    ref_len_series = pd.read_table(altref_file, header=0)
    ref_len_series.index = ref_len_series['#CHROM']
    ref_len_series = ref_len_series['END']

    return ref_len_series


def annotate_variant_info(variant_table, ref_len_series, flank):
    """
    Annotate variant info with locations reads will be extracted from.

    :param variant_table: Variant info table.
    :param ref_len_series: Series of contig sizes.
    :param flank: Number of bases from variant breakpoints.

    :return: `variant_table` with additional fields.
    """
    
    # Annotate variant info with flank locations
    variant_table['FLANK_L_REF'] = variant_table['POS'] - flank
    variant_table['FLANK_L_REF'] = variant_table['FLANK_L_REF'].apply(lambda pos: pos if pos > 0 else 0)

    variant_table['FLANK_R_REF'] = variant_table['END'] + flank
    variant_table['FLANK_R_REF'] = variant_table.apply(lambda row: min(row['FLANK_R_REF'], ref_len_series[row['#CHROM']]), axis=1)

    variant_table['FLANK_L_CTG'] = variant_table['CONTIG_START'] - flank
    variant_table['FLANK_L_CTG'] = variant_table['FLANK_L_CTG'].apply(lambda pos: pos if pos > 0 else 0)

    variant_table['FLANK_R_CTG'] = variant_table['CONTIG_END'] + flank
    variant_table['FLANK_R_CTG'] = variant_table.apply(lambda row: min(row['FLANK_R_CTG'], ref_len_series[row['CONTIG']]), axis=1)

    # Annotate with the midpoint of the variant sequence
    variant_table['VAR_CONTIG'] = variant_table.apply(lambda row: row['#CHROM'] if row['SVTYPE'] == 'DEL' else row['CONTIG'], axis=1)
    variant_table['VAR_MIDPOINT'] = variant_table.apply(
        lambda row:
            (row['POS'] + row['END']) / 2 if row['SVTYPE'] == 'DEL' else (row['CONTIG_START'] + row['CONTIG_END']) / 2,
        axis=1)

    variant_table['VAR_MIDPOINT'] = variant_table['VAR_MIDPOINT'].astype(np.int64)

    return variant_table


# Main
if __name__ == '__main__':

    # Get arguments
    arg_parser = argparse.ArgumentParser(description='Get insert size deltas on the reference over the SV breakpoints.')

    arg_parser.add_argument('bam', help='BAM file of short read alignments.')

    arg_parser.add_argument('bed', help='SV info BED file with columns "#CHROM", "POS", "END", "SVTYPE", "CONTIG", '
                                        '"CONTIG_START", and "CONTIG_END", including a header line.')

    arg_parser.add_argument('alt_info', help='BED file of contigs in the reference.')

    arg_parser.add_argument('out', help='Output file.')

    arg_parser.add_argument('--out_stats',
                            help='Output depth distribution statistics.')

    arg_parser.add_argument('--mapq', type=int, default=20,
                            help='Minimum mapping quality of aligned reads.')

    arg_parser.add_argument('--flank', type=int, default=100,
                            help='Number of reference bases on each side of the SV for flanking regions.')

    arg_parser.add_argument('--ref', nargs='?',
                            default=None, help='Reference for records are aligned against.')

    args = arg_parser.parse_args()

    # Check arguments
    if not os.path.isfile(args.bam):
        raise RuntimeError('Input BAM file does not exist or is not a regular file: {}'.format(args.bam))

    if args.mapq < 0:
        raise RuntimeError('Mapping quality is negative: {}'.format(args.mapq))

    if args.flank < 0:
        raise RuntimeError('Flank is negative: {}'.format(args.flank))

    args.out = args.out.strip()

    if not args.out:
        raise RuntimeError('Output file name is empty.')

    # Get variant info
    df_bed = pd.read_table(args.bed, header=0)

    # Get reference chromosome sizes
    ref_len = get_ref_contig_sizes(args.alt_info)

    # Annotate variant info with locations reads are extracted from
    df_bed = annotate_variant_info(df_bed, ref_len, args.flank)

    # Count reads over variant midpoint
    df_bed['DP_N_VAR'] =\
        get_read_depth(df_bed.loc[:, ['VAR_CONTIG', 'VAR_MIDPOINT']], args.bam, args.mapq, ref_filename=args.ref)

    # Count reads over reference flank
    df_bed['DP_N_PROX_REF'] =\
        get_read_depth(df_bed.loc[:, ['#CHROM', 'FLANK_L_REF', 'FLANK_R_REF']], args.bam, args.mapq, ref_filename=args.ref)

    # Count reads over contig flank
    df_bed['DP_N_PROX_CTG'] =\
        get_read_depth(df_bed.loc[:, ['CONTIG', 'FLANK_L_CTG', 'FLANK_R_CTG']], args.bam, args.mapq, ref_filename=args.ref)

    # Get global stats
    ref_mean = np.mean(df_bed['DP_N_PROX_REF'])
    ref_sd = np.std(df_bed['DP_N_PROX_REF'])

    if ref_mean == 0:
        raise RuntimeError('Cannot compute global depth stats: Global mean of proximal reference breakpoint depths is 0')

    # Combine total depths
    df_bed['DP_N_VAR_PROX_REF'] = df_bed['DP_N_VAR'] + df_bed['DP_N_PROX_REF']
    df_bed['DP_N_VAR_PROX_CTG'] = df_bed['DP_N_VAR'] + df_bed['DP_N_PROX_CTG']

    # Set relative ratios
    df_bed['DP_VAR_REF'] = df_bed.apply(
        lambda row: row['DP_N_VAR'] / row['DP_N_VAR_PROX_REF'] if row['DP_N_VAR_PROX_REF'] > 0 else 0,
        axis=1
    )

    df_bed['DP_VAR_CTG'] = df_bed.apply(
        lambda row: row['DP_N_VAR'] / row['DP_N_VAR_PROX_CTG'] if row['DP_N_VAR_PROX_CTG'] > 0 else 0,
        axis=1
    )

    df_bed['DP_VAR_GLOBAL'] = df_bed['DP_N_VAR'] / ref_mean

    # Write
    df_features = df_bed.loc[
             :, ('INDEX', 'DP_VAR_REF', 'DP_VAR_CTG', 'DP_VAR_GLOBAL', 'DP_N_VAR', 'DP_N_PROX_REF', 'DP_N_PROX_CTG')
    ]

    df_features.to_csv(
        args.out, sep='\t', index=False
    )

    # Write stats
    if args.out_stats:
        with open(args.out_stats, 'w') as stats_out:
            stats_out.write('ref_mean\t{:.6f}\n'.format(ref_mean))
            stats_out.write('ref_sd\t{:.6f}\n'.format(ref_sd))
