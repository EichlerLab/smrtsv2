#!/usr/bin/env python3

"""
Get read depth over breakpoints.

An SV has 3 breakpoints. Two are on one contig and one is on another.
* DEL: 2 breakpoints on reference, 1 on the local assembly.
* INS: 1 breakpoint on reference, 2 on the local assembly.

For the purposes of this script, the two breakpoints on the same contig are bp_l (left) and bp_r (right), and the
breakpoint on a single contig is bp_s (single).

Breakpoints depths over all 3 are calculated by finding all reads that overlap each breakpoint, choosing the best
alignment for each read (CIGAR string and edit distance), and normalizing over bp_l and bp_r ((bp_l + bp_r) / 2).
"""

import argparse
import gc
import collections
import os
import pandas as pd
import pysam

from scipy.stats import binom

# CIGAR operations
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8

# Set of CIGAR operations that indicate mis-matched alignments (ignoring operations bwa would not write).
CIGAR_NO_ALIGN = {BAM_CINS, BAM_CDEL, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CDIFF}

CIGAR_CLIPPED = {BAM_CSOFT_CLIP, BAM_CHARD_CLIP}

# Number of SVs to process before resetting pysam (close and re-open file). Avoids a memory leak in pysam.
PYSAM_RESET_INTERVAL = 1000


class AlignRecord:
    """
    One alignment record with required fields extracted from the pysam segment.

    Fields:
    * qname: Query name.
    * pos: Position.
    * mapq: Mapping quality.
    * align_distance: Distance calculated by summing CIGAR string entries of non-matching bases.
    * edit_distance: Number of mismatched bases in the aligned regions of the sequence.
    * read_count: Number of reads this record represents.
    """

    def __init__(self, segment):
        """
        Create a new alignment record.

        :param segment: Pysam segment (one BAM entry).
        """

        # Get basic information
        self.qname = segment.query_name

        self.pos = segment.reference_start

        self.mapq = segment.mapping_quality

        # Count number of clipped bases on each end
        self.clip_l = 0
        self.clip_r = 0

        index = 0

        while segment.cigartuples[index][0] in CIGAR_CLIPPED:
            self.clip_l += segment.cigartuples[index][1]
            index += 1

        index = len(segment.cigartuples) - 1

        while segment.cigartuples[index][0] in CIGAR_CLIPPED:
            self.clip_r += segment.cigartuples[index][1]
            index -= 1

        self.clip_lr = self.clip_l + self.clip_r

        # Get CIGAR distance
        self.align_distance = 0

        for cigar_entry in segment.cigartuples:
            if cigar_entry[0] in CIGAR_NO_ALIGN:
                self.align_distance += cigar_entry[1]

        # Get edit distance (NM tag)
        try:
            self.edit_distance = segment.get_tag('NM')

        except KeyError:
            self.edit_distance = None

        # Read count: Number of breakpoints this record represents. If the same read is mapped over two ends of
        # an event, count it twice initially so that it normalizes properly.
        self.read_count = 1

    def compare(self, other):
        """
        Compare the distance from this record to the reference to the distance of `other` from the reference. The
        mapping quality is analyzed first and returned if they are not the same. The CIGAR operations are analyzed
        second, and if the number of bases differing is equal, then the NM tag is analyzed (if it exists
        in both records).

        The difference between the mapping quality, CIGAR distance, or the edit distance (see above) is returned
        as `self - other` for the first comparison that is non-zero.

        :param other: Other record to compare.

        :return: A negative number if this record is closer to the reference, a positive number if the other record
            is closer to the reference, and 0 if both records align equally well to the reference by CIGAR and NM tag
            edit distance.
        """

        # Compare mapping quality
        record_distance = self.mapq - other.mapq

        if record_distance != 0:
            return record_distance

        # Compare CIGAR difference
        record_distance = self.align_distance - other.align_distance

        if record_distance != 0:
            return record_distance

        # Compare edit distance (NM tag) if set
        if self.edit_distance is not None and other.edit_distance is not None:
            record_distance = self.edit_distance - other.edit_distance

        return record_distance

    def __repr__(self):
        """
        Get a string representation of this object.

        :return: String representation of this object.
        """
        return 'AlignRecord[qname={qname}, pos={pos}, dist_cigar={align_distance}, dist_nm={edit_distance}]'.format(
            **self.__dict__
        )


def get_best_alignment(record_list):
    """
    Get the best alignment from a list of alignments. The best alignment has the lowest distance to the reference. If
    more than one alignment has the minimum distance, then the first one in the list is returned.

    :param record_list: List of records.

    :return: Best record in `record_list`.
    """

    best_alignment = record_list[0]

    for alignment in record_list[1:]:
        if alignment.compare(best_alignment) < 0:
            best_alignment = alignment

    return best_alignment


def merge_breakpoint_records(records_l, records_r):
    """
    Merge records found on both the left and the right breakpoints. If a read or its paired end overlaps both
    breakpoints, then choose the lesser of the two alignments to represent the read (lesser by furthest distance
    from the reference or contig it belonged to).

    :param records_l: Left breakpoint records.
    :param records_r: Right breakpoint records.

    :return: Merged left and right breakpoints.
    """

    records_lr = dict()

    # Get records in each set
    qname_l = set(records_l.keys())
    qname_r = set(records_r.keys())

    qname_lr = qname_l.intersection(qname_r)

    # Merge records with the same name from the left and right breakpoints. Choose the lesser of the two alignments.
    for record_name in qname_lr:
        if records_l[record_name].compare(records_r[record_name]) >= 0:
            records_lr[record_name] = records_l[record_name]
        else:
            records_lr[record_name] = records_r[record_name]

        records_lr[record_name].read_count += 1  # Count both records

    # Add records unique to each set
    for record_name in (qname_l - qname_lr):
        records_lr[record_name] = records_l[record_name]

    for record_name in (qname_r - qname_lr):
        records_lr[record_name] = records_r[record_name]

    # Return the merged set
    return records_lr


def choose_best_records(records_lr, records_s):
    """
    Get the best records for each set. For all records in both sets, the best is chosen. All records unique to either
    set is also selected and returned.

    :param records_lr: Records across the left and right breakpoints.
    :param records_s: Records acroos the single breakpoint.

    :return: A tuple of two elements: best records in left-right, and best records in single.
    """

    # Init lists
    best_lr = list()
    best_s = list()

    # Get sets of record names
    qname_lr = set(records_lr.keys())
    qname_s = set(records_s.keys())

    qname_all = qname_lr.intersection(qname_s)

    # Iterate over shared alignments
    for record_name in qname_all:
        diff = records_lr[record_name].compare(records_s[record_name])

        if diff < 0:
            best_lr.append(records_lr[record_name])
        elif diff > 0:
            best_s.append(records_s[record_name])

        # Records in both sets with an equal quality are discarded

    # Add records unique to each
    for record_name in qname_lr - qname_all:
        best_lr.append(records_lr[record_name])

    for record_name in qname_s - qname_all:
        best_s.append(records_s[record_name])

    # Return lists of records
    return best_lr, best_s


def get_clipped_freq(records_l, records_r, records_s, minclip=4):
    """
    Get the number of reads that are clipped over the breakpoints in a meaningful way. For the left and right
    breakpoints, look for clipping in the SV. For the single breakpoint, look for clipping on either side of the read.

    :param records_l: Records over the left breakpoint.
    :param records_r: Records over the right breakpoint.
    :param records_s: Records over the single breakpoint.
    :param minclip: Minimum number of bases clipped to count the read as clipped.

    :return: A tuple of two elements: Proportion of reads clipped over the left and right breakpoints, and proportion
        of reads clipped over the single breakpoint (in that order).
    """

    if len(records_l) > 0 or len(records_r) > 0:
        clipped_lr = (
            len([True for key in records_l if records_l[key].clip_r >= minclip]) +
            len([True for key in records_r if records_r[key].clip_l >= minclip])
        ) / (len(records_l) + len(records_r))
    else:
        clipped_lr = 0

    if len(records_s) > 0:
        clipped_s = (
            len([True for key in records_s if records_s[key].clip_l >= minclip or records_s[key].clip_r >= minclip]) /
            len(records_s)
        )
    else:
        clipped_s = 0

    return clipped_lr, clipped_s


def get_bp_depth(sv_record, bam_file, minclip=4, mapq=20):
    """
    Get the average depth over reference and alternate contigs.

    :param sv_record: SV record from the input BED file.
    :param bam_file: Pysam opened BAM file.
    :param minclip: The minimum number of bases clipped on one end to count that end of the read as clipped.
    :param mapq: Minimum mapping quality of reads to be counted.

    :return: A 5-element tuple: Depth over the reference, depth over the alternate contig, number of clipped reads over
        over the reference reads counted toward the read depth, number of clipped reads over the alternate contig reads
        counted toward the read depth, and proportion of reads clipped on the primary contig.
    """

    # Get breakpoint locations
    # * bp_l: Left breakpoint (contig if INS, primary if DEL)
    # * bp_r: Right breakpoint (contig if INS, primary if DEL)
    # * bp_s: Single breakpoing (primary if INS (point of insertion), contig if DEL)
    if sv_record['SVTYPE'] == 'INS':

        is_ins = True

        bp_lr_chr = sv_record['CONTIG']
        bp_l = sv_record['CONTIG_START']
        bp_r = sv_record['CONTIG_END']

        bp_s_chr = sv_record['#CHROM']
        bp_s = sv_record['POS']

    elif sv_record['SVTYPE'] == 'DEL':

        is_ins = False

        bp_lr_chr = sv_record['#CHROM']
        bp_l = sv_record['POS']
        bp_r = sv_record['END']

        bp_s_chr = sv_record['CONTIG']
        bp_s = sv_record['CONTIG_START']

    else:
        raise RuntimeError(
            'Cannot process unknown variant type {} (must be "INS" or "DEL")'.format(sv_record['SVTYPE'])
        )

    # Get records over each breakpoint
    records_l = collections.defaultdict(list)
    records_r = collections.defaultdict(list)
    records_s = collections.defaultdict(list)

    for record in bam_file.fetch(str(bp_lr_chr), bp_l, bp_l + 1):
        if record.mapping_quality >= mapq and record.is_proper_pair:
            records_l[record.query_name].append(AlignRecord(record))

    for record in bam_file.fetch(str(bp_lr_chr), bp_r, bp_r + 1):
        if record.mapping_quality >= mapq and record.is_proper_pair:
            records_r[record.query_name].append(AlignRecord(record))

    for record in bam_file.fetch(str(bp_s_chr), bp_s, bp_s + 1):
        if record.mapping_quality >= mapq and record.is_proper_pair:
            records_s[record.query_name].append(AlignRecord(record))

    # Get best representation for each read over each breakpoint
    for key in records_l.keys():
        records_l[key] = get_best_alignment(records_l[key])

    for key in records_r.keys():
        records_r[key] = get_best_alignment(records_r[key])

    for key in records_s.keys():
        records_s[key] = get_best_alignment(records_s[key])

    # Get clipped count
    ref_clipped_lr, ref_clipped_s = get_clipped_freq(records_l, records_r, records_s, minclip)

    # Merge records in l and r breakpoints.
    records_lr = merge_breakpoint_records(records_l, records_r)

    # Get the best records after comparing alignments over the breakpoints
    records_lr, records_s = choose_best_records(records_lr, records_s)

    clip_lr = sum([record.read_count for record in records_lr if record.clip_lr >= minclip]) / 2
    clip_s = float(sum([record.read_count for record in records_s if record.clip_lr >= minclip]))

    # Return the depths over the breakpoints
    if is_ins:
        return (
            float(sum(record.read_count for record in records_s)),
            sum([record.read_count for record in records_lr]) / 2,
            clip_s,
            clip_lr,
            ref_clipped_s
        )

    else:
        return (
            sum(record.read_count for record in records_lr) / 2,
            float(sum([record.read_count for record in records_s])),
            clip_lr,
            clip_s,
            ref_clipped_lr
        )


# Main
if __name__ == '__main__':

    # Get arguments
    arg_parser = argparse.ArgumentParser(description='Calculate read depths over breakpoints for structural variants.')

    arg_parser.add_argument('bam', help='BAM file of short read alignments.')

    arg_parser.add_argument('bed',
                            help='SV info BED file with columns "#CHROM", "POS", "END", "SVTYPE", "CONTIG", '
                                 '"CONTIG_START", and "CONTIG_END", including a header line.')

    arg_parser.add_argument('out',
                            help='Output file.')

    arg_parser.add_argument('--mapq', type=int, default=20,
                            help='Minimum mapping quality of aligned reads.')

    arg_parser.add_argument('--minclip', type=int, default=4,
                            help='Minimum number of clipped bases on one end to count a read as clipped on that end.')

    arg_parser.add_argument('--force', '-f', action='store_true',
                            help='Force overwrite the output file if it exists.')

    arg_parser.add_argument('--phom', type=float, default=0.05,
                            help='Probability of a read aligning better to the alternate breakpoint(s) assuming the '
                                 'genotype is homozygous-reference (no actual alt allele present). The probability of '
                                 'homozygous-alternate is calculated by subtracting this value from 1.')

    arg_parser.add_argument('--phet', type=float, default=0.5,
                            help='Probability of a read aligning better to the alternate breakpoint(s) assuming the '
                                 'genotype is heterozygous.')

    args = arg_parser.parse_args()

    # Check arguments
    if not os.path.isfile(args.bam):
        raise RuntimeError('Input BAM file does not exist or is not a regular file: {}'.format(args.bam))

    if args.mapq < 0:
        raise RuntimeError('Mapping quality is negative: {}'.format(args.mapq))

    if args.minclip < 1:
        raise RuntimeError('Minimum clipping threshold must be at least 1: {}'.format(args.minclip))

    if args.out is not None:
        args.out = args.out.strip()

        if not args.out:
            raise RuntimeError('Output file name is empty.')

        if not args.force and os.path.exists(args.out):
            raise RuntimeError('Output file exists and --force was not set: {}'.format(args.out))

    if args.phet > 1 or args.phet < 0:
        raise RuntimeError('ALT-probability given HET must be in [0, 1]: {}'.format(args.phet))

    if args.phom > 1 or args.phom < 0:
        raise RuntimeError('ALT-probability given HOM-REF must be in [0, 1]: {}'.format(args.phom))

    # Save arguments that are reused often
    phet = args.phet
    phomref = args.phom
    phomalt = 1 - phomref

    # Get variant info
    df_bed = pd.read_table(args.bed, header=0)

    # Open files and process
    bam_file_in = None  # Init within df_bed iterator loop

    with open(args.out, 'w') as out_file:

        # Write header
        out_file.write(
            'INDEX\tBP_REF_COUNT\tBP_ALT_COUNT\t'
            'BP_REF_REL\tBP_ALT_REL\t'
            'BP_PROB_HOMREF\tBP_PROB_HET\tBP_PROB_HOMALT\t'
            'BP_REF_CLIP\tBP_ALT_CLIP\tBP_PRIMARY_CLIP\n'
        )

        # Iterate over SV calls
        for index in range(df_bed.shape[0]):

            # Close and re-open pysam at set intervals (memory-leak work-around)
            if index % PYSAM_RESET_INTERVAL == 0:
                if bam_file_in is not None:
                    bam_file_in.close()

                gc.collect()

                bam_file_in = pysam.AlignmentFile(args.bam, 'r')

            # Get record
            sv_rec = df_bed.iloc[index]

            # Get counts
            ref_count, alt_count, ref_clip, alt_clip, primary_clip = \
                get_bp_depth(sv_rec, bam_file_in, args.minclip, args.mapq)

            # Normalize to ratios of depth and get binomial probability of HOM-REF and HET given read depth
            total_depth = ref_count + alt_count

            if total_depth > 0:
                ref_rel = ref_count / total_depth
                alt_rel = alt_count / total_depth

                prob_homref = binom(total_depth, phomref).pmf(alt_count)
                prob_het = binom(total_depth, phet).pmf(alt_count)
                prob_homalt = binom(total_depth, phomalt).pmf(alt_count)

            else:
                ref_rel = 0
                alt_rel = 0

                prob_homref = 0
                prob_het = 0
                prob_homalt = 0

            # Normalize clipping counts
            if ref_count > 0:
                ref_clip /= ref_count
            else:
                ref_clip = 0

            if alt_count > 0:
                alt_clip /= alt_count
            else:
                alt_clip = 0

            # Update Series
            sv_rec['BP_REF_COUNT'] = ref_count
            sv_rec['BP_ALT_COUNT'] = alt_count
            sv_rec['BP_REF_REL'] = ref_rel
            sv_rec['BP_ALT_REL'] = alt_rel
            sv_rec['BP_PROB_HOMREF'] = prob_homref
            sv_rec['BP_PROB_HET'] = prob_het
            sv_rec['BP_PROB_HOMALT'] = prob_homalt
            sv_rec['BP_REF_CLIP'] = ref_clip
            sv_rec['BP_ALT_CLIP'] = alt_clip
            sv_rec['BP_PRIMARY_CLIP'] = primary_clip

            # Write
            out_file.write(
                '{INDEX}\t{BP_REF_COUNT}\t{BP_ALT_COUNT}\t'
                '{BP_REF_REL}\t{BP_ALT_REL}\t'
                '{BP_PROB_HOMREF}\t{BP_PROB_HET}\t{BP_PROB_HOMALT}\t'
                '{BP_REF_CLIP}\t{BP_ALT_CLIP}\t{BP_PRIMARY_CLIP}\n'.format(**sv_rec)
            )

    # Gracefully close input file
    if bam_file_in is not None:
        bam_file_in.close()
