import argparse
from collections import defaultdict
import logging
import pysam

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

BAM_CMATCH = 0


def calculate_differences_in_alignment(alignment):
    # Calculate indel and softclipped bases from the sum of lengths for CIGAR
    # options that aren't matches.
    indels_and_softclips = sum([cigar[1] for cigar in alignment.cigar if cigar[0] != BAM_CMATCH])

    # Calculate mismatches from the NM tags.
    mismatches = alignment.get_tag("NM")

    differences = indels_and_softclips + mismatches

    logger.debug("Found %i differences in alignment %s", differences, alignment)
    return differences


def get_best_alignment(alignments_by_read_name, filtered_reads):
    """
    For a given dictionary of lists of alignments indexed by read name and
    number (e.g., "read_name/1"), find the alignment with the fewest differences
    (mismatches, indels, and softclips). If all alignments for a given read have the
    same number of differences, don't return any alignment for that read since its
    mapping quality would have been 0 by traditional BWA MEM alignment rules.

    Example duplicates:

    HJ2HJCCXX160108:5:2220:13311:66496      83      40
    HJ2HJCCXX160108:5:2220:13311:66496      2131    40
    HJ2HJCCXX160108:5:2220:13311:66496      163     60

    Should be reduced to:

    HJ2HJCCXX160108:5:2220:13311:66496      83      40
    HJ2HJCCXX160108:5:2220:13311:66496      163     60
    """
    for name, alignments in alignments_by_read_name.iteritems():
        if len(alignments) > 1:
            # Calculate differences between read and the reference for each
            # alignment.
            alignments_by_differences = defaultdict(list)
            for alignment in alignments:
                differences = calculate_differences_in_alignment(alignment)
                alignments_by_differences[differences].append(alignment)

            # Find the alignments(s) with the fewest differences.
            min_differences = min(alignments_by_differences.keys())

            # If there is more than one alignment with the fewest differences,
            # don't emit any alignments to be consistent with a mapping quality
            # of 0. If there is only one best alignment, emit that.
            if len(alignments_by_differences[min_differences]) == 1:
                filtered_reads.write(alignments_by_differences[min_differences][0])
                logger.debug("Found best alignment for %s", name)
            else:
                logger.debug("Skipping equally scored alignments for %s", name)
        else:
            # Emit the only read in the set.
            filtered_reads.write(alignments[0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("reads")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)

    reads = pysam.AlignmentFile(args.reads, "rb")
    filtered_reads = pysam.AlignmentFile("-", "wbu", template=reads)

    current_read_name = None
    alignments_by_read_name = defaultdict(list)

    reads_processed = 0
    for read in reads:
        reads_processed += 1
        if read.query_name != current_read_name:
            if current_read_name is not None:
                get_best_alignment(alignments_by_read_name, filtered_reads)

            current_read_name = read.query_name
            alignments_by_read_name = defaultdict(list)

        # Index reads by query name and read number to avoid treating reads 1
        # and 2 of the same pair as duplicates.
        read_name = "%s/%s" % (read.query_name, read.is_read1 and 1 or 2)
        alignments_by_read_name[read_name].append(read)
    else:
        # Get best reads from the last set in the for loop.
        get_best_alignment(alignments_by_read_name, filtered_reads)

    logger.debug("Processed %i reads for filtering", reads_processed)
