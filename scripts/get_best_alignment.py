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


def calculate_differences_in_read(read):
    # Calculate indel and softclipped bases from the sum of lengths for CIGAR
    # options that aren't matches.
    indels_and_softclips = sum([cigar[1] for cigar in read.cigar if cigar[0] != BAM_CMATCH])

    # Calculate mismatches from the NM tags.
    mismatches = read.get_tag("NM")

    differences = indels_and_softclips + mismatches

    logger.debug("Found %i differences in read %s", differences, read)
    return differences


def get_best_reads(alignments_by_read_name, filtered_reads):
    """
    Example duplicates:

    HJ2HJCCXX160108:5:2220:13311:66496      83      40
    HJ2HJCCXX160108:5:2220:13311:66496      2131    40
    HJ2HJCCXX160108:5:2220:13311:66496      163     60

    Should be reduced to:

    HJ2HJCCXX160108:5:2220:13311:66496      83      40
    HJ2HJCCXX160108:5:2220:13311:66496      163     60
    """
    for name, reads in alignments_by_read_name.iteritems():
        if len(reads) > 1:
            # Sort reads by total differences from reference in ascending order.
            reads = sorted(reads, key=lambda read: calculate_differences_in_read(read))

        # Either print the only read or print the best read.
        filtered_reads.write(reads[0])


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

    for read in reads:
        if read.query_name != current_read_name:
            if current_read_name is not None:
                get_best_reads(alignments_by_read_name, filtered_reads)

            current_read_name = read.query_name
            alignments_by_read_name = defaultdict(list)

        # Index reads by query name and read number to avoid treating reads 1
        # and 2 of the same pair as duplicates.
        read_name = "%s/%s" % (read.query_name, read.is_read1 and 1 or 2)
        alignments_by_read_name[read_name].append(read)
    else:
        # Get best reads from the last set in the for loop.
        get_best_reads(alignments_by_read_name, filtered_reads)
