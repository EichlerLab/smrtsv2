"""
Genotype a set of SV calls in a given set of BAMs.
"""
import argparse
from collections import defaultdict
import intervaltree
from joblib import Parallel, delayed
import logging
import numpy as np
import operator
import pybedtools
import pysam
from scipy import stats

np.random.seed(1)

# create logger
logger = logging.getLogger("calculate_depth")
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


# Input file looks like this with SV coordinates in local assembly at columns 6-8.
# chr1    350793  350794  insertion       40      chr1-350756-362784|utg7180000000000|merged      37      77
CHROMOSOME=0
START=1
END=2
EVENT_TYPE=3
EVENT_LENGTH=4
CONTIG_NAME=8
CONTIG_START=9
CONTIG_END=10

# CHROMOSOME=0
# START=1
# END=2
# EVENT_TYPE=3
# EVENT_LENGTH=4
# CONTIG_NAME=8
# CONTIG_START=9
# CONTIG_END=10

BAM_CMATCH = 0
BAM_CSOFT_CLIP = 4
DEVIATIONS_FOR_THRESHOLD = 1.5
# Allow at most 2 mismatches in a 100 bp read.
MAX_ERROR_RATE = 0.02
MAX_LIKELIHOOD = 100
MAX_DEPTH = 100
MIN_DEPTH = 5


def expand_and_merge_intervals(breakpoints, chromosome_size, slop_size):
    chromosome = breakpoints[0][0]
    intervals = intervaltree.IntervalTree()

    # Add slop to both sides of the breakpoint.
    for breakpoint in breakpoints:
        intervals.add(intervaltree.Interval(max(0, breakpoint[1] - slop_size), min(breakpoint[2] + slop_size, chromosome_size)))

    # Merge overlapping intervals.
    intervals.merge_overlaps()

    # Return a list of interval tuples.
    return [(chromosome,) + interval[:2] for interval in intervals]


def has_perfect_mapping(read):
    """
    Returns True if the given pysam read maps perfectly.

    A "perfect" mapping is either a mapping quality greater than zero and no
    more than 2% mismatches or a full-length alignment of the read without any
    mismatches. Also consider soft-clipping and other non-match alignment
    features in the CIGAR string by summing all non-match bases and comparing
    against max errors allowed.
    """
    return (not read.is_unmapped and read.mapq > 20)


def spans_region(read, region):
    """
    Returns True if the given pysam read spans the given pybedtools.Interval,
    ``region``.
    """
    return has_perfect_mapping(read) and read.reference_start <= region[1] and read.reference_end >= region[2]


def has_gaps_in_region(read, region):
    """
    Returns True if the given pysam read spans the given pybedtools.Interval,
    ``region``.
    """
    # If the given read has gaps in its alignment to the reference inside the
    # given interval (more than one block inside the SV event itself), there are
    # gaps inside the SV.
    tree = intervaltree.IntervalTree()
    for block in read.get_blocks():
        tree[block[0]:block[1]] = block

    return len(tree[region[1]:region[2]]) > 1


def pair_spans_regions(read_pair, regions):
    """
    Returns True if the given pysam reads spans the given pybedtools.Interval
    list, ``regions``. Read pairs can span the regions even if either read
    overlaps the region breakpoints.
    """
    return len(read_pair) == 2 and read_pair[0].reference_start < regions[0][1] and read_pair[1].reference_end > regions[-1][2]


def pair_has_proper_orientation(read_pair):
    """
    Returns True if the given pysam reads map towards each other (i.e., read
    with lowest coordinate is in forward orientation and other read is in
    reverse).
    """
    return len(read_pair) == 2 and not read_pair[0].is_reverse and read_pair[1].is_reverse


def soft_clips_at_breakpoint(read, region):
    """
    Returns True if the given pysam read maps up to the edge of the given
    pybedtools.Interval, ``region`` and soft clips at that edge.

    Two cases are possible:

    1. The end of the read soft clips such that the read's reference end stops
    at the breakpoint start and the end of the read is softclipped.

    2. The beginning of the read soft clips such that the read's reference start
    begins at the breakpoint end and the beginning of the read is softclipped.

    In both cases, the read should be mapped in its best location in the
    reference ("perfect mapping").
    """
    return (
        (read.reference_end == region[1] and read.cigar[-1][0] == BAM_CSOFT_CLIP) or
        (read.reference_start == region[2] + 1 and read.cigar[0][0] == BAM_CSOFT_CLIP)
    )


def maps_outside_regions(read, regions):
    """
    Returns True if the given pysam read maps outside the range of the given
    regions.
    """
    return (
        (read.reference_start < regions[0][1] and read.reference_end < regions[0][1]) or
        (read.reference_start > regions[-1][2] and read.reference_end > regions[-1][2])
    )


def is_proper_pair(read_pair, lower_insert_size_threshold, upper_insert_size_threshold):
    """
    Returns True if both reads in the given pair map perfectly, map within the
    expected insert size, and map with the proper orientation.
    """
    return (
        lower_insert_size_threshold <= np.abs(read_pair[0].isize) <= upper_insert_size_threshold and
        pair_has_proper_orientation(read_pair) and
        any(map(has_perfect_mapping, read_pair))
    )


def get_insert_sizes_for_region(bam, region):
    """
    Return insert sizes reads in the given pysam.AlignmentFile, ``bam``, and the
    given pybedtools.Interval, ``region``.
    """
    return [read.tlen for read in bam.fetch(region[0], int(region[1]), int(region[2]))
            if has_perfect_mapping(read) and not read.mate_is_unmapped and read.tlen >= 0 and read.tlen <= 1000]


def get_depth_for_region(bam_fields, regions, breakpoints, region_type="control"):
    """
    Return read depth for concordant reads in the given dictionary
    ``bam_fields`` with a pysam.AlignmentFile in the field "bam" and the given
    pybedtools.BedTool of regions, ``regions``.

    The depth calculated across the given regions depends on the type of region
    which is "control" by default. Other acceptable values for the region type
    are "insertion" or "deletion". If the region type is not "control" then the
    coordinates given by ``breakpoints`` are used to evaluate concordant and
    discordant read pairs.
    """
    bam = bam_fields["file"]
    lower_insert_threshold = bam_fields["lower_insert_threshold"]
    upper_insert_threshold = bam_fields["upper_insert_threshold"]

    # Convert breakpoints to a list for easier indexing.
    breakpoints = list(breakpoints)

    # If the event is an insertion, calculate the midpoint of the insertion
    # using the given breakpoints.
    if region_type == "insertion":
        chromosome = breakpoints[0][0]
        # Midpoint is ((end of second breakpoint) - (start of first breakpoint) / 2) + start of first breakpoint
        midpoint = int((breakpoints[1][2] - breakpoints[0][1]) / 2.0) + breakpoints[0][1]
        breakpoints = ((chromosome, midpoint, midpoint + 1),)

    # Get all distinct reads in given regions.
    reads = set()
    for region in regions:
        for read in bam.fetch(region[0], region[1], region[2]):
            # Exclude secondary and supplementary alignments (flags 0x100 and 0x800).
            if not read.is_secondary and (read.flag & 0x800 == 0):
                reads.add(read)

    # Group reads by read name.
    reads_by_name = defaultdict(list)
    for read in reads:
        reads_by_name[read.qname].append(read)

    logger.debug("Found %i potential read pairs for %s", len(reads_by_name), str(region).strip())

    concordant_pairs = []

    for read_name, read_pair in reads_by_name.iteritems():
        # Sort reads by position to get reads "1" and "2" in expected order.
        read_pair = sorted(read_pair, key=operator.attrgetter("pos"))

        is_concordant = False
        if region_type == "insertion":
            # Reads are in a proper pair and either of the reads maps across
            # the insertion midpoint.
            if is_proper_pair(read_pair, lower_insert_threshold, upper_insert_threshold) and any([spans_region(read, breakpoint) for read in read_pair for breakpoint in breakpoints]):
                is_concordant = True
        elif region_type == "deletion":
            # Reads are in a proper pair and either of the reads maps across the
            # deletion breakpoint.
            if is_proper_pair(read_pair, lower_insert_threshold, upper_insert_threshold) and any([spans_region(read, breakpoint) for read in read_pair for breakpoint in breakpoints]):
                is_concordant = True
        else:
            # Control regions look for all properly paired reads with proper
            # orientation.
            if is_proper_pair(read_pair, lower_insert_threshold, upper_insert_threshold):
                is_concordant = True

        if is_concordant:
            concordant_pairs.append(read_pair)

    return (concordant_pairs, len(reads_by_name.keys()))


def get_min_depth_for_region(bam_fields, breakpoints, region_type="control"):
    """
    Return read depth for concordant reads in the given dictionary
    ``bam_fields`` with a pysam.AlignmentFile in the field "bam" and the given
    pybedtools.BedTool of regions, ``regions``.

    The depth calculated across the given regions depends on the type of region
    which is "control" by default. Other acceptable values for the region type
    are "insertion" or "deletion". If the region type is not "control" then the
    coordinates given by ``breakpoints`` are used to evaluate concordant and
    discordant read pairs.
    """
    # Convert breakpoints to a list for easier indexing.
    region = "%s:%s-%s" % list(breakpoints)[0]
    region_size = breakpoints[0][2] - breakpoints[0][1]
    logger.debug("Min depth for region %s (size: %s)", region, region_size)
    depths = [int(line.rstrip().split("\t")[-1])
              for line in pysam.depth("-Q", "20", "-q", "20", "-r", region, bam_fields["filename"])]

    # If there are not depths for every base of the given region, fill those
    # missing values with enough additional zeroes to equal the region size.
    if len(depths) < region_size:
        depths = depths + [0] * (region_size - len(depths))

    logger.debug("Found %i depths for region %s: %s", len(depths), region, depths)
    logger.debug("Median depths: %s", np.median(depths))
    logger.debug("Standard deviation depths: %s", np.std(depths))
    logger.debug("Standard error depths: %s", stats.sem(depths))
    high_quality_depth = max(int(np.round(np.median(depths) - np.std(depths))), 0)
    # median_total_depth = int(np.median([int(line.rstrip().split("\t")[-1])
    #                                     for line in pysam.depth("-r", region, bam_fields["filename"])]))

    return high_quality_depth


def genotype_call_with_read_pair(concordant, discordant, std_depth):
    """
    Genotype call based on total depth of concordant and discordant reads a la
    Hormozdiari et al. 2010 (Genome Research).
    """
    #homozygous_deletion_threshold = std_depth
    homozygous_deletion_threshold = 5

    if concordant < homozygous_deletion_threshold and discordant < homozygous_deletion_threshold:
        genotype = "./."
        discordant_genotype_likelihood = np.power(2, discordant) / np.power(2, homozygous_deletion_threshold)
        concordant_genotype_likelihood = np.power(2, concordant) / np.power(2, homozygous_deletion_threshold)
        shared_likelihood = np.sqrt(np.square(concordant_genotype_likelihood) + np.square(discordant_genotype_likelihood))
        genotype_likelihood = np.floor(-10 * np.log10(shared_likelihood))
    else:
        expected_discordant_lower_bound = concordant * 0.25
        expected_discordant_upper_bound = concordant * 4

        if discordant < expected_discordant_lower_bound:
            genotype = "1/1"
            # Calculate likelihood of homozygous alternate genotype as
            # Phred-scaled proportion of distance between the observed
            # discordant depth and expected depth for the corresponding
            # concordant depth.
            genotype_likelihood = np.floor(-10 * np.log10(np.power(2, discordant) / np.power(2, expected_discordant_lower_bound)))
        elif expected_discordant_lower_bound <= discordant < expected_discordant_upper_bound:
            # Calculate likelihood of heterozygous genotype as Phred-scaled
            # proportion of distance between the observed discordant depth and
            # expected depth for the corresponding concordant depth.
            genotype_ratio = np.power(2, np.abs(discordant - expected_discordant_upper_bound)) / np.power(2, expected_discordant_upper_bound)
            genotype_likelihood = np.floor(-10 * np.log10(1 - min(genotype_ratio, 1)))
            genotype = "1/0"
        else:
            # Calculate likelihood of homozygous "reference" genotype as
            # Phred-scaled proportion of distance between the observed
            # discordant depth and expected depth for the corresponding
            # concordant depth.
            genotype_ratio = np.power(2, np.abs(discordant - expected_discordant_upper_bound)) / np.power(2, expected_discordant_upper_bound)
            genotype_likelihood = np.floor(-10 * np.log10(1 - min(genotype_ratio, 1)))
            genotype = "0/0"

    return genotype, genotype_likelihood


def get_depth_for_sv_call(sv_call, bams_by_name, chromosome_sizes):
    SLOP_FOR_BREAKPOINTS = 50

    if sv_call[EVENT_TYPE] == "deletion":
        breakpoint_intervals = [
            (sv_call[CONTIG_NAME], int(sv_call[CONTIG_START]), int(sv_call[CONTIG_END]) + SLOP_FOR_BREAKPOINTS),
        ]
        reference_intervals = [
            (sv_call[CHROMOSOME], int(sv_call[START]), int(sv_call[END]))
        ]
        reference_call_type = "insertion"
    else:
        breakpoint_intervals = [
            (sv_call[CONTIG_NAME], int(sv_call[CONTIG_START]), int(sv_call[CONTIG_END]))
        ]
        reference_intervals = [
            (sv_call[CHROMOSOME], int(sv_call[START]), int(sv_call[START]) + SLOP_FOR_BREAKPOINTS),
        ]
        reference_call_type = "deletion"

    logger.debug("Breakpoint intervals: %s", breakpoint_intervals)

    for bam_name, bam in bams_by_name.iteritems():
        breakpoint_concordant_depth = get_min_depth_for_region(bam, breakpoint_intervals, sv_call[EVENT_TYPE])
        breakpoint_discordant_depth = get_min_depth_for_region(bam, reference_intervals, reference_call_type)
        logger.debug("Found concordant depth for %s: %s", sv_call[EVENT_TYPE], breakpoint_concordant_depth)
        logger.debug("Found discordant depth for %s: %s", sv_call[EVENT_TYPE], breakpoint_discordant_depth)

        print("\t".join(map(str, (
            bam["sample"],
            sv_call[CHROMOSOME],
            sv_call[START],
            sv_call[END],
            sv_call[EVENT_TYPE],
            sv_call[CONTIG_NAME],
            sv_call[CONTIG_START],
            sv_call[CONTIG_END],
            breakpoint_concordant_depth,
            breakpoint_discordant_depth
        ))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sv_calls", help="BED file of SV calls with reference coordinates in columns 1-3, SV type in 4, length in 5, and contig coordinates in 6-8")
    parser.add_argument("bams", nargs="+", help="one or more BAMs per sample to genotype; read groups must be present in BAM")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--nproc", type=int, default=1, help="number of processes to run in parallel")
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        verbosity = 5
    else:
        verbosity = 0

    chromosome_sizes = None
    bams_by_name = defaultdict(dict)
    for bam in args.bams:
        bam_file = pysam.AlignmentFile(bam, "rb")

        # Prepare chromosome lengths based on one BAM's header.
        if chromosome_sizes is None:
            chromosome_sizes = dict(zip(bam_file.references, [(0, length) for length in bam_file.lengths]))
            logger.debug("Got %i chromosome sizes", len(chromosome_sizes))

        bams_by_name[bam]["file"] = None
        bams_by_name[bam]["filename"] = bam
        bams_by_name[bam]["sample"] = bam_file.header["RG"][0]["SM"]
        bam_file.close()

    columns = ("sample", "chr", "start", "end", "sv_call", "contig", "contig_start", "contig_end", "concordant", "discordant", "total_concordant", "total_discordant")
    print("\t".join(columns))

    sv_calls = pybedtools.BedTool(args.sv_calls)
    # for sv_call in sv_calls:
    #     get_depth_for_sv_call(sv_call, bams_by_name, chromosome_sizes)

    Parallel(n_jobs=args.nproc, verbose=verbosity)(delayed(get_depth_for_sv_call)(sv_call, bams_by_name, chromosome_sizes) for sv_call in sv_calls)
