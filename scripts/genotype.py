"""
Genotype a set of SV calls in a given set of BAMs.
"""
import argparse
from collections import defaultdict
from joblib import Parallel, delayed
import logging
import numpy as np
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


def get_min_depth_for_region(bam_fields, breakpoints):
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

    return high_quality_depth


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
        breakpoint_concordant_depth = get_min_depth_for_region(bam, breakpoint_intervals)
        breakpoint_discordant_depth = get_min_depth_for_region(bam, reference_intervals)
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
