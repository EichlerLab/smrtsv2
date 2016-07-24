"""
Genotype a set of SV calls in a given set of BAMs.
"""
import argparse
from collections import defaultdict
import logging
import numpy as np
import pybedtools
import pysam

from get_best_alignment import get_best_alignments, get_depth_by_reference_and_position

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


def get_depth_for_regions(bam_fields, alt_breakpoints, ref_breakpoints, min_mapping_quality, min_base_quality):
    """
    Return median minus standard deviation of read depth across variant
    breakpoints in the given alternate and reference haplotypes.
    """
    # Convert breakpoints to a list for easier indexing.
    alt_regions = ["%s:%s-%s" % alt_breakpoint for alt_breakpoint in alt_breakpoints]
    ref_regions = ["%s:%s-%s" % ref_breakpoint for ref_breakpoint in ref_breakpoints]
    # alt_region_size = alt_breakpoints[2] - alt_breakpoints[1]
    # ref_region_size = ref_breakpoints[2] - ref_breakpoints[1]
    # logger.debug("Min depth for alt region %s (size: %s)", alt_region, alt_region_size)
    # logger.debug("Min depth for ref region %s (size: %s)", ref_region, ref_region_size)

    best_alignments = get_best_alignments(bam_fields["file"], ref_regions + alt_regions, quality=min_mapping_quality)
    depth_by_reference_and_position = get_depth_by_reference_and_position(best_alignments, bam_fields["file"], min_base_quality)

    alt_depths = [depth_by_reference_and_position.get(alt_breakpoint[0], {}).get(i, 0)
                  for alt_breakpoint in alt_breakpoints
                  for i in xrange(alt_breakpoint[1], alt_breakpoint[2] + 1)]
    ref_depths = [depth_by_reference_and_position.get(ref_breakpoint[0], {}).get(i, 0)
                  for ref_breakpoint in ref_breakpoints
                  for i in xrange(ref_breakpoint[1], ref_breakpoint[2] + 1)]

    logger.debug("Found %i depths for regions %s: %s", len(alt_depths), alt_regions, alt_depths)
    logger.debug("Found %i depths for regions %s: %s", len(ref_depths), ref_regions, ref_depths)
    logger.debug("Median alt depths: %s", np.median(alt_depths))
    logger.debug("Standard deviation alt depths: %s", np.std(alt_depths))
    logger.debug("Median ref depths: %s", np.median(ref_depths))
    logger.debug("Standard deviation ref depths: %s", np.std(ref_depths))
    high_quality_alt_depth = max(int(np.round(np.median(alt_depths) - np.std(alt_depths))), 0)
    high_quality_ref_depth = max(int(np.round(np.median(ref_depths) - np.std(ref_depths))), 0)

    return high_quality_alt_depth, high_quality_ref_depth


def add_slop_to_breakpoint(breakpoint, slop):
    return (breakpoint[0], max(0, breakpoint[1] - slop), breakpoint[2] + slop)


def get_depth_for_sv_call(sv_call, bams_by_name, chromosome_sizes, min_mapping_quality, min_base_quality, slop_for_breakpoints):
    # Deletions are single-base events in the alternate haplotype and multiple-base events in the reference.
    # Insertions are multiple-base events in the alternate haplotype and single-base events in the reference.
    if sv_call[EVENT_TYPE] == "deletion":
        breakpoint_intervals = ((sv_call[CONTIG_NAME], int(sv_call[CONTIG_START]), int(sv_call[CONTIG_START])),)
        reference_intervals = ((sv_call[CHROMOSOME], int(sv_call[START]), int(sv_call[START])),
                               (sv_call[CHROMOSOME], int(sv_call[END]), int(sv_call[END])))
        reference_call_type = "insertion"
    else:
        breakpoint_intervals = ((sv_call[CONTIG_NAME], int(sv_call[CONTIG_START]), int(sv_call[CONTIG_START])),
                                (sv_call[CONTIG_NAME], int(sv_call[CONTIG_END]), int(sv_call[CONTIG_END])))
        reference_intervals = ((sv_call[CHROMOSOME], int(sv_call[START]), int(sv_call[START])),)
        reference_call_type = "deletion"

    # Add slop to breakpoints.
    breakpoint_intervals = [add_slop_to_breakpoint(breakpoint, slop_for_breakpoints) for breakpoint in breakpoint_intervals]
    reference_intervals = [add_slop_to_breakpoint(breakpoint, slop_for_breakpoints) for breakpoint in reference_intervals]

    logger.debug("Alternate intervals: %s", breakpoint_intervals)
    logger.debug("Reference intervals: %s", reference_intervals)

    for bam_name, bam in bams_by_name.iteritems():
        breakpoint_concordant_depth, breakpoint_discordant_depth = get_depth_for_regions(bam, breakpoint_intervals, reference_intervals, min_mapping_quality, min_base_quality)
        logger.debug("Found concordant/discordant depth for %s: %s / %s", sv_call[EVENT_TYPE], breakpoint_concordant_depth, breakpoint_discordant_depth)

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
    parser.add_argument("--min_mapping_quality", type=int, default=20, help="minimum mapping quality of alignments to consider for genotyping")
    parser.add_argument("--min_base_quality", type=int, default=20, help="minimum base quality to consider for calculating read depth per haplotype")
    parser.add_argument("--slop_for_breakpoints", type=int, default=25, help="number of bases on either side of an insertion breakpoint to calculate read depth across")
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

        bams_by_name[bam]["file"] = bam_file
        bams_by_name[bam]["filename"] = bam
        bams_by_name[bam]["sample"] = bam_file.header["RG"][0]["SM"]

    columns = ("sample", "chr", "start", "end", "sv_call", "contig", "contig_start", "contig_end", "concordant", "discordant")
    print("\t".join(columns))

    sv_calls = pybedtools.BedTool(args.sv_calls)
    for sv_call in sv_calls:
        get_depth_for_sv_call(sv_call, bams_by_name, chromosome_sizes, args.min_mapping_quality, args.min_base_quality, args.slop_for_breakpoints)
