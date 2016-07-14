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
import subprocess
import tempfile

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


def get_depth_for_regions(bam_fields, alt_breakpoints, ref_breakpoints):
    """
    Return median minus standard deviation of read depth across variant
    breakpoints in the given alternate and reference haplotypes.
    """
    # Convert breakpoints to a list for easier indexing.
    alt_region = "%s:%s-%s" % alt_breakpoints
    ref_region = "%s:%s-%s" % ref_breakpoints
    alt_region_size = alt_breakpoints[2] - alt_breakpoints[1]
    ref_region_size = ref_breakpoints[2] - ref_breakpoints[1]
    logger.debug("Min depth for alt region %s (size: %s)", alt_region, alt_region_size)
    logger.debug("Min depth for ref region %s (size: %s)", ref_region, ref_region_size)

    arguments = {
        "BAM": bam_fields["filename"],
        "REF_REGION": ref_region,
        "ALT_REGION": alt_region,
        "TMPDIR": tempfile.gettempdir()
    }
    command = "samtools view -bu %(BAM)s '%(REF_REGION)s' '%(ALT_REGION)s' | samtools sort -l 0 -n -O bam -T %(TMPDIR)s/reads - | python ~jlhudd/src/smrtsv/scripts/get_best_alignment.py /dev/stdin | samtools sort -l 0 -O bam -T %(TMPDIR)s/filtered_reads - | samtools depth -q 20 -Q 20 -" % arguments
    logger.debug("Running command: %s", command)
    result = subprocess.check_output(command, shell=True)
    alt_depths = []
    ref_depths = []
    for line in result.splitlines():
        contig, position, depth = line.split("\t")
        position = int(position)
        depth = int(depth)

        if contig == ref_breakpoints[0] and position >= ref_breakpoints[1] and position <= ref_breakpoints[2]:
            ref_depths.append(depth)
        elif contig == alt_breakpoints[0] and position >= alt_breakpoints[1] and position <= alt_breakpoints[2]:
            alt_depths.append(depth)

    # If there are not depths for every base of the given region, fill those
    # missing values with enough additional zeroes to equal the region size.
    if len(alt_depths) < alt_region_size:
        alt_depths = alt_depths + [0] * (alt_region_size - len(alt_depths))

    if len(ref_depths) < ref_region_size:
        ref_depths = ref_depths + [0] * (ref_region_size - len(ref_depths))

    logger.debug("Found %i depths for region %s: %s", len(alt_depths), alt_region, alt_depths)
    logger.debug("Found %i depths for region %s: %s", len(ref_depths), ref_region, ref_depths)
    logger.debug("Median alt depths: %s", np.median(alt_depths))
    logger.debug("Standard deviation alt depths: %s", np.std(alt_depths))
    logger.debug("Median ref depths: %s", np.median(ref_depths))
    logger.debug("Standard deviation ref depths: %s", np.std(ref_depths))
    high_quality_alt_depth = max(int(np.round(np.median(alt_depths) - np.std(alt_depths))), 0)
    high_quality_ref_depth = max(int(np.round(np.median(ref_depths) - np.std(ref_depths))), 0)

    return high_quality_alt_depth, high_quality_ref_depth


def get_depth_for_sv_call(sv_call, bams_by_name, chromosome_sizes):
    SLOP_FOR_BREAKPOINTS = 25

    if sv_call[EVENT_TYPE] == "deletion":
        breakpoint_intervals = (sv_call[CONTIG_NAME], max(0, int(sv_call[CONTIG_START]) - SLOP_FOR_BREAKPOINTS), int(sv_call[CONTIG_END]) + SLOP_FOR_BREAKPOINTS)
        reference_intervals = (sv_call[CHROMOSOME], int(sv_call[START]), int(sv_call[END]))
        reference_call_type = "insertion"
    else:
        breakpoint_intervals = (sv_call[CONTIG_NAME], int(sv_call[CONTIG_START]), int(sv_call[CONTIG_END]))
        reference_intervals = (sv_call[CHROMOSOME], max(0, int(sv_call[START]) - SLOP_FOR_BREAKPOINTS), int(sv_call[START]) + SLOP_FOR_BREAKPOINTS)
        reference_call_type = "deletion"

    logger.debug("Breakpoint intervals: %s", breakpoint_intervals)

    for bam_name, bam in bams_by_name.iteritems():
        breakpoint_concordant_depth, breakpoint_discordant_depth = get_depth_for_regions(bam, breakpoint_intervals, reference_intervals)
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

    columns = ("sample", "chr", "start", "end", "sv_call", "contig", "contig_start", "contig_end", "concordant", "discordant")
    print("\t".join(columns))

    sv_calls = pybedtools.BedTool(args.sv_calls)
    # for sv_call in sv_calls:
    #     get_depth_for_sv_call(sv_call, bams_by_name, chromosome_sizes)

    Parallel(n_jobs=args.nproc, verbose=verbosity)(delayed(get_depth_for_sv_call)(sv_call, bams_by_name, chromosome_sizes) for sv_call in sv_calls)
