import argparse
from collections import defaultdict
import logging
import numpy as np
import pysam
import sys

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


def get_depth_by_reference_and_position(alignments, bam, base_quality):
    # Index number of reads supporting a given position in a given reference.
    depth_by_reference_and_position = {}

    for alignment in alignments:
        # Index depth by human readable reference name for reuse by other tools
        # that don't have access to the samtools reference id lookup table.
        ref = bam.getrname(alignment.reference_id)

        # Create an empty counter by position for the current reference contig
        # if one doesn't exist.
        if ref not in depth_by_reference_and_position:
            depth_by_reference_and_position[ref] = defaultdict(int)

        # Filter aligned base pairs between query and reference sequences to
        # exclude gaps in the query. These gaps should not count toward the
        # overall depth across a given reference position and they cause the
        # number of aligned pairs to exceed the number of bases in the query.
        pairs = np.array([pair for pair in alignment.aligned_pairs if pair[0] is not None])

        # Consider only aligned pairs where the query's base qualities are
        # greater than or equal to the minimum base quality threshold.
        for base in pairs[np.array(alignment.query_qualities) >= base_quality][:,1]:
            if base is not None:
                depth_by_reference_and_position[ref][base] += 1

    return depth_by_reference_and_position


def calculate_differences_in_alignment(alignment):
    # Calculate indel and softclipped bases from the sum of lengths for CIGAR
    # options that aren't matches.
    indels_and_softclips = sum([cigar[1] for cigar in alignment.cigar if cigar[0] != BAM_CMATCH])

    # Calculate mismatches from the NM tags.
    mismatches = alignment.get_tag("NM")

    differences = indels_and_softclips + mismatches

    logger.debug("Found %i differences in alignment %s", differences, alignment)
    return differences


def get_best_alignment(alignments_by_read_name):
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
    filtered_alignments = []
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
                filtered_alignments.append(alignments_by_differences[min_differences][0])
                logger.debug("Found best alignment for %s", name)
            else:
                logger.debug("Skipping equally scored alignments for %s", name)
        else:
            # Emit the only read in the set.
            filtered_alignments.append(alignments[0])

    return filtered_alignments


def get_best_alignments(bam, regions, quality):
    # Load alignments for requested regions into memory and sort by read name.
    # Do our own fancy splitting of regions by colon and then by hyphen to allow
    # unpleasant local assembly contig names to be parsed correctly.
    alignments_from_regions = [read
                               for region in regions
                               for read in bam.fetch(*([region.split(":")[0]] + map(int, region.split(":")[1].split("-"))))
                               if read.mapping_quality >= quality]
    logger.debug("Loaded %i reads from BAM from %i regions", len(alignments_from_regions), len(regions))

    # Sort alignments by read name.
    logger.debug("Sorting alignments by read name")
    alignments_from_regions = sorted(alignments_from_regions, key=lambda alignment: alignment.query_name)

    current_read_name = None
    alignments_by_read_name = defaultdict(list)

    logger.debug("Filtering alignments")
    alignments_processed = 0
    filtered_alignments_from_regions = []
    for alignment in alignments_from_regions:
        alignments_processed += 1
        if alignment.query_name != current_read_name:
            if current_read_name is not None:
                filtered_alignments_from_regions.extend(get_best_alignment(alignments_by_read_name))

            current_read_name = alignment.query_name
            alignments_by_read_name = defaultdict(list)

        # Index alignments by query name and read number to avoid treating readsn 1
        # and 2 of the same pair as duplicates.
        read_name = "%s/%s" % (alignment.query_name, alignment.is_read1 and 1 or 2)
        alignments_by_read_name[read_name].append(alignment)
    else:
        # Get best alignments from the last set in the for loop.
        filtered_alignments_from_regions.extend(get_best_alignment(alignments_by_read_name))

    # Sort alignments by chromosome and start position.
    logger.debug("Sorting alignments by chromosome and start position")
    filtered_alignments_from_regions = sorted(
        filtered_alignments_from_regions,
        key=lambda alignment: (bam.getrname(alignment.reference_id), alignment.reference_start)
    )

    logger.debug("Processed %i alignments and retained %i", alignments_processed, len(filtered_alignments_from_regions))

    # Return filtered alignments .
    return filtered_alignments_from_regions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignments")
    parser.add_argument("output")
    parser.add_argument("regions", nargs="+", help="one or more samtools-style regions (e.g., chr1:1000-2000) to filter")
    parser.add_argument("--quality", type=int, default=20, help="minimum mapping quality for reads to output")
    parser.add_argument("--sam", action="store_true", help="output SAM for debugging")
    parser.add_argument("--uncompressed", action="store_true", help="output uncompressed BAM for piping")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    if args.sam:
        write_mode = "w"
    elif args.uncompressed:
        write_mode = "wbu"
    else:
        write_mode = "wb"

    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)

    with pysam.AlignmentFile(args.alignments, "rb") as bam:
        best_alignments = get_best_alignments(bam, args.regions, args.quality)

        with pysam.AlignmentFile(args.output, write_mode, template=alignments) as filtered_alignments:
            for alignment in best_alignments:
                filtered_alignments.write(alignment)
