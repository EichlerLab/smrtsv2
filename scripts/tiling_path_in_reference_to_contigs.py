import argparse
from joblib import Parallel, delayed
import pybedtools
import pysam


def get_path_for_interval(interval, bam_name):
    bam = pysam.AlignmentFile(bam_name, "rb")
    reads = [read for read in bam.fetch(interval.chrom, interval.start, interval.end)]
    query = interval.name
    query_read = None

    for read in reads:
        if read.query_name == query:
            query_read = read
            break

    query_start = None
    query_end = None
    if query_read is not None:
        for pair in query_read.get_aligned_pairs():
            if query_start is None and pair[0] is not None and pair[1] is not None and pair[1] >= interval.start:
                query_start = pair[0]

            if query_start is not None and pair[0] is not None and pair[1] is not None:
                query_end = pair[0]

            if pair[1] is not None and pair[1] >= interval.end:
                break

    if query_start is not None and query_end is not None:
        print("\t".join(map(str, (interval.chrom, interval.start, interval.end, query, query_start, query_end))))

    bam.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="BAM of local assemblies aligned to a reference")
    parser.add_argument("tiling_path", help="BED file with tiling path across the reference built from local assembly alignments")
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    bed = pybedtools.BedTool(args.tiling_path)

    if args.debug:
        verbosity = 5
    else:
        verbosity = 0

    Parallel(n_jobs=args.nproc, verbose=verbosity)(delayed(get_path_for_interval)(interval, args.bam) for interval in bed)
