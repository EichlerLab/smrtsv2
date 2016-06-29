import argparse
import pysam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bam", help="BAM to filter")
    parser.add_argument("queries_to_keep", help="a text file of query/read names to keep in the output BAM")
    parser.add_argument("output_bam", help="filtered BAM")
    args = parser.parse_args()

    input_bam = pysam.AlignmentFile(args.input_bam, "rb")
    output_bam = pysam.AlignmentFile(args.output_bam, "wb", template=input_bam)

    with open(args.queries_to_keep, "r") as fh:
        queries = set([query.strip() for query in fh])

    for alignment in input_bam:
        if alignment.query_name in queries:
            output_bam.write(alignment)

    input_bam.close()
    output_bam.close()
