import argparse
import pysam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq")
    parser.add_argument("--proportion_of_Ns_allowed", type=float, default=0.5)
    args = parser.parse_args()

    fh = pysam.FastqFile(args.fastq)
    for record in fh:
        if record.sequence.count("N") < args.proportion_of_Ns_allowed * len(record.sequence):
            print "@%s" % record.name
            print record.sequence
            print "+"
            print record.quality

    fh.close()
