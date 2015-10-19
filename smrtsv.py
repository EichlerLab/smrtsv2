import argparse
import subprocess
import sys

CLUSTER_SETTINGS = '" -q all.q -V -cwd -e ./log -o ./log {params.sge_opts} -w n -S /bin/bash"'
CLUSTER_FLAG = ("--drmaa", CLUSTER_SETTINGS, "-w", "30")

def _build_prefix(args):
    prefix = ["snakemake", "-pq"]
    if args.distribute:
        prefix.extend(CLUSTER_FLAG)

    return tuple(prefix)

def index(args):
    prefix = _build_prefix(args)
    args = prefix + ("prepare_reference", "--config", "reference=%s" % args.reference)
    args = " ".join(args)
    print args
    return subprocess.call(args, shell=True)

def align(args):
    print "Align reads"

def call(args):
    print "Call variants"

def genotype(args):
    print "Genotype"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--distribute", action="store_true", help="Distribute analysis to Grid Engine-style cluster")
    subparsers = parser.add_subparsers()

    # Index a reference for use by BLASR.
    parser_index = subparsers.add_parser("index", help="index a reference sequence for use by BLASR")
    parser_index.add_argument("reference", help="FASTA file of reference to index")
    parser_index.set_defaults(func=index)

    # Align PacBio reads to an indexed reference with BLASR.
    parser_align = subparsers.add_parser("align", help="align PacBio reads to an indexed reference with BLASR")
    parser_align.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_align.add_argument("reads", help="text file with one absolute path to a PacBio reads file (.bax.h5) per line")
    parser_align.add_argument("--alignments", help="text file with one absolute path to a BLASR alignments file (.bam) per line", default="alignments.fofn")
    parser_align.add_argument("--alignments_dir", help="absolute path of directory for BLASR alignment files", default="alignments")
    parser_align.set_defaults(func=align)

    # Call SVs and indels from BLASR alignments.
    parser_caller = subparsers.add_parser("call", help="call SVs and indels by local assembly of BLASR-aligned reads")
    parser_caller.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_caller.add_argument("alignments", help="text file with one absolute path to a BLASR alignments file (.bam) per line")
    parser_caller.add_argument("variants", help="VCF of variants called by local assembly alignments")
    parser_caller.set_defaults(func=call)

    # Genotype SVs with Illumina reads.
    parser_genotyper = subparsers.add_parser("genotype", help="Genotype SVs with Illumina reads")
    parser_genotyper.add_argument("variants", help="VCF of SMRT SV variants to genotype")
    parser_genotyper.add_argument("genotyped_variants", help="VCF of SMRT SV variant genotypes for the given sample-level BAMs")
    parser_genotyper.add_argument("samples", nargs="+", help="one or more sample-level BAMs to genotype for the given variants")
    parser_genotyper.set_defaults(func=genotype)

    args = parser.parse_args()
    return_code = args.func(args)
    sys.exit(return_code)
