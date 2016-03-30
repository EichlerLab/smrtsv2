import argparse
import logging
import subprocess
import sys
import os

logging.basicConfig(filename="smrtsv.log", level=logging.DEBUG)

CLUSTER_SETTINGS = '"qsub -sync y -q all.q -V -cwd -e ./log -o ./log {params.sge_opts} -S /bin/bash"'
CLUSTER_FLAG = ("--cluster-sync", CLUSTER_SETTINGS, "-w", "120")

def _get_dist_dir():
    dirname, filename = os.path.split(os.path.abspath(__file__))
    return dirname

def _build_prefix(args):
    prefix = ["snakemake", "--snakefile", os.path.join(_get_dist_dir(), "Snakefile"), "-pq", "-j", str(args.jobs)]
    if args.dryrun:
        prefix.append("-n")

    if args.distribute:
        prefix.extend(CLUSTER_FLAG)

    return tuple(prefix)

def index(args):
    prefix = _build_prefix(args)
    command = prefix + ("prepare_reference", "--config", "reference=%s" % args.reference)
    return subprocess.call(" ".join(command), shell=True)

def align(args):
    prefix = _build_prefix(args)
    command = prefix + (
        "align_reads",
        "--config",
        "reference=%s" % args.reference,
        "reads=%s" % args.reads,
        "alignments=%s" % args.alignments,
        "alignments_dir=%s" % args.alignments_dir,
        "batches=%s" % args.batches,
        "threads=%s" % args.threads
    )
    return subprocess.call(" ".join(command), shell=True)

def detect(args):
    """
    Detect SVs from signatures in read alignments.
    """
    # Find candidate regions in alignments.
    sys.stdout.write("Searching for candidate regions\n")
    prefix = _build_prefix(args)

    command = prefix + (
        "get_regions",
        "--config",
        "reference=%s" % args.reference,
        "alignments=%s" % args.alignments
    )

    if args.exclude:
        command = command + ("regions_to_exclude=%s" % args.exclude,)

    if args.candidates:
        command = command + ("candidates=%s" % args.candidates,)

    return subprocess.call(" ".join(command), shell=True)

def assemble(args):
    """
    Assemble candidate regions from raw reads aligned to regions.
    """
    # Generate local assemblies across the genome.
    sys.stdout.write("Starting local assemblies\n")
    prefix = _build_prefix(args)

    command = prefix + (
        "collect_assembly_alignments",
        "--config",
        "reference=%s" % args.reference,
        "alignments=%s" % args.alignments,
        "reads=%s" % args.reads
    )

    if args.regions:
        command = command + ("regions_to_assemble=%s" % args.regions,)

    if args.assembly_alignments:
        command = command + ("assembly_alignments=%s" % args.assembly_alignments,)

    logging.debug("Assembly command: %s", " ".join(command))
    return subprocess.call(" ".join(command), shell=True)

def call(args):
    # Call SVs, indels, and inversions.
    sys.stdout.write("Calling variants\n")
    prefix = _build_prefix(args)
    command = prefix + (
        "call_variants",
        "--config",
        "reference=%s" % args.reference,
        "alignments=%s" % args.alignments,
        "local_assembly_alignments=%s" % args.assembly_alignments,
        "variants=%s" % args.variants
    )
    return_code = subprocess.call(" ".join(command), shell=True)

    if return_code != 0:
        sys.stderr.write("Failed to call variants\n")
        return return_code

def run(args):
    # Detect SVs.
    return_code = detect(args)
    if return_code != 0:
        sys.stderr.write("Failed to identify candidate regions\n")
        return return_code

    # Run local assemblies.
    return_code = assemble(args)
    if return_code != 0:
        sys.stderr.write("Failed to generate local assemblies\n")
        return return_code

    # Call SVs, indels, and inversions.
    return_code = call(args)
    if return_code != 0:
        sys.stderr.write("Failed to call variants\n")
        return return_code

def genotype(args):
    print("Genotype")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dryrun", "-n", action="store_true", help="Print commands that will run without running them")
    parser.add_argument("--distribute", action="store_true", help="Distribute analysis to Grid Engine-style cluster")
    parser.add_argument("--jobs", help="number of jobs to run simultaneously", type=int, default=1)
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
    parser_align.add_argument("--batches", help="number of batches to split input reads into such that there will be one BAM output file per batch", type=int, default=1)
    parser_align.add_argument("--threads", help="number of threads to use for each BLASR alignment job", type=int, default=1)
    parser_align.set_defaults(func=align)

    # Detect SV signatures in BLASR alignments and build sliding windows to assemble.
    parser_detector = subparsers.add_parser("detect", help="detect SV signatures in BLASR-aligned reads")
    parser_detector.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_detector.add_argument("alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line")
    parser_detector.add_argument("candidates", help="BED file of candidates detected in read alignments")
    parser_detector.add_argument("--exclude", help="BED file of regions to exclude from local assembly (e.g., heterochromatic sequences, etc.)")
    parser_detector.set_defaults(func=detect)

    # Assemble candidate regions and align assemblies back to the reference.
    parser_assembler = subparsers.add_parser("assemble", help="assemble candidate regions and align assemblies back to the reference")
    parser_assembler.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_assembler.add_argument("reads", help="text file with one absolute path to a PacBio reads file (.bax.h5) per line")
    parser_assembler.add_argument("alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line")
    parser_assembler.add_argument("regions", help="BED file of regions to assemble from raw read alignments")
    parser_assembler.add_argument("assembly_alignments", help="BAM file with BLASR alignments of local assemblies against the reference")
    parser_assembler.set_defaults(func=assemble)

    # Call SVs and indels from BLASR alignments of local assemblies.
    parser_caller = subparsers.add_parser("call", help="call SVs and indels by BLASR alignments of local or whole genome assemblies")
    parser_caller.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_caller.add_argument("alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line")
    parser_caller.add_argument("assembly_alignments", help="BAM file with BLASR alignments of local assemblies against the reference")
    parser_caller.add_argument("variants", help="VCF of variants called by local assembly alignments")
    parser_caller.set_defaults(func=call)

    # Call SVs and indels from BLASR alignments of raw reads.
    parser_runner = subparsers.add_parser("run", help="call SVs and indels by local assembly of BLASR-aligned reads")
    parser_runner.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_runner.add_argument("reads", help="text file with one absolute path to a PacBio reads file (.bax.h5) per line")
    parser_runner.add_argument("alignments", help="text file with one absolute path to a BLASR alignments file (.bam) per line")
    parser_runner.add_argument("variants", help="VCF of variants called by local assembly alignments")
    parser_runner.add_argument("--exclude", help="BED file of regions to exclude from local assembly (e.g., heterochromatic sequences, etc.)")
    parser_runner.set_defaults(func=call)

    # Genotype SVs with Illumina reads.
    parser_genotyper = subparsers.add_parser("genotype", help="Genotype SVs with Illumina reads")
    parser_genotyper.add_argument("variants", help="VCF of SMRT SV variants to genotype")
    parser_genotyper.add_argument("genotyped_variants", help="VCF of SMRT SV variant genotypes for the given sample-level BAMs")
    parser_genotyper.add_argument("samples", nargs="+", help="one or more sample-level BAMs to genotype for the given variants")
    parser_genotyper.set_defaults(func=genotype)

    args = parser.parse_args()
    return_code = args.func(args)
    sys.exit(return_code)
