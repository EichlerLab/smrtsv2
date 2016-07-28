#!/bin/env python
import argparse
import logging
import subprocess
import sys
import os

# Set logging
logging.basicConfig(filename="smrtsv.log", level=logging.DEBUG)

# Set cluster parameters
CLUSTER_SETTINGS = ' -V -cwd -e ./log -o ./log {params.sge_opts} -w n -S /bin/bash'
CLUSTER_FLAG = ("--drmaa", CLUSTER_SETTINGS, "-w", "60")

# Setup environment for executing commands
PROCESS_ENV = os.environ.copy()

# Prepend to PROCESS_ENV["PATH"]
INSTALL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

INSTALL_PATH = [  # List of paths relative to INSTALL_DIR to be added to the environment $PATH
    "bin",
    "dist/miniconda/envs/python2/bin",
    "dist/miniconda/envs/python3/bin",
    "dist/miniconda/bin",
    "dist/celera/wgs-8.3rc2/Linux-amd64/bin/",
    "dist/amos-3.1.0/bin",
    "canu/Linux-amd64/bin"
]

PROCESS_ENV_PATH = ":".join([os.path.join(INSTALL_DIR, THIS_PATH) for THIS_PATH in INSTALL_PATH])

if "PATH" in PROCESS_ENV:
    PROCESS_ENV["PATH"] = PROCESS_ENV_PATH + ":" + PROCESS_ENV["PATH"]
else:
    PROCESS_ENV["PATH"] = PROCESS_ENV_PATH

# Prepend to PROCESS_ENV["LD_LIBRARY_PATH"]
INSTALL_LD_PATH = [
    "dist/hdf5/lib"
]

PROCESS_ENV_LD_PATH = ":".join([os.path.join(INSTALL_DIR, THIS_PATH) for THIS_PATH in INSTALL_LD_PATH])

if "LD_LIBRARY_PATH" in PROCESS_ENV:
    PROCESS_ENV["LD_LIBRARY_PATH"] = PROCESS_ENV_LD_PATH + ":" + PROCESS_ENV["LD_LIBRARY_PATH"]
else:
    PROCESS_ENV["LD_LIBRARY_PATH"] = PROCESS_ENV_LD_PATH

os.environ["LD_LIBRARY_PATH"] = PROCESS_ENV["LD_LIBRARY_PATH"]


# Function definitions
def _get_dist_dir():
    dirname, filename = os.path.split(os.path.abspath(__file__))
    return dirname

# def _build_prefix(args):
#     prefix = ["snakemake", "-T", "--rerun-incomplete", "--snakefile", os.path.join(os.path.dirname(_get_dist_dir()), "Snakefile"), "-j", str(args.jobs)]
#     if args.dryrun:
#         prefix.append("-n")
#
#     if args.distribute:
#         prefix.extend(CLUSTER_FLAG)
#
#     return tuple(prefix)

def _run_cmd(args):
    """
    Run a command with the proper environment set.

    :param args: A tuple of arguments starting with the command name.

    :return: Return code or -1 if the process did not complete.
    """
    sys.stdout.flush()

    p = subprocess.Popen(args, env=PROCESS_ENV)

    p.wait()

    ret_code = p.returncode

    return ret_code if ret_code is not None else -1

def _run_snake_target(args, *cmd):
    """
    Run a snakemake target.

    :param args: Arguments processed from the command line.
    :param cmd: The command to run as a tuple starting with the name of the snakemake target.

    :return: Return code from snakemake.
    """

    # Setup snakemake command
    prefix = [
        "snakemake",
        "-T",
        "--rerun-incomplete",
        "--snakefile", os.path.join(os.path.dirname(_get_dist_dir()), "Snakefile"),
        "-j", str(args.jobs)
    ]

    if args.dryrun:
        prefix.append("-n")

    if args.distribute:
        prefix.extend(CLUSTER_FLAG)

    # Append command
    prefix.extend(cmd)

    # Report (verbose)
    if args.verbose:
        print("Running snakemake command: %s" % " ".join(prefix))

    # Run snakemake command
    return _run_cmd(prefix)

def index(args):
    return _run_snake_target(
        args,
        "prepare_reference",
        "--config",
        "reference=%s" % args.reference
    )

def align(args):
    return _run_snake_target(
        args,
        "align_reads",
        "--config",
        "reference=%s" % args.reference,
        "reads=%s" % args.reads,
        "alignments=%s" % args.alignments,
        "alignments_dir=%s" % args.alignments_dir,
        "batches=%s" % args.batches,
        "threads=%s" % args.threads,
        "tmp_dir=%s" % args.tmpdir
    )

def detect(args):
    """
    Detect SVs from signatures in read alignments.
    """
    # Find candidate regions in alignments.
    sys.stdout.write("Searching for candidate regions\n")

    command = (
        "get_regions",
        "--config",
        "reference=%s" % args.reference,
        "alignments=%s" % args.alignments,
        "assembly_window_size=%s" % args.assembly_window_size,
        "assembly_window_slide=%s" % args.assembly_window_slide
    )

    if args.exclude:
        command = command + ("regions_to_exclude=%s" % args.exclude,)

    if args.candidates:
        command = command + ("candidates=%s" % args.candidates,)

    return _run_snake_target(args, *command)

def assemble(args):
    """
    Assemble candidate regions from raw reads aligned to regions.
    """
    # Generate local assemblies across the genome.
    sys.stdout.write("Starting local assemblies\n")

    base_command = (
        "collect_assembly_alignments",
        "--config",
        "reference=%s" % args.reference,
        "alignments=%s" % args.alignments,
        "reads=%s" % args.reads,
        "tmp_dir=%s" % args.tmpdir,
        "ld_path=%s" % PROCESS_ENV["LD_LIBRARY_PATH"],
        "path=%s" % PROCESS_ENV["PATH"]
    )

    if args.regions:
        # For each contig/chromosome in the regions file, submit a separate
        # Snakemake command. To do so, first split regions to assemble into one
        # file per contig in a temporary directory.
        tmpdir = os.path.join(os.getcwd(), "regions_by_contig")

        rebuild_regions_by_contig = False
        if not args.dryrun and (not os.path.exists(tmpdir) or args.rebuild_regions):
            rebuild_regions_by_contig = True

        if rebuild_regions_by_contig:
            try:
                os.mkdir(tmpdir)
            except OSError:
                pass

        previous_contig = None
        with open(args.regions, "r") as fh:
            contigs = set()
            for line in fh:
                contig = line.strip().split()[0]

                if previous_contig != contig:
                    if previous_contig is not None and rebuild_regions_by_contig:
                        contig_file.close()

                    previous_contig = contig
                    contigs.add(contig)

                    if rebuild_regions_by_contig:
                        contig_file = open(os.path.join(tmpdir, "%s.bed" % contig), "w")

                if rebuild_regions_by_contig:
                    contig_file.write(line)

        if rebuild_regions_by_contig:
            contig_file.close()

        # Assemble regions per contig creating a single merged BAM for each contig.
        local_assembly_basename = os.path.basename(args.assembly_alignments)
        local_assemblies = set()

        return_code = 0

        for contig in contigs:
            contig_local_assemblies = os.path.join("local_assemblies", local_assembly_basename.replace(".bam", ".%s.bam" % contig))
            local_assemblies.add(contig_local_assemblies)

            if os.path.exists(contig_local_assemblies):
                sys.stdout.write("Local assemblies already exist for %s\n" % contig)
                continue

            command = base_command + ("regions_to_assemble=%s" % os.path.join(tmpdir, "%s.bed" % contig),)
            command = command + ("assembly_alignments=%s" % contig_local_assemblies,)
            sys.stdout.write("Starting local assemblies for %s\n" % contig)
            logging.debug("Assembly command: %s", " ".join(command))

            return_code = _run_snake_target(args, *command)

            if return_code != 0:
                break

        # If the last command executed successfully, try to merge all local
        # assemblies per contig into a single file.
        if not args.dryrun and return_code == 0:
            return_code = _run_cmd(["samtools", "merge", args.assembly_alignments] + list(local_assemblies))

            if return_code == 0:
                return_code = _run_cmd(["samtools", "index", args.assembly_alignments])

        # Return the last return code.
        return return_code
    else:
        if args.assembly_alignments:
            command = base_command + ("assembly_alignments=%s" % args.assembly_alignments,)

            logging.debug("Assembly command: %s", " ".join(command))
            return _run_cmd(command)

def call(args):
    # Call SVs, indels, and inversions.
    sys.stdout.write("Calling variants\n")

    return_code = _run_snake_target(
        args,
        "call_variants",
        "--config",
        "reference=%s" % args.reference,
        "alignments=%s" % args.alignments,
        "local_assembly_alignments=%s" % args.assembly_alignments,
        "variants=%s" % args.variants,
        "species=\"%s\"" % args.species,
        "sample=\"%s\"" % args.sample
    )

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


# Main
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--dryrun", "-n", action="store_true", help="Print commands that will run without running them")
    parser.add_argument("--distribute", action="store_true", help="Distribute analysis to Grid Engine-style cluster")
    parser.add_argument("--jobs", help="number of jobs to run simultaneously", type=int, default=1)
    parser.add_argument("--tmpdir", help="temporary directory to use for distributed jobs", default="/var/tmp")
    parser.add_argument("--verbose", "-v", help="print extra runtime information", action="store_true")
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
    parser_detector.add_argument("--assembly_window_size", type=int, help="size of reference window for local assemblies", default=60000)
    parser_detector.add_argument("--assembly_window_slide", type=int, help="size of reference window slide for local assemblies", default=30000)
    parser_detector.set_defaults(func=detect)

    # Assemble candidate regions and align assemblies back to the reference.
    parser_assembler = subparsers.add_parser("assemble", help="assemble candidate regions and align assemblies back to the reference")
    parser_assembler.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_assembler.add_argument("reads", help="text file with one absolute path to a PacBio reads file (.bax.h5) per line")
    parser_assembler.add_argument("alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line")
    parser_assembler.add_argument("regions", help="BED file of regions to assemble from raw read alignments")
    parser_assembler.add_argument("assembly_alignments", help="BAM file with BLASR alignments of local assemblies against the reference")
    parser_assembler.add_argument("--rebuild_regions", action="store_true", help="rebuild subset of regions to assemble")
    parser_assembler.set_defaults(func=assemble)

    # Call SVs and indels from BLASR alignments of local assemblies.
    parser_caller = subparsers.add_parser("call", help="call SVs and indels by BLASR alignments of local or whole genome assemblies")
    parser_caller.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_caller.add_argument("alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line")
    parser_caller.add_argument("assembly_alignments", help="BAM file with BLASR alignments of local assemblies against the reference")
    parser_caller.add_argument("variants", help="VCF of variants called by local assembly alignments")
    parser_caller.add_argument("--sample", help="Sample name to use in final variant calls", default="UnnamedSample")
    parser_caller.add_argument("--species", help="Common or scientific species name to pass to RepeatMasker", default="human")
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

    if args.verbose:
        print("PATH=%s" % PROCESS_ENV["PATH"])
        print("LD_LIBRARY_PATH=%s" % PROCESS_ENV["LD_LIBRARY_PATH"])
        sys.stdout.flush()

    # Make a log directory for grid-engine-style error logs if commands are
    # being distributed in non-dryrun mode.
    if args.distribute and not args.dryrun and not os.path.isdir("log"):
        os.mkdir("log")

    return_code = args.func(args)
    sys.exit(return_code)
