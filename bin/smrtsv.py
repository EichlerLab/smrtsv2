#!/bin/env python
import argparse
import logging
import subprocess
import sys
import os
import re

# Set logging
logging.basicConfig(filename="smrtsv.log", level=logging.DEBUG)

# Set cluster parameters
CLUSTER_SETTINGS = ' -V -cwd -e ./log -o ./log {cluster.params} -w n -S /bin/bash'
CLUSTER_FLAG = ("--drmaa", CLUSTER_SETTINGS, "-w", "120")

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

# Prepend to PERL5LIB
INSTALL_PERL5LIB = [
    "dist/pm"
]

PROCESS_PERL5LIB = ":".join([os.path.join(INSTALL_DIR, THIS_PATH) for THIS_PATH in INSTALL_PERL5LIB])

if "PERL5LIB" in PROCESS_ENV:
    PROCESS_ENV["PERL5LIB"] = PROCESS_PERL5LIB + ":" + PROCESS_ENV["PERL5LIB"]
else:
    PROCESS_ENV["PERL5LIB"] = PROCESS_PERL5LIB

os.environ["PERL5LIB"] = PROCESS_ENV["PERL5LIB"]


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
    # Use the user-defined cluster config path if one is given. Otherwise, use
    # an empty config that comes with the SMRT-SV distribution.
    if args.cluster_config is not None:
        cluster_config_path = args.cluster_config
    else:
        cluster_config_path = os.path.join(os.path.dirname(_get_dist_dir()), "cluster.template.json")

    # Setup snakemake command
    prefix = [
        "snakemake",
        "-T",
        "--rerun-incomplete",
        "--cluster-config", cluster_config_path,
        "--snakefile", os.path.join(os.path.dirname(_get_dist_dir()), "Snakefile"),
        "-j", str(args.jobs)
    ]

    if args.dryrun:
        prefix.append("-n")

    if args.distribute:
        prefix.extend(CLUSTER_FLAG)

    # Append command
    prefix.extend(cmd)

    # Append path and ld_path
    prefix.extend([
        "ld_path=%s" % PROCESS_ENV["LD_LIBRARY_PATH"],
        "path=%s" % PROCESS_ENV["PATH"],
        "perl5lib=%s" % PROCESS_ENV["PERL5LIB"]
    ])

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
        "tmp_dir=%s" % args.tmpdir,
        "alignment_parameters=\"%s\"" % args.alignment_parameters
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
        "assembly_window_slide=%s" % args.assembly_window_slide,
        "min_length=%s" % args.min_length,
        "min_support=%s" % args.min_support,
        "max_support=%s" % args.max_support,
        "min_coverage=%s" % args.min_coverage,
        "max_coverage=%s" % args.max_coverage,
        "min_hardstop_support=%s" % args.min_hardstop_support,
        "max_candidate_length=%s" % args.max_candidate_length
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
        "asm_alignment_parameters=\"%s\"" % args.asm_alignment_parameters,
        "mapping_quality=\"%s\"" % args.mapping_quality,
        "minutes_to_delay_jobs=\"%s\"" % args.minutes_to_delay_jobs,
        "assembly_log=\"%s\"" % args.assembly_log
    )

    if args.candidates:
        # For each contig/chromosome in the candidates file, submit a separate
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
        contig_file = None

        with open(args.candidates, "r") as fh:
            contigs = set()
            for line in fh:
                contig = line.strip().split()[0]

                if previous_contig != contig:
                    if previous_contig is not None and rebuild_regions_by_contig:
                        contig_file.close()
                        contig_file = None

                    previous_contig = contig
                    contigs.add(contig)

                    if rebuild_regions_by_contig:
                        contig_file = open(os.path.join(tmpdir, "%s.bed" % contig), "w")

                if rebuild_regions_by_contig:
                    contig_file.write(line)

        if rebuild_regions_by_contig:
            if contig_file is not None:
                contig_file.close()

        # Assemble regions per contig creating a single merged BAM for each contig.
        local_assembly_basename = os.path.basename(args.assembly_alignments)
        local_assemblies = set()

        return_code = 0

        contig_count = 0  # Number of contigs assemblies were generated for

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

            contig_count += 1

            if return_code != 0:
                break

        # If the last command executed successfully, try to merge all local
        # assemblies per contig into a single file. Only build if at least one set of local assemblies was performed
        # or the local assemblies file does not exist.
        if not args.dryrun and return_code == 0 and (contig_count > 0 or not os.path.exists(args.assembly_alignments)):
            if len(local_assemblies) > 1:
                return_code = _run_cmd(["samtools", "merge", args.assembly_alignments] + list(local_assemblies))
            else:
                return_code = _run_cmd(["samtools", "view", "-b", "-o", args.assembly_alignments] + list(local_assemblies))

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

    # Get default jobs
    if "jobs" in args:
        default_jobs = args.jobs
    else:
        default_jobs = 1

    # Get the number of jobs for each step
    job_step = re.split("\\s*[,;:]\\s*", args.runjobs.strip())  # Split into array
    job_step = [job_step[i] if len(job_step) > i else '' for i in range(4)]  # Extend to length 4

    # Convert each number of jobs to integers
    for i in range(4):
        if job_step[i] != '':
            try:
                job_step[i] = int(job_step[i])
            except ValueError:
                sys.stderr.write("Invalid number of jobs for step %d: Must be an integer: \"%s\"\n" % ((i + 1), job_step[i]))
                return 1
        else:
            job_step[i] = default_jobs

    # Report the number of jobs for each task
    if args.verbose and args.distribute:
        print("Jobs per task:")
        print("\t*    Align: %s" % job_step[0])
        print("\t*   Detect: %s" % job_step[1])
        print("\t* Assemble: %s" % job_step[2])
        print("\t*     Call: %s" % job_step[3])

    # Build reference indices
    return_code = index(args)
    if return_code != 0:
        sys.stderr.write("Failed to index reference\n")
        return return_code

    # Align
    args.jobs = job_step[0]

    return_code = align(args)
    if return_code != 0:
        sys.stderr.write("Failed to align reads\n")
        return return_code

    # Detect SVs.
    args.jobs = job_step[1]

    return_code = detect(args)
    if return_code != 0:
        sys.stderr.write("Failed to identify candidate regions\n")
        return return_code

    # Run local assemblies.
    args.jobs = job_step[2]

    return_code = assemble(args)
    if return_code != 0:
        sys.stderr.write("Failed to generate local assemblies\n")
        return return_code

    # Call SVs, indels, and inversions.
    args.jobs = job_step[3]

    return_code = call(args)
    if return_code != 0:
        sys.stderr.write("Failed to call variants\n")
        return return_code

    return 0

def genotype(args):
    # Genotype SVs.
    sys.stdout.write("Genotyping SVs\n")

    return_code = _run_snake_target(
        args,
        "convert_genotypes_to_vcf",
        "--config",
        "genotyper_config=%s" % args.genotyper_config,
        "genotyped_variants=%s" % args.genotyped_variants,
        "threads=%s" % args.threads
    )

    if return_code != 0:
        sys.stderr.write("Failed to genotype SVs\n")

    return return_code


# Main
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--dryrun", "-n", action="store_true", help="Print commands that will run without running them")
    parser.add_argument("--distribute", action="store_true", help="Distribute analysis to Grid Engine-style cluster")
    parser.add_argument("--jobs", help="number of jobs to run simultaneously", type=int, default=1)
    parser.add_argument("--tmpdir", help="temporary directory to use for distributed jobs. Keyword \"default\" tells Python to pick a directory.", default="default")
    parser.add_argument("--verbose", "-v", help="print extra runtime information", action="store_true")
    parser.add_argument("--cluster_config", help="JSON/YAML file specifying cluster configuration parameters to pass to Snakemake's --cluster-config option")
    parser.add_argument("--drmaalib", help="For jobs that are distributed, this is the location to the DRMAA library (libdrmaa.so) installed with Grid Engine. Use this to set DRMAA_LIBRARY_PATH in the environment for pipelined commands. If DRMAA_LIBRARY_PATH is already set in the environment when calling this program, this option is not required.")
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
    parser_align.add_argument("--alignment_parameters", help="BLASR parameters to use to align raw reads", default="-bestn 2 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 50")
    parser_align.set_defaults(func=align)

    # Detect SV signatures in BLASR alignments and build sliding windows to assemble.
    parser_detector = subparsers.add_parser("detect", help="detect SV signatures in BLASR-aligned reads")
    parser_detector.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_detector.add_argument("alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line")
    parser_detector.add_argument("candidates", help="BED file of candidates detected in read alignments")
    parser_detector.add_argument("--exclude", help="BED file of regions to exclude from local assembly (e.g., heterochromatic sequences, etc.)")
    parser_detector.add_argument("--assembly_window_size", type=int, help="size of reference window for local assemblies", default=60000)
    parser_detector.add_argument("--assembly_window_slide", type=int, help="size of reference window slide for local assemblies", default=20000)
    parser_detector.add_argument("--min_length", type=int, help="minimum length required for SV candidates", default=50)
    parser_detector.add_argument("--min_support", type=int, help="minimum number of supporting reads required to flag a region as an SV candidate", default=5)
    parser_detector.add_argument("--max_support", type=int, help="maximum number of supporting reads allowed to flag a region as an SV candidate", default=100)
    parser_detector.add_argument("--min_coverage", type=int, help="minimum number of total reads required to flag a region as an SV candidate", default=5)
    parser_detector.add_argument("--max_coverage", type=int, help="maximum number of total reads allowed to flag a region as an SV candidate", default=100),
    parser_detector.add_argument("--min_hardstop_support", type=int, help="minimum number of reads with hardstops required to flag a region as an SV candidate", default=11)
    parser_detector.add_argument("--max_candidate_length", type=int, help="maximum length allowed for an SV candidate region", default=60000)
    parser_detector.set_defaults(func=detect)

    # Assemble candidate regions and align assemblies back to the reference.
    parser_assembler = subparsers.add_parser("assemble", help="assemble candidate regions and align assemblies back to the reference")
    parser_assembler.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_assembler.add_argument("reads", help="text file with one absolute path to a PacBio reads file (.bax.h5) per line")
    parser_assembler.add_argument("alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line")
    parser_assembler.add_argument("candidates", help="BED file of regions to assemble from raw read alignments")
    parser_assembler.add_argument("assembly_alignments", help="BAM file with BLASR alignments of local assemblies against the reference")
    parser_assembler.add_argument("--rebuild_regions", action="store_true", help="rebuild subset of regions to assemble")
    parser_assembler.add_argument("--asm_alignment_parameters", help="BLASR parameters to use to align local assemblies", default="-affineAlign -affineOpen 8 -affineExtend 0 -bestn 1 -maxMatch 30 -sdpTupleSize 13")
    parser_assembler.add_argument("--mapping_quality", type=int, help="minimum mapping quality of raw reads to use for local assembly", default=30)
    parser_assembler.add_argument("--minutes_to_delay_jobs", type=int, help="maximum number of minutes to delay local assembly jobs to limit simultaneous I/O on shared storage", default=1)
    parser_assembler.add_argument("--assembly_log", help="name of log file for local assemblies", default="assembly.log")
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

    # Run: Call SVs and indels from BLASR alignments of raw reads.
    parser_runner = subparsers.add_parser("run", help="call SVs and indels by local assembly of BLASR-aligned reads")
    parser_runner.add_argument("reference", help="FASTA file of indexed reference with .ctab and .sa in the same directory")
    parser_runner.add_argument("reads", help="text file with one absolute path to a PacBio reads file (.bax.h5) per line")
    parser_runner.add_argument("--variants", help="VCF of variants called by local assembly alignments", default="variants.vcf")
    parser_runner.add_argument("--alignments", help="text file with one absolute path to a BLASR raw reads alignments file (.bam) per line", default="alignments.fofn")
    parser_runner.add_argument("--alignments_dir", help="absolute path of directory for BLASR alignment files", default="alignments")
    parser_runner.add_argument("--candidates", help="BED file of candidates detected in read alignments", default="candidates.bed")
    parser_runner.add_argument("--assembly_alignments", help="BAM file with BLASR alignments of local assemblies against the reference", default="local_assembly_alignments.bam")
    parser_runner.add_argument("--batches", help="number of batches to split input reads into such that there will be one BAM output file per batch", type=int, default=1)
    parser_runner.add_argument("--threads", help="number of threads to use for each BLASR alignment job", type=int, default=1)
    parser_runner.add_argument("--exclude", help="BED file of regions to exclude from local assembly (e.g., heterochromatic sequences, etc.)")
    parser_runner.add_argument("--assembly_window_size", type=int, help="size of reference window for local assemblies", default=60000)
    parser_runner.add_argument("--assembly_window_slide", type=int, help="size of reference window slide for local assemblies", default=20000)
    parser_runner.add_argument("--min_length", type=int, help="minimum length required for SV candidates", default=50)
    parser_runner.add_argument("--min_support", type=int, help="minimum number of supporting reads required to flag a region as an SV candidate", default=5)
    parser_runner.add_argument("--max_support", type=int, help="maximum number of supporting reads allowed to flag a region as an SV candidate", default=100)
    parser_runner.add_argument("--min_coverage", type=int, help="minimum number of total reads required to flag a region as an SV candidate", default=5)
    parser_runner.add_argument("--max_coverage", type=int, help="maximum number of total reads allowed to flag a region as an SV candidate", default=100),
    parser_runner.add_argument("--rebuild_regions", action="store_true", help="rebuild subset of regions to assemble")
    parser_runner.add_argument("--sample", help="Sample name to use in final variant calls", default="UnnamedSample")
    parser_runner.add_argument("--species", help="Common or scientific species name to pass to RepeatMasker", default="human")
    parser_runner.add_argument("--runjobs", help="A comma-separated list of jobs for each step: align, detect, assemble, and call (in that order). A missing number uses the value set by --jobs (or 1 if --jobs was not set).", default="")
    parser_runner.add_argument("--alignment_parameters", help="BLASR parameters to use to align raw reads", default="-bestn 2 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 50")
    parser_runner.add_argument("--asm_alignment_parameters", help="BLASR parameters to use to align local assemblies", default="-affineAlign -affineOpen 8 -affineExtend 0 -bestn 1 -maxMatch 30 -sdpTupleSize 13")
    parser_runner.add_argument("--mapping_quality", type=int, help="minimum mapping quality of raw reads to use for local assembly", default=30)
    parser_runner.add_argument("--minutes_to_delay_jobs", type=int, help="maximum number of minutes to delay local assembly jobs to limit simultaneous I/O on shared storage", default=1)
    parser_runner.add_argument("--assembly_log", help="name of log file for local assemblies", default="assembly.log")
    parser_runner.add_argument("--min_hardstop_support", type=int, help="minimum number of reads with hardstops required to flag a region as an SV candidate", default=11)
    parser_runner.add_argument("--max_candidate_length", type=int, help="maximum length allowed for an SV candidate region", default=60000)
    parser_runner.set_defaults(func=run)

    # Genotype SVs with Illumina reads.
    parser_genotyper = subparsers.add_parser("genotype", help="Genotype SVs with Illumina reads")
    parser_genotyper.add_argument("genotyper_config", help="JSON configuration file with SV reference paths, samples to genotype as BAMs, and their corresponding references")
    parser_genotyper.add_argument("genotyped_variants", help="VCF of SMRT SV variant genotypes for the given sample-level BAMs")
    parser_genotyper.add_argument("--threads", help="number of threads to use for each BWA MEM alignment job", type=int, default=1)
    parser_genotyper.set_defaults(func=genotype)

    args = parser.parse_args()

    # Set DRMAA library path
    if args.drmaalib is not None:
        PROCESS_ENV["DRMAA_LIBRARY_PATH"] = args.drmaalib
    elif args.distribute and "DRMAA_LIBRARY_PATH" not in PROCESS_ENV:
        sys.stderr.write("WARNING: --distribute is set, but DRMAA_LIBRARY_PATH is not set in the environment or via the --drmaalib option: Searching only in Python's library path for libdrmaa.so\n")

    # Report paths if verbose
    if args.verbose:

        # Print python version
        print('Python version: {0}'.format(re.sub('\s*\n\s*', ' - ', sys.version)))

        # Print environment
        print("PATH:")
        for PATH_ELEMENT in PROCESS_ENV["PATH"].split(":"):
            print("\t* %s" % PATH_ELEMENT)

        print("LD_LIBRARY_PATH:")
        for PATH_ELEMENT in PROCESS_ENV["LD_LIBRARY_PATH"].split(":"):
            print("\t* %s" % PATH_ELEMENT)

        print("PERL5LIB:")
        for PATH_ELEMENT in PROCESS_ENV["PERL5LIB"].split(":"):
            print("\t* %s" % PATH_ELEMENT)

        if "DRMAA_LIBRARY_PATH" in PROCESS_ENV:
            print("DRMAA_LIBRARY_PATH: %s" % PROCESS_ENV["DRMAA_LIBRARY_PATH"])
        else:
            print("DRMAA_LIBRARY_PATH: <NOT_SET>\n\t* Not required unless --distribute is set")

        # Print arguments
        print("Arguments:")
        for key in sorted(vars(args).keys()):
            print('\t* %s = %s' % (key, getattr(args, key)))

        # Flush output
        sys.stdout.flush()

    # Make a log directory for grid-engine-style error logs if commands are
    # being distributed in non-dryrun mode.
    if args.distribute and not args.dryrun and not os.path.isdir("log"):
        os.mkdir("log")

    # Run target command
    return_code = args.func(args)
    sys.exit(return_code)
