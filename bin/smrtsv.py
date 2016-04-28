#!/bin/env python
import argparse
import logging
import subprocess
import sys
import os

logging.basicConfig(filename="smrtsv.log", level=logging.DEBUG)
CLUSTER_SETTINGS = '" -V -cwd -e ./log -o ./log {params.sge_opts} -w n -S /bin/bash"'
CLUSTER_FLAG = ("--drmaa", CLUSTER_SETTINGS, "-w", "60")

def _get_dist_dir():
    dirname, filename = os.path.split(os.path.abspath(__file__))
    return dirname

def _build_prefix(args):
    prefix = ["snakemake", "--snakefile", os.path.join(os.path.dirname(_get_dist_dir()), "Snakefile"), "-j", str(args.jobs)]
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
        "threads=%s" % args.threads,
        "tmp_dir=%s" % args.tmpdir
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

    base_command = prefix + (
        "collect_assembly_alignments",
        "--config",
        "reference=%s" % args.reference,
        "alignments=%s" % args.alignments,
        "reads=%s" % args.reads,
        "tmp_dir=%s" % args.tmpdir
    )

    if args.regions:
        # For each contig/chromosome in the regions file, submit a separate
        # Snakemake command. To do so, first split regions to assemble into one
        # file per contig in a temporary directory.
        tmpdir = os.path.join(os.getcwd(), "regions_by_contig")

        rebuild_regions_by_contig = False
        if not args.dryrun and (not os.path.exists(tmpdir) or os.stat(args.regions).st_mtime > os.stat(tmpdir).st_mtime):
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
            command = base_command + ("regions_to_assemble=%s" % os.path.join(tmpdir, "%s.bed" % contig),)
            command = command + ("assembly_alignments=%s" % contig_local_assemblies,)
            sys.stdout.write("Starting local assemblies for %s\n" % contig)
            logging.debug("Assembly command: %s", " ".join(command))
            return_code = subprocess.call(" ".join(command), shell=True)
            if return_code != 0:
               break

        # If the last command executed successfully, try to merge all local
        # assemblies per contig into a single file.
        if not args.dryrun and return_code == 0:
            return_code = subprocess.call(" ".join (["samtools", "merge", args.assembly_alignments] + list(local_assemblies)), shell=True)

            if return_code == 0:
                return_code = subprocess.call(["samtools", "index", args.assembly_alignments])

        # Return the last return code.
        return return_code
    else:
        if args.assembly_alignments:
            command = base_command + ("assembly_alignments=%s" % args.assembly_alignments,)

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
        "variants=%s" % args.variants,
        "species=%s" % args.species,
        "sample=%s" % args.sample
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
    parser.add_argument("--tmpdir", help="temporary directory to use for distributed jobs", default="/var/tmp")
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
    parser_caller.add_argument("--sample", help="Sample name to use in final variant calls", default="UnnamedSample")
    parser_caller.add_argument("--species", help="Species name to use for repeat masking", default="Homo sapiens")
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

    # Make a log directory for grid-engine-style error logs if commands are
    # being distributed in non-dryrun mode.
    if args.distribute and not args.dryrun and not os.path.isdir("log"):
        os.mkdir("log")

    return_code = args.func(args)
    sys.exit(return_code)
