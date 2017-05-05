#!/bin/env python
import argparse
import logging
import sys
import os
import re

from smrtsvlib.args import args_dict
from smrtsvlib import smrtsvrunner

# Set logging
logging.basicConfig(filename='smrtsv.log', level=logging.DEBUG)

# Set cluster parameters
CLUSTER_SETTINGS = ' -V -cwd -e ./log -o ./log {cluster.params} -w n -S /bin/bash'
CLUSTER_FLAG = ('--drmaa', CLUSTER_SETTINGS, '-w', '120')

# Get SMRTSV base directory
SMRTSV_DIR = os.path.dirname(os.path.abspath(__file__))

# Modify environment and get a copy
PROCESS_ENV = smrtsvrunner.get_env(SMRTSV_DIR)


def index(args):
    
    print('Indexing reference')
    
    return smrtsvrunner.run_snake_target(
        args,
        'prepare_reference',
        '--config',
        'reference={}'.format(args.reference)
    )


def align(args):

    print('Aligning sequence reads')

    return smrtsvrunner.run_snake_target(
        args, PROCESS_ENV, SMRTSV_DIR, CLUSTER_FLAG,
        'align_reads',
        '--config',
        'reference={}'.format(args.reference),
        'reads={}'.format(args.reads),
        'alignments={}'.format(args.alignments),
        'alignments_dir={}'.format(args.alignments_dir),
        'batches={}'.format(args.batches),
        'threads={}'.format(args.threads),
        'tmp_dir={}'.format(args.tmpdir),
        'alignment_parameters="{}"'.format(args.alignment_parameters)
    )


def detect(args):
    """
    Detect SVs from signatures in read alignments.
    """

    print('Searching for variant signatures')

    command = (
        'get_regions',
        '--config',
        'reference={}'.format(args.reference),
        'alignments={}'.format(args.alignments),
        'assembly_window_size={}'.format(args.assembly_window_size),
        'assembly_window_slide={}'.format(args.assembly_window_slide),
        'min_length={}'.format(args.min_length),
        'min_support={}'.format(args.min_support),
        'max_support={}'.format(args.max_support),
        'min_coverage={}'.format(args.min_coverage),
        'max_coverage={}'.format(args.max_coverage),
        'min_hardstop_support={}'.format(args.min_hardstop_support),
        'max_candidate_length={}'.format(args.max_candidate_length)
    )

    if args.exclude:
        command = command + ('regions_to_exclude={}'.format(args.exclude),)

    if args.candidates:
        command = command + ('candidates={}'.format(args.candidates),)

    return smrtsvrunner.run_snake_target(args, PROCESS_ENV, SMRTSV_DIR, CLUSTER_FLAG, *command)


def assemble(args):
    """
    Assemble candidate regions from raw reads aligned to regions.
    """

    print('Starting local assemblies')

    base_command = (
        'collect_assembly_alignments',
        '--config',
        'reference={}'.format(args.reference),
        'alignments={}'.format(args.alignments),
        'reads={}'.format(args.reads),
        'tmp_dir={}'.format(args.tmpdir),
        'asm_alignment_parameters="{}"'.format(args.asm_alignment_parameters),
        'mapping_quality="{}"'.format(args.mapping_quality),
        'minutes_to_delay_jobs="{}"'.format(args.minutes_to_delay_jobs),
        'assembly_log="{}"'.format(args.assembly_log)
    )

    if args.candidates:
        # For each contig/chromosome in the candidates file, submit a separate
        # Snakemake command. To do so, first split regions to assemble into one
        # file per contig in a temporary directory.
        tmpdir = os.path.join(os.getcwd(), 'regions_by_contig')

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
                        contig_file = open(os.path.join(tmpdir, '{}.bed'.format(contig)), 'w')

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
            contig_local_assemblies = os.path.join(
                'local_assemblies',
                local_assembly_basename.replace('.bam', '.{}.bam'.format(contig))
            )

            local_assemblies.add(contig_local_assemblies)

            if os.path.exists(contig_local_assemblies):
                sys.stdout.write('Local assemblies already exist for %s\n' % contig)
                continue

            command = base_command + ('regions_to_assemble=%s' % os.path.join(tmpdir, '%s.bed' % contig),)
            command = command + ('assembly_alignments={}'.format(contig_local_assemblies),)
            sys.stdout.write('Starting local assemblies for {}\n'.format(contig))
            logging.debug('Assembly command: %s', ' '.join(command))

            return_code = smrtsvrunner.run_snake_target(args, PROCESS_ENV, SMRTSV_DIR, CLUSTER_FLAG, *command)

            contig_count += 1

            if return_code != 0:
                break

        # If the last command executed successfully, try to merge all local
        # assemblies per contig into a single file. Only build if at least one set of local assemblies was performed
        # or the local assemblies file does not exist.
        if not args.dryrun and return_code == 0 and (contig_count > 0 or not os.path.exists(args.assembly_alignments)):
            if len(local_assemblies) > 1:
                return_code = smrtsvrunner.run_cmd(['samtools', 'merge', args.assembly_alignments] +
                                                   list(local_assemblies), PROCESS_ENV)
            else:
                return_code = smrtsvrunner.run_cmd(['samtools', 'view', '-b', '-o', args.assembly_alignments] +
                                                   list(local_assemblies), PROCESS_ENV)

            if return_code == 0:
                return_code = smrtsvrunner.run_cmd(['samtools', 'index', args.assembly_alignments], PROCESS_ENV)

        # Return the last return code.
        return return_code
    else:
        if args.assembly_alignments:
            command = base_command + ('assembly_alignments={}'.format(args.assembly_alignments),)

            logging.debug('Assembly command: %s', ' '.join(command))
            return smrtsvrunner.run_cmd(command, PROCESS_ENV)


def call(args):
    # Call SVs, indels, and inversions.
    sys.stdout.write("Calling variants\n")

    return_code = smrtsvrunner.run_snake_target(
        args, PROCESS_ENV, SMRTSV_DIR, CLUSTER_FLAG,
        'call_variants',
        '--config',
        'reference={}'.format(args.reference),
        'alignments={}'.format(args.alignments),
        'local_assembly_alignments={}'.format(args.assembly_alignments),
        'variants={}'.format(args.variants),
        'species="{}"'.format(args.species),
        'sample="{}""'.format(args.sample)
    )

    if return_code != 0:
        sys.stderr.write('Failed to call variants\n')

    return return_code


def run(args):
    """
    Run all steps of the variant calling pipeline from raw reads.

    :param args: Command arguments.

    :return: 0 if the command ran successfully, and a non-zero code otherwise.
    """

    # Get default jobs
    if 'jobs' in args:
        default_jobs = args.jobs
    else:
        default_jobs = 1

    # Get the number of jobs for each step
    job_step = re.split('\\s*[,;:]\\s*', args.runjobs.strip())  # Split into array
    job_step = [job_step[i] if len(job_step) > i else '' for i in range(4)]  # Extend to length 4

    # Convert each number of jobs to integers
    for i in range(4):

        if job_step[i] != '':
            try:
                job_step[i] = int(job_step[i])

            except ValueError:
                sys.stderr.write(
                    'Invalid number of jobs for step {}: Must be an integer: "{}"\n'.format((i + 1), job_step[i])
                )

                return 1
        else:
            job_step[i] = default_jobs

    # Report the number of jobs for each task
    if args.verbose and args.distribute:
        print('Jobs per task:')
        print('\t*    Align: {}'.format(job_step[0]))
        print('\t*   Detect: {}'.format(job_step[1]))
        print('\t* Assemble: {}'.format(job_step[2]))
        print('\t*     Call: {}'.format(job_step[3]))

    # Build reference indices
    return_code = index(args)
    if return_code != 0:
        sys.stderr.write('Failed to index reference\n')
        return return_code

    # Align
    args.jobs = job_step[0]

    return_code = align(args)
    if return_code != 0:
        sys.stderr.write('Failed to align reads\n')
        return return_code

    # Detect SVs.
    args.jobs = job_step[1]

    return_code = detect(args)
    if return_code != 0:
        sys.stderr.write('Failed to identify candidate regions\n')
        return return_code

    # Run local assemblies.
    args.jobs = job_step[2]

    return_code = assemble(args)
    if return_code != 0:
        sys.stderr.write('Failed to generate local assemblies\n')
        return return_code

    # Call SVs, indels, and inversions.
    args.jobs = job_step[3]

    return_code = call(args)
    if return_code != 0:
        sys.stderr.write('Failed to call variants\n')
        return return_code

    return 0


def genotype(args):
    """
    Genotype structural variants.

    :param args: Command arguments.

    :return: 0 if the command ran successfully, and a non-zero code otherwise.
    """

    print('Genotyping SVs')

    return_code = smrtsvrunner.run_snake_target(
        args, PROCESS_ENV, SMRTSV_DIR, CLUSTER_FLAG,
        'convert_genotypes_to_vcf',
        '--config',
        'genotyper_config={}'.format(args.genotyper_config),
        'genotyped_variants={}'.format(args.genotyped_variants),
        'threads={}'.format(args.threads)
    )

    if return_code != 0:
        sys.stderr.write('Failed to genotype SVs\n')

    return return_code


# Main
if __name__ == '__main__':

    # Setup Parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--dryrun', '-n', action='store_true',
                        help='Print commands that will run without running them')

    parser.add_argument('--distribute', action='store_true',
                        help='Distribute analysis to Grid Engine-style cluster')

    parser.add_argument('--jobs', type=int, default=1,
                        help='number of jobs to run simultaneously')

    parser.add_argument('--tmpdir', default='default',
                        help='temporary directory to use for distributed jobs. Keyword "default" tells Python to pick '
                             'a directory.'
                        )

    parser.add_argument('--verbose', '-v', action='store_true',
                        help='print extra runtime information')

    parser.add_argument('--cluster_config',
                        help='JSON/YAML file specifying cluster configuration parameters to pass to Snakemake\'s'
                             ' --cluster-config option'
                        )

    parser.add_argument('--drmaalib',
                        help='For jobs that are distributed, this is the location to the DRMAA library (libdrmaa.so) '
                             'installed with Grid Engine. If DRMAA_LIBRARY_PATH is already set in the environment, '
                             'then this option is not required.'
                        )
    subparsers = parser.add_subparsers()

    # SMRTSV command: Index reference
    parser_index = subparsers.add_parser('index', help='Index a reference for BLASR')
    parser_index.add_argument('reference', **args_dict['reference'])
    parser_index.set_defaults(func=index)

    # SMRTSV command: Align PacBio reads
    parser_align = subparsers.add_parser('align', help='Align reads.')
    parser_align.add_argument('reference', **args_dict['reference'])
    parser_align.add_argument('reads', **args_dict['reads'])
    parser_align.add_argument('--alignments', **args_dict['alignments'])
    parser_align.add_argument('--alignments_dir', **args_dict['alignments_dir'])
    parser_align.add_argument('--batches', **args_dict['batches'])
    parser_align.add_argument('--threads', **args_dict['threads'])
    parser_align.add_argument('--alignment_parameters', **args_dict['alignment_parameters'])
    parser_align.set_defaults(func=align)

    # SMRTSV command: Detect windows
    parser_detector = subparsers.add_parser('detect', help='Detect SV signatures in aligned reads')
    parser_detector.add_argument('reference', **args_dict['reference'])
    parser_detector.add_argument('alignments', **args_dict['alignments'])
    parser_detector.add_argument('candidates', **args_dict['candidates'])
    parser_detector.add_argument('--exclude', **args_dict['exclude'])
    parser_detector.add_argument('--assembly_window_size', **args_dict['assembly_window_size'])
    parser_detector.add_argument('--assembly_window_slide', **args_dict['assembly_window_slide'])
    parser_detector.add_argument('--min_length', **args_dict['min_length'])
    parser_detector.add_argument('--min_support', **args_dict['min_support'])
    parser_detector.add_argument('--max_support', **args_dict['max_support'])
    parser_detector.add_argument('--min_coverage', **args_dict['min_coverage'])
    parser_detector.add_argument('--max_coverage', **args_dict['max_coverage'])
    parser_detector.add_argument('--min_hardstop_support', **args_dict['min_hardstop_support'])
    parser_detector.add_argument('--max_candidate_length', **args_dict['max_candidate_length'])
    parser_detector.set_defaults(func=detect)

    # SMRTSV command: Assemble regions
    parser_assembler = subparsers.add_parser('assemble', help='Assemble candidate regions.')
    parser_assembler.add_argument('reference', **args_dict['reference'])
    parser_assembler.add_argument('reads', **args_dict['reads'])
    parser_assembler.add_argument('alignments', **args_dict['alignments'])
    parser_assembler.add_argument('candidates', **args_dict['candidates'])
    parser_assembler.add_argument('assembly_alignments', **args_dict['assembly_alignments'])
    parser_assembler.add_argument('--rebuild_regions', **args_dict['rebuild_regions'])
    parser_assembler.add_argument('--asm_alignment_parameters', **args_dict['asm_alignment_parameters'])
    parser_assembler.add_argument('--mapping_quality', **args_dict['mapping_quality'])
    parser_assembler.add_argument('--minutes_to_delay_jobs', **args_dict['minutes_to_delay_jobs'])
    parser_assembler.add_argument('--assembly_log', **args_dict['assembly_log'])
    parser_assembler.set_defaults(func=assemble)

    # SMRTSV command: Call variants
    parser_caller = subparsers.add_parser('call', help='Call variants from assemblies.')
    parser_caller.add_argument('reference', **args_dict['reference'])
    parser_caller.add_argument('alignments', **args_dict['alignments'])
    parser_caller.add_argument('assembly_alignments', **args_dict['assembly_alignments'])
    parser_caller.add_argument('variants', **args_dict['variants'])
    parser_caller.add_argument('--sample', **args_dict['sample'])
    parser_caller.add_argument('--species', **args_dict['species'])
    parser_caller.set_defaults(func=call)

    # SMRTSV command: Run all steps of the variant calling pipeline from raw reads
    parser_runner = subparsers.add_parser('run', help='call SVs and indels by local assembly of BLASR-aligned reads')
    parser_runner.add_argument('reference', **args_dict['reference'])
    parser_runner.add_argument('reads', **args_dict['reads'])
    parser_runner.add_argument('--variants', **args_dict['variants'])
    parser_runner.add_argument('--alignments', **args_dict['alignments'])
    parser_runner.add_argument('--alignments_dir', **args_dict['alignments_dir'])
    parser_runner.add_argument('--candidates', **args_dict['candidates'])
    parser_runner.add_argument('--assembly_alignments', **args_dict['assembly_alignments'])
    parser_runner.add_argument('--batches', **args_dict['batches'])
    parser_runner.add_argument('--threads', **args_dict['threads'])
    parser_runner.add_argument('--exclude', **args_dict['exclude'])
    parser_runner.add_argument('--assembly_window_size', **args_dict['assembly_window_size'])
    parser_runner.add_argument('--assembly_window_slide', **args_dict['assembly_window_slide'])
    parser_runner.add_argument('--min_length', **args_dict['min_length'])
    parser_runner.add_argument('--min_support', **args_dict['min_support'])
    parser_runner.add_argument('--max_support', **args_dict['max_support'])
    parser_runner.add_argument('--min_coverage', **args_dict['min_coverage'])
    parser_runner.add_argument('--max_coverage', **args_dict['max_coverage']),
    parser_runner.add_argument('--rebuild_regions', **args_dict['rebuild_regions'])
    parser_runner.add_argument('--sample', **args_dict['sample'])
    parser_runner.add_argument('--species', **args_dict['species'])
    parser_runner.add_argument('--runjobs', **args_dict['runjobs'])
    parser_runner.add_argument('--alignment_parameters', **args_dict['alignment_parameters'])
    parser_runner.add_argument('--asm_alignment_parameters', **args_dict['asm_alignment_parameters'])
    parser_runner.add_argument('--mapping_quality', **args_dict['mapping_quality'])
    parser_runner.add_argument('--minutes_to_delay_jobs', **args_dict['minutes_to_delay_jobs'])
    parser_runner.add_argument('--assembly_log', **args_dict['assembly_log'])
    parser_runner.add_argument('--min_hardstop_support', **args_dict['min_hardstop_support'])
    parser_runner.add_argument('--max_candidate_length', **args_dict['max_candidate_length'])
    parser_runner.set_defaults(func=run)

    # SMRTSV command: Genotype SMRTSV SV calls with Illumina reads
    parser_genotyper = subparsers.add_parser('genotype', help='Genotype SVs with Illumina reads')
    parser_genotyper.add_argument('genotyper_config', **args_dict['genotyper_config'])
    parser_genotyper.add_argument('genotyped_variants', **args_dict['genotyped_variants'])
    parser_genotyper.add_argument('--threads', **args_dict['threads'])
    parser_genotyper.set_defaults(func=genotype)

    cmd_args = parser.parse_args()

    # Set DRMAA library path
    if cmd_args.drmaalib is not None:
        PROCESS_ENV['DRMAA_LIBRARY_PATH'] = cmd_args.drmaalib

    elif cmd_args.distribute and 'DRMAA_LIBRARY_PATH' not in PROCESS_ENV:
        sys.stderr.write(
            'WARNING: --distribute is set, but DRMAA_LIBRARY_PATH is not set in the environment or via the '
            '--drmaalib option: Searching only in Python\'s library path for libdrmaa.so\n'
        )

    # Report paths if verbose
    if cmd_args.verbose:

        # Print python version
        print('Python version: {0}'.format(re.sub('\s*\n\s*', ' - ', sys.version)))

        # Print environment
        print('PATH:')
        for PATH_ELEMENT in PROCESS_ENV['PATH'].split(':'):
            print('\t* {}'.format(PATH_ELEMENT))

        print('LD_LIBRARY_PATH:')
        for PATH_ELEMENT in PROCESS_ENV['LD_LIBRARY_PATH'].split(':'):
            print('\t* {}'.format(PATH_ELEMENT))

        print('PERL5LIB:')
        for PATH_ELEMENT in PROCESS_ENV['PERL5LIB'].split(':'):
            print('\t* {}'.format(PATH_ELEMENT))

        if 'DRMAA_LIBRARY_PATH' in PROCESS_ENV:
            print('DRMAA_LIBRARY_PATH: {}'.format(PROCESS_ENV['DRMAA_LIBRARY_PATH']))
        else:
            print('DRMAA_LIBRARY_PATH: <NOT_SET>\n\t* Not required unless --distribute is set')

        # Print arguments
        print('Arguments:')
        for key in sorted(vars(cmd_args).keys()):
            print('\t* {} = {}'.format(key, getattr(cmd_args, key)))

        # Flush output
        sys.stdout.flush()

    # Make log directory for distributed jobs
    if cmd_args.distribute and not cmd_args.dryrun and not os.path.isdir('log'):
        os.mkdir('log')

    # Run target command
    cmd_return_code = cmd_args.func(cmd_args)
    sys.exit(cmd_return_code)
