#!/bin/env python3

import argparse
import sys
import os
import re

from smrtsvlib.args import args_dict
from smrtsvlib.args import get_arg
from smrtsvlib import smrtsvrunner


# Get SMRTSV base directory
SMRTSV_DIR = os.path.dirname(os.path.abspath(__file__))

# Modify environment and get a copy
PROCESS_ENV = smrtsvrunner.get_env(SMRTSV_DIR)


def index(args):
    
    print('Preparing reference')
    
    return smrtsvrunner.run_snake_target(
        'rules/reference.snakefile', args, PROCESS_ENV, SMRTSV_DIR,
        (
            'ref_all',
            '--config',
            'reference={}'.format(args.reference)
        ),
        cmd_log='reference/log'
    )


def align(args):

    print('Aligning sequence reads')

    return smrtsvrunner.run_snake_target(
        'rules/align.snakefile', args, PROCESS_ENV, SMRTSV_DIR,
        (
            'aln_run',
            '--config',
            'reads={}'.format(args.reads),
            'batches={}'.format(args.batches),
            'threads={}'.format(args.threads),
            'alignment_parameters="{}"'.format(args.alignment_parameters)
        ),
        cmd_log='align/log'
    )


def detect(args):
    """
    Detect SVs from signatures in read alignments.
    """

    print('Detecting variant signatures')

    command = (
        'detect_group_merge_regions',
        '--config',
        'assembly_window_size={}'.format(args.assembly_window_size),
        'assembly_window_slide={}'.format(args.assembly_window_slide),
        'mapping_quality={}'.format(args.mapping_quality),
        'min_length={}'.format(args.min_length),
        'min_support={}'.format(args.min_support),
        'max_support={}'.format(args.max_support),
        'min_coverage={}'.format(args.min_coverage),
        'max_coverage={}'.format(args.max_coverage),
        'min_hardstop_support={}'.format(args.min_hardstop_support),
        'max_candidate_length={}'.format(args.max_candidate_length),
        'candidate_group_size={}'.format(args.candidate_group_size)
    )

    if args.exclude:
        command = command + ('exclude={}'.format(args.exclude),)

    return smrtsvrunner.run_snake_target('rules/detect.snakefile', args, PROCESS_ENV, SMRTSV_DIR, command, cmd_log='detect/log')


def assemble(args):
    """
    Assemble candidate regions from raw reads aligned to regions.
    """

    print('Running local assemblies')

    # Quote alignment params
    args.asm_alignment_parameters = '"{}"'.format(args.asm_alignment_parameters)

    # Run
    command = (
        'asm_merge_groups',
        '--config',
        'asm_alignment_parameters={}'.format(args.asm_alignment_parameters),
        'mapping_quality={}'.format(args.mapping_quality),
        'asm_cpu={}'.format(args.asm_cpu),
        'asm_mem={}'.format(args.asm_mem),
        'asm_parallel={}'.format(args.asm_parallel),
        'asm_polish={}'.format(args.asm_polish),
        'no_rm_temp={}'.format(args.nt)
    )

    return smrtsvrunner.run_snake_target('rules/assemble.snakefile', args, PROCESS_ENV, SMRTSV_DIR, command, cmd_log='assemble/log')


def call(args):
    # Call SVs, indels, and inversions.
    sys.stdout.write("Calling variants\n")

    return_code = smrtsvrunner.run_snake_target(
        'rules/call.snakefile', args, PROCESS_ENV, SMRTSV_DIR,
        (
            'call_variant_vcf',
            '--config',
            'variants={}'.format(args.variants),
            'species="{}"'.format(args.species),
            'sample="{}"'.format(args.sample),
            'rmsk={}'.format(str(args.rmsk).lower())
        ),
        cmd_log='call/log'
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
        'rules/genotype.snakefile', args, PROCESS_ENV, SMRTSV_DIR,
        (
            'gt_vcf_write',
            '--config',
            'genotyper_config={}'.format(args.genotyper_config),
            'genotyped_variants={}'.format(args.genotyped_variants),
            'gt_mapq={}'.format(args.gt_mapq),
            'gt_map_cpu={}'.format(args.gt_map_cpu),
            'gt_map_mem={}'.format(args.gt_map_mem),
            'gt_map_disk={}'.format(args.gt_map_disk),
            'gt_map_time={}'.format(args.gt_map_time),
            'gt_keep_temp={}'.format(args.gt_keep_temp)
         )
    )

    if return_code != 0:
        sys.stderr.write('Failed to genotype SVs\n')

    return return_code


# Main
if __name__ == '__main__':

    # Setup Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--dryrun', **args_dict['dryrun'])
    parser.add_argument('--distribute', **args_dict['distribute'])
    parser.add_argument('--jobs', **args_dict['jobs'])
    parser.add_argument('--tempdir', **args_dict['tempdir'])
    parser.add_argument('--verbose', '-v', **args_dict['verbose'])
    parser.add_argument('--cluster-config', dest='cluster_config', **args_dict['cluster_config'])
    parser.add_argument('--cluster-params', dest='cluster_params', **args_dict['cluster_params'])
    parser.add_argument('--drmaalib', **args_dict['drmaalib'])
    parser.add_argument('--job-prefix', **args_dict['job_prefix'])
    parser.add_argument('--keep-going', '-k', dest='keep_going', **args_dict['keep_going'])
    parser.add_argument('--nt', **args_dict['nt'])
    parser.add_argument('--log', **args_dict['log'])
    parser.add_argument('--wait-time', dest='wait_time', **args_dict['wait_time'])
    subparsers = parser.add_subparsers()

    # SMRTSV command: Index reference
    parser_index = subparsers.add_parser('index', help='Index a reference for BLASR')
    parser_index.add_argument('reference', **args_dict['reference'])
    parser_index.set_defaults(func=index)

    # SMRTSV command: Align PacBio reads
    parser_align = subparsers.add_parser('align', help='Align reads.')
    parser_align.add_argument('reads', **args_dict['reads'])
    parser_align.add_argument('--batches', **args_dict['batches'])
    parser_align.add_argument('--threads', **args_dict['threads'])
    parser_align.add_argument('--alignment-parameters', dest='alignment_parameters', **args_dict['alignment_parameters'])
    parser_align.set_defaults(func=align)

    # SMRTSV command: Detect windows
    parser_detector = subparsers.add_parser('detect', help='Detect SV signatures in aligned reads')
    parser_detector.add_argument('--mapping-quality', dest='mapping_quality', **args_dict['mapping_quality'])
    parser_detector.add_argument('--exclude', **args_dict['exclude'])
    parser_detector.add_argument('--assembly-window-size', dest='assembly_window_size', **args_dict['assembly_window_size'])
    parser_detector.add_argument('--assembly-window-slide', dest='assembly_window_slide', **args_dict['assembly_window_slide'])
    parser_detector.add_argument('--min-length', dest='min_length', **args_dict['min_length'])
    parser_detector.add_argument('--min-support', dest='min_support', **args_dict['min_support'])
    parser_detector.add_argument('--max-support', dest='max_support', **args_dict['max_support'])
    parser_detector.add_argument('--min-coverage', dest='min_coverage', **args_dict['min_coverage'])
    parser_detector.add_argument('--max-coverage', dest='max_coverage', **args_dict['max_coverage'])
    parser_detector.add_argument('--min-hardstop-support', dest='min_hardstop_support', **args_dict['min_hardstop_support'])
    parser_detector.add_argument('--max-candidate-length', dest='max_candidate_length', **args_dict['max_candidate_length'])
    parser_detector.add_argument('--candidate-group-size', dest='candidate_group_size', **args_dict['candidate_group_size'])
    parser_detector.set_defaults(func=detect)

    # SMRTSV command: Assemble regions
    parser_assembler = subparsers.add_parser('assemble', help='Assemble candidate regions.')
    parser_assembler.add_argument('--asm-alignment-parameters', dest='asm_alignment_parameters', **args_dict['asm_alignment_parameters'])
    parser_assembler.add_argument('--mapping-quality', dest='mapping_quality', **args_dict['mapping_quality'])
    parser_assembler.add_argument('--asm-cpu', dest='asm_cpu', **args_dict['asm_cpu'])
    parser_assembler.add_argument('--asm-mem', dest='asm_mem', **args_dict['asm_mem'])
    parser_assembler.add_argument('--asm-polish', dest='asm_polish', **args_dict['asm_polish'])
    parser_assembler.add_argument('--asm-group-rt', dest='asm_group_rt', **args_dict['asm_group_rt'])
    parser_assembler.add_argument('--asm-rt', dest='asm_rt', **args_dict['asm_rt'])
    parser_assembler.add_argument('--asm-parallel', dest='asm_parallel', **args_dict['asm_parallel'])
    parser_assembler.set_defaults(func=assemble)

    # SMRTSV command: Call variants
    parser_caller = subparsers.add_parser('call', help='Call variants from assemblies.')
    parser_caller.add_argument('variants', **args_dict['variants'])
    parser_caller.add_argument('--sample', **args_dict['sample'])
    parser_caller.add_argument('--species', **args_dict['species'])
    parser_caller.add_argument('--rmsk', **args_dict['rmsk'])
    parser_caller.set_defaults(func=call)

    # SMRTSV command: Run all steps of the variant calling pipeline from raw reads
    parser_runner = subparsers.add_parser('run', help='call SVs and indels by local assembly of BLASR-aligned reads')
    parser_runner.add_argument('reference', **args_dict['reference'])
    parser_runner.add_argument('reads', **args_dict['reads'])
    parser_runner.add_argument('--variants', **args_dict['variants'])
    parser_runner.add_argument('--batches', **args_dict['batches'])
    parser_runner.add_argument('--threads', **args_dict['threads'])
    parser_runner.add_argument('--exclude', **args_dict['exclude'])
    parser_runner.add_argument('--assembly-window-size', dest='assembly_window_size', **args_dict['assembly_window_size'])
    parser_runner.add_argument('--assembly-window-slide', dest='assembly_window_slide', **args_dict['assembly_window_slide'])
    parser_runner.add_argument('--asm-cpu', dest='asm_cpu', **args_dict['asm_cpu'])
    parser_runner.add_argument('--asm-mem', dest='asm_mem', **args_dict['asm_mem'])
    parser_runner.add_argument('--asm-polish', dest='asm_polish', **args_dict['asm_polish'])
    parser_runner.add_argument('--min-length', dest='min_length', **args_dict['min_length'])
    parser_runner.add_argument('--min-support', dest='min_support', **args_dict['min_support'])
    parser_runner.add_argument('--max-support', dest='max_support', **args_dict['max_support'])
    parser_runner.add_argument('--min-coverage', dest='min_coverage', **args_dict['min_coverage'])
    parser_runner.add_argument('--max-coverage', dest='max_coverage', **args_dict['max_coverage']),
    parser_runner.add_argument('--sample', **args_dict['sample'])
    parser_runner.add_argument('--species', **args_dict['species'])
    parser_runner.add_argument('--runjobs', **args_dict['runjobs'])
    parser_runner.add_argument('--alignment-parameters', dest='alignment_parameters', **args_dict['alignment_parameters'])
    parser_runner.add_argument('--asm-alignment-parameters', dest='asm_alignment_parameters', **args_dict['asm_alignment_parameters'])
    parser_runner.add_argument('--mapping-quality', dest='mapping_quality', **args_dict['mapping_quality'])
    parser_runner.add_argument('--min-hardstop-support', dest='min_hardstop_support', **args_dict['min_hardstop_support'])
    parser_runner.add_argument('--max-candidate-length', dest='max_candidate_length', **args_dict['max_candidate_length'])
    parser_runner.add_argument('--candidate-group-size', dest='candidate_group_size', **args_dict['candidate_group_size'])
    parser_runner.add_argument('--asm-group-rt', dest='asm_group_rt', **args_dict['asm_group_rt'])
    parser_runner.add_argument('--asm-rt', dest='asm_rt', **args_dict['asm_rt'])
    parser_runner.add_argument('--asm-parallel', dest='asm_parallel', **args_dict['asm_parallel'])
    parser_runner.add_argument('--rmsk', **args_dict['rmsk'])
    parser_runner.set_defaults(func=run)

    # SMRTSV command: Genotype SMRTSV SV calls with Illumina reads
    parser_genotyper = subparsers.add_parser('genotype', help='Genotype SVs with Illumina reads')
    parser_genotyper.add_argument('genotyper_config', **args_dict['genotyper_config'])
    parser_genotyper.add_argument('genotyped_variants', **args_dict['genotyped_variants'])
    parser_genotyper.add_argument('--gt-mapq', '--mapq', dest='gt_mapq', **args_dict['gt_mapq'])
    parser_genotyper.add_argument('--gt-map-cpu', dest='gt_map_cpu', **args_dict['gt_map_cpu'])
    parser_genotyper.add_argument('--gt-map-mem', dest='gt_map_mem', **args_dict['gt_map_mem'])
    parser_genotyper.add_argument('--gt-map-disk', dest='gt_map_disk', **args_dict['gt_map_disk'])
    parser_genotyper.add_argument('--gt-map-time', dest='gt_map_time', **args_dict['gt_map_time'])
    parser_genotyper.add_argument('--gt-keep-temp', dest='gt_keep_temp', **args_dict['gt_keep_temp'])
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

    # Check target
    if 'func' not in cmd_args:
        print('Unrecognized SMRT-SV command or command-line options could not be parsed. Try "-h" for help.', file=sys.stderr)
        sys.exit(1)

    # Run target command
    cmd_return_code = cmd_args.func(cmd_args)
    sys.exit(cmd_return_code)
