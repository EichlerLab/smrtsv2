"""
Coordinates assembly of each group of regions and merges results from each group.
"""

import pandas as pd
import shutil
import subprocess
import socket

if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'include.snakefile'

from smrtsvlib import smrtsvrunner
from smrtsvlib import smrtsvutil


###################
### Definitions ###
###################

def _get_group_bams(wildcards):
    """
    Return a list of all BAM files produced by running local assemblies on each group. There is one BAM file per
    group, and each are output by a call to rule "asm_assemble_group".

    :param wildcards: Ignored

    :return: A list of BAM files output by local assemblies.
    """
    return [
        'assemble/group/{}/contig.bam'.format(group_id)
        for group_id in
            pd.read_table('detect/candidate_groups.bed', header=0, usecols=('GROUP_ID', ), squeeze=True).tolist()
    ]


#############
### Rules ###
#############

# asm_merge_assembled_groups
#
# Merge local assemblies into one BAM.
rule asm_merge_assembled_groups:
    input:
        bam=_get_group_bams
    output:
        bam='assemble/local_assemblies.bam',
        bai='assemble/local_assemblies.bam.bai'
    log:
        'assemble/log/local_assemblies.merge.log'
    run:

        # Merge BAMs from each group into one
        if len(input.bam) > 1:
            shell("""samtools merge {output.bam} {input.bam} >{log} 2>&1""")
        else:
            # Copy if there was only one group
            shell("""cp {input.bam} {output.bam}""")

        # Index
        shell("""samtools index {output.bam}""")

# asm_assemble_group
#
# Assemble one group of regions.
rule asm_assemble_group:
    input:
        align_fofn='align/alignments.fofn',
        bed_grp='detect/candidate_groups.bed',
        bed_can='detect/candidates.bed',
        ref_fa='reference/ref.fasta',
        ref_fai='reference/ref.fasta.fai',
        ref_sa='reference/ref.fasta.sa',
        ref_ctab='reference/ref.fasta.ctab'
    output:
        bam='assemble/group/{group_id}/contig.bam',
        bai='assemble/group/{group_id}/contig.bam.bai'
    params:
        mapq=get_config_param('mapping_quality'),
        align_params=get_config_param('asm_alignment_parameters'),
        threads=get_config_param('asm_cpu'),  # Parses into cluster params
        mem=get_config_param('asm_mem'),      # Parses into cluster params
        asm_polish=get_config_param('asm_polish'),
        no_rm_temp=get_config_param('no_rm_temp'),
    log:
        contig_group='assemble/group/{group_id}/contig_group.log'
    benchmark:
        'assemble/group/{group_id}/contig_group_bm.log'
    run:
        # Note: Runs assemble_group snakemake with resources for "threads". This allows single-core steps to be
        # run in parallel, but runs multi-core steps to one at a time.

        # Set assemble_temp (will be deleted if not None)
        assemble_temp = None

        try:

            # Create temporary directories
            assemble_temp = os.path.join(TEMP_DIR, 'asm_group_{}'.format(wildcards.group_id))

            os.makedirs(assemble_temp, exist_ok=True)

            # Create log directory
            log_dir = os.path.abspath(os.path.join(os.path.dirname(log.contig_group), 'log'))
            os.makedirs(log_dir, exist_ok=True)

            # Setup sub-Snake command
            command = [
                'merge_group_contigs',
                '-f',
                '--jobs',
                str(params.threads)
            ]

            if smrtsvutil.as_bool(params.no_rm_temp):
                command += ['--nt']

            command += [
                '--config',
                'contig_out={}'.format(os.path.abspath(output.bam)),
                'log_dir={}'.format(log_dir),
                'mapping_quality={}'.format(params.mapq),
                'asm_alignment_parameters={}'.format(params.align_params),  # String is already quoted to prevent Snakemake from trying to interpret the command options
                'asm_polish={}'.format(params.asm_polish),
                'group_id={}'.format(wildcards.group_id),
                'align_fofn={}'.format(os.path.abspath(input.align_fofn)),
                'bed_groups={}'.format(os.path.abspath(input.bed_grp)),
                'bed_candidates={}'.format(os.path.abspath(input.bed_can)),
                'threads={:d}'.format(params.threads),
                'ref_fa={}'.format(os.path.abspath(input.ref_fa))
            ]

            # Clear locks
            #shell("""rm -f {working_dir}/.snakemake/locks/*""")

            # Run assembly
            with open(log.contig_group, 'w') as log_file:
                log_file.write('Assembling group: {}\n'.format(wildcards.group_id))
                log_file.write('Assemble temp: {}\n'.format(assemble_temp))
                log_file.write('Hostname: {}\n'.format(socket.gethostname()))
                log_file.write('Threads: {}\n'.format(params.threads))
                log_file.write('Memory: {}\n'.format(params.mem))
                log_file.write('Polish: {}\n'.format(params.asm_polish))
                log_file.write('MAPQ: {}\n\n'.format(params.mapq))
                log_file.flush()

                return_code = smrtsvrunner.run_snake_target(
                    'rules/assemble_group.snakefile', None, PROCESS_ENV, SMRTSV_DIR, command,
                    stdout=log_file, stderr=subprocess.STDOUT, cwd=assemble_temp,
                    resources=['threads={:d}'.format(params.threads)]
                )

            if return_code != 0:
                raise RuntimeError('Failed to assemble group {}: See log {}'.format(wildcards.group_id, log.contig_group))

        finally:

            # Clean up temp directory
            if assemble_temp is not None and not params.no_rm_temp:
                shutil.rmtree(assemble_temp)
