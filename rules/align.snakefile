"""
Align sequence reads in batches and produce an FOFN file listing each BAM file.
"""

import os

import numpy as np

if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'include.snakefile'

localrules: aln_run, aln_assign_batches


###################
### Definitions ###
###################

# List of batch numbers
ALIGN_BATCH_LIST = list(range(int(get_config_param('batches'))))


#############
### Rules ###
#############

# aln_run
#
# Run all alignments in batches and write and FOFN file of each BAM (one per batch).
rule aln_run:
    input:
        bam=expand('align/bam/{batch_id}.bam', batch_id=ALIGN_BATCH_LIST),
        bai=expand('align/bam/{batch_id}.bam.bai', batch_id=ALIGN_BATCH_LIST),
    output:
        fofn='align/alignments.fofn'
    run:
        with open(output.fofn, 'w') as out_file:
            out_file.write('\n'.join([os.path.abspath(bam_file) for bam_file in input.bam]))
            out_file.write('\n')

# aln_align_batch
#
# Align one batch of reads to the reference.
rule aln_align_batch:
    input:
        reads='align/batches/{batch_id}.fofn',
        ref_fa='reference/ref.fasta',
        ref_fai='reference/ref.fasta.fai',
        ref_sa='reference/ref.fasta.sa'
    output:
        bam=protected('align/bam/{batch_id}.bam'),
        bai=protected('align/bam/{batch_id}.bam.bai'),
        unaligned=protected('align/bam/unaligned/{batch_id}.fa.gz')
    params:
        threads=get_config_param('threads'),
        samtools_threads="1",
        samtools_memory="4G",
        bwlimit="30000",
        align_params=config.get('alignment_parameters', '').strip('"')
    log:
        'align/bam/log/{batch_id}.log'
    run:

        # Get unaligned FASTA file name before compression
        unaligned_file_name = output.unaligned[:-3]  # Removed .gz extension

        # Detect ctab
        ctab_file_name = '{}.ctab'.format(input.ref_fa)

        if os.path.isfile(ctab_file_name):
            ctab_op = '--ctab {input.ref_ctab} '.format(ctab_file_name)
        else:
            ctab_op = ''

        # Map
        if os.path.getsize(input.reads) > 0:

            shell(
                """ALIGN_BATCH_TEMP={TEMP_DIR}/aln_align_batch_{wildcards.batch_id}; """
                """mkdir -p ${{ALIGN_BATCH_TEMP}}; """
                """echo "Aligning batch: {wildcards.batch_id}" >{log}; """
                """echo "Temp dir: ${{ALIGN_BATCH_TEMP}}" >> {log}; """
                """blasr {input.reads} {input.ref_fa} """
                    """--unaligned {unaligned_file_name} """
                    """--out ${{ALIGN_BATCH_TEMP}}/batch_out.bam """
                    """--sam """
                    """--sa {input.ref_sa} """
                    """{ctab_op}"""
                    """--nproc {params.threads} """
                    """--clipping subread """
                    """{params.align_params} """
                    """>>{log} 2>&1; """
                """echo "Sorting..." >>{log}; """
                """samtools sort """
                    """-@ {params.samtools_threads} """
                    """-m {params.samtools_memory} """
                    """-O bam """
                    """-T ${{ALIGN_BATCH_TEMP}}/{wildcards.batch_id} """
                    """-o {output.bam} """
                    """${{ALIGN_BATCH_TEMP}}/batch_out.bam """
                    """>>{log} 2>&1; """
                """echo "Indexing..." >>{log} 2>&1; """
                """samtools index {output.bam} >>{log} 2>&1; """
                """echo "Compressing unaligned reads..." >>{log} 2>&1; """
                """gzip {unaligned_file_name} >>{log} 2>&1; """
                """echo "Cleaning temp \\"${{ALIGN_BATCH_TEMP}}\\"..." >>{log}; """
                """rm -rf ${{ALIGN_BATCH_TEMP}}"""
                """echo "Done aligning batch {wildcards.batch_id}" >>{log}; """
            )

        else:
            # Empty FOFN, generate empty BAM
            shell(
                """echo "Empty FOFN for batch {wildcards.batch_id}." >{log}; """
                """echo "Creating empty BAM with headers..." >>{log}; """
                """echo -e "@HD\tVN:1.3\tSO:coordinate" | """
                """samtools view -b > {output.bam}; """
                """samtools index {output.bam}; """
                """touch {unaligned_file_name}; """
                """gzip {unaligned_file_name}; """
                """echo "Done" >>{log}; """
            )

# aln_assign_batches
#
# Assign input reads to batches for alignment.
rule aln_assign_batches:
    input:
        reads=get_config_param('reads')
    output:
        fofn=expand('align/batches/{batch_id}.fofn', batch_id=ALIGN_BATCH_LIST)
    run:

        # Check input
        if not input.reads:
            raise RuntimeError('Missing input reads (.fofn) file in parameter "reads"')

        if not input.reads.endswith('.fofn'):
            raise RuntimeError('Input file must end with ".fofn" and contain a paths to sequence files (one per line)')

        if not os.path.isfile(input.reads):
            raise RuntimeError('Input file does not exist or is not a regular file: {}'.format(input.reads))

        # Read input file
        with open(input.reads, 'r') as fh:
            input_files = [os.path.realpath(line.rstrip()) for line in fh]

        input_files = [file_name for file_name in input_files if file_name]  # Remove blank lines

        # Check input file list
        for input_file in input_files:

            # File exists
            if not os.path.isfile(input_file):
                raise RuntimeError('Input file does not exist or is not a regular file: {}'.format(input_file))

            # Must be a subreads BAM file
            if not input_file.lower().endswith('.subreads.bam'):

                if input_file.lower().endswith('.bax.h5'):
                    raise RuntimeError('File in input FOFN is a BAX file: Please use dep/bin/bax2bam to convert and include only ".subreads.bam" files in the input FOFN: {}'.format(input_file))

                if input_file.lower().endswith('.scraps.bam'):
                    raise RuntimeError('File in input FOFN is a scraps BAM file: Please include only ".subreads.bam" BAM files in the input FOFN: {}'.format(input_file))

                raise RuntimeError('Unrecognized file type in input FOFN: Expect a list of only ".subreads.bam" files: {}'.format(input_file))

        # Get the number of files in each batch
        files_per_batch = int(len(input_files) / len(output.fofn))

        files_in_batch = [files_per_batch for i in output.fofn]

        for i in range(len(input_files) % len(output.fofn)):
            files_in_batch[i] += 1

        batch_indices = list(np.cumsum([0] + files_in_batch))

        # Separate into batches
        for batch_id in range(len(output.fofn)):
            with open(output.fofn[batch_id], 'w') as out_file:

                batch_input_file_list = input_files[batch_indices[batch_id]:batch_indices[batch_id + 1]]

                if batch_input_file_list:
                    out_file.write('\n'.join(batch_input_file_list))
                    out_file.write('\n')
