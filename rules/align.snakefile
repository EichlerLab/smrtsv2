"""
Align sequence reads in batches and produce an FOFN file listing each BAM file.
"""

import os

import numpy as np

include: 'include.snakefile'

localrules: aln_run, aln_assign_batches


####################
### Declarations ###
####################

# Run parameters
READS = config.get('reads')
BATCH_COUNT = int(config.get('batches', '1'))
THREADS = int(config.get('threads', '1'))
ALIGN_PARAMS = config.get("alignment_parameters", "").strip('"')

# List of batches
BATCHES = list(range(BATCH_COUNT))


#############
### Rules ###
#############

# aln_run
#
# Run all alignments in batches and write and FOFN file of each BAM (one per batch).
rule aln_run:
    input:
        bam=expand('align/bam/{batch_id}.bam', batch_id=BATCHES),
        bai=expand('align/bam/{batch_id}.bam.bai', batch_id=BATCHES),
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
        reads='align/batches/fofn/{batch_id}.fofn',
        ref_fa='reference/ref.fasta',
        ref_fai='reference/ref.fasta.fai',
        ref_sa='reference/ref.fasta.sa',
        ref_ctab='reference/ref.fasta.ctab'
    output:
        bam=protected('align/bam/{batch_id}.bam'),
        bai=protected('align/bam/{batch_id}.bam.bai')
    params:
        threads=str(THREADS),
        samtools_threads="1",
        samtools_memory="4G",
        bwlimit="30000",
        align_params=ALIGN_PARAMS
    log:
        'align/batches/fofn/log/{batch_id}.log'
    run:

        if os.path.getsize(input.reads) > 0:

            shell(
                """ALIGN_BATCH_TEMP={TEMP_DIR}/aln_align_batch_{wildcards.batch_id}; """
                """mkdir -p ${{ALIGN_BATCH_TEMP}}; """
                """echo "Aligning batch {wildcards.batch_id}..." >{log}; """
                """blasr {input.reads} {input.ref_fa} """
                    """--unaligned /dev/null """
                    """--out ${{ALIGN_BATCH_TEMP}}/batch_out.bam """
                    """--sam """
                    """--sa {input.ref_sa} """
                    """--ctab {input.ref_ctab} """
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
                """samtools index {output.bam}; """
                """echo "Cleaning temp \\"${{ALIGN_BATCH_TEMP}}\\"..." >>{log}; """
                """rm -rf ${{ALIGN_BATCH_TEMP}}"""
                """echo "Indexing..." >>{log} 2>&1; """
                """echo "Done aligning batch {wildcards.batch_id}" >>{log}; """
            )

        else:
            # Empty FOFN, generate empty BAM
            shell(
                """echo "Empty FOFN for batch {wildcards.batch_id}." >{log}; """
                """echo "Creating empty BAM with headers..." >>{log}; """
                """echo -e "@HD\tVN:1.3\tSO:coordinate" | """
                """samtools view -b > {output.bam}; """
                """samtools index {output.bam}"""
                """echo "Done" >>{log}; """
            )

# aln_assign_batches
#
# Assign input reads to batches for alignment.
rule aln_assign_batches:
    input:
        reads=config.get('reads', '')
    output:
        fofn=expand('align/batches/fofn/{batch_id}.fofn', batch_id=BATCHES)
    run:

        # Check input
        if not input.reads:
            raise RuntimeError('Missing input reads (.fofn) file')

        if not input.reads.endswith('.fofn'):
            raise RuntimeError('Input file must end with ".fofn" and contain a paths to sequence files (one per line)')

        if not os.path.isfile(input.reads):
            raise RuntimeError('Input file does not exist or is not a regular file: %s' % input.reads)

        # Read input file
        with open(input.reads, 'r') as fh:
            input_files = [os.path.realpath(line.rstrip()) for line in fh]

        input_files = [file_name for file_name in input_files if file_name]  # Remove blank lines

        # Check input file list
        for input_file in input_files:
            if not os.path.isfile(input_file):
                raise RuntimeError('Input file does not exist or is not a regular file: {}'.format(input_file))

        # Get the number of files in each batch
        files_per_batch = int(len(input_files) / len(output.fofn))

        files_in_batch = [files_per_batch for i in output.fofn]

        for i in range(len(input_files) % len(output.fofn)):
            files_in_batch[i] += 1

        batch_indices = list(np.cumsum([0] + files_in_batch))

        # Separate into batches
        for batch_id in range(len(output.fofn)):
            with open(output.fofn[batch_id], 'w') as out_file:
                out_file.write('\n'.join(input_files[batch_indices[batch_id]:batch_indices[batch_id + 1]]))
                out_file.write('\n')
