"""
Rules for preparing the reference.
"""

import os

localrules: ref_run, ref_set_fa


# ref_run
#
# Create all reference files that do not exist.
rule ref_run:
    input:
        fa='reference/ref.fasta',
        fai='reference/ref.fasta.fai',
        sa='reference/ref.fasta.sa'

# ref_make_sa
#
# Bulid index suffix array.
rule ref_make_sa:
    input:
        fa='reference/ref.fasta'
    output:
        sa='reference/ref.fasta.sa'
    log:
        'reference/log/ref_make_sa.log'
    shell:
        """sawriter {output} {input} >{log} 2>&1"""

# ref_make_fai
#
# Make FASTA index.
rule ref_make_fai:
    input:
        fa='reference/ref.fasta'
    output:
        fai='reference/ref.fasta.fai'
    log:
        'reference/log/ref_make_fai.log'
    shell:
        """samtools faidx {input.fa} >{log} 2>&1"""

# ref_set_fa
#
# Link the reference FASTA to the reference directory.
rule ref_set_fa:
    output:
        fa='reference/ref.fasta'
    run:

        reference = config.get('reference', None)

        # Check reference
        if reference is None:
            raise RuntimeError('Cannot find reference sequence: Required configuration parameter "reference" is not set')

        reference = reference.strip()

        if not reference:
            raise RuntimeError('Cannot find reference sequence: Required configuration parameter "reference" is empty')

        # Get absolute path
        reference = os.path.abspath(reference)

        shell('ln -sf {} {}'.format(reference, output.fa))
