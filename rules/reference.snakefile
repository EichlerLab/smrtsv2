"""
Rules for preparing the reference.
"""

import os

include: 'include.snakefile'

localrules: ref_run, ref_set_fa


#
# Parameters
#

REFERENCE = config.get('reference', None)
LINK_INDEX = config.get('link_index', True)


#
# Rules
#

# ref_run
#
# Create all reference files that do not exist.
rule ref_run:
    input:
        ref_fa='reference/ref.fasta',
        ref_fai='reference/ref.fasta.fai',
        ref_ctab='reference/ref.fasta.ctab',
        ref_sa='reference/ref.fasta.sa'

# ref_make_ctab
#
# Count k-mers in the reference.
rule ref_make_ctab:
    input:
        ref_fa='reference/ref.fasta'
    output:
        ref_ctab='reference/ref.fasta.ctab'
    run:

        reference_ctab = REFERENCE + '.ctab'

        if os.path.isfile(reference_ctab) and LINK_INDEX:
            shell("""ln -sf {} {}""".format(reference_ctab, output.ref_ctab))
        else:
            shell("""printTupleCountTable {input.ref_fa} > {output.ref_ctab}""")

# ref_make_sa
#
# Bulid index suffix array.
rule ref_make_sa:
    input:
        ref_fa='reference/ref.fasta'
    output:
        ref_sa='reference/ref.fasta.sa'
    log:
        'reference/log/ref_make_sa.log'
    run:

        reference_sa = REFERENCE + '.sa'

        if os.path.isfile(reference_sa) and LINK_INDEX:
            shell("""ln -sf {} {}""".format(reference_sa, output.ref_sa))
        else:
            shell("""sawriter {output} {input} >{log} 2>&1""")

# ref_make_fai
#
# Make FASTA index.
rule ref_make_fai:
    input:
        ref_fa='reference/ref.fasta'
    output:
        ref_fai='reference/ref.fasta.fai'
    log:
        'reference/log/ref_make_fai.log'
    run:

        REFERENCE_FAI = REFERENCE + '.fai'

        if os.path.isfile(REFERENCE_FAI) and LINK_INDEX:
            shell("""ln -sf {} {}""".format(REFERENCE_FAI, output.ref_fai))
        else:
            shell("""samtools faidx {input.fa} >{log} 2>&1""")

# ref_set_fa
#
# Link the reference FASTA to the reference directory.
rule ref_set_fa:
    output:
        ref_fa='reference/ref.fasta'
    run:

        reference = REFERENCE

        # Check reference
        if reference is None:
            raise RuntimeError('Cannot find reference sequence: Required configuration parameter "reference" is not set')

        reference = reference.strip()

        if not reference:
            raise RuntimeError('Cannot find reference sequence: Required configuration parameter "reference" is empty')

        # Get absolute path
        reference = os.path.abspath(reference)

        shell('ln -sf {} {}'.format(reference, output.ref_fa))
