"""
Rules for preparing the reference.
"""

import os

include: 'include.snakefile'

localrules: ref_run, ref_set_fa


####################
### Declarations ###
####################

REFERENCE = config.get('reference', None)
LINK_INDEX = str(config.get('link_index', True))


#############
### Rules ###
#############

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
    shell:
        """if [ "{LINK_INDEX}" = "True" -a -f {REFERENCE}.ctab ]; then """
            """ln -sf {REFERENCE}.ctab {output.ref_ctab}; """
        """else """
            """printTupleCountTable {input.ref_fa} > {output.ref_ctab}; """
            """chmod a-w {output.ref_ctab}; """
        """fi"""

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
    shell:
        """if [ "{LINK_INDEX}" = "True" -a -f {REFERENCE}.sa ]; then """
            """ln -sf {REFERENCE}.sa {output.ref_sa}; """
        """else """
            """sawriter {output.ref_sa} {input.ref_fa} >{log} 2>&1; """
            """chmod a-w {output.ref_sa}; """
        """fi"""

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
    shell:
        """if [ "{LINK_INDEX}" = "True" -a -f {REFERENCE}.fai ]; then """
            """ln -sf {REFERENCE}.fai {output.ref_fai}; """
        """else """
            """samtools faidx {input.ref_fa} >{log} 2>&1; """
            """chmod a-w {output.ref_fai}; """
        """fi"""

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
