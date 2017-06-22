"""
Run assemblies for all windows in a group.
"""

import os
import pandas as pd
import shutil

if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'include.snakefile'


###################
### Definitions ###
###################

### Get parameters ###

# Params
MAPQ = config['mapping_quality']
ALN_PARAMS=config['asm_alignment_parameters'].strip('"')
GROUP_ID = config['group_id']
ASSEMBLE_TEMP = config['assemble_temp']

# Input files
ALIGN_FOFN = config['align_fofn']
BED_GROUPS = config['bed_groups']
BED_CANDIDATES = config['bed_candidates']


### Utility Functions ###

def _asm_temp(file_name):
    """
    Get a path relative to the temporary directory `ASSEMBLE_TEMP`.

    :param temp_file: Relative path to 'ASSEMBLE_TEMP`.

    :return: Full path in `ASSEMBLE_TEMP`.
    """
    return os.path.join(ASSEMBLE_TEMP, file_name)


### Candidate Regions ###

# Get region
DF_CANDIDATES = pd.read_table(BED_CANDIDATES, header=0, index_col='ID')

if GROUP_ID not in set(DF_CANDIDATES['GROUP_ID']):
    raise RuntimeError('Group ID {} is not in the candidates file {}'.format(GROUP_ID, BED_CANDIDATES))

DF_CANDIDATES = DF_CANDIDATES.loc[DF_CANDIDATES['GROUP_ID'] == GROUP_ID, :]
DF_CANDIDATES['POS1'] = DF_CANDIDATES['POS'] + 1

del(DF_CANDIDATES['GROUP_ID'])


### Group Series ###

GROUP = pd.read_table(BED_GROUPS, header=0, index_col='GROUP_ID')
GROUP = GROUP.loc[GROUP_ID].copy()
GROUP['POS1'] = GROUP['POS'] + 1


#############
### Rules ###
#############

# assemble_reads
#
# Assemble sequence reads.
rule assemble_reads:
    input:
        fasta=_asm_temp('region/{region_id}/reads/reads.fasta'),
    output:
        fasta=_asm_temp('region/{region_id}/asm/contigs.fasta')
    params:
        threads='4',
        read_length='1000',
        partitions='50',
        max_runtime='10m'
    log:
        'canu.log'
    run:

        # Setup output locations and flag
        assembly_exists = False

        # Calculate genome size
        df_region = DF_CANDIDATES.loc[wildcards.region_id]
        genome_size = df_region['END'] - df_region['POS']

        # Get output directory
        canu_dir = os.path.join(os.path.dirname(output.fasta), 'canu')
        shutil.rmtree(canu_dir, ignore_errors=True)

        # Run assembly
        try:
            shell(
                """timeout {params.max_runtime} canu """
                """-pacbio-raw """
                """{input.fasta} """
                """genomeSize={genome_size} """
                """-d {canu_dir} """
                """-p asm """
                """useGrid=false """
                """corMhapSensitivity=high """
                """corMinCoverage=2 """
                """errorRate=0.035 """
                """>{log}"""
            )
        except:
            print('Assembly crashed on region: {wildcards.region_id}')


        # Find assembly and copy to output.
        if os.path.exists(assembly_output) and os.stat(assembly_output).st_size > 0:
            shell(
                """cat {assembly_output} > {output.fasta}; """
                """echo -e "{REGION}\tassembly_exists" >> %s""" % config['log']
            )
            assembly_exists = True

        elif os.path.exists(unitig_output) and os.stat(unitig_output).st_size > 0:
            shell(
                """cat {unitig_output} > {output.fasta}; """ +
                """echo -e "{REGION}\tunitig_assembly_exists" >> %s""" % config['log']
            )
            assembly_exists = True

        else:
            shell("""echo -e "{REGION}\tno_assembly_exists" >> %s""" % config['log'])

        # Create an empty assembly for failed regions.
        if not assembly_exists:
            shell("echo -e '>{REGION}\nN' > {output}")


# convert_reads_to_fasta
#
# Get FASTA and FASTQ file of extracted sequence reads.
rule asm_group_reads_to_fasta:
    input:
        bam=_asm_temp('region/{region_id}/reads/reads.bam')
    output:
        fasta=_asm_temp('region/{region_id}/reads/reads.fasta'),
        fastq=_asm_temp('region/{region_id}/reads/reads.fastq')
    shell:
        """echo "### Entering: asm_group_get_region_bam"; """
        """echo "Writing FASTA: {output.fasta}"; """
        """echo "Writing FASTQ: {output.fastq}"; """
        """{SMRTSV_DIR}/scripts/assemble/BamToFasta.py {input.bam} {output.fasta} --fakename --fastq {output.fastq}"""

# asm_group_get_region_bam
#
# Get reads for one region from the alignment cache (align/reads.bam).
rule asm_group_get_region_bam:
    input:
        bam=_asm_temp('group/reads.bam'),
        bai=_asm_temp('group/reads.bam.bai')
    output:
        bam=temp(_asm_temp('region/{region_id}/reads/reads.bam')),
        bai=temp(_asm_temp('region/{region_id}/reads/reads.bam.bai'))
    run:

        # Get region
        if not wildcards.region_id in DF_CANDIDATES.index:
            raise RuntimeError('Region ID {} is not in the candidates file {}'.format(wildcards.region_id, BED_CANDIDATES))

        region = '{#CHROM}:{POS1}-{END}'.format(**DF_CANDIDATES.loc[wildcards.region_id])

        # Extract reads
        print('Extracting over region: {}'.format(region))

        shell(
            """samtools view -hb {input.bam} {region} """
            """>{output.bam}; """
            """samtools index {output.bam}"""
        )

# asm_group_get_reads
#
# Cache reads for all regions within this group. When distributed, this makes the pipeline extract reads once from
# shared storage, and each region fetches reads from local storage.
rule asm_group_get_reads:
    input:
        fofn=ALIGN_FOFN
    output:
        bam=_asm_temp('group/reads.bam'),
        bai=_asm_temp('group/reads.bam.bai')
    params:
        mapq=MAPQ,
        group_id=GROUP_ID
    run:

        # Read input alignment batches into a list
        bam_file_list = list()

        with open(input.fofn, 'r') as in_file:
            for line in in_file:

                line = line.strip()

                if not line:
                    continue

                bam_file_list.append(line)

        group_region = '{#CHROM}:{POS1}-{END}'.format(**GROUP)

        print('Extracting over region: {}'.format(group_region))

        # Extract reads
        BATCH_TEMP = _asm_temp('group/batch_temp')
        os.makedirs(BATCH_TEMP, exist_ok=True)

        try:
            for bam_file in bam_file_list:
                batch_index = os.path.basename(bam_file).rstrip('.bam')

                shell(
                    """echo "Extracting reads from batch {batch_index}..."; """
                    """samtools view -hb -q {params.mapq} {bam_file} {group_region} """
                    """>{BATCH_TEMP}/{batch_index}.bam"""
                )

            shell(
                """echo "Merging batches..."; """
                """samtools merge {output.bam} {BATCH_TEMP}/*.bam; """
                """samtools index {output.bam}; """
            )

        finally:
            shutil.rmtree(BATCH_TEMP, ignore_errors=True)
