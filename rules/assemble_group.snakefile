"""
Run assemblies for all windows in a group.
"""

import os
import pandas as pd
import shutil
import sys

from Bio import SeqIO

if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'include.snakefile'


###################
### Definitions ###
###################

### Get parameters ###

# Output
CONTIG_BAM = config['contig_out']
CONTIG_BAM_INDEX = CONTIG_BAM + '.bai'

LOG_DIR = config['log_dir']

# Params
MAPQ = config['mapping_quality']
ALN_PARAMS=config['asm_alignment_parameters'].strip('"')
GROUP_ID = config['group_id']
THREADS = config['threads']
POLISH_METHOD = config['asm_polish']
REF_FA = config['ref_fa']

# Input files
ALIGN_FOFN = config['align_fofn']
BED_GROUPS = config['bed_groups']
BED_CANDIDATES = config['bed_candidates']


### Candidate Regions ###

# Get region
DF_CANDIDATES = pd.read_table(BED_CANDIDATES, header=0, index_col='ID')

if GROUP_ID not in set(DF_CANDIDATES['GROUP_ID']):
    raise RuntimeError('Group ID {} is not in the candidates file {}'.format(GROUP_ID, BED_CANDIDATES))

DF_CANDIDATES = DF_CANDIDATES.loc[DF_CANDIDATES['GROUP_ID'] == GROUP_ID, :].copy()
DF_CANDIDATES['POS1'] = DF_CANDIDATES['POS'] + 1

del(DF_CANDIDATES['GROUP_ID'])


### Group Series ###

GROUP = pd.read_table(BED_GROUPS, header=0, index_col='GROUP_ID')
GROUP = GROUP.loc[GROUP_ID].copy()
GROUP['POS1'] = GROUP['POS'] + 1


### Init Paths ###
os.makedirs(LOG_DIR, exist_ok=True)


#############
### Rules ###
#############

# merge_group_contigs
#
# Get final contig for each assembly in this group.
rule merge_group_contigs:
    input:
        bam=expand('region/{region_id}/asm/contig_fixup.bam', region_id=DF_CANDIDATES.index)
    output:
        bam=CONTIG_BAM,
        bai=CONTIG_BAM_INDEX
    resources:
        threads=1
    shell:
        """samtools merge {output.bam} {input.bam}; """
        """samtools index {output.bam}"""

# assemble_align_fixup
#
# Translate alignment from region back to reference and prepend each contig name with
# the region to ensure they are globally unique within each sample.
rule assemble_align_fixup:
    input:
        sam='region/{region_id}/asm/contig.sam',
        align_fofn=ALIGN_FOFN,
        ref_fa=REF_FA,
        group_bam='group/reads.bam'
    output:
        bam=temp('region/{region_id}/asm/contig_fixup.bam')
    resources:
        threads=1
    run:

        bam_temp='temp/region/{}/asm/assemble_align_fixup'.format(wildcards.region_id)
        bam_usort='temp/region/{}/asm/assemble_align_fixup_usort.bam'.format(wildcards.region_id)
        bam_tlen='temp/region/{}/asm/assemble_align_fixup_tlen.bam'.format(wildcards.region_id)

        if os.stat(input.sam).st_size > 0:
            shell(
                """mkdir -p $(dirname {bam_temp}); """
                """alignfixup """
                    """-i {input.sam} """
                    """-d $(head -n 1 {input.align_fofn}) """
                    """-r {REF_FA} -g {wildcards.region_id} """
                    """-o {bam_usort}; """
                """tlenadd -i {bam_usort} -o {bam_tlen} -r {input.ref_fa}; """
                """rm {bam_usort}; """
                """{{ """
                    """samtools view -H {bam_tlen} | grep -Ev '^@(RG|PG)'; """
                    """samtools view {bam_tlen}; """
                """}} | """
                """bamleftalign -f {input.ref_fa} | """
                """samtools sort -T {bam_temp} -o {output.bam}; """
                """rm {bam_tlen}"""
            )
        else:
            # Output bam with headers only
            shell(
                """samtools view -H {input.group_bam} | """
                """grep -E '^@(HD|SQ)' | """
                """samtools view -o {output.bam};"""
            )

# assemble_align_ref_region
#
# Align contigs to the reference region.
rule assemble_align_ref_region:
    input:
        ref='region/{region_id}/asm/ref_region.fasta',
        contig='region/{region_id}/asm/contigs_named.fasta'
    output:
        sam=temp('region/{region_id}/asm/contig.sam')
    params:
        threads=THREADS
    resources:
        threads=THREADS
    run:

        if os.stat(input.contig).st_size > 0:
            shell(
                """blasr """
                    """{input.contig} {input.ref} """
                    """--sam """
                    """--unaligned /dev/null """
                    """--out {output.sam} """
                    """--clipping subread """
                    """--nproc {params.threads} """
                    """{ALN_PARAMS};"""
            )

        else:
            open(output.sam, 'w').close()  # Touch and/or clear file

# assemble_get_ref_region
#
# Get reference region.
rule assemble_get_ref_region:
    input:
        fasta=REF_FA,
        fasta_contig='region/{region_id}/asm/contigs.fasta'
    output:
        fasta=temp('region/{region_id}/asm/ref_region.fasta')
    resources:
        threads=1
    run:

        if os.stat(input.fasta_contig).st_size > 0:
            # Get region
            if not wildcards.region_id in DF_CANDIDATES.index:
                raise RuntimeError('Region ID {} is not in the candidates file {}'.format(wildcards.region_id, BED_CANDIDATES))

            region = '{#CHROM}:{POS1}-{END}'.format(**DF_CANDIDATES.loc[wildcards.region_id])

            # Extract region
            if os.stat(input.fasta).st_size > 0:
                shell("""samtools faidx {input.fasta} {region} >{output.fasta}""")
            else:
                open(output.fasta, 'w').close()  # Touch and/or clear file

        else:
            open(output.fasta, 'w').close()

# assemble_set_pb_seq_name
#
# BLASR SAM output is designed for PacBio SAM/BAM input. To align a FASTA file, the sequence names must follow
# the PacBio convention "movie/zmw/start_end". For this pipeline, the movie name will be the name of the
# polished contig, the ZMW will start at 1 and be incremented, the start will be 0, and the end will be the
# the length of the sequence.
rule assemble_set_pb_seq_name:
    input:
        fasta='region/{region_id}/asm/contigs_polished.fasta'
    output:
        fasta=temp('region/{region_id}/asm/contigs_named.fasta')
    resources:
        threads=1
    run:

        if os.stat(input.fasta).st_size > 0:
            zmw_id = 0

            seq_list = list()

            # Read and set name
            with open(input.fasta, 'r') as in_file:
                with open(output.fasta, 'w') as out_file:
                    for record in SeqIO.parse(in_file, 'fasta'):
                        zmw_id += 1

                        # movie/zmw/start_end
                        record.id = '{}/{}/0_{}'.format(record.name, zmw_id, len(record.seq))
                        record.description = ''

                        seq_list.append(record)

            # Write
            with open(output.fasta, 'w') as out_file:
                SeqIO.write(seq_list, out_file, 'fasta')

        else:
            open(output.fasta, 'w').close()

# assemble_polish
#
# Polish assembly.
rule assemble_polish:
    input:
        fasta='region/{region_id}/asm/contigs.fasta',
        fai='region/{region_id}/asm/contigs.fasta.fai',
        bam='region/{region_id}/asm/contig_aligned_reads.bam',
        pbi='region/{region_id}/asm/contig_aligned_reads.bam.pbi'
    output:
        fasta=temp('region/{region_id}/asm/contigs_polished.fasta')
    params:
        algorithm=POLISH_METHOD,
        threads=THREADS
    resources:
        threads=THREADS
    run:

        if os.stat(input.fasta).st_size > 0:
            shell(
                """variantCaller """
                    """--referenceFilename {input.fasta} """
                    """{input.bam} """
                    """-o {output.fasta} """
                    """-j {params.threads} """
                    """--algorithm={params.algorithm}; """
            )

        else:  # Touch and/or clear file
            open(output.fasta, 'w').close()

# assemble_align_org
#
# Align original reads back to assembly
rule assemble_align_org:
    input:
        fasta='region/{region_id}/asm/contigs.fasta',
        bam='region/{region_id}/reads/reads.bam'
    output:
        bam=temp('region/{region_id}/asm/contig_aligned_reads.bam'),
        pbi=temp('region/{region_id}/asm/contig_aligned_reads.bam.pbi')
    params:
        threads=THREADS
    resources:
        threads=THREADS
    run:

        # Make temporary file location (usorted BAM)
        bam_usort='temp/region/{}/reads/align_reads.bam'.format(wildcards.region_id)
        bam_stemp='temp/region/{}/reads/align_reads_slice'.format(wildcards.region_id)

        # Map reads
        if os.stat(input.fasta).st_size > 0:
            os.makedirs(os.path.dirname(bam_usort), exist_ok=True)

            shell(
                """blasr """
                """{input.bam} {input.fasta} """
                """--bam --bestn 1 """
                """--unaligned /dev/null """
                """--out {bam_usort} """
                """--nproc {params.threads}; """
                """samtools sort -O bam -T {bam_stemp} -o {output.bam} {bam_usort}; """
                """pbindex {output.bam}; """
                """rm {bam_usort}"""
            )

        else:
            open(output.bam, 'w').close()  # Touch and/or clear file
            open(output.pbi, 'w').close()  # Touch and/or clear file

# assemble_reads
#
# Assemble sequence reads.
rule assemble_reads:
    input:
        fasta='region/{region_id}/reads/reads.fasta'
    output:
        fasta=temp('region/{region_id}/asm/contigs.fasta'),
        fai=temp('region/{region_id}/asm/contigs.fasta.fai'),
        preads=temp('region/{region_id}/asm/corrected_reads.fastq.gz')  # Corrected reads
    params:
        threads=THREADS,
        read_length='1000',
        partitions='50',
        max_runtime='30m'
    resources:
        threads=THREADS
    log:
        os.path.join(LOG_DIR, '{region_id}.log')
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
        shell(
            """set +e; """
            """timeout {params.max_runtime} canu """
                """-pacbio-raw """
                """{input.fasta} """
                """genomeSize={genome_size} """
                """-d {canu_dir} """
                """-p asm """
                """useGrid=false """
                """corMhapSensitivity=high """
                """corMinCoverage=2 """
                """correctedErrorRate=0.045 """  # Was errorRate=0.035
                """-pacbio-raw """
                """>{log} 2>&1; """
            """RET_CODE=$?; """
            """if [ ${{RET_CODE}} -ne 0 ]; then """
                """if [ ${{RET_CODE}} -eq 124 ]; then """
                    """echo "Assembly timeout ({params.max_runtime})"; """
                """else """
                    """echo "Assembly error: Return code = ${{RET_CODE}}"; """
                """fi; """
            """fi; """

        )

        # Copy unitigs and corrected reads
        contig_file_name = 'region/{region_id}/asm/canu/asm.unitigs.fasta'.format(**wildcards)
        preads_file_name = 'region/{region_id}/asm/canu/asm.correctedReads.fasta.gz'.format(**wildcards)

        if os.path.exists(contig_file_name) and os.stat(contig_file_name).st_size > 0:
            shutil.copyfile(contig_file_name, output.fasta)
            shutil.copyfile(preads_file_name, output.preads)

            shell("""samtools faidx {output.fasta}""")

        else:  # Touch and/or clear files
            open(output.fasta, 'w').close()
            open(output.preads, 'w').close()
            open(output.fai, 'w').close()

        # Clean assembly directory
        shutil.rmtree(canu_dir, ignore_errors=True)

# convert_reads_to_fasta
#
# Get FASTA and FASTQ file of extracted sequence reads.
rule asm_group_reads_to_fasta:
    input:
        bam='region/{region_id}/reads/reads.bam',
        bai='region/{region_id}/reads/reads.bam.bai'
    output:
        fasta=temp('region/{region_id}/reads/reads.fasta'),
        fastq=temp('region/{region_id}/reads/reads.fastq')
    resources:
        threads=1
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
        bam='group/reads.bam',
        bai='group/reads.bam.bai'
    output:
        bam=temp('region/{region_id}/reads/reads.bam'),
        bai=temp('region/{region_id}/reads/reads.bam.bai')
    resources:
        threads=1
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
        bam='group/reads.bam',
        bai='group/reads.bam.bai'
    params:
        mapq=MAPQ,
        group_id=GROUP_ID
    resources:
        threads=1
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
        batch_temp = 'group/batch_temp'
        os.makedirs(batch_temp, exist_ok=True)

        try:
            for bam_file in bam_file_list:
                batch_index = os.path.basename(bam_file).rstrip('.bam')

                shell(
                    """echo "Extracting reads from batch {batch_index}..."; """
                    """samtools view -hb -q {params.mapq} {bam_file} {group_region} """
                    """>{batch_temp}/{batch_index}.bam"""
                )

            shell(
                """echo "Merging batches..."; """
                """samtools merge {output.bam} {batch_temp}/*.bam; """
                """samtools index {output.bam}; """
            )

        finally:
            shutil.rmtree(batch_temp, ignore_errors=True)
