"""
Coordinates assembly of each group of regions and merges results from each group.
"""

import pandas as pd
import shutil
import subprocess

from smrtsvlib import smrtsvrunner

if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'include.snakefile'


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

        # TODO: See old rule collect_assembly_alignments for commands to left-align, etc

# asm_assemble_group
#
# Assemble one group of regions.
rule asm_assemble_group:
    input:
        align_fofn='align/alignments.fofn',
        bed_grp='detect/candidate_groups.bed',
        bed_can='detect/candidates.bed'
    output:
        bam='assemble/group/{group_id}/contig.bam',
        bai='assemble/group/{group_id}/contig.bam.bai'
    params:
        mapq=get_config_param('mapping_quality'),
        align_params=get_config_param('asm_alignment_parameters')
    log:
        'assemble/group/{group_id}/contig.log'
    run:

        # Set assemble_temp (will be deleted if not None)
        assemble_temp = None

        try:

            # Create temporary directories
            assemble_temp = tempfile.mkdtemp(
                prefix=os.path.join(
                    TEMP_DIR,
                    'asm_assemble_group_{}_'.format(wildcards.group_id)
                )
            )

            # Create log directory
            log_dir = os.path.abspath('assemble/group/{group_id}/log'.format(wildcards.group_id))
            os.makedirs(log_dir, exist_ok=True)

            # Setup sub-Snake command
            command = (
                'get_contig',
                '-f',
                '--config',
                'contig_bam={}'.format(os.path.abspath(output.bam)),
                'log_dir={}'.format(log_dir),
                'mapping_quality={}'.format(params.mapq),
                'asm_alignment_parameters={}'.format(params.align_params),  # String is already quoted to prevent Snakemake from trying to interpret the command options
                'group_id={}'.format(wildcards.group_id),
                'align_fofn={}'.format(os.path.abspath(input.align_fofn)),
                'bed_groups={}'.format(os.path.abspath(input.bed_grp)),
                'bed_candidates={}'.format(os.path.abspath(input.bed_can))
            )

            # Clear locks
            shell("""rm -f {working_dir}/.snakemake/locks/*""")

            # Run assembly
            with open(log, 'w') as log_file:
                smrtsvrunner.run_snake_target(
                    'rules/assemble_group.snakefile', None, PROCESS_ENV, SMRTSV_DIR, command,
                    stdout=log_file, stderr=subprocess.STDOUT, cwd=assemble_temp
                )

        finally:

            # Clean up temp directory
            if assemble_temp is not None:
                shutil.rmtree(assemble_temp)




#######################################################################################################
#######################################################################################################
#######                                          OLD CODE                                        ######
#######################################################################################################
#######################################################################################################
##
## Load regions to assemble.
#REGIONS_TO_ASSEMBLE = config.get("regions_to_assemble", "filtered_assembly_candidates_with_coverage.bed")
#MAPPING_QUALITY = get_config_param('mapping_quality')
#ASM_ALIGNMENT_PARAMETERS = get_config_param('asm_alignment_parameters')
##
#
####################
#### Definitions ###
####################
##
#def _get_assembly_alignments(wildcards):
#    if os.path.exists(REGIONS_TO_ASSEMBLE) and not os.path.exists(LOCAL_ASSEMBLY_ALIGNMENTS):
#        with open(REGIONS_TO_ASSEMBLE, "r") as fh:
#            LIST_OF_REGIONS_TO_ASSEMBLE = ["-".join(line.rstrip().split("\t")[:3]) for line in fh if line.rstrip()]
##
#        return ["mhap_assembly/{chromosome}/{region}/consensus_reference_alignment.sam".format(chromosome=region.split("-")[0], region=region) for region in LIST_OF_REGIONS_TO_ASSEMBLE]
#    else:
#        return []
##
##############
#### Rules ###
##############
##
##
## Collect assemblies.
##
##
## collect_assembly_alignments
##
## Collect all assembly alignments.
#rule collect_assembly_alignments:
#    input:
#        alignments=_get_assembly_alignments,
#        chromosome_lengths=CHROMOSOME_LENGTHS,
#        reference=config["reference"],
#        regions=REGIONS_TO_ASSEMBLE
#    output:
#        LOCAL_ASSEMBLY_ALIGNMENTS
#    run:
#        list_filename = output[0].replace("bam", "list.txt")
#        base_filename = os.path.basename(list_filename)
#        with open(list_filename, "w") as oh:
#            for i in input.alignments:
#                oh.write("%s\n" % i)
##
#        shell("""mkdir -p {TMP_DIR}; while read file; do sed 's/\/0_[0-9]\+//' $file; done < %s | samtools view -Sbu -t {input.chromosome_lengths} - | bamleftalign -f {input.reference} | samtools sort -O bam -T {TMP_DIR}/%s -o {output}""" % (list_filename, base_filename))
##
## assemble_region
##
## Assemble one region. Extracts reads from the region, assembles, polishes, and re-aligns the contig back to the
## reference sequence.
#rule assemble_region:
#    input:
#        alignments=ALIGNMENTS
#    output:
#        sam=touch(ASSEMBLY_DIR + '/{chromosome}/{region}/consensus_reference_alignment.sam')
#    log:
#        ASSEMBLY_DIR + '/{chromosome}/{region}/mhap.log'
#    params:
#        threads='4'
#    shell:
#        # Setup Log
#        """JOB_LOG=$(readlink -f {log}); """
##
#        # Random sleep up to 60 seconds; prevents overwhelming IO starting many jobs
#        """sleep $(shuf -i 0-60 -n 1); """
##
#        # Setup paths
#        """export ANALYSIS_DIR=`pwd`; """
#        """export ALIGNMENTS_PATH=`readlink -f {input.alignments}`; """
#        """export REFERENCE_PATH=`readlink -f {REFERENCE}`; """
#        """export INPUT_READS_PATH=`readlink -f {INPUT_READS}`; """
#        """export LOG_FILE=`readlink -f {LOG_FILE}`; """
##
#        # Get locations of output files
#        """ALIGN_FILE={TMP_DIR}/{wildcards.region}/consensus_reference_alignment.sam; """
#        """ALIGN_LOG_FILE={TMP_DIR}/{wildcards.region}/assembly.log; """
##
#        # Make temp directory
#        """mkdir -p {TMP_DIR}/{wildcards.region}; """
#        """echo "Temp directory: {TMP_DIR}/{wildcards.region}" > $JOB_LOG; """
##
#        # Enter temp directory
#        """echo -n "Entering temp directory (popd): " >> $JOB_LOG; """
#        """pushd {TMP_DIR}/{wildcards.region} >> $JOB_LOG 2>&1; """
##
#        # Copy config.json
#        """if [[ -e "$ANALYSIS_DIR/config.json" ]]; then rsync $ANALYSIS_DIR/config.json {TMP_DIR}/{wildcards.region}/; fi; """
##
#        # Run assembly
#        """echo "Running snakemake" >> $JOB_LOG; """
#        """snakemake -j {params.threads} -s {SNAKEMAKE_DIR}/rules/assemble_group.snakefile """
#            """--config alignments=$ALIGNMENTS_PATH """
#            """reference=$REFERENCE_PATH """
#            """reads=$INPUT_READS_PATH """
#            """region={wildcards.region} """
#            """log=$LOG_FILE """
#            """mapping_quality={MAPPING_QUALITY} """
#            """asm_alignment_parameters=\\"{ASM_ALIGNMENT_PARAMETERS}\\" """  # Quote BLASR parameters to prevent Snakemake from interpreting them
#            """>> $JOB_LOG 2>&1; """
#        """echo "Snakemake complete" >> $JOB_LOG; """
#        """popd; """
##
#        # Copy alignment and alignment log
#        """echo -n "$ALIGN_FILE exists: " >>$JOB_LOG; """
#            """if [[ -e $ALIGN_FILE ]]; then echo "True" >>$JOB_LOG; else echo "False" >>$JOB_LOG; fi; """
#        """rsync -a $ALIGN_FILE $ALIGN_LOG_FILE $(dirname {output.sam})/ 2>>$JOB_LOG; """
##
#        # Cleanup
#        """rm -rf {TMP_DIR}/{wildcards.region}; """
