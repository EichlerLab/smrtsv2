"""
Rules for local assembly of genomic regions.
"""
import csv
import os

#
# Inputs:
#  1. Text file with a list of absolute paths for BAMs with reads for assembly
#  2. BED file with a list of regions to assemble
#
# For signature-based SV calling, the list of regions is based on signatures of
# SVs. For the tiled-based SV calling, the regions are sliding windows across
# the entire genome.

# Load regions to assemble.
INPUT_READS = config.get("reads", "input.fofn")
REGIONS_TO_ASSEMBLE = config.get("regions_to_assemble", "filtered_assembly_candidates_with_coverage.bed")
LOCAL_ASSEMBLY_ALIGNMENTS = config.get("assembly_alignments", "local_assembly_alignments.bam")
MAPPING_QUALITY = config.get("mapping_quality")
ASM_ALIGNMENT_PARAMETERS = config.get("asm_alignment_parameters")

#
# Define internal constants.
#
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
ASSEMBLY_DIR = "mhap_assembly"
LOG_FILE = config.get("assembly_log", "assembly.log")
REFERENCE = config["reference"]

# User-defined file of alignments with one absolute path to a BAM per line.
ALIGNMENTS = config.get("alignments", "")

#
# Define helper functions.
#
def _get_assembly_alignments(wildcards):
    if os.path.exists(REGIONS_TO_ASSEMBLE) and not os.path.exists(LOCAL_ASSEMBLY_ALIGNMENTS):
        with open(REGIONS_TO_ASSEMBLE, "r") as fh:
            LIST_OF_REGIONS_TO_ASSEMBLE = ["-".join(line.rstrip().split("\t")[:3]) for line in fh if line.rstrip()]

        return ["mhap_assembly/{chromosome}/{region}/consensus_reference_alignment.sam".format(chromosome=region.split("-")[0], region=region) for region in LIST_OF_REGIONS_TO_ASSEMBLE]
    else:
        return []

#
# Define rules.
#

#
# Collect assemblies.
#

# collect_assembly_alignments
#
# Collect all assembly alignments.
rule collect_assembly_alignments:
    input:
        alignments=_get_assembly_alignments,
        chromosome_lengths=CHROMOSOME_LENGTHS,
        reference=config["reference"],
        regions=REGIONS_TO_ASSEMBLE
    output:
        LOCAL_ASSEMBLY_ALIGNMENTS
    run:
        list_filename = output[0].replace("bam", "list.txt")
        base_filename = os.path.basename(list_filename)
        with open(list_filename, "w") as oh:
            for i in input.alignments:
                oh.write("%s\n" % i)

        shell("""mkdir -p {TMP_DIR}; while read file; do sed 's/\/0_[0-9]\+//' $file; done < %s | samtools view -Sbu -t {input.chromosome_lengths} - | bamleftalign -f {input.reference} | samtools sort -O bam -T {TMP_DIR}/%s -o {output}""" % (list_filename, base_filename))

# assemble_region
#
# Assemble one region. Extracts reads from the region, assembles, polishes, and re-aligns the contig back to the
# reference sequence.
rule assemble_region:
    input:
        alignments=ALIGNMENTS
    output:
        sam=touch(ASSEMBLY_DIR + '/{chromosome}/{region}/consensus_reference_alignment.sam')
    log:
        ASSEMBLY_DIR + '/{chromosome}/{region}/mhap.log'
    params:
        threads='4'
    shell:
        # Setup Log
        """JOB_LOG=$(readlink -f {log}); """

        # Random sleep up to 60 seconds; prevents overwhelming IO starting many jobs
        """sleep $(shuf -i 0-60 -n 1); """

        # Setup paths
        """export ANALYSIS_DIR=`pwd`; """
        """export ALIGNMENTS_PATH=`readlink -f {input.alignments}`; """
        """export REFERENCE_PATH=`readlink -f {REFERENCE}`; """
        """export INPUT_READS_PATH=`readlink -f {INPUT_READS}`; """
        """export LOG_FILE=`readlink -f {LOG_FILE}`; """

        # Get locations of output files
        """ALIGN_FILE={TMP_DIR}/{wildcards.region}/consensus_reference_alignment.sam; """
        """ALIGN_LOG_FILE={TMP_DIR}/{wildcards.region}/assembly.log; """

        # Make temp directory
        """mkdir -p {TMP_DIR}/{wildcards.region}; """
        """echo "Temp directory: {TMP_DIR}/{wildcards.region}" > $JOB_LOG; """

        # Enter temp directory
        """echo -n "Entering temp directory (popd): " >> $JOB_LOG; """
        """pushd {TMP_DIR}/{wildcards.region} >> $JOB_LOG 2>&1; """

        # Copy config.json
        """if [[ -e "$ANALYSIS_DIR/config.json" ]]; then rsync $ANALYSIS_DIR/config.json {TMP_DIR}/{wildcards.region}/; fi; """

        # Run assembly
        """echo "Running snakemake" >> $JOB_LOG; """
        """snakemake -j {params.threads} -s {SNAKEMAKE_DIR}/rules/assemble_group.snakefile """
            """--config alignments=$ALIGNMENTS_PATH """
            """reference=$REFERENCE_PATH """
            """reads=$INPUT_READS_PATH """
            """region={wildcards.region} """
            """log=$LOG_FILE """
            """mapping_quality={MAPPING_QUALITY} """
            """asm_alignment_parameters=\\"{ASM_ALIGNMENT_PARAMETERS}\\" """  # Quote BLASR parameters to prevent Snakemake from interpreting them
            """>> $JOB_LOG 2>&1; """
        """echo "Snakemake complete" >> $JOB_LOG; """
        """popd; """

        # Copy alignment and alignment log
        """echo -n "$ALIGN_FILE exists: " >>$JOB_LOG; """
            """if [[ -e $ALIGN_FILE ]]; then echo "True" >>$JOB_LOG; else echo "False" >>$JOB_LOG; fi; """
        """rsync -a $ALIGN_FILE $ALIGN_LOG_FILE $(dirname {output.sam})/ 2>>$JOB_LOG; """

        # Cleanup
        """rm -rf {TMP_DIR}/{wildcards.region}; """
