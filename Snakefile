"""
Structural variant caller for PacBio reads.

See also: https://github.com/EichlerLab/pacbio_variant_caller
"""
import math
import os

#
# Define internal constants.
#
BLASR_BIN = "/net/eichler/vol20/projects/pacbio/nobackups/users/jlhudd/blasr_jlhudd/alignment/bin/blasr"
CWD = os.getcwd()

#
# Load user variables.
#
configfile: "config.json"
TMP_DIR = config["tmp_dir"]

#
# Define rules
#

# Create a list of BAM files for downstream analysis.
rule collect_alignments:
    input: dynamic("alignments/{batch_id}.bam")
    output: "alignments.txt"
    params: sge_opts=""
    shell: "echo {input} > {output}"

# Sync input reads and reference assembly to local disk, align reads, sort
# output, and write final BAM to shared disk.
rule align_reads:
    input: reads="batched_reads/{batch_id}.fofn", reference=config["reference"]["assembly"], suffix=config["reference"]["suffix_array"], ctab=config["reference"]["ctab"]
    output: "alignments/{batch_id}.bam"
    params: sge_opts="-l disk_free=70G -l mfree=3G -pe serial 12 -N align_batch_{batch_id}", threads="8", samtools_threads="4", samtools_memory="4G"
    shell:
        "mkdir -p {TMP_DIR};"
        "cd {TMP_DIR};"
        "{BLASR_BIN} {CWD}/{input.reads} {input.reference} -out /dev/stdout -sam -sa {input.suffix} -ctab {input.ctab} -nproc {params.threads} -bestn 2 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 50 -clipping subread | samtools sort -@ {params.samtools_threads} -m {params.samtools_memory} -O bam -T {wildcards.batch_id} -o {wildcards.batch_id}.bam -;"
        "rsync --bwlimit=20000 --remove-source-files -W {wildcards.batch_id}.bam {CWD}/{output}"

# Divide input reads into batches for alignment.
rule assign_batches:
    input: config["input"]["reads"]
    output: dynamic("batched_reads/{batch_id}.fofn")
    params: sge_opts=""
    run:
        shell("mkdir -p batched_reads")

        with open(input[0], "r") as fh:
            input_files = [line for line in fh]

        # Total batches is either the number of batches requested or number of
        # input files (when fewer files exist than batches requested).
        total_input_files = len(input_files)
        files_per_batch = config["alignment"]["files_per_batch"]

        for files_processed in range(total_input_files):
            batch_id = math.floor(float(files_processed) / files_per_batch)
            current_output = open("batched_reads/%s.fofn" % batch_id, "a")
            current_output.write(input_files[files_processed])
            current_output.close()
