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
    input: dynamic("gaps_in_aligned_reads/{batch_id}.bed")
    output: "alignments.txt"
    params: sge_opts=""
    shell: "echo {input} > {output}"

rule find_gaps_in_aligned_reads:
    input: alignments="alignments/{batch_id}.bam", reference=config["reference"]["assembly"]
    output: "gaps_in_aligned_reads/{batch_id}.bed"
    params: sge_opts=""
    shell:
        "samtools view -h -F 0x4 {input.alignments} "
            "| python scripts/PrintGaps.py {input.reference} /dev/stdin --tsd 10 --condense 20 "
            "| python scripts/rmdup.py /dev/stdin /dev/stdout "
            "| sort -k 1,1 -k 2,2n > {output}"

# Sync input reads and reference assembly to local disk, align reads, sort
# output, and write final BAM to shared disk.
rule align_reads:
    input: reads="batched_reads/{batch_id}.fofn", reference=config["reference"]["assembly"], suffix=config["reference"]["suffix_array"], ctab=config["reference"]["ctab"]
    output: "alignments/{batch_id}.bam"
    params: sge_opts="-l disk_free=70G -l mfree=3G -pe serial 12 -N align_batch_{batch_id}", threads="8", samtools_threads="4", samtools_memory="4G"
    shell:
        "mkdir -p {TMP_DIR}/{wildcards.batch_id};"
        "cd {TMP_DIR}/{wildcards.batch_id};"
        "rsync --bwlimit=20000 -LW --no-relative --files-from={CWD}/{input.reads} / .;"
        """find ./ -name "*.bax.h5" > input.fofn;"""
        "{BLASR_BIN} input.fofn {input.reference} -unaligned /dev/null -out /dev/stdout -sam -sa {input.suffix} -ctab {input.ctab} -nproc {params.threads} -bestn 2 -maxAnchorsPerPosition 100 -advanceExactMatches 10 -affineAlign -affineOpen 100 -affineExtend 0 -insertion 5 -deletion 5 -extend -maxExtendDropoff 50 -clipping subread | samtools sort -@ {params.samtools_threads} -m {params.samtools_memory} -O bam -T {wildcards.batch_id} -o {wildcards.batch_id}.bam -;"
        "samtools index {wildcards.batch_id}.bam;"
        "rsync --bwlimit=20000 --remove-source-files -W {wildcards.batch_id}.bam* {CWD}/`dirname {output}`/;"
        "rm -rf {TMP_DIR}/{wildcards.batch_id}"

# Divide input reads into batches for alignment.
rule assign_batches:
    input: config["input"]["reads"]
    output: dynamic("batched_reads/{batch_id}.fofn")
    params: sge_opts=""
    run:
        output_dir = os.path.dirname(output[0])
        shell("mkdir -p %s" % output_dir)

        with open(input[0], "r") as fh:
            input_files = [os.path.realpath(line.rstrip()) for line in fh]

        sizes_by_file = dict([(fofn, os.path.getsize(fofn))
                              for fofn in input_files if os.path.exists(fofn)])

        # Assign batches based on the total size of the files per batch. This
        # should produce roughly equally-sized output files.
        size_per_batch = float(config["alignment"]["size_per_batch"])
        total_file_sizes = 0

        for fofn, fofn_size in sizes_by_file.items():
            batch_id = math.floor(total_file_sizes / size_per_batch)
            current_output = open("%s/%s.fofn" % (output_dir, batch_id), "a")
            current_output.write("%s\n" % fofn)
            current_output.close()
            total_file_sizes += fofn_size
