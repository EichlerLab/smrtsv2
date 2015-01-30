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

# Merge filtered candidates.
rule merge_filtered_candidates:
    input: "filtered_candidates_for_{event_type}.bed"
    output: "merged_filtered_candidates_for_{event_type}.bed"
    params: sge_opts="", merge_distance="1"
    shell: "bedtools merge -i {input} -d {params.merge_distance} > {output}"

# Filter candidates by support and coverage.
rule filter_candidates:
    input: "coverage_and_merged_support_for_{event_type}.bed"
    output: "filtered_candidates_for_{event_type}.bed"
    params: sge_opts="", min_length="50", min_coverage="5", max_coverage="100"
    shell: "awk '$3 - $2 >= {params.min_length} && $9 >= {params.min_coverage} && $9 <= {params.max_coverage}' {input} > {output}"
    #shell: "python scripts/FilterSupport.py --Mc {params.max_coverage} {input} > {output}"

# Annotate merged gap support with alignment coverage.
rule annotate_coverage_of_merged_gap_support:
    input: support="merged_support_for_{event_type}.bed", coverage="coverage.bed"
    output: "coverage_and_merged_support_for_{event_type}.bed"
    params: sge_opts=""
    shell: "bedtools intersect -a {input.support} -b {input.coverage} -sorted -wao | groupBy -i stdin -g 1,2,3,4,5,6,7,8 -c 12 -o mean > {output}"

# Merge gap support for each type of event.
rule merge_gap_support_from_aligned_reads:
    input: dynamic("aligned_reads_{{event_type}}/{batch_id}.bed")
    output: "merged_support_for_{event_type}.bed"
    params: sge_opts=""
    shell: "sort -k 1,1 -k 2,2n -m {input} | python scripts/MergeGapSupport.py > {output}"

# Classify insertions and deletions into their own output files.
rule classify_gaps_in_aligned_reads:
    input: "gaps_in_aligned_reads/{batch_id}.bed"
    output: "aligned_reads_{event_type}/{batch_id}.bed"
    params: sge_opts=""
    shell: """awk '$4 == "{wildcards.event_type}"' {input} > {output}"""

# Parse CIGAR string of aligned reads for insertions and deletions.
rule find_gaps_in_aligned_reads:
    input: alignments="alignments/{batch_id}.bam", reference=config["reference"]["assembly"]
    output: "gaps_in_aligned_reads/{batch_id}.bed"
    params: sge_opts="", mapping_quality_threshold=str(config["alignment"]["mapping_quality"])
    shell:
        "samtools view -h -q {params.mapping_quality_threshold} {input.alignments} "
            "| python scripts/PrintGaps.py {input.reference} /dev/stdin --tsd 0 --condense 20 "
            "| sort -k 1,1 -k 2,2n > {output}"

# Calculate coverage from all alignments.
rule calculate_coverage:
    input: dynamic("alignments/{batch_id}.bam")
    output: "coverage.bed"
    params: sge_opts="-l mfree=8G"
    run:
        bam_inputs = " -in ".join(input)
        command = "~mchaisso/projects/mcst/coverage {output} -in %s" % bam_inputs
        print(command)
        shell(command)

# Create a list of BAM files for downstream analysis.
rule collect_alignment_summaries:
    input: dynamic("alignment_lengths/{batch_id}.tab")
    output: "alignment_lengths.tab"
    params: sge_opts=""
    shell: "cat {input} > {output}"

# Summarize alignments by length.
rule get_alignment_lengths:
    input: "alignments/{batch_id}.bam"
    output: "alignment_lengths/{batch_id}.tab"
    params: sge_opts=""
    shell: """samtools view {input} | awk 'OFS="\\t" {{ if ($3 == "*") {{ print "unmapped",length($10) }} else {{ print "mapped",$9 }} }}' > {output}"""

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
        shell("rm -rf %s; mkdir -p %s" % (output_dir, output_dir))

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
