"""
Rules to identify SV candidates from read alignments confirmation by for local
assembly.
"""

localrules: detect_get_regions

include: 'include.snakefile'

import os


####################
### Declarations ###
####################

# Get the region filter file
REGIONS_TO_EXCLUDE = config.get('regions_to_exclude', None)

if REGIONS_TO_EXCLUDE == 'None':
    REGIONS_TO_EXCLUDE = None


# Get a list of alignment batches
BATCHES = list()

with open('align/alignments.fofn', 'r') as in_file:
    for line in in_file:
        line = line.strip()

        if not line:
            continue

        BATCHES.append(os.path.basename(line).rstrip('.bam'))


#############
### Rules ###
#############

#
# Plots
#


# plot_candidate_summary
#
# Plot candidate summary.
rule plot_candidate_summary:
    input:
        tab='detect/gaps/reads/plot/candidate_summary.tab'
    output:
        lengths='detect/gaps/plot/sv_candidate_lengths.pdf',
        support='detect/gaps/plot/sv_candidate_support.pdf'
    shell:
        """Rscript {SMRTSV_DIR}/scripts/detect/plot_SV_candidate_summary.R {input.tab} {output.lengths} {output.support}"""

# detect_gaps_ref_make_candidates_summary
#
# Summarize filtered candidates by event attributes.
rule detect_gaps_ref_make_candidates_summary:
    input:
        bed=expand('detect/gaps/gaps_with_cov_{svtype}_filtered.bed', svtype=INSDEL)
    output:
        tab='detect/gaps/reads/plot/candidate_summary.tab'
    shell:
        """awk 'OFS="\\t" {{ if (NR == 1) {{ print "event_type","mean_length","support" }} print $4,$5,$6 }}' {input.bed} """
        """>{output.tab}"""


#
# Get assembly candidate regions
#

# detect_get_regions
#
# Make BED of candidate regions.
rule detect_get_regions:
    input:
        bed='detect/candidates/assembly_candidates_with_coverage.bed'
    output:
        bed='detect/candidates.bed'
    params:
        min_coverage=str(config.get('min_coverage')),
        max_coverage=str(config.get('max_coverage')),
        max_length=str(config.get('max_candidate_length'))
    run:
        if REGIONS_TO_EXCLUDE is not None:
            shell(
                """bedtools intersect -a {input.bed} -b {REGIONS_TO_EXCLUDE} -wa -v | """
                """awk '$4 >= {params.min_coverage} && $4 <= {params.max_coverage} && $3 - $2 <= {params.max_length}' """
                """>{output.bed}"""
            )
        else:
            shell(
                """awk '$4 >= {params.min_coverage} && $4 <= {params.max_coverage} && $3 - $2 <= {params.max_length}' {input} """
                """> {output}"""
            )

# detect_candidates_annotate_coverage
#
# Annotate assembly candidates with coverage.
rule detect_candidates_annotate_coverage:
    input:
        bed_cand='detect/candidates/assembly_candidates_with_windows.bed',
        bed_cov='detect/coverage/coverage.bed'
    output:
        bed='detect/candidates/assembly_candidates_with_coverage.bed'
    shell:
        """bedtools intersect -a {input.bed_cand} -b {input.bed_cov} -sorted -wao | """
        """awk 'OFS="\\t" {{ if ($7 == ".") {{ $7 = 0 }} print }}' | """
        """groupBy -i stdin -g 1,2,3 -c 7 -o mean """
        """>{output.bed}"""

# detect_candidates_merge_windows
#
# Merge filtered candidates with tiled windows.
rule detect_candidates_merge_windows:
    input:
        bed_detect='detect/candidates/assembly_candidates.bed',
        bed_window='detect/windows/windows.bed'
    output:
        bed='detect/candidates/assembly_candidates_with_windows.bed'
    shell:
        """sort -k 1,1 -k 2,2n -m {input.bed_detect} {input.bed_window} """
        """>{output.bed}"""

# detect_candidates_merge
#
# Merge filtered gaps and hardstops.
rule detect_candidates_merge:
    input:
        bed_gap='detect/gaps/gaps.bed',
        bed_hardstop='detect/stops/hardstops_merged.bed',
        ref_fai='reference/ref.fasta.fai'
    output:
        bed='detect/candidates/assembly_candidates.bed'
    params:
        merge_distance='500',
        slop='10000'
    shell:
        """set -o pipefail; """
        """cut -f 1-3 {input.bed_gap} {input.bed_hardstop} | """
        """sort -k 1,1 -k 2,2n | """
        """bedtools merge -i stdin -d {params.merge_distance} | """
        """bedtools slop -i stdin -g {input.ref_fai} -b {params.slop} """
        """>{output}"""


#
# Windows for tiled assemblies.
#

# detect_windows_create
#
# Create a sliding windows over the genome reference.
rule detect_windows_create:
    input:
        ref_fai='reference/ref.fasta.fai'
    output:
        bed='detect/windows/windows.bed'
    params:
        window=str(config.get('assembly_window_size')),
        slide=str(config.get('assembly_window_slide'))
    shell:
        """bedtools makewindows -g {input.ref_fai} -w {params.window} -s {params.slide} | """
        """sort -k 1,1 -k 2,2n """
        """>{output.bed}"""


#
# Find gaps in reference (hardstops)
#

# detect_stops_filter_merge
#
# Filter hardstops to exclude reference gaps and small SV candidates. Then merge
# adjacent bins that pass filters.
rule detect_stops_filter_merge:
    input:
        bed_hardstop='detect/stops/hardstops/hardstop_breakpoints_binned.bed',
        bed_refgap='detect/stops/hardstops/ref_gap.bed',
        bed_gaps='detect/gaps/candidates.bed'
    output:
        bed='detect/stops/hardstops_merged.bed'
    params:
        min_support=str(config.get('min_hardstop_support'))
    shell:
        """awk '$4 > {params.min_support}' {input.bed_hardstop} | """
        """bedtools window -w 1000 -a stdin -b {input.bed_refgap} -v | """
        """bedtools window -w 1000 -a stdin -b {input.bed_gaps} -v | """
        """bedtools merge -i stdin -d 1 """
        """>{output.bed}"""

# detect_stops_ref_gaps
#
# Find gap bases in the reference assembly to exclude from hardstop collection.
rule detect_stops_ref_gaps:
    input:
        ref_fa='reference/ref.fasta'
    output:
        bed='detect/stops/hardstops/ref_gap.bed'
    shell:
        """python {SMRTSV_DIR}/scripts/detect/find_fasta_gaps.py {input.ref_fa} """
        """>{output.bed}"""

# detect_stops_count_per_bin
#
# Count hardstops per bin over the reference.
rule detect_stops_count_per_bin:
    input:
        bed_stops='detect/stops/hardstops/hardstop_breakpoints.bed',
        bed_bins='detect/stops/hardstop_bins.bed'
    output:
        bed='detect/stops/hardstops/hardstop_breakpoints_binned.bed'
    shell:
        """bedtools intersect -a {input.bed_bins} -b {input.bed_stops} -sorted -c """
        """>{output.bed}"""

# detect_stops_find_breakpoints
#
# Create hardstop breakpoints from hardstop locations (either left, right, or
# both). If a read is clipped on the left, use its start position as the
# breakpoint. If it is clipped on the right, use its end position. If a read is
# clipped from both sides, print two separate breakpoints using these same rules
# for left and right breakpoints.
rule detect_stops_find_breakpoints:
    input:
        bed='detect/stops/hardstops/hardstop.bed'
    output:
        bed='detect/stops/hardstops/hardstop_breakpoints.bed'
    shell:
        """awk 'OFS="\\t" {{\n"""
        """    if ($6 == "left") {{ print $1,$2,$2 + 1,$4,$7,"left" }}\n"""
        """    else if ($6 == "right") {{ print $1,$3 - 1,$3,$4,$8,"right" }}\n"""
        """    else if ($6 == "both") {{ print $1,$2,$2 + 1,$4,$7,"left"; print $1,$3 - 1,$3,$4,$8,"right" }}\n"""
        """}}' {input.bed} | """
        """sort -k 1,1 -k 2,2n """
        """>{output.bed}"""

# detect_stops_merge_batches
#
# Collect gaps in one command
rule detect_stops_merge_batches:
    input:
        bed=expand('detect/stops/hardstops/batches/hardstops_{batch_id}.bed', batch_id=BATCHES)
    output:
        bed='detect/stops/hardstops/hardstop.bed'
    shell:
        """sort -k 1,1 -k 2,2n -m {input.bed} """
        """>{output.bed}"""

# detect_stops_find_hardstops
#
# Parse CIGAR string of aligned reads in one batch for clipped alignments.
rule detect_stops_find_hardstops:
    input:
        bam='align/bam/{batch_id}.bam'
    output:
        bed='detect/stops/hardstops/batches/hardstops_{batch_id}.bed'
    params:
        mapq=str(config.get('mapping_quality')),
        min_clipping='500'
    shell:
        """{SNAKEMAKE_DIR}/scripts/mcst/hardstop {input.bam} {params.mapq} {params.min_clipping} {output.bam}; """
        """sort -k 1,1 -k 2,2n -o {output.bed} {output.bed}"""

# detect_stops_make_bins
#
# Make 500 bp windows across the reference.
rule detect_stops_make_bins:
    input:
        ref_fai='reference/ref.fasta.fai'
    output:
        bed='detect/stops/hardstop_bins.bed'
    params:
        bin_size='500'
    shell:
        """bedtools makewindows -g {input.ref_fai} -w {params.bin_size} | """
        """sort -k 1,1 -k 2,2n """
        """>{output.bed}"""


#
# Find gaps in reads (ins/del in alignments)
#

# detect_gaps_make_candidates
#
# Summarize filtered candidates by event attributes.
rule detect_gaps_make_candidates:
    input:
        bed=expand('detect/gaps/gaps_with_cov_{svtype}_filtered.bed', svtype=INSDEL)
    output:
        bed='detect/gaps/gaps.bed'
    shell:
        """sort -k 1,1 -k 2,2n -m {input.bed} | """
        """cut -f 1-4,6 """
        """>{output.bed}"""

# detect_gaps_filter
#
# Filter candidates by support and coverage.
rule detect_gaps_filter:
    input:
        bed='detect/gaps/gaps_with_cov_{svtype}.bed'
    output:
        bed='detect/gaps/gaps_with_cov_{svtype}_filtered.bed',
        temp=temp('detect/gaps/gaps_with_cov_{svtype}_filtered.temp')
    params:
        min_support=str(config.get('min_support')),
        max_support=str(config.get('max_support')),
        min_length=str(config.get('min_length')),
        min_coverage=str(config.get('min_coverage')),
        max_coverage=str(config.get('max_coverage'))
    run:

        # Filter
        shell(
            """awk '"""
                """$4 >= {params.min_length} && """
                """$5 >= {params.min_support} && """
                """$5 <= {params.max_support} && """
                """$9 >= {params.min_coverage} && """
                """$9 <= {params.max_coverage}"""
            """' {input.bed} """
            """>{output.temp}"""
        )

        # Merge regions
        if os.stat(output.temp).st_size > 0:
            shell("""bedtools merge -i {output.temp} -d 1 -c 6,4,5 -o distinct,mean,mean > {output.bed}""")
        else:
            shell("""touch {output.bed}""")

# detect_gaps_annotate_coverage
#
# Annotate merged gap support with alignment coverage.
rule detect_gaps_annotate_coverage:
    input:
        bed_gap='detect/gaps/gaps_{svtype}.bed',
        bed_cov='detect/coverage/coverage.bed'
    output:
        bed='detect/gaps/gaps_with_cov_{svtype}.bed'
    shell:
        """bedtools intersect -a {input.bed_gap} -b {input.bed_cov} -sorted -wao | """
        """awk 'OFS="\\t" {{ if ($13 == ".") {{ $13 = 0 }} print }}' | """
        """cut -f 1-6,8- | """
        """groupBy -i stdin -g 1,2,3,4,5,6,7,8 -c 12 -o mean """
        """>{output}"""

# detect_gaps_merge_batches
#
# Merge gap support for each variant type over all alignment batches.
rule detect_gaps_merge_batches:
    input:
        bed=expand('detect/gaps/batch/{batch_id}_{{svtype}}.bed', batch_id=BATCHES)
    output:
        bed='detect/gaps/gaps_{svtype}.bed'
    shell:
        """sort -k 1,1 -k 2,2n -m {input.bed} | """
        """python {SMRTSV_DIR}/scripts/PrintGapSupport.py /dev/stdin /dev/stdout | """
        """sort -k 1,1 -k 2,2n -k 3,3n -k 4,4n -k 5,5n -k 6,6 -k 7,7 -k 8,8 -k 9,9 """
        """>{output.bed}"""

# detect_gaps_by_svtype
#
# Separate gaps by variant type.
rule detect_gaps_by_svtype:
    input:
        bed='detect/gaps/batch/{batch_id}_ALL.bed'
    output:
        bed='detect/gaps/batch/{batch_id}_{svtype}.bed'
    shell:
        """awk '$4 == "{wildcards.svtype}"' {input.bed} | """
        """sort -k 1,1 -k 2,2n > {output.bed}"""

# detect_gaps_search
#
# Parse CIGAR string of aligned reads for insertions and deletions.
rule detect_gaps_search:
    input:
        bam='align/bam/{batch_id}.bam'
        ref_fa='reference/ref.fasta',
        ref_fai='reference/ref.fasta.fai',
    output:
        bed=temp('detect/gaps/batch/{batch_id}_ALL.bed')
    log:
        'detect/gaps/log/{batch_id}.log'
    params:
        mapq=str(config.get('mapping_quality'))
    shell:
        """samtools view -F 0x4 -q {params.mapq} {input.bam} | """
        """python {SMRTSV_DIR}/scripts/PrintGaps.py {input.ref_fa} /dev/stdin --tsd 0 --condense 20 """
        """>{output.bed} 2> {log}"""


#
# Get alignment coverage
#

# detect_merge_coverage
#
# Collect coverages from each alignment batch.
rule detect_merge_coverage:
    input:
        bed=expand('detect/coverage/batch/{batch_id}.bed', batch_id=BATCHES)
    output:
        bed='detect/coverage/coverage.bed'
    shell:
        """paste {input.bed} | """
        """awk 'OFS="\\t" {{ sum = 0; for (i = 4; i <= NF; i++) {{ if (i % 4 == 0) {{ sum += $i }} }} print $1,$2,$3,sum }}' | """
        """sort -k 1,1 -k 2,2n > {output.bed}"""

# detect_coverage_per_batch
#
# Calculate coverage from each batch.
rule detect_coverage_per_batch:
    input:
        bam='align/bam/{batch_id}.bam'
    output:
        bed='detect/coverage/batch/{batch_id}.bed'
    shell:
        """{SMRTSV_DIR}/scripts/mcst/coverage {output.bed} -in {input.bam}"""
