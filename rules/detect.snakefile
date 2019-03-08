"""
Rules to identify SV candidates from read alignments confirmation by for local
assembly.
"""

import os
import pandas as pd
import numpy as np

localrules: detect_get_regions

if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'include.snakefile'


###################
### Definitions ###
###################

# Set list of alignment batches
DETECT_BATCH_LIST = list()

if os.path.isfile('align/alignments.fofn'):
    with open('align/alignments.fofn', 'r') as in_file:
        for line in in_file:
            line = line.strip()

            if not line or line.startswith('#'):
                continue

            DETECT_BATCH_LIST.append(int(os.path.basename(line).rstrip('.bam')))

CANDIDATE_BED_DTYPES = {
    '#CHROM': np.object,
    'POS': np.int64,
    'END': np.int64,
    'ID': np.object,
    'DEPTH': np.float32,
    'GROUP_ID': np.object
}

#############
### Rules ###
#############

#
# Group candidate regions
#

# detect_group_merge_regions
#
# Each candidate region within a group is run as one unit. Reads are extracted from the reference for the group,
# the reads are cached locally, and each candidate assembly will extract reads from that cache. Candidate regions
# are assigned to the best group, but they extend past the ends of the groups. This rule generates a file of full
# regions for each group from the first base of the left-most candidate window to the last base of the right-most
# candidate window.
rule detect_group_merge_regions:
    input:
        bed_can='detect/candidates.bed',
        bed_win='detect/group/regions/windows.bed'
    output:
        bed='detect/candidate_groups.bed'
    run:

        # Read candidates and windows
        df_can = pd.read_table(
            input.bed_can,
            header=0,
            dtype=CANDIDATE_BED_DTYPES,
            usecols=('#CHROM', 'POS', 'END', 'GROUP_ID')
        )

        df_can['#CHROM'] = df_can['#CHROM'].apply(str)

        df_win = pd.read_table(input.bed_win, header=0)

        # Filter windows by those with candidate regions
        group_set = set(df_can['GROUP_ID'].tolist())

        df_win = df_win.loc[df_win['GROUP_ID'].apply(lambda set_name: set_name in group_set)]

        # Merge
        df = pd.concat([df_can, df_win])

        # Group and merge
        df_group = df.groupby('GROUP_ID')

        df_range = pd.concat([
            df_group['#CHROM'].apply(lambda vals: vals.iloc[0]),  # Group is on one chr
            df_group['POS'].min(),
            df_group['END'].max(),
        ], axis=1)

        df_range['GROUP_ID'] = df_range.index

        # Add actual size after merging the candidate regions with groups
        df_range['SIZE'] = df_range.apply(lambda row: row['END'] - row['POS'], axis=1)

        # Write
        df_range.to_csv(output.bed, sep='\t', header=True, index=False)


# detect_group_merge_candidates
#
# Assign each candidate region to a window containing many regions.
rule detect_group_assign_candidates:
    input:
        bed_can='detect/candidates/candidates.bed',
        bed_win='detect/group/regions/windows.bed',
        ref_sizes='reference/ref.fasta.sizes'
    output:
        bed='detect/candidates.bed'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tID\tDEPTH\tGROUP_ID" > {output.bed}; """
        """bedtools closest -a {input.bed_can} -b {input.bed_win} -t first -g {input.ref_sizes} -sorted | """
        """awk -vOFS="\\t" '{{print $1, $2, $3, $4, $5, $9}}' """
        """>>{output.bed}"""

# detect_group_make_windows
#
# Make windows candidate regions will be grouped into. To reduce redundant IO on a cluster, candidate regions will
# grouped into these windows, reads over the window will be queried stored on a compute node, and a local assembly
# will be performed for each candidate within the group.
rule detect_group_make_windows:
    input:
        ref_sizes='reference/ref.fasta.sizes'
    output:
        bed_win='detect/group/regions/windows.bed'
    params:
        group_size=get_config_param('candidate_group_size')
    shell:
        """echo -e "#CHROM\tPOS\tEND\tGROUP_ID" > {output.bed_win}; """
        """bedtools makewindows -g {input.ref_sizes} -w {params.group_size} | """
            """awk -vOFS="\t" '{{print $0, "gr-" $1 "-" $2 "-" ($3 - $2)}}' """
            """>>{output.bed_win}; """


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
        bed='detect/candidates/candidates.bed'
    params:
        exclude_regions=get_config_param('exclude', as_is=True),
        min_coverage=get_config_param('min_coverage'),
        max_coverage=get_config_param('max_coverage'),
        max_length=get_config_param('max_candidate_length')
    run:

        # Filter
        if params.exclude_regions is not None:
            shell(
                """bedtools intersect -a {input.bed} -b {params.exclude_regions} -wa -v | """
                """awk '$4 >= {params.min_coverage} && $4 <= {params.max_coverage} && $3 - $2 <= {params.max_length}' """
                """>{output.bed}.temp"""
            )
        else:
            shell(
                """awk '$4 >= {params.min_coverage} && $4 <= {params.max_coverage} && $3 - $2 <= {params.max_length}' {input} """
                """>{output.bed}.temp"""
            )

        # Add header and id
        shell(
            """awk -vOFS="\\t" '"""
                """BEGIN {{print "#CHROM", "POS", "END", "ID", "DEPTH"}} """
                """{{print $1, $2, $3, $1 "-" $2 "-" ($3 - $2), $4}}"""
            """' {output.bed}.temp """
            """>{output.bed}"""
        )

        # Remove temp
        os.remove(output.bed + '.temp')


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
        bed_hardstop='detect/stops/hardstops.bed',
        ref_sizes='reference/ref.fasta.sizes'
    output:
        bed='detect/candidates/assembly_candidates.bed'
    params:
        merge_distance='500',
        slop='10000'
    shell:
        """cut -f 1-3 {input.bed_gap} {input.bed_hardstop} | """
        """sort -k 1,1 -k 2,2n | """
        """bedtools merge -i stdin -d {params.merge_distance} | """
        """bedtools slop -i stdin -g {input.ref_sizes} -b {params.slop} """
        """>{output}"""


#
# Windows for tiled assemblies.
#

# detect_windows_create
#
# Create a sliding windows over the genome reference.
rule detect_windows_create:
    input:
        ref_sizes='reference/ref.fasta.sizes'
    output:
        bed='detect/windows/windows.bed'
    params:
        window=get_config_param('assembly_window_size'),
        slide=get_config_param('assembly_window_slide')
    shell:
        """bedtools makewindows -g {input.ref_sizes} -w {params.window} -s {params.slide} | """
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
        bed_gaps='detect/gaps/gaps.bed'
    output:
        bed='detect/stops/hardstops.bed'
    params:
        min_support=get_config_param('min_hardstop_support')
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
        """python2 -s {SMRTSV_DIR}/scripts/detect/find_fasta_gaps.py {input.ref_fa} """
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
        bed=expand('detect/stops/hardstops/batches/{batch_id}.bed', batch_id=DETECT_BATCH_LIST)
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
        bed='detect/stops/hardstops/batches/{batch_id}.bed'
    params:
        mapq=get_config_param('mapping_quality'),
        min_clipping='500'
    log:
        'detect/stops/hardstops/batches/log/{batch_id}.log'
    shell:
        """{SMRTSV_DIR}/scripts/mcst/hardstop {input.bam} {params.mapq} {params.min_clipping} {output.bed} """
        """>{log} 2>&1; """
        """sort -k 1,1 -k 2,2n -o {output.bed} {output.bed}"""

# detect_stops_make_bins
#
# Make 500 bp windows across the reference.
rule detect_stops_make_bins:
    input:
        ref_sizes='reference/ref.fasta.sizes'
    output:
        bed='detect/stops/hardstop_bins.bed'
    params:
        bin_size='500'
    shell:
        """bedtools makewindows -g {input.ref_sizes} -w {params.bin_size} | """
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
        bed_ins='detect/gaps/gaps_ins.bed',
        bed_del='detect/gaps/gaps_del.bed'
    output:
        bed='detect/gaps/gaps.bed'
    shell:
        """sort -k 1,1 -k 2,2n -m {input.bed_ins} {input.bed_del} | """
        """cut -f 1-4,6 """
        """>{output.bed}"""

# detect_gaps_filter
#
# Filter candidates by support and coverage.
rule detect_gaps_filter:
    input:
        bed='detect/gaps/gaps_with_cov_{svtype}.bed'
    output:
        bed='detect/gaps/gaps_{svtype,ins|del}.bed',
        temp=temp('detect/gaps/gaps_{svtype,ins|del}.temp')
    params:
        min_support=get_config_param('min_support'),
        max_support=get_config_param('max_support'),
        min_length=get_config_param('min_length'),
        min_coverage=get_config_param('min_coverage'),
        max_coverage=get_config_param('max_coverage')
    run:

        # Filter
        shell(
            """awk -vOFS="\\t" '("""
                """$1 ~ /^#/ || ("""
                """$4 >= {params.min_length} && """
                """$5 >= {params.min_support} && """
                """$5 <= {params.max_support} && """
                """$9 >= {params.min_coverage} && """
                """$9 <= {params.max_coverage}"""
                """)) {{print}}"""
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
# Annotate gap support with alignment coverage.
rule detect_gaps_annotate_coverage:
    input:
        bed_gap='detect/gaps/gaps_no_cov_{svtype}.bed',
        bed_cov='detect/coverage/coverage.bed'
    output:
        bed='detect/gaps/gaps_with_cov_{svtype,ins|del}.bed'
    shell:
        """bedtools intersect -a {input.bed_gap} -b {input.bed_cov} -sorted -wao | """
        """awk 'OFS="\\t" {{ if ($13 == ".") {{ $13 = 0 }} print }}' | """
        """cut -f 1-6,8- | """
        """bedtools groupby -i stdin -g 1,2,3,4,5,6,7,8 -c 12 -o mean """
        """>{output.bed}"""

# detect_gaps_merge_batches
#
# Merge gaps and annotate the number of reads for each variant type.
rule detect_gaps_merge_batches:
    input:
        bed=expand('detect/gaps/batch/{batch_id}_{{svtype}}.bed', batch_id=DETECT_BATCH_LIST)
    output:
        bed=temp('detect/gaps/gaps_no_cov_{svtype,ins|del}.bed')
    shell:
        """sort -k 1,1 -k 2,2n -m {input.bed} | """
        """python2 -s {SMRTSV_DIR}/scripts/detect/PrintGapSupport.py /dev/stdin /dev/stdout | """
        """sort -k 1,1 -k 2,2n -k 3,3n -k 4,4n -k 5,5n -k 6,6 -k 7,7 -k 8,8 -k 9,9 """
        """>{output.bed}"""

# detect_gaps_by_svtype
#
# Separate gaps by variant type.
rule detect_gaps_by_svtype:
    input:
        bed='detect/gaps/batch/{batch_id}_all.bed'
    output:
        bed=temp('detect/gaps/batch/{batch_id}_{svtype,ins|del}.bed')
    run:
        sv_type = wildcards.svtype.upper()

        shell(
            """awk '($1 ~ /^#/ || $4 == "{sv_type}")' {input.bed} | """
            """sort -k 1,1 -k 2,2n """
            """>{output.bed}"""
        )

# detect_gaps_search
#
# Parse CIGAR string of aligned reads for insertions and deletions.
rule detect_gaps_search:
    input:
        bam='align/bam/{batch_id}.bam',
        ref_fa='reference/ref.fasta',
        ref_fai='reference/ref.fasta.fai',
    output:
        bed='detect/gaps/batch/{batch_id}_all.bed'
    log:
        'detect/gaps/batch/log/{batch_id}.log'
    params:
        mapq=get_config_param('mapping_quality')
    shell:
        """samtools view -F 0x4 -q {params.mapq} {input.bam} | """
        """python2 -s {SMRTSV_DIR}/scripts/PrintGaps.py {input.ref_fa} /dev/stdin --condense 20 """
        """>{output.bed} 2>{log}"""


#
# Get alignment coverage
#

# detect_merge_coverage
#
# Coalesce alignment coverage from each batch.
rule detect_merge_coverage:
    input:
        bed=expand('detect/coverage/batch/{batch_id}.bed', batch_id=DETECT_BATCH_LIST)
    output:
        bed='detect/coverage/coverage.bed'
    shell:
        """echo -e "#CHROM\tPOS\tEND\tDEPTH" > {output.bed}; """
        """paste {input.bed} | """
        """awk 'OFS="\\t" {{ sum = 0; for (i = 4; i <= NF; i++) {{ if (i % 4 == 0) {{ sum += $i }} }} print $1,$2,$3,sum }}' | """
        """sort -k 1,1 -k 2,2n >> {output.bed}"""

# detect_coverage_per_batch
#
# Calculate alignment coverage over the reference on each alignment batch.
rule detect_coverage_per_batch:
    input:
        bam='align/bam/{batch_id}.bam'
    output:
        bed='detect/coverage/batch/{batch_id}.bed'
    log:
        'detect/coverage/batch/log/{batch_id}.log'
    run:

        if os.path.getsize(input.bam) > 0:
            shell(
                """{SMRTSV_DIR}/scripts/mcst/coverage {output.bed} -in {input.bam} >{log} 2>&1"""
            )
        else:
            with open(output.bed, 'w'):
                pass
