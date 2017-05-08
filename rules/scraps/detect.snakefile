#
# Inaccessible regions
#

# identify_inaccessible_regions
#
#
rule identify_inaccessible_regions:
    input:
        coverage='coverage.bed',
        excluded_regions='regions_to_exclude_as_inaccessible.bed'
    output:
        'inaccessible_regions.bed'
    params:
        max_support=str(config.get('max_inaccessible_support'))
    shell:
        """awk '$4 <= {params.max_support}' {input.coverage} | """
        """bedtools merge -i stdin -d 1000 | """
        """bedtools intersect -a stdin -b {input.excluded_regions} -v """
        """>{output}"""

# collect_regions_to_exclude_as_inaccessible_regions
#
#
rule collect_regions_to_exclude_as_inaccessible_regions:
    input:
        'merged_hardstops_per_bin.bed',
        'gaps_in_reference_assembly.bed',
        'filtered_candidates.tab'
    output:
        'regions_to_exclude_as_inaccessible.bed'
    shell:
        """cut -f 1-3 {input} | """
        """sort -k 1,1 -k 2,2n | """
        """bedtools merge -i stdin -d 1 """
        """>{output}"""
