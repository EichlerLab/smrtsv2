"""
Optional plots for the detect stage.
"""

if not 'INCLUDE_SNAKEFILE' in globals():
    include: '../include.snakefile'


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
