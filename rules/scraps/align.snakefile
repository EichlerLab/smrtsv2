# Plot alignment lengths.
rule plot_alignment_summaries:
    input: "alignment_lengths.tab"
    output: "alignment_lengths.pdf"
    shell: "Rscript plot_read_lengths_by_type.R {input} {output}"


# Collect summary of aligned read lengths.
rule collect_alignment_summaries:
    input: expand("alignment_lengths/{batch_id}.tab", batch_id=BATCHES)
    output: "alignment_lengths.tab"
    shell: """awk 'OFS="\\t" {{ if (NR == 1) {{ print "alignment_status","subread_length","aligned_length" }} print }}' {input} > {output}"""

# Summarize alignments by length.
# TODO: get subread lengths at the same time as the aligned lengths.
rule get_subread_and_alignment_lengths:
    input: "%s/{batch_id}.bam" % ALIGNMENTS_DIR
    output: "alignment_lengths/{batch_id}.tab"
    params: mapping_quality_threshold=config.get("mapping_quality")
    shell:
        "mkdir -p {TMP_DIR}; "
        """samtools view {input} | awk 'OFS="\\t" {{ if ($3 == "*" || $5 >= {params.mapping_quality_threshold}) {{ num_of_pieces = split($1, pieces, "/"); num_of_coords = split(pieces[3], coords, "_"); subread_length = coords[2] - coords[1]; if ($3 == "*") {{ print "unmapped",subread_length,length($10) }} else if ($5 >= 30) {{ print "mapped",subread_length,$9 }} }} }}' > {TMP_DIR}/lengths.`basename {output}`; """
        "rsync --remove-source-files {TMP_DIR}/lengths.`basename {output}` {output}; "
