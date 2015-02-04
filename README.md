# PacBio Variant Caller

Structural variant caller for PacBio reads based on methods from [Chaisson et
al. 2014](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature13907.html).

## Prepare a reference

Select a reference assembly to use for variant calling. Prepare a suffix array
and ctab file for this assembly with the following commands.

    code

## Configure variant calling parameters.

Modify ''config.json''.

## Call variants

The following command will align PacBio reads, parse alignments to identify SV
candidates, and produce a list of candidate regions for local assembly using no
more than 20 CPUs at any given time.

    snakemake --cluster "qsub {params.sge_opts}" -j 20 assembly_candidates.bed
