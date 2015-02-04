# PacBio Variant Caller

Structural variant caller for PacBio reads based on methods from [Chaisson et
al. 2014](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature13907.html).

## Install dependencies

The PacBio variant caller has the following dependencies:

  - anaconda (>= 2.1.0)
  - bamtools (>= 2.3.0)
  - bedtools (>= 2.21.0)
  - samtools (>= 1.1)
  - snakemake (>= 3.2.1)

If modules are already installed for these tools, source the configuration file.

    . config.sh

## Prepare a reference

Select a reference assembly to use for variant calling. Prepare a suffix array
and ctab file for this assembly with the following commands.

    export PACBIO_DIR=/net/eichler/vol20/projects/pacbio/opt/smrtanalysis/current/analysis/bin
	$(PACBIO_DIR)/printTupleCountTable ucsc.hg38.no_alts.fasta > ucsc.hg38.no_alts.fasta.ctab
	$(PACBIO_DIR)/sawriter ucsc.hg38.no_alts.fasta.sa ucsc.hg38.no_alts.fasta

## Create a manifest of input reads

Create a manifest of input reads (.bax.h5 files) in a text file with one
absolute path per file per line.

## Configure variant calling parameters.

Modify ''config.json'' to include the correct paths to the reference assembly's
files and input reads manifest. Also indicated how many bytes of input reads to
include per alignment batch. The default is 30GB.

## Align reads and call variants

The following command will align PacBio reads, parse alignments to identify SV
candidates, and produce a list of candidate regions for local assembly using no
more than 20 CPUs at any given time.

    snakemake --cluster "qsub {params.sge_opts}" -j 20 assembly_candidates.bed
