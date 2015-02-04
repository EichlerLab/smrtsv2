# PacBio Variant Caller

Structural variant caller for PacBio reads based on methods from [Chaisson et
al. 2014](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature13907.html).

## Install dependencies

The PacBio variant caller has the following dependencies:

  - anaconda (>= 2.1.0)
  - bamtools (>= 2.3.0)
  - bedtools (>= 2.21.0)
  - R (>= 3.1.0)
  - samtools (>= 1.1)
  - snakemake (>= 3.2.1)

If modules are already installed for these tools, source the configuration file.

```bash
. config.sh
```

## Build binaries

Build binaries for custom tools that detect variants in SAM alignments.

```bash
pushd scripts/mcst; make; popd
```

## Prepare a reference

Select a reference assembly to use for variant calling. Prepare a suffix array
and ctab file for this assembly with the following commands.

```bash
export PACBIO_DIR=/net/eichler/vol20/projects/pacbio/opt/smrtanalysis/current/analysis/bin
${PACBIO_DIR}/printTupleCountTable ucsc.hg38.no_alts.fasta > ucsc.hg38.no_alts.fasta.ctab
${PACBIO_DIR}/sawriter ucsc.hg38.no_alts.fasta.sa ucsc.hg38.no_alts.fasta
```

## Create a manifest of input reads

Create a manifest of input reads (.bax.h5 files) in a text file with one
absolute path per file per line. For example, create a file named `input.fofn`
with the following lines.

```
/nfs/shared/CHM1/PacBio/m130928_232712_42213_c100518541910000001823079209281310_s1_p0.1.bax.h5
/nfs/shared/CHM1/PacBio/m130928_232712_42213_c100518541910000001823079209281310_s1_p0.2.bax.h5
/nfs/shared/CHM1/PacBio/m130928_232712_42213_c100518541910000001823079209281310_s1_p0.3.bax.h5
```

## Configure variant calling parameters.

Modify `config.json` to include the correct paths to the reference assembly's
files and input reads manifest. Also indicated how many bytes of input reads to
include per alignment batch. The default is 30GB.

## Align reads and call variants

The following command will align PacBio reads, parse alignments to identify SV
candidates, and produce a list of candidate regions for local assembly using no
more than 20 CPUs at any given time.

```bash
snakemake --cluster "qsub {params.sge_opts}" -j 20 assembly_candidates.bed
```
