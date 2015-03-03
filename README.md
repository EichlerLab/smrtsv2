# PacBio Variant Caller

Structural variant caller for PacBio reads based on methods from [Chaisson et
al. 2014](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature13907.html).

## Install dependencies

The PacBio variant caller has the following dependencies:

  - anaconda (>= 2.1.0)
  - bamtools (>= 2.3.0)
  - bedtools (>= 2.21.0)
  - BioPython (>= 1.6.5)
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

Select a reference assembly to use for variant calling and acquire the FASTA
sequence for that reference. Prepare a samtools FASTA index, suffix array, and
tuple count table (ctab) for this sequence with the following command. This step
creates the appropriate configuration file entries for the reference in the file
``reference_config.json``. Add these lines to your existing config file based on
the template in this repository.

```bash
snakemake prepare_reference --config reference_fasta=/path/to/reference.fasta
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

## Align reads

Align PacBio reads to the prepared reference.

```bash
snakemake align_reads
```

## Identify SV candidate regions

Parse alignments to identify SV candidates and produce a list of candidate
regions for local assembly using no more than 20 CPUs at any given time. All
output will be written to the current working directory. The `-w` flag instructs
Snakemake to wait 20 seconds for output files when there is excessive latency on
the output filesystem.

```bash
snakemake --cluster "qsub {params.sge_opts}" -w 20 -j 20 assembly_candidates.bed
```

The caller assumes that the analysis will need to start with raw input reads
prior to alignment against the reference. The complete pipeline can be visualize
through Snakemake's DAG as follows.

![alt text](https://raw.githubusercontent.com/EichlerLab/pacbio_variant_caller/master/pipeline.png?token=AAFNfD2eMQivs1n82iTL6pD7Mo5TorUuks5U27iIwA%3D%3D "Snakemake DAG for the SV caller pipeline")
