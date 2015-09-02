# SMRT SV

Structural variant caller for PacBio reads based on methods from [Chaisson et
al. 2014](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature13907.html).

## Install dependencies

The PacBio variant caller has the following Python dependencies:

  - anaconda (>= 2.1.0)
  - BioPython (>= 1.6.5)
  - intervaltree (>= 2.1.0)
  - snakemake (>= 3.2.1)

Additionally, SMRT SV relies on these external tools:

  - [PacBio SMRT analysis](http://www.pacb.com/devnet/)
  - [BLASR](https://github.com/EichlerLab/blasr) (>= 1.MC.rc42)
  - [PBcR](http://wgs-assembler.sourceforge.net/wiki/index.php/PBcR) (Celera Assembler with MHAP overlapper)
  - bedtools (>= 2.23.0)
  - freebayes (>= 0.9.14)
  - perl (>= 5.14.2)
  - R (>= 3.1.0)
  - RepeatMasker (>= 3.3.0)
  - samtools (>= 1.1)

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
snakemake prepare_reference --config reference=/path/to/reference.fasta
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
include per alignment batch. The default is ~200GB.

## Align reads

Align PacBio reads to the prepared reference. You must specify a file with paths
to PacBio reads and a reference to align the reads to. By default, the reference
is assumed to have a suffix array and ctab file withe same name as the reference
FASTA plus the corresponding extension (e.g., "reference.fasta" has a suffix
array of "reference.fasta.sa" and ctab of "reference.fasta.ctab").

```bash
snakemake align_reads --config reads=input.fofn reference=/path/to/ucsc.hg38.no_alts.fasta
```

By default, all alignments are stored in a single BAM. Alignments can be divided
into any number of batches with the `batches` configuration parameter. The file
containing the final alignment paths can be specified with the `alignments`
parameter. The following example shows the use of both parameters.

```bash
snakemake align_reads --config reads=input.fofn reference=/path/to/ucsc.hg38.no_alts.fasta \
    batches=4 alignments=custom_alignments.fofn alignments_dir=custom_alignments
```

Note that if you already have alignments of PacBio reads to your reference in BAM format, you can skip this step and create an `alignments.fofn` file for the subsequent steps.

```bash
find /path/to/your/alignments -name "*.bam" > alignments.fofn
```

### Alignment parameters

| Parameter | Definition |
| --------- | ---------- |
| reads | a text file of absolute paths to PacBio reads in HDF5 format (i.e., .bax.h5 files) with one path per line |
| reference | a FASTA sequence to align reads to with a .sa and .ctab file in the same directory |
| batches (default: 1) | number of batches to split input reads into such that there will be one BAM output file per batch |
| alignments (default: "alignments.fofn") | name of output file with list of absolute paths to BAM output files |
| alignments_dir (default: "alignments") | name of directory where BAM output files will be written |
| threads (default: 1) | number of threads to use for each BLASR alignment job |

## Identify SV candidate regions

Parse alignments to identify SV candidates and produce a list of candidate
regions for local assembly.

```bash
snakemake get_regions --config alignments=alignments.fofn reference=reference.fasta
```

### Candidate region parameters

| Parameter | Definition |
| --------- | ---------- |
| alignments | a text file of absolute paths to PacBio reads alignments in BAM format |
| reference | a FASTA sequence to align local assemblies to with a .sa and .ctab file in the same directory |
| regions_to_exclude *(optional)* | a BED file of regions to exclude from local assembly (e.g., heterochromatic sequences, gaps, etc.) |

## Assemble candidate regions

Assemble SV candidate regions with MHAP/Celera. The final assemblies and their
alignments against the reference are in `local_assembly_alignments.bam`. After
producing local assemblies, call SVs from assemblies based on gaps in their
alignments back to the reference. The final output is in `sv_calls.bed`.

```bash
snakemake call_variants --config alignments=alignments.fofn
```

### Assembly parameters

| Parameter | Definition |
| --------- | ---------- |
| alignments | a text file of absolute paths to PacBio reads alignments in BAM format |
| reference | a FASTA sequence to align local assemblies to with a .sa and .ctab file in the same directory |
| regions_to_assemble | a BED file of regions to assemble for variant detection |
| log (default: assembly.log) | file to log assembly results to |

## Distributing analyses on a cluster

SMRT SV is designed to run on a grid engine-style cluster like UGE with 1 or
more CPUs. For example, to run the SV signature detection method with 20 cores,
run `snakemake` with the following parameters.

```bash
snakemake --cluster "qsub {params.sge_opts}" -j 20
```

If the DRMAA library is available, use this to get better control over job
execution on the cluster.

```bash
snakemake --drmaa " -q all.q -V -cwd -e ./log -o ./log {params.sge_opts} -w n -S /bin/bash" -j 20 -w 30 align_reads
```
