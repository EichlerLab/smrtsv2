# SMRT-SV

Structural variant (SV) and indel caller for PacBio reads based on methods from
[Chaisson et
al. 2014](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature13907.html).

## What's new in SMRT-SV

SMRT-SV provides an official software package for tools described in [Chaisson et
al. 2014](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature13907.html) and adds several key features including the following.

  * Unified variant calling user interface with built-in cluster compute support
  * Small indel calling (2-49 bp)
  * Improved inversion calling (`screenInversions`)
  * Quality metric for SV calls based on number of local assemblies supporting each call
  * Higher sensitivity for SV calls using tiled local assemblies across the entire genome instead of "signature" regions
  * Genotyping of SVs with Illumina paired-end reads from WGS samples

## Installation

SMRT-SV requires git, Python (2.6.6 or later) and Perl (5.10.1 or later) for
installation.

SMRT-SV has been tested on CentOS 6.8 and should work with most Linux-style
distributions.

### Get the code

Clone the repository into your desired installation directory and build SMRT-SV
dependencies.

```bash
mkdir /usr/local/smrtsv
cd /usr/local/smrtsv
git clone --recursive https://github.com/EichlerLab/pacbio_variant_caller.git .
make
```

Note that some dependencies (e.g., RepeatMasker) require hardcoded paths to this
installation directory. If you need to move SMRT-SV to another directory, it is
easier to change to that directory, clone the repository, and rebuild the
dependencies there.

### Test installation

Add the installation directory to your path.

```bash
export PATH=/usr/local/smrtsv:$PATH
```

Print SMRT-SV help to confirm installation.

```bash
smrtsv.py --help
```

Alternately, run `smrtsv.py` directly from the installation directory.

```bash
/usr/local/smrtsv/bin/smrtsv.py --help
```

## Configure distributed environment

SMRT-SV uses DRMAA to submit jobs to a grid-engine-style cluster. To enable the `--distribute` option of SMRT SV, add the following line to your `.bash_profile` with the correct path to the DRMAA library for your cluster.

```bash
export DRMAA_LIBRARY_PATH=/opt/uge/lib/lx-amd64/libdrmaa.so.1.0
```

Alternately, provide the path to your DRMAA library with the SMRT-SV
`--drmaalib` option.

Additionally, you may need to configure resource requirements depending on your
cluster and PacBio data. Use the `--cluster_config` option when running SMRT-SV
to pass a JSON file that specifies [Snakemake-style cluster
parameters](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-cluster-configuration). An
example configuration used to run SMRT-SV with human genomes on the Eichler lab
cluster is provided in this repository in the file `cluster.eichler.json`.

## Tutorial for variant calling

The following tutorial shows how to call structural variants and indels in
yeast.

### Download PacBio reads

```bash
# List of AWS-hosted files from PacBio including raw reads and an HGAP assembly.
wget https://gist.githubusercontent.com/pb-jchin/6359919/raw/9c172c7ff7cbc0193ce89e715215ce912f3f30e6/gistfile1.txt

# Keep only .xml, .bas.h5, and .bax.h5 files.
sed '/fasta/d;/fastq/d;/celera/d;/HGAP/d' gistfile1.txt > gistfile1.keep.txt

# Download data into a raw reads directory.
mkdir -p raw_reads
cd raw_reads
for f in `cat ../gistfile1.keep.txt`; do wget --force-directories $f; done

# Create a list of reads for analysis.
cd ..
find ./raw_reads -name "*.bax.h5" -exec readlink -f {} \; > reads.fofn
```

### Prepare the reference assembly

Download the reference assembly (sacCer3) from UCSC.

```bash
mkdir -p reference
cd reference
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
```

Unpack the reference tarball and concatenate individual chromosome files into a
single reference FASTA file.

```bash
tar zxvf chromFa.tar.gz
cat *.fa > sacCer3.fasta
rm -f *.fa *.gz
cd ..
```

Prepare the reference sequence for alignment with PacBio reads. This step
produces suffix array and ctab files used by BLASR to speed up alignments.

```bash
smrtsv.py index reference/sacCer3.fasta
```

### Align reads to the reference

Align reads to the reference with BLASR.

```bash
smrtsv.py align reference/sacCer3.fasta reads.fofn
```

### Find signatures of variants in raw reads

Find candidate regions to search for SVs based on SV signatures.

```bash
smrtsv.py detect reference/sacCer3.fasta alignments.fofn candidates.bed
```

### Assemble regions

Assemble local regions of the genome that have SV signatures or tile across the
genome.

```bash
smrtsv.py assemble reference/sacCer3.fasta reads.fofn alignments.fofn candidates.bed local_assembly_alignments.bam
```

### Call variants

Call variants by aligning tiled local assemblies back to the
reference. Optionally, specify the sample name for annotation of the final VCF
file and a species name (common or scientific as supported by
[RepeatMasker](http://www.repeatmasker.org/)) for repeat masking of structural
variants.

```bash
smrtsv.py call reference/sacCer3.fasta alignments.fofn local_assembly_alignments.bam variants.vcf --sample UCSF_Yeast9464 --species yeast
```

## Genotyping

After discovery of SVs with SMRT-SV, use SMRT Genotyper to determine whether
those SVs are present in one or more Illumina-sequenced samples. The genotyper
provides homozygous reference, heterozygous, and homozygous alternate genotypes
for each SV when 5 or more reads are present at any of the SV breakpoints.

To run the genotyper, first prepare a configuration file that looks like the
following example.

```json
{
    "homozygous_binomial_probability": 0.95,
    "heterozygous_binomial_probability": 0.5,
    "sample_manifest": "/home/jlhudd/samples.tab",
    "local_assembly_alignments": "/home/jlhudd/CHM1_local_assembly_alignments.bam",
    "sv_calls": "/home/jlhudd/CHM1_variants.vcf.gz",
    "sv_reference": "/home/jlhudd/ucsc.hg38.no_alts.fasta",
    "sv_reference_lengths": "/home/jlhudd/ucsc.hg38.no_alts.fasta.fai",
    "bam_reference": {
        "human_1kg_v37": "/home/jlhudd/human_1kg_v37.fasta",
        "hg38": "/home/jlhudd/ucsc.hg38.no_alts.fasta",
    },
    "default_bam_reference": "human_1kg_v37",
    "sample_bam_reference": {
    },
    "samples": {
        "CHM1": "/home/jlhudd/CHM1_illumina_reads.bam",
        "CHM13": "/home/jlhudd/CHM13_illumina_reads.bam"
    }
}
```

The parameters in this JSON file are described in the table below.

Parameter | Description
--------- | -----------
homozygous_binomial_probability | the probability to use in the binomial probability calculation for the heterozygous genotype state
heterozygous_binomial_probability | the probability to use in the binomial probability calculation for the heterozygous genotype state
sample_manifest | a headered tab-delimited manifest with "sample" and "sex" columns for each sample being genotyped where the sample name must match the sample name in the corresponding BAM's read group
local_assembly_alignments | the absolute path to a BAM file containing BLASR alignments of local assemblies to the SV reference
sv_calls | a VCF of variants including SVs (insertions and deletions >=50 bp)
sv_reference | the absolute path to the FASTA for the reference used to call SVs
sv_reference_lengths | the absolute path to the FASTA index (.fai) for the reference used to call SVs or chromInfo.txt file
bam_reference | a dictionary of reference names and absolute paths to their corresponding FASTA sequence and BWA index
default_bam_reference | the name of the reference to use by default when one isn't specified for a sample
sample_bam_reference | a dictionary of sample names and their corresponding reference names if they differ from the default reference
samples | a dictionary of sample names and absolute paths to BAMs containing paired-end Illumina sequences for each sample

Finally, genotype SVs using the configuration file and specifying the name of
the final compressed VCF with genotypes.

```bash
smrtsv.py genotype genotyper.config.json genotypes.vcf.gz
```

Note that the genotyper assumes that:

 1. input genomes are in BAM format with alignments generated by BWA MEM against an existing reference assembly
 2. BAMs have the sample tag ("SM") defined in the read group
