# How to call SVs in chr20 of GRCh38 with CHM1 P5/C3 data

## Prepare the reference assembly

Prepare a configuration file using the given template.

```bash
cp config.template.json config.json
```

Download and unpack the sequence for chr20 of the GRCh38 reference assembly from
[UCSC's Genome Browser](http://genome.ucsc.edu/).

```bash
mkdir -p reference
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz -O reference/chr20.fa.gz
gunzip reference/chr20.fa.gz
```

Prepare the reference sequence for alignment with PacBio reads. This step
produces suffix array and ctab files used by BLASR to speed up alignments.

```bash
snakemake prepare_reference --config reference=reference/chr20.fa
```

## Download PacBio reads for CHM1 (P5/C3 chemistry)

Download PacBio reads for CHM1 from the [Sequence Read Archive
(SRA)](http://www.ncbi.nlm.nih.gov/sra/). For this tutorial, we will use a small
subset of the [full 54-fold
coverage](http://www.ncbi.nlm.nih.gov/Traces/sra/?study=SRP043558).

```bash
mkdir -p data

wget http://sra-download.ncbi.nlm.nih.gov/srapub_files/SRR1304331_SRR1304331_hdf5.tgz -P data
wget http://sra-download.ncbi.nlm.nih.gov/srapub_files/SRR1304332_SRR1304332_hdf5.tgz -P data
wget http://sra-download.ncbi.nlm.nih.gov/srapub_files/SRR1304333_SRR1304333_hdf5.tgz -P data
```

Unpack reads.

```bash
tar zxvf data/*.tgz
```

Create a list of input reads for analysis.

```bash
find data/ -name "*.bax.h5" -exec readlink -f {} \; > input.fofn
```

## Align reads to the reference

Align reads to the reference with BLASR.

```bash
snakemake align_reads --config reads=input.fofn reference=reference/chr20.fa alignments=alignments.fofn
```
