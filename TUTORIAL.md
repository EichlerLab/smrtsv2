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
python smrtsv.py index reference/chr20.fa
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
python smrtsv.py align reference/chr20.fa input.fofn
```

## Find signatures of variants in raw reads

Find candidate regions to search for SVs based on SV signatures.

```bash
python smrtsv.py detect reference/chr20.fa alignments.fofn candidates.bed
```

## Assemble regions

Assemble local regions of the genome that have SV signatures or tile across the
genome.

```bash
python smrtsv.py assemble reference/chr20.fa input.fofn alignments.fofn candidates.bed local_assembly_alignments.bam
```

## Call variants

Call variants by creating tiles of local assemblies across the reference and
aligning assemblies back to the reference.

```bash
python smrtsv.py call reference/chr20.fa alignments.fofn local_assembly_alignments.bam variants.vcf
```
