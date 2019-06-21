# SMRT-SV Genotyper

The SMRT-SV genotyper takes a VCF of structural variant (SV) calls and some number of sample alignment files (BAM/CRAM)
containing short-read Illumina data. In each sample, it searches for evidence of each SV calls and attempts to assign a
genotype (0/0, 0/1, or 1/1).

Using the reference FASTA SVs were called against, the genotyper builds an ALT-reference by taking the primary contigs
from the SV reference and adding local-assembly contigs from SMRT-SV as ALT sequences. BWA-MEM is then used to remap all
short-reads to this ALT reference. Therefore, each SV has both a primary (from primary contigs) and alt (from SV
contigs) representation, and so reads associated with the presence or the absense of an SV have a place to map. This
greatly improves sensitivity over genotyping methods that use the reference without SV ALTs, which is especially true
for insertions. However, this comes at the high cost of remapping all sequence reads to the ALT reference.

## Setup

The SMRT-SV genotyper requires two input files, a genotyper config file and a table of short-read samples. These files
may be called anything, but in the following examples, we will use `genotyper.json` and `samples.tab`.

We recommend running the genotype in its own directory separate from other project files. Do no run in the SMRT-SV
install directory. To setup, change to a clean location with no files, create the configuration and sample table in
that location, and run the genotyper from there.

### Genotyper config

Example configuration file:
```
{
  "sample_manifest": "samples.tab",
  "sv_reference": "/path/to/hg38.no_alt.fa",
  "sv_calls": "/path/to/sv_calls.vcf.gz",
  "sv_contigs": "/path/to/assemble/local_assemblies.bam",
  "model": "30x-4",
  "min_call_depth": 8
}
```

* sample_manifest: Location of the sample table (see next section)
* sv_reference: Reference FASTA SVs were called against
* sv_calls: VCF of SV calls.
* sv_contigs: A BAM of SV-containing contigs
* model: Genotyping model
* min_call_depth: Minimum breakpoint alignment depth over both reference and alternate breakpoints needed to make a call
  * SVs that do not have this many reads over breakpoints are given a no-call genotype (`./.`) 

The VCF must be output by SMRT-SV or be similar to SMRT-SV. It requires INFO fields `CONTIG`, `CONTIG_START`, and
`CONTIG_END`. These are the name of the local assembly containing the SV (`CONTIG`) and BED-like coordinates of the
SV breakpoints on that contig (`CONTIG_START` and `CONTIG_END`). These contigs must be found in the BAM referenced by
`sv_contigs`.

Currently, the reference file must be an uncompressed FASTA. We are expecting to add support for gzipped and bgzipped
FASTA files, but that will currently not work correctly.

The genotyping model is the name of a directory in `files/gtmodel`. Several models were trained at different read depths
and with no-call cut-off criteria. For high-coverage data, we used `30x-4` with a `min_call_depth` of `8`. For
low-coverage 1000 Genomes data, we used `8x2` with a `min_call_depth` of `2`.

### Sample table

The sample table contains each sample to be genotyped, sample sex, and paths to their alignment files (BAM or CRAM).

Example sample table:
```
SAMPLE	SEX	DATA
sample1	F	/path/to/sample1.cram
sample2	M	/path/to/sample2.cram
sample3	M	/path/to/sample3.cram
sample4	U	/path/to/sample4.cram
```

Columns:
* SAMPLE: Name of the sample. These become sample columns in the output VCF.
* SEX: Male (M), Female (F), or unknown (U). Used to correct genotypes on sex chromosome.
* DATA: Path to an alignment file (BAM or CRAM)

To support aligned CRAM files, an optional `REF` column may be specified to point to the reference each sample is
to. If the `REF` column is present and the value is not `NA` for a sample, then that reference is used to read the
alignment file. If the reference has been seen by samtools or htslib before, it may be cached in
`~/.cache/hts-ref`, and in this case, the `REF` field may not be required for to read these CRAM files.

## Run

To run, call `smrtsv` from the install directory and run the `genotype` module. It will need the 

`/path/to/smrtsv genotype <genotyper.json> `

Example:
```
SMRTSV_DIR=/path/to/smrtsv2

${SMRTSV_DIR}/smrtsv --jobs 10 genotype genotyper.json variants.vcf.gz
```


Example distributed:
```
SMRTSV_DIR=/path/to/smrtsv2

${SMRTSV_DIR}/smrtsv --cluster-config=${SMRTSV_DIR}/cluster.eichler.json --drmaalib /path/to/libdrmaa.so.1.0 --distribute --jobs 20 genotype genotyper.json variants.vcf.gz
```

## Results

All VCF fields up to the `FILTER` column are copied from the input VCF. `INFO` is altered after building the VCF using
`vcffixup` in `vcflib`.

The `FORMAT` of each sample is `GT:GQ:GL:DPR:DPA`

* GT: Genotype.
  * ./.: No call. Insufficient read depth over breakpoints to make a genotype call.
  * 0/0: Homozygous reference (absent).
  * 0/1: Heterozygous (present in one copy).
  * 1/1: Homozygous (present in two copies).
* GQ: Genotype quality computed as a PHRED-scaled score of genotype density (see below).
* GL: Genotype likelihood based on the relative genotype density predicted by the genotype model (see below).
* DPR: Read read depth over reference breakpoints.
* DPA: Read depth over alternate breakpoints on SV contig.


Note: For male samples, `0/1` will never be observed on chromosome X. However, X is present in only one copy, so
allele number calculations should adjust for this.

### GQ and GL FORMAT fields

Given a set of features (breakpoint depth, split reads, intsert size, etc), the relative density of each genotype call
(0/0, 0/1, and 1/1) is computed by the machine learning model, and the sum of these three densities adds up to 1. This
is encoded in the GL field as a comma-separated list for 0/0, 0/1, and 1/1, respectively. For example,
`0.0441,0.8809,0.0750` is a 0/1 call with a relative density of 0.8809 (compared to 0.0441 for 0/0 and 0.0750 for 1/1).

The genotype quality is a Phred-scaled score of the maximum likelihood. It is computed as `10 * -log10(1 - D)`, where
`D` is relative density of the accepted genotype call (0.8809 in the example above).
