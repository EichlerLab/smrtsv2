# SMRT-SV Genotyper

The SMRT-SV genotyper takes a VCF of structural variant (SV) calls and some number of sample alignment files (BAM/CRAM).
In each sample, it searches for evidence of each SV calls and attempts to assign a genotype (0/0, 0/1, or 1/1).

Using the reference FASTA SVs were called against, the genotyper builds an ALT-reference by taking the primary contigs
from the SV reference and adding local-assembly contigs from SMRT-SV as ALT sequences. BWA-MEM is then used to remap all
reads to this ALT reference. Therefore, each SV has both a primary (from primary contigs) and alt (from SV contigs)
representation, and so reads associated with the presence or the absense of an SV have a place to map. This greatly
improves sensitivity over genotyping methods that use the reference without SV ALTs, which is especially true for
insertions. However, this comes at the high cost of remapping all sequence reads to the ALT reference.

## Setup

The SMRT-SV genotyper requires two input files, a genotyper config file and a table of samples. These files may be
called anything, but in the following examples, we will use `genotyper.json` and `samples.tab'.

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
