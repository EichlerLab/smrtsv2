# SMRT-SV long read structural variant caller

SMRT-SV calls structural variants (SVs) using long reads (PacBio RS II or Sequel). The input is an FOFN file that lists
each input file (BAX H5 or BAM) and a genome reference. Reads are aligned to the reference and local assemblies are
done in kilobase-scale overlapping windows across the genome. Assemblies are aligned back to the reference, and
structural variants are called from them.

Because of this approach, SMRT-SV yields both SV calls and contigs containing the SV with know breakpoints on the
contig. These data feed into powerful analysis tools that use the reference and the contigs, such as the SMRT-SV
genotyper, which is also included in this package.


## Short-read genotyper

Variant calls discovered with long reads using the SMRT-SV method can be genotyped in short-reads with the SMRT-SV
genotyper, which is also included in this package. See [GENOTYPE.md](:doc:/GENOTYPE.md) for information about
configuring and running the SMRT-SV genotyper.


## Pipeline process

SMRT-SV is run in several steps:

1) Align
2) Detect SV signatures
3) Assemble
4) Call variants

### Align

All reads are aligned to the reference. Alignments are done in batches, which is controlled by the `--batches` parameter
(default = 20). Batching is useful for distributing alignment workload over a cluster. Alignments are done in the
`alignment` subdirectory, and `align/alignments.fofn` contains a list of aligned BAM files.

### Detect SV signatures

Overlapping windows are tiled across the genome. By default, 60 kbp windows are tiled each offset by 20 kbp.
The pipeline also searches for signatures of structural variation, and it places additional windows around those. The
output from this step goes into the `detect` subdirectory.

### Assemble

For each window, reads are pulled down and assembled with Canu. The assemblies are polished (Arrow or Quiver), and
then aligned back to the reference. This step takes the most time and compute resources. Output is in the `assemble`
subdirectory with `assemble/local_assemblies.bam` containing all assemblies.

### Call

From the contigs, variants are called from the alignments. The output from this step is in `call`, and the final
output is written to the root of the working directory (default = `variantsn.vcf.gz`)



## Build dependencies

In the install directory, go into the `dep` directory and run `make`. This will build a set of programs SMRT-SV relies
on. This must be completed before attempting to run SMRT-SV.

A `RepeatMasker` step is run by the pipeline, but it does not install this dependency. Your environment should have a
working version of RepeatMasker accessible through `PATH`.


## SMRT-SV structure


### Steps and rules

SMRT-SV is made of several steps (`run`, `genotype`, `align`, see above), and each step is implemented by one or more rules. Each rule performs a
simple task, such as aligning a set of reads or generating a data table for variant calling. Snakemake coordinates these
rules to accomplish the step.

### Command-line structure

SMRT-SV has two sets of command-line options, global and per-step. All global options must appear before the step name
and all per-step options must appear after the step name.

Global options are not specific to the pipeline step. For example, `--tempdir` or `--distribute`. To see help for the
global options, run `smrtsv -h`.

Per-step options are specific to the step being executed. For example, `align` and `run` have the option `--batches`,
which tells SMRT-SV how many batches to split input reads into before aligning. Since the `--batches` option does not
make sense for downstream analyses or the genotyper, those steps don't have it. To see the options for a specific step,
run `smrtsv <step> -h` where `<step>` is the name of the step (e.g. `run` or `genotype`).

For example, to specify a temp directory (global option) and the number of batches (per-step option), the command would
look something like `smrtsv --temp-dir /path/to/temp run --batches 15 ...`.

### Run script

In the SMRT-SV directory, there is a `smrtsv` and a `smrtsv.py`. Be sure to run `smrtsv`. It will setup the correct
runtime environment for SMRT-SV (using the `dep` directory) and run it.


### Environment

SMRTSV only requires that a version of `RepeatMasker` be runnable and in the `PATH`. It will ignore almost everthing
else in the environment and use the packages in `dep`. Some distributed systems clear `LD_LIBRARY_PATH` when a job
is run on a compute node, but SMRT-SV will set it back to `LD_LIBRARY_PATH` from calling environment.


## Running SMRT-SV

### Run

Change to a clean working directory that has nothing in it. The SMRT-SV pipeline should be installed in another
directory. Setup your environment so all dependencies are available. All output will go into this working directory.
The command below assumes a variable `SMRTSV_DIR` is defined to point to the root of the SMRT-SV install directory
containing `smrtsv.py`.

To run all steps:

`${SMRTSV_DIR}/smrtsv run --batches 20 --threads 8 path/to/ref.fasta path/to/reads.fofn`

The `run` command will execute the whole pipeline, but each step (see above) can be run on its own by replacing `run`
with the name of the step. Commands later in the pipeline expect to find output from previous commands (e.g. `detect`
depends on output from `align`), so these steps still need to be run in order. Note that each command has its own set of
options. See `smrtsv -h` for help, or use `-h` on a command to see its options.

### Run distributed

SMRT-SV can be run on a local machine or distributed over a cluster. Distributing is managed by Snakemake. When
it is distributed, see `cluster.eichler.json` for expected resource usage for each rule. The cluster configuration
file will probably need to be modified for your environment.


### Temporary directories

Several steps in SMRT-SV will use temporary directories to store files. SMRT-SV uses Python's `tempfile.gettempdir()`
function to choose a temporary location (See https://docs.python.org/3.7/library/tempfile.html#tempfile.gettempdir).
By default, it uses environment variables `TMPDIR`, `TEMP`, and `TMP` to choose a temporary location. If they are not
defined, then in falls back to `/tmp`. The temporary directory can be explicitly set by using the SMRT-SV command-line
option `--tempdir` option (e.g. `--tempdir /path/to/temp`), which will ignore these defaults and use the specified
location. Note that `--tempdir` must appear on the command-line before the SMRT-SV command name because it is a
global option. If the temporary directory does not exist, SMRT-SV will create it the first time it is accessed.

Each rule within SMRT-SV creates a unique path within this temporary directory so that one rule will not overwrite
temporary files from another rule. Files in the temporary directory are never shared among rules, so once a rule
completes or SMRT-SV stops, the temporary directory can be cleared without losing progress. Each rule should clear
it's temporary directory unless the rule is stopped before it reaches the end.

When SMRT-SV is distributed over a cluster, the temporary directory should be set to a location on the compute node
running the job. This saves distributed filesystem from unnecessary I/O. If there is a fast disk available, such as
a solid-state drive (SSD), then the pipeline can be sped up by placing a temp directory on that storage.


## Filtering SVs

SMRT-SV outputs variants with varying levels of confidence, and this includes false calls with weak support.

The best way to filter SMRT-SV calls is by using the `CONTIG_SUPPORT` value in the `INFO` field of the output VCF. By
selecting SVs with a contig support value of at least 2, most false calls can be eliminated.

The `QUAL` field is based on a phred-scaled ratio of contig support and total contig depth at the SV locus
(`CONTIG_SUPPORT / CONTIG_DEPTH`), which is then further scaled by the `CONTIG_DEPTH` to give more weight to SVs where
there are many contigs. A maximum value of `100` is imposed. This score has limited value as a quality measurement.

### Mixed haplotype assemblies and HETs

SMRT-SV aligns variants to the reference and pulls them down in windows to do assemblies. These come from taking
sliding-windows across the genome (default: 60 kbp widows each offset by 20 kbp) and from adding additional windows
around signatures of SVs seen in the sequence read alignments. This attempts to build multiple assemblies over each
SV for quality control.

When pulling reads down and performing an assembly, both haplotypes are mixed (assuming it is a diploid sample). For
heterozygous SVs, the assembly could yeild either the SV-containing haplotype or the reference haplotype. By adding
additional windows around SV signature, it increases the chances that the SV haplotype is represented in one or more
assemblies.

As a result of heterozygosity and mixed-haplotype assemblies, many heterozygous calls go undetected. It takes a phasing
approach or an approach based only on read alignments to recover these.

For a full treatment of this phenomenon using SMRT-SV, see http://genome.cshlp.org/lookup/doi/10.1101/gr.214007.116

For a newer method that employs phasing, see https://www.biorxiv.org/content/early/2018/06/13/193144
