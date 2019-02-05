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


## Running SMRT-SV

### Build dependencies

In the install directory, go into the `dep` directory and run `make`. This will build a set of programs SMRT-SV relies
on.

A `RepeatMasker` step is run by the pipeline, but it does not install this dependency. Your environment should have a
working version of RepeatMasker accessible through `PATH`.


### Run

Change to a clean working directory that has nothing in it. The SMRT-SV pipeline should be installed in another
directory. Setup your environment so all dependencies are available. All output will go into this working directory.
The command below assumes a variable `SMRTSV_DIR` is defined to point to the root of the SMRT-SV install directory
containing `smrtsv.py`.

To run all steps:

`${SMRTSV_DIR}/smrtsv.py run --batches 20 --threads 8 path/to/ref.fasta path/to/reads.fofn`

The `run` command will execute the whole pipeline, but each step (see above) can be run on its own by replacing `run`
with the name of the step. Commands later in the pipeline expect to find output from previous commands (e.g. `detect`
depends on output from `align`), so these steps still need to be run in order. Note that each command has its own set of
options. See `smrtsv.py -h` for help, or use `-h` on a command to see its options.

### Run distributed

SMRT-SV can be run on a local machine or distributed over a cluster. Distributing is managed by Snakemake. When
it is distributed, see `cluster.eichler.json` for expected resource usage for each rule. The cluster configuration
file will probably need to be modified for your environment.
