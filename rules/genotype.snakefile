import json
import numpy as np
import os
import pandas as pd
import pysam
import shutil
import tempfile
import subprocess

from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'include.snakefile'

from smrtsvlib import genotype
from smrtsvlib import ml
from smrtsvlib import smrtsvrunner
from smrtsvlib import smrtsvutil

localrules: gt_vcf_write


###################
### Definitions ###
###################

### Read Genotyper Config File ###

CONFIG_FILE = config['genotyper_config']

with open(CONFIG_FILE, 'r') as CONFIG_FILE_FH:
    CONFIG_GT = json.load(CONFIG_FILE_FH)


### Find genotyping model ###

# Scaler: Scales feature matrix for model prediction.
GT_SCALER = CONFIG_GT.get('model', {}).get('scaler', None)

if GT_SCALER is None:
    GT_SCALER = os.path.join(SMRTSV_DIR, 'files/gtmodel/scaler.pkl')

# Predictor: Predicts genotype from scaled features.
GT_PREDICTOR = CONFIG_GT.get('model', {}).get('predictor', None)

if GT_PREDICTOR is None:
    GT_PREDICTOR = os.path.join(SMRTSV_DIR, 'files/gtmodel/predictor.pkl')

### Find bwa-postalt.js ###

POSTALT_PATH = None

for path in PROCESS_ENV['PATH'].split(':'):
    check_path = os.path.join(path, 'bwa-postalt.js')

    if os.path.isfile(check_path):
        POSTALT_PATH = check_path
        break

if POSTALT_PATH is None:
    raise RuntimeError('Cannot find "bwa-postalt.js" in PATH (part of bwakit)')


### Other Parameters ###

FINAL_GENOTYPES = config['genotyped_variants']
SAMPLES = sorted(CONFIG_GT['samples'].keys())
SVMAP_REF_WINDOW = 5000
SVMAP_CONTIG_WINDOW = 500
MIN_CALL_DEPTH = CONFIG_GT.get('min_call_depth', 4)

KEEP_TEMP = smrtsvutil.as_bool(config.get('gt_keep_temp', False))


###  Get sample manifest ###
SAMPLE_MANIFEST = CONFIG_GT.get('sample_manifest', None)

if SAMPLE_MANIFEST is None:
    if not smrtsvutil.as_bool(CONFIG_GT.get('default_sample_manifest', False)):
        raise RuntimeError('Configuration file "{}" is missing entry "sample_manifest" (only allowed if "default_sample_manifest" is set to "True")'.format(CONFIG_FILE))


### Utility Functions ###

def _get_sv_regions_for_sample(wildcards):
    """
    Get the BAM file specifying which regions reads are extracted from. The input BAM for each sample is aligned to
    a reference sequence, which may not be the same as the reference SVs were called against. It may also be different
    among the samples. This rule finds the BED file containing regions reads should be extracted from for each sample.

    :param wildcards: Rule wildcards. Must contain attribute `sample`.
    """

    if 'sample_bam_reference' in CONFIG_GT and CONFIG_GT["sample_bam_reference"].get(wildcards.sample, None) is not None:
        reference = CONFIG_GT['sample_bam_reference'].get(wildcards.sample)
    else:
        reference = CONFIG_GT['default_bam_reference']

    return 'svmap/mapping/{}/pseudoreads.bed'.format(reference)

def _get_reference_for_sample(wildcards):
    """
    Get the reference FASTA file for a sample.

    :param wildcards: Rule wildcards. Must contain attribute `sample`.
    """

    if 'sample_bam_reference' in CONFIG_GT and CONFIG_GT["sample_bam_reference"].get(wildcards.sample, None) is not None:
        reference = CONFIG_GT['sample_bam_reference'].get(wildcards.sample)
    else:
        reference = CONFIG_GT['default_bam_reference']

    return CONFIG_GT['bam_reference'][reference]


#############
### Rules ###
#############


#
# Format VCF output
#

# gt_vcf_write
#
# Write final VCF.
rule gt_vcf_write:
    input:
       vcf='gt/temp/variants_fixup.vcf'
    output:
        vcf=FINAL_GENOTYPES
    run:

        if output.vcf.lower().endswith('.vcf.gz'):
            shell("""bgzip -c {input.vcf} > {output.vcf}; tabix {output.vcf}""")

        elif output.vcf.lower().endswith('.vcf'):
            shell("""cp {input.vcf} {output.vcf}""")

        else:
            raise ValueError('Unknown output file format: Expected output file with ".vcf" or ".vcf.gz" extension')

# gt_vcf_fixup
#
# Annotate INFO fields with genotype information.
rule gt_vcf_fixup:
    input:
        vcf='gt/temp/variants_genotyped.vcf'
    output:
        vcf=temp('gt/temp/variants_fixup.vcf')
    shell:
        """vcffixup {input.vcf} > {output.vcf}"""

# gt_vcf_merge
#
# Merge genotype calls with the original VCF file.
rule gt_vcf_merge:
    input:
        vcf='sv_calls/sv_calls.vcf.gz',
        genotype=expand('samples/{sample}/genotype.tab', sample=SAMPLES)
    output:
        vcf=temp('gt/temp/variants_genotyped.vcf')
    run:

        # Get VCF header lines and body dataframe (body contains all columns fields before samples)
        header_line_list = genotype.vcf_header_lines(input.vcf)
        vcf_body = genotype.vcf_table(input.vcf)

        # Get sexes from the sample manifest
        sexes = genotype.preprocess_manifest(SAMPLE_MANIFEST, SAMPLES)

        # Merge samples into VCF
        vcf_df = pd.concat(
            [vcf_body] + [genotype.get_sample_column('samples/{}/genotype.tab'.format(sample), sample, sexes[sample]) for sample in SAMPLES],
            axis=1
        )

        # Write headers and VCF
        with open(output.vcf, 'w') as out_file:
            out_file.write(''.join(header_line_list))
            vcf_df.to_csv(out_file, sep='\t', index=False)

#
# Call genotypes (apply learned model to features)
#

# gt_call_predict
#
# Predict genotype
rule gt_call_predict:
    input:
        tab='samples/{sample}/gt_features.tab',
        predictor=GT_PREDICTOR,
        scaler=GT_SCALER
    output:
        tab='samples/{sample}/genotype.tab'
    params:
        min_depth=MIN_CALL_DEPTH
    run:

        # Read features
        features = pd.read_table(input.tab, header=0)

        # Annotate no-calls
        features['CALLABLE'] = features.apply(
            lambda row: True if row['BP_REF_COUNT'] + row['BP_ALT_COUNT'] >= params.min_depth else False,
            axis=1
        )

        # Predict genotype and estimate density
        model = ml.GtModel(input.predictor, input.scaler)

        genotype, class_density = model.genotype_and_density(features)

        features['CLASS'] = genotype

        features = pd.concat([features, class_density], axis=1)

        # Write
        features.to_csv(output.tab, sep='\t', index=False)

# gt_call_sample_merge
#
# Merge genotyping information tables. The resulting table will be used to make genotype calls.
rule gt_call_sample_merge:
    input:
        sv_bed='sv_calls/sv_calls.bed',
        bp_tab='samples/{sample}/temp/breakpoint_depth.tab',
        insert_tab='samples/{sample}/temp/insert_delta.tab',
        depth_tab='samples/{sample}/temp/depth_delta.tab'
    output:
        tab='samples/{sample}/gt_features.tab'
    run:

        # Read input tables (SVs, breakpoint  depths, and insert size deltas)
        df_sv = pd.read_table(input.sv_bed, header=0)
        df_sv.index = df_sv['INDEX']
        del(df_sv['CONTIG'])
        del(df_sv['CONTIG_START'])
        del(df_sv['CONTIG_END'])

        df_bp = pd.read_table(input.bp_tab, header=0)
        df_bp.index = df_bp['INDEX']
        del(df_bp['INDEX'])

        df_ins = pd.read_table(input.insert_tab, header=0)
        df_ins.index = df_ins['INDEX']
        del(df_ins['INDEX'])

        df_depth = pd.read_table(input.depth_tab, header=0)
        df_depth.index = df_depth['INDEX']
        del(df_depth['INDEX'])

        # Merge
        df = pd.concat(
            [df_sv, df_bp, df_ins, df_depth],
            axis=1
        )

        # Write
        df.to_csv(output.tab, sep='\t', index=False)

# gt_call_sample_read_depth
#
# Get read depths in and around altered bases.
rule gt_call_sample_read_depth:
    input:
        bed='sv_calls/sv_calls.bed',
        bam='samples/{sample}/alignments.bam',
        alt_info='altref/alt_info.bed'
    output:
        tab=temp('samples/{sample}/temp/depth_delta.tab'),
        stats='samples/{sample}/depth_delta_stats.tab'
    params:
        mapq=get_config_param('gt_mapq'),
        flank=100
    shell:
        """python {SMRTSV_DIR}/scripts/genotype/GetReadDepthDiff.py """
            """{input.bam} {input.bed} {input.alt_info} {output.tab} """
            """--out_stats {output.stats} """
            """--mapq {params.mapq} """
            """--flank {params.flank}"""

# gt_call_sample_insert_delta
#
# For each variant, determine if the expected insert size of alignments over the altered reference changed.
rule gt_call_sample_insert_delta:
    input:
        bed='sv_calls/sv_calls.bed',
        bam='samples/{sample}/alignments.bam'
    output:
        tab=temp('samples/{sample}/temp/insert_delta.tab'),
        stats='samples/{sample}/insert_size_stats.tab'
    params:
        mapq=get_config_param('gt_mapq'),
        ref_flank=int(5e3),
        sample_size=int(5e6),
        size_limit=1500,
        z_threshold=1.5
    shell:
        """python {SMRTSV_DIR}/scripts/genotype/GetInsertSizeDelta.py """
            """-f """
            """--mapq {params.mapq} """
            """--sample_size {params.sample_size} """
            """--size_limit {params.size_limit} """
            """--ref_flank {params.ref_flank} """
            """--z_threshold {params.z_threshold} """
            """--out_stats {output.stats} """
            """{input.bam} {input.bed} {output.tab}"""

# gt_call_sample_breakpoint_depth
#
# Get breakpoints
rule gt_call_sample_breakpoint_depth:
    input:
        bed='sv_calls/sv_calls.bed',
        bam='samples/{sample}/alignments.bam'
    output:
        tab=temp('samples/{sample}/temp/breakpoint_depth.tab')
    params:
        mapq=get_config_param('gt_mapq'),
    shell:
        """python {SMRTSV_DIR}/scripts/genotype/GetBreakpointReadDepth.py """
            """-f """
            """--mapq {params.mapq} """
            """{input.bam} {input.bed} {output.tab}"""

#
# Map reads.
#

# gt_map_sample_reads
#
# Extract reads from the original BAM guided. Mapped read locations discovered by mapping pseudoreads to the sample's
# reference are extracted along with unmapped reads.
rule gt_map_sample_reads:
    input:
        sample_bam=lambda wildcards: CONFIG_GT['samples'][wildcards.sample],
        sample_ref=_get_reference_for_sample,
        sample_regions=_get_sv_regions_for_sample,
        sv_ref='altref/ref.fasta',
        sv_ref_bwt='altref/ref.fasta.bwt',
        sv_ref_alts='altref/ref.fasta.alt',
        sv_ref_alt_info='altref/alt_info.bed'
    output:
        bam='samples/{sample}/alignments.bam',
        bai='samples/{sample}/alignments.bam.bai'
    params:
        mapq=get_config_param('gt_mapq'),
        threads=get_config_param('gt_map_cpu'),  # Parses into cluster params
        mem=get_config_param('gt_map_mem')       # Parses into cluster params
    benchmark:
        'samples/{sample}/bm/alignments.txt'
    log:
        align='samples/{sample}/alignments.log',
        map='samples/{sample}/primary_mapping.log'
    run:

        # Set mapping_temp (will be deleted if not None)
        mapping_temp = None

        try:

            # Create temporary directories
            mapping_temp = tempfile.mkdtemp(
                prefix=os.path.join(
                    TEMP_DIR,
                    'map_sample_{}_'.format(wildcards.sample)
                )
            )

            # Setup sub-Snake command
            command = (
                'gt_map_postalt_merge',
                '-f',
                '--jobs', str(params.threads),
                '--config',
                'sample={}'.format(wildcards.sample),
                'sample_bam={}'.format(os.path.abspath(input.sample_bam)),
                'sample_ref={}'.format(os.path.abspath(input.sample_ref)),
                'sample_regions={}'.format(os.path.abspath(input.sample_regions)),
                'sv_ref={}'.format(os.path.abspath(input.sv_ref)),
                'sv_ref_alt={}'.format(os.path.abspath(input.sv_ref_alts)),
                'sv_ref_alt_info={}'.format(os.path.abspath(input.sv_ref_alt_info)),
                'output_bam={}'.format(os.path.abspath(output.bam)),
                'primary_map_log={}'.format(os.path.abspath(log.map)),
                'mapq={}'.format(params.mapq),
                'threads={}'.format(params.threads),
                'smrtsv_dir={}'.format(SMRTSV_DIR),
                'postalt_path={}'.format(POSTALT_PATH)
            )

            # Run mapping step
            with open(log.align, 'w') as log_file:
                return_code = smrtsvrunner.run_snake_target(
                    'rules/genotype_map.snakefile', None, PROCESS_ENV, SMRTSV_DIR, command,
                    stdout=log_file, stderr=subprocess.STDOUT, cwd=mapping_temp
                )

            # Check return code
            if return_code != 0:
                raise RuntimeError(
                    'Alignment step failed with return code {}: See {} for errors.'.format(return_code, log.align)
                )

        finally:

            # Clean up temp directory
            if mapping_temp is not None and not KEEP_TEMP:
                shutil.rmtree(mapping_temp)


#
# Prepare local assembly sequences and references.
#

# gt_altref_alt_info_bed
#
# Make a BED file for each chromosme and ALT contig (records cover the whole sequence) with the last column identifying
# which primary assembly contig it belongs to. For primary assembly records in this BED file, the last column is the
# contig name itself.
rule gt_altref_alt_info_bed:
    input:
        sv_ref_fai=CONFIG_GT['sv_reference_fai'],
        altref_fai='altref/ref.fasta.fai',
        contig_chr='temp/altref/alt_to_primary.tab'
    output:
        bed='altref/alt_info.bed'
    run:

        # Read reference FAI (Series of lengths keyed by chromosome name)
        seq_len = pd.read_table(
            input.altref_fai,
            header=None,
            names=('CHROM', 'LEN', 'OFFSET', 'LINE_LEN', 'LINE_BYTES'),
            usecols=('CHROM', 'LEN'),
            index_col='CHROM',
            squeeze=True
        )

        # Get a list of primary contigs (from reference sequence only, no alt contigs)
        primary_seqs = set(pd.read_table(input.sv_ref_fai, header=None, usecols=(0, ))[0].tolist())

        # Get the alternates file (a SAM file of ALT contigs aligned to a primary contig)
        alt_to_primary = pd.read_table(input.contig_chr, header=None, index_col=0, squeeze=True)

        # Make BED of all contigs as a dataframe
        df_bed = pd.DataFrame(
            np.concatenate([[chr, 0, seq_len[chr]] for chr in seq_len.index]).reshape(-1, 3),
            columns=('#CHROM', 'POS', 'END')
        )

        # Add a column indicating whether this is a primary contig (True) or ALT (False)
        df_bed['IS_PRIMARY'] = df_bed.apply(
            lambda row: row['#CHROM'] in primary_seqs,
            axis=1
        )

        # Add a column of the primary contig name. If the contig is primary, then add the name of itself for that
        # record. If the contig is alt, then add the name of the primary contig it is an alternate of.
        df_bed['PRIMARY'] = df_bed.apply(
            lambda row: row['#CHROM'] if row['IS_PRIMARY'] else alt_to_primary[row['#CHROM']],
            axis=1
        )

        # Write
        df_bed.to_csv(output.bed, sep='\t', header=True, index=False)

# gt_altref_alt_info_lengths
#
# Get a table of contigs (col 1) and the chromosome they map to (col 2).
rule gt_altref_alt_contig_to_chr:
    input:
        altref_alt='altref/ref.fasta.alt'
    output:
        tab=temp('temp/altref/alt_to_primary.tab')
    shell:
        """awk -vOFS="\t" '($1 !~ /^@/) {{print $1, $3}}' {input.altref_alt} >{output.tab}"""

# gt_altref_make_alts
#
# The alternates file is a SAM file of the alternate contigs and how they map to the primary assembly. This
# file tells BWA which contigs are alternates and bwa-postalts.js uses it to fill the alternate alignments and
# adjust alignment scores.
rule gt_altref_make_alts:
    input:
        sam='contigs/contigs.sam'
    output:
        alt='altref/ref.fasta.alt'
    shell:
        """python {SMRTSV_DIR}/scripts/genotype/SamToAlt.py {input.sam} {output.alt}"""

# gt_altref_index
#
# Index the reference with contigs as alternate sequences.
rule gt_altref_index:
    input:
        fasta='altref/ref.fasta'
    output:
        bwt='altref/ref.fasta.bwt'
    params:
        block_size=1000000000
    log:
        'altref/log/sv_altref_index.log'
    shell:
        """bwa index -b {params.block_size} {input.fasta} >{log} 2>&1"""


# gt_altref_prepare
#
# Create a reference with contigs as alternate sequences.
rule gt_altref_prepare:
    input:
        ref=CONFIG_GT['sv_reference'],
        contig='contigs/contigs.fasta'
    output:
        fasta='altref/ref.fasta',
        fai='altref/ref.fasta.fai'
    shell:
        """cat {input.ref} {input.contig} >{output.fasta}; """
        """samtools faidx {output.fasta}; """


#
# Find locations where SV associated reads might have mapped in the input BAM file. This is accomplished
# by making psedoreads from SV associated regions (breakpoints and inserted sequence) and mapping them against
# the reference sequence of the aligned BAM file.
#

# gt_svmap_locate_regions
#
# Locate regions where the pseudoreads mapped.
rule gt_svmap_locate_regions:
    input:
        bam='svmap/mapping/{reference}/pseudoreads.bam',
        ref=lambda wildcards: CONFIG_GT['bam_reference'][wildcards.reference]
    output:
        bed='svmap/mapping/{reference}/pseudoreads.bed'
    benchmark:
        'svmap/mapping/{reference}/bm/gt_svmap_locate_regions.txt'
    params:
        slop=SVMAP_REF_WINDOW
    shell:
        """bedtools bamtobed -i {input.bam} | """
        """bedtools slop -i stdin -g {input.ref}.fai -b {params.slop} | """
        """sort -k 1,1 -k 2,2n | """
        """bedtools merge -i stdin -d 0 """
        """>{output.bed}"""

# gt_svmap_map_pseudoreads
#
# Map pseudoreads to the reference.
rule gt_svmap_map_pseudoreads:
    input:
        ref=lambda wildcards: CONFIG_GT['bam_reference'][wildcards.reference],
        fasta='svmap/pseudoreads.fasta.gz'
    output:
        bam='svmap/mapping/{reference}/pseudoreads.bam'
    params:
        threads=12
    log:
        'svmap/mapping/{reference}/pseudoreads.log'
    benchmark:
        'svmap/mapping/{reference}/bm/gt_svmap_map_pseudoreads.txt'
    run:

        mapping_temp = None

        try:
            # Create temporary directory
            mapping_temp = tempfile.mkdtemp(
                prefix=os.path.join(
                    TEMP_DIR,
                    'map_pseudo_{}_'.format(wildcards.reference)
                )
            )

            shell(
                """echo "Temp directory: {mapping_temp}" >{log}; """
                """bwa mem -t {params.threads} {input.ref} {input.fasta} 2>>{log} | """
                """samtools sort -o {output.bam} -O bam -T {mapping_temp} 2>>{log}"""
            )

        finally:
            # Remove temp
            if mapping_temp is not None:
                shutil.rmtree(mapping_temp)

# gt_svmap_compress_fragmented_fasta
#
# Compress fragmented FASTA.
rule gt_svmap_compress_fragmented_fasta:
    input:
        fasta='svmap/pseudoreads.fasta'
    output:
        fasta='svmap/pseudoreads.fasta.gz'
    shell:
        """bgzip {input.fasta}"""

# gt_svmap_fragment_fasta
#
# Create pseudoreads by taking all SV associated sequences (breakpoints on the reference and contigs) and
# breaking them into 250 bp sequences (overlapping by 125 bp).
rule gt_svmap_fragment_fasta:
    input:
        fasta='svmap/regions.fasta'
    output:
        fasta=temp('svmap/pseudoreads.fasta')
    params:
        window="250",
        slide="125"
    shell:
        """python {SMRTSV_DIR}/scripts/genotype/FragmentFastaRecords.py {input.fasta} {output.fasta} {params.window} --slide {params.slide}"""

# gt_svmap_merge_fasta
#
# Merge FASTA sequences from inserted sequence and breakpoints on the reference.
rule gt_svmap_merge_fasta:
    input:
        fasta_contig='svmap/regions/contig_insdel.fasta',
        fasta_ref='svmap/regions/ref_insdel.fasta'
    output:
        fasta=temp('svmap/regions.fasta')
    shell:
        """cat {input.fasta_contig} {input.fasta_ref} """
        """>{output.fasta}"""

# gt_svmap_contig_sequences
#
# Get sequences around SVs over contigs.
rule gt_svmap_contig_sequences:
    input:
        sv_bed='svmap/regions/contig_insdel.bed',
        contigs='contigs/contigs.fasta'
    output:
        fasta=temp('svmap/regions/contig_insdel.fasta')
    shell:
        """bedtools getfasta -fi {input.contigs} -bed {input.sv_bed} -fo {output.fasta}"""

# gt_svmap_ref_sequences
#
# Get sequences around SVs.
rule gt_svmap_sequences:
    input:
        sv_bed='svmap/regions/ref_insdel.bed',
        ref=CONFIG_GT['sv_reference']
    output:
        fasta=temp('svmap/regions/ref_insdel.fasta')
    shell:
        """bedtools getfasta -fi {input.ref} -bed {input.sv_bed} -fo {output.fasta}"""

# gt_svmap_contig_positions
#
# Get SV positions over contigs.
rule gt_svmap_contig_positions:
    input:
        sv_bed='sv_calls/sv_calls.bed',
        contig_fai='contigs/contigs.fasta.fai'
    output:
        bed='svmap/regions/contig_insdel.bed'
    params:
        slop=SVMAP_CONTIG_WINDOW
    shell:
        """awk -vOFS="\\t" '(NR > 1) {{print $8, $9, $10}}' {input.sv_bed} | """
        """bedtools slop -i stdin -g {input.contig_fai} -b {params.slop} | """
        """sort -k 1,1 -k 2,2n | """
        """bedtools merge -i stdin -d 0 """
        """>{output.bed}"""

# gt_svmap_merge_ref_positions
#
# Merge positions from insertions and deletions.
rule gt_svmap_merge_ref_positions:
    input:
        bed_ins='svmap/regions/ref_ins.bed',
        bed_del='svmap/regions/ref_del.bed'
    output:
        bed='svmap/regions/ref_insdel.bed'
    shell:
        """sort -k 1,1 -k 2,2n -m {input.bed_ins} {input.bed_del} | """
        """bedtools merge -i stdin -d 0 """
        """>{output.bed}"""

# gt_svmap_ref_ins_positions
#
# Get locations of insertions.
rule gt_svmap_ref_ins_positions:
    input:
        sv_bed='sv_calls/sv_calls.bed',
        ref_fai=CONFIG_GT['sv_reference_fai']
    output:
        bed='svmap/regions/ref_ins.bed'
    params:
        slop=SVMAP_CONTIG_WINDOW
    shell:
        """awk '$6 == "INS"' {input.sv_bed} | """
        """awk 'OFS="\\t" {{ print $1,$2,$2+1 }}' | """
        """bedtools slop -i stdin -g {input.ref_fai} -b {params.slop} | """
        """sort -k 1,1 -k 2,2n | """
        """bedtools merge -i stdin -d 0 """
        """>{output.bed}"""

# gt_svmap_ref_del_positions
#
# Get locations of deletions.
rule gt_svmap_ref_del_positions:
    input:
        sv_bed='sv_calls/sv_calls.bed',
        ref_fai=CONFIG_GT['sv_reference_fai']
    output:
        bed='svmap/regions/ref_del.bed'
    params:
        slop=SVMAP_CONTIG_WINDOW
    shell:
        """awk '$6 == "DEL"' {input.sv_bed} | """
        """cut -f 1-3 | """
        """bedtools slop -i stdin -g {input.ref_fai} -b {params.slop} | """
        """bedtools merge -i stdin -d 0 """
        """>{output.bed}"""


#
# Get Contigs
#

# gt_contig_bam_to_fasta
#
# Get a FASTA file of local assemblies.
rule gt_contig_bam_to_fasta:
    input:
        sam='contigs/contigs.sam'
    output:
        fasta='contigs/contigs.fasta',
        fai='contigs/contigs.fasta.fai'
    log:
        'contigs/log/gt_contig_bam_to_fasta.log'
    run:

        seq_list = list()

        # Get SAM Records as FASTA records
        with open(input.sam, 'r') as in_file:
            for line in in_file:

                line = line.strip()

                if not line or line.startswith('@'):
                    continue

                tok = line.split('\t')

                seq_list.append(SeqRecord.SeqRecord(Seq.Seq(tok[9]), id=tok[0], description=''))

        # Write FASTA Records
        with open(output.fasta, 'w') as out_file:
            SeqIO.write(seq_list, out_file, 'fasta')

        # Index FASTA
        shell("""samtools faidx {output.fasta}""")


# gt_contig_filter
#
# Retrieve all contigs with at least one structural variant.
rule gt_contig_filter:
    input:
        contigs=CONFIG_GT['sv_contigs'],
        list='contigs/contig_list.txt'
    output:
        sam='contigs/contigs.sam',
        sizes='contigs/contigs.sizes'
    shell:
        """python {SMRTSV_DIR}/scripts/genotype/FilterContigs.py {input.contigs} {input.list} {output.sam} {output.sizes}"""

# gt_contig_list
#
# Get a list of assembled contigs with variant calls.
rule gt_contig_list:
    input:
        vcf='sv_calls/sv_calls.vcf.gz'
    output:
        txt='contigs/contig_list.txt'
    run:
        with pysam.VariantFile(input.vcf) as vcf:
            with open(output.txt, 'w') as oh:
                for contig in sorted(set([record.info['CONTIG'] for record in vcf])):
                    oh.write('%s\n' % contig)


#
# Get data from VCF
#

# gt_sv_to_bed
#
# Get BED file of variants.
rule gt_sv_to_bed:
    input:
        vcf='sv_calls/sv_calls.vcf.gz'
    output:
        vcf='sv_calls/sv_calls.bed'
    run:
        record_count = 0

        with pysam.VariantFile(input.vcf) as vcf:
            with open(output[0], "w") as oh:

                # Write header
                oh.write('#CHROM\tPOS\tEND\tID\tINDEX\tSVTYPE\tSVLEN\tCONTIG\tCONTIG_START\tCONTIG_END\n')

                # Process records
                for record in vcf:
                    record_count += 1

                    # Get length and type
                    sv_length = abs(record.info['SVLEN'])
                    sv_type= record.info['SVTYPE'].upper()

                    # Set end
                    if sv_type == 'INS':
                        sv_end = record.start + 1
                    elif sv_type == 'DEL':
                        sv_end = record.start + sv_length
                    else:
                        raise RuntimeError('Unrecognized SVTYPE {} in record {}'.format(sv_type, record_count))

                    # Check for CONTIG, CONTIG_START, and CONTIG_END
                    if not 'CONTIG' in record.info:
                        raise RuntimeError('Missing CONTIG INFO in record {}'.format(record_count))

                    if not 'CONTIG_START' in record.info:
                        raise RuntimeError('Missing CONTIG_START INFO in record {}'.format(record_count))

                    if not 'CONTIG_END' in record.info:
                        raise RuntimeError('Missing CONTIG_END INFO in record {}'.format(record_count))

                    # Check ID
                    if record.id is not None:
                        record_id = record.id.strip()

                        if not record_id:
                            record_id = '.'
                    else:
                        record_id = '.'

                    # Write
                    oh.write("%s\n" % "\t".join(map(str, (
                        record.chrom,
                        record.start,
                        sv_end,
                        record_id,
                        record_count - 1,
                        sv_type,
                        sv_length,
                        record.info['CONTIG'],
                        record.info['CONTIG_START'],
                        record.info['CONTIG_END']
                    ))))

# gt_sv_filter_vcf_for_sv
#
# Extract structural variants from variant calls.
rule gt_sv_filter_vcf_for_sv:
    input:
        vcf=CONFIG_GT['sv_calls']
    output:
        vcf='sv_calls/sv_calls.vcf.gz',
        tbi='sv_calls/sv_calls.vcf.gz.tbi'
    shell:
        """bcftools filter -O z -i "(SVTYPE='INS' | SVTYPE='DEL') & (SVLEN >= 50 | SVLEN <= -50)" {input.vcf} """
        """>{output.vcf}; """
        """tabix -p vcf {output.vcf}"""
