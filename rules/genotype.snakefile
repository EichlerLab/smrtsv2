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

### Default Values ###

GT_DEFAULT_MODEL_NAME = '30x-4'
GT_DEFAULT_MIN_CALL_DEPTH = 8


### Read Genotyper Config File ###

# Get config file
CONFIG_FILE = config.get('genotyper_config', None)

if CONFIG_FILE is None:
    raise RuntimeError('No genotyper config file specified (e.g. "snakemake ... --config genotyper_config=path/to/config.json")')

if not os.path.isfile(CONFIG_FILE):
    raise RuntimeError('Missing genotyper config file: {}'.format(CONFIG_FILE))

# Load config
with open(CONFIG_FILE, 'r') as CONFIG_FILE_FH:
    CONFIG_GT = json.load(CONFIG_FILE_FH)


### Find genotyping model ###

# Get model name
GT_MODEL_NAME = CONFIG_GT.get('model', None)

if GT_MODEL_NAME is None:
    GT_MODEL_NAME = GT_DEFAULT_MODEL_NAME
    GT_MODEL_DEFAULT = True
else:
    GT_MODEL_DEFAULT = False

# Find path to model files
GT_MODEL_PATH = os.path.join(GT_MODEL_NAME, 'model')

if not os.path.isdir(GT_MODEL_PATH):
    GT_MODEL_PATH = os.path.join(SMRTSV_DIR, 'files', 'gtmodel', GT_MODEL_NAME, 'model')

    if not os.path.isdir(GT_MODEL_PATH):
        if GT_MODEL_DEFAULT:
            raise RuntimeError('Cannot locate default genotyper model "{}": {} (no such directory)'.format(GT_MODEL_NAME, GT_MODEL_PATH))
        else:
            raise RuntimeError('Cannot locate genotyper model "{}": {} (no such directory) and {} (no such directory)'.format(GT_MODEL_NAME, GT_MODEL_PATH, os.path.join(GT_MODEL_NAME, 'model')))

GT_MODEL_PATH = os.path.abspath(GT_MODEL_PATH)

# Verify model files
GT_SCALER = os.path.join(GT_MODEL_PATH, 'scaler.pkl')
GT_PREDICTOR = os.path.join(GT_MODEL_PATH, 'predictor.pkl')

if not os.path.isfile(GT_SCALER):
    raise RuntimeError('Missing scaler "scaler.pkl" for model {}: {}'.format(GT_MODEL_NAME, GT_SCALER))

if not os.path.isfile(GT_PREDICTOR):
    raise RuntimeError('Missing predictor "predictor.pkl" for model {}: {}'.format(GT_MODEL_NAME, GT_PREDICTOR))


### Get paths to variants, contigs, and reference ###

# Contigs
SV_CONTIGS = CONFIG_GT.get('sv_contigs', None)

if SV_CONTIGS is None:
    raise ValueError('Genotyper config file is missing element "sv_contigs"')

if not os.path.isfile(SV_CONTIGS):
    raise ValueError('Genotyper config contains a path to a missing file for element "sv_contigs": {}'.format(SV_CONTIGS))

# Calls
SV_CALLS = CONFIG_GT.get('sv_calls', None)

if SV_CALLS is None:
    raise ValueError('Genotyper config file is missing element "sv_calls"')

if not os.path.isfile(SV_CALLS):
    raise ValueError('Genotyper config contains a path to a missing file for element "sv_calls": {}'.format(SV_CALLS))

# Reference
SV_REFERENCE = CONFIG_GT.get('sv_reference', None)

if SV_REFERENCE is None:
    raise ValueError('Genotyper config file is missing element "sv_reference"')

if not os.path.isfile(SV_REFERENCE):
    raise ValueError('Genotyper config contains a path to a missing file for element "sv_reference": {}'.format(SV_REFERENCE))

SV_REFERENCE_FAI = '{}.fai'.format(SV_REFERENCE)

if not os.path.isfile(SV_REFERENCE_FAI):
    raise ValueError('Genotyper config element "sv_reference" points to a reference with no index (.fai file): {}'.format(SV_REFERENCE))


### Find bwa-postalt.js ###

POSTALT_PATH = None

for path in PROCESS_ENV['PATH'].split(':'):
    check_path = os.path.join(path, 'bwa-postalt.js')

    if os.path.isfile(check_path):
        POSTALT_PATH = check_path
        break

if POSTALT_PATH is None:
    raise RuntimeError('Cannot find "bwa-postalt.js" in PATH (part of bwakit)')


### Get sample manifest ###

SAMPLE_MANIFEST_FILE = CONFIG_GT.get('sample_manifest', None)

if SAMPLE_MANIFEST_FILE is None:
    raise RuntimeError('Configuration file "{}" is missing entry "sample_manifest"'.format(CONFIG_FILE))

if not os.path.isfile(SAMPLE_MANIFEST_FILE):
    raise RuntimeError('Cannot find sample manifest file: {}'.format(SAMPLE_MANIFEST_FILE))

SAMPLE_TABLE = genotype.preprocess_manifest(SAMPLE_MANIFEST_FILE)

SAMPLES = sorted(list(SAMPLE_TABLE.index))


### Other Parameters ###

FINAL_GENOTYPES = config.get('genotyped_variants', 'variants_gt.vcf.gz')

if os.path.exists(FINAL_GENOTYPES):
    raise RuntimeError('Aborting: Output file {} exists'.format(FINAL_GENOTYPES))

MIN_CALL_DEPTH = CONFIG_GT.get('min_call_depth', GT_DEFAULT_MIN_CALL_DEPTH)

KEEP_TEMP = smrtsvutil.as_bool(config.get('gt_keep_temp', False))



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
        genotype=expand('samples/{sample}/temp/vcf_column.tab', sample=SAMPLES)
    output:
        vcf=temp('gt/temp/variants_genotyped.vcf')
    run:

        # Get VCF header lines and body dataframe (body contains all columns fields before samples)
        header_line_list = genotype.vcf_header_lines(input.vcf)
        vcf_body = genotype.vcf_table(input.vcf)

        # Merge samples into VCF
        vcf_df = pd.concat(
            [vcf_body] + [pd.read_csv('samples/{}/temp/vcf_column.tab'.format(sample), sep='\t') for sample in SAMPLES],
            axis=1
        )

        # Write headers and VCF
        with open(output.vcf, 'w') as out_file:
            out_file.write(''.join(header_line_list))
            vcf_df.to_csv(out_file, sep='\t', index=False)

# gt_vcf_get_sample_column
#
# Get VCF column
rule gt_vcf_get_sample_column:
    input:
        tab='samples/{sample}/genotype.tab'
    output:
        tab=temp('samples/{sample}/temp/vcf_column.tab')
    run:

        # Get VCF column and write
        genotype.get_sample_column(
            input.tab, wildcards.sample, SAMPLE_TABLE.loc[wildcards.sample, 'SEX']
        ).to_csv(output.tab, sep='\t', index=False, header=True)


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
        features = pd.read_csv(input.tab, sep='\t', header=0)

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
        tab=temp('samples/{sample}/gt_features.tab')
    run:

        # Read input tables (SVs, breakpoint  depths, and insert size deltas)
        df_sv = pd.read_csv(input.sv_bed, sep='\t', header=0)
        df_sv.index = df_sv['INDEX']
        del(df_sv['CONTIG'])
        del(df_sv['CONTIG_START'])
        del(df_sv['CONTIG_END'])

        df_bp = pd.read_csv(input.bp_tab, sep='\t', header=0)
        df_bp.index = df_bp['INDEX']
        del(df_bp['INDEX'])

        df_ins = pd.read_csv(input.insert_tab, sep='\t', header=0)
        df_ins.index = df_ins['INDEX']
        del(df_ins['INDEX'])

        df_depth = pd.read_csv(input.depth_tab, sep='\t', header=0)
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
        bam='samples/{sample}/alignments.cram',
        alt_info='altref/alt_info.bed'
    output:
        tab=temp('samples/{sample}/temp/depth_delta.tab'),
        stats='samples/{sample}/depth_delta_stats.tab'
    params:
        mapq=get_config_param('gt_mapq'),
        flank=100
    shell:
        """python3 -s {SMRTSV_DIR}/scripts/genotype/GetReadDepthDiff.py """
            """{input.bam} {input.bed} {input.alt_info} {output.tab} """
            """--out_stats {output.stats} """
            """--mapq {params.mapq} """
            """--flank {params.flank}"""

# gt_call_sample_insert_delta
#
# For each variant, determine if the expected insert size of alignments over the altered reference changed.
rule gt_call_sample_insert_delta:
    input:
        sv_ref='altref/ref.fasta',
        bed='sv_calls/sv_calls.bed',
        bam='samples/{sample}/alignments.cram'
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
        """python3 -s {SMRTSV_DIR}/scripts/genotype/GetInsertSizeDelta.py """
            """-f """
            """--mapq {params.mapq} """
            """--sample_size {params.sample_size} """
            """--size_limit {params.size_limit} """
            """--ref_flank {params.ref_flank} """
            """--z_threshold {params.z_threshold} """
            """--out_stats {output.stats} """
            """--ref {input.sv_ref} """
            """{input.bam} {input.bed} {output.tab}"""

# gt_call_sample_breakpoint_depth
#
# Get breakpoints
rule gt_call_sample_breakpoint_depth:
    input:
        bed='sv_calls/sv_calls.bed',
        bam='samples/{sample}/alignments.cram'
    output:
        tab=temp('samples/{sample}/temp/breakpoint_depth.tab')
    params:
        mapq=get_config_param('gt_mapq'),
    shell:
        """python3 -s {SMRTSV_DIR}/scripts/genotype/GetBreakpointReadDepth.py """
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
        aln=lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'DATA'],
        sv_ref='altref/ref.fasta',
        sv_ref_bwt='altref/ref.fasta.bwt',
        sv_ref_alts='altref/ref.fasta.alt',
        sv_ref_alt_info='altref/alt_info.bed',
        map_regions_bed='altref/map_regions.bed'
    output:
        aln='samples/{sample}/alignments.cram',
        aln_index='samples/{sample}/alignments.cram.crai'
    params:
        mapq=get_config_param('gt_mapq'),
        threads=get_config_param('gt_map_cpu'),    # Parses into cluster params
        mem=get_config_param('gt_map_mem'),        # Parses into cluster params
        rt=get_config_param('gt_map_time'),        # Parses into cluster params
        disk_free=get_config_param('gt_map_disk')  # Parses into cluster params
    log:
        align='samples/{sample}/alignments.log'
    run:

        # Get alignment reference
        aln_ref = SAMPLE_TABLE.loc[wildcards.sample, 'REF']

        if aln_ref != 'NA':
            aln_ref = os.path.abspath(aln_ref)

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
                'gt_map_align_sample_reads',
                '-f',
                '--jobs', str(params.threads),
                '--config',
                'sample={}'.format(wildcards.sample),
                'sample_aln={}'.format(os.path.abspath(input.aln)),
                'sample_aln_ref={}'.format(aln_ref),
                'sv_ref={}'.format(os.path.abspath(input.sv_ref)),
                'sv_ref_alt={}'.format(os.path.abspath(input.sv_ref_alts)),
                'sv_ref_alt_info={}'.format(os.path.abspath(input.sv_ref_alt_info)),
                'map_regions_bed={}'.format(os.path.abspath(input.map_regions_bed)),
                'output_aln={}'.format(os.path.abspath(output.aln)),
                'output_aln_index={}'.format(os.path.abspath(output.aln_index)),
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

# gt_altref_map_regions_merge_bed
#
# Merge regions around SVs. Reads mapping outside these regions will be dropped.
#
# Note: params.ref_flank must be larger than params.ref_flank in rule gt_call_sample_insert_delta
rule gt_altref_map_regions_merge_bed:
    input:
        bed='temp/altref/map_regions.bed',
        fai='altref/ref.fasta.fai'
    output:
        bed='altref/map_regions.bed'
    params:
        ref_flank=int(5.2e3)
    shell:
        """bedtools slop -i {input.bed} -g {input.fai} -b {params.ref_flank} | """
        """sort -k1,1 -k2,2n | """
        """bedtools merge """
        """> {output.bed}"""

# gt_altref_map_regions_bed
#
# Make a BED file of primary and SV-contig regions.
rule gt_altref_map_regions_bed:
    input:
        bed='sv_calls/sv_calls.bed'
    output:
        bed=temp('temp/altref/map_regions.bed')
    run:

        # Read breakpoints on primary contigs
        df_primary = pd.read_csv(input.bed, sep='\t', usecols=('#CHROM', 'POS', 'END'))

        # Read breakpoints on SV contigs
        df_sv = pd.read_csv(input.bed, sep='\t', usecols=('CONTIG', 'CONTIG_START', 'CONTIG_END'))
        df_sv.columns = ('#CHROM', 'POS', 'END')

        # Merge and write
        pd.concat([df_primary, df_sv], axis=0).to_csv(output.bed, sep='\t', index=False)


# gt_altref_alt_info_bed
#
# Make a BED file for each chromosme and ALT contig (records cover the whole sequence) with the last column identifying
# which primary assembly contig it belongs to. For primary assembly records in this BED file, the last column is the
# contig name itself.
rule gt_altref_alt_info_bed:
    input:
        sv_ref_fai=SV_REFERENCE_FAI,
        altref_fai='altref/ref.fasta.fai',
        contig_chr='temp/altref/alt_to_primary.tab'
    output:
        bed='altref/alt_info.bed'
    run:

        # Read reference FAI (Series of lengths keyed by chromosome name)
        seq_len = pd.read_csv(
            input.altref_fai,
            sep='\t',
            header=None,
            names=('CHROM', 'LEN', 'OFFSET', 'LINE_LEN', 'LINE_BYTES'),
            usecols=('CHROM', 'LEN'),
            index_col='CHROM',
            squeeze=True
        )

        # Get a list of primary contigs (from reference sequence only, no alt contigs)
        primary_seqs = set(pd.read_csv(input.sv_ref_fai, sep='\t', header=None, usecols=(0, ))[0].tolist())

        # Get the alternates file (a SAM file of ALT contigs aligned to a primary contig)
        alt_to_primary = pd.read_csv(input.contig_chr, sep='\t', header=None, index_col=0, squeeze=True)

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
        """python3 -s {SMRTSV_DIR}/scripts/genotype/SamToAlt.py {input.sam} {output.alt}"""

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
        ref=SV_REFERENCE,
        contig='contigs/contigs.fasta'
    output:
        fasta='altref/ref.fasta',
        fai='altref/ref.fasta.fai'
    shell:
        """cat {input.ref} {input.contig} >{output.fasta}; """
        """samtools faidx {output.fasta}; """


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
        contigs=SV_CONTIGS,
        list='contigs/contig_list.txt'
    output:
        sam='contigs/contigs.sam',
        sizes='contigs/contigs.sizes'
    shell:
        """python3 -s {SMRTSV_DIR}/scripts/genotype/FilterContigs.py {input.contigs} {input.list} {output.sam} {output.sizes} {SV_REFERENCE}"""

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
            with open(output.txt, 'w') as out_file:
                for contig in sorted(set([record.info['CONTIG'] for record in vcf])):
                    out_file.write('%s\n' % contig)


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
        bed='sv_calls/sv_calls.bed'
    run:
        record_count = 0

        with pysam.VariantFile(input.vcf) as vcf:
            with open(output.bed, 'w') as out_file:

                # Write header
                out_file.write('#CHROM\tPOS\tEND\tID\tINDEX\tSVTYPE\tSVLEN\tCONTIG\tCONTIG_START\tCONTIG_END\n')

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
                    out_file.write("%s\n" % "\t".join(map(str, (
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
        vcf=SV_CALLS
    output:
        vcf='sv_calls/sv_calls.vcf.gz',
        tbi='sv_calls/sv_calls.vcf.gz.tbi'
    shell:
        """bcftools filter -O z -i "(SVTYPE='INS' | SVTYPE='DEL') & (SVLEN >= 50 | SVLEN <= -50)" {input.vcf} """
        """>{output.vcf}; """
        """tabix -p vcf {output.vcf}"""
