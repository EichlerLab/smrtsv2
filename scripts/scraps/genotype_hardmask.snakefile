
# gt_contig_bam
#
# Replace BAM query sequences with their masked version.
rule gt_contig_bam:
    input:
        bam='temp/contigs/contigs.bam',
        fasta='contigs/contigs.fasta',
        fai='contigs/contigs.fasta.fai'
    output:
        bam='contigs/contigs.bam',
        bai='contigs/contigs.bam.bai'
    run:

        # Read from input BAM, change query sequence to the masked version, and write to output BAM
        with pysam.AlignmentFile(input.bam, 'r') as in_bam:
            with pysam.FastaFile(input.fasta) as in_fa:
                with pysam.AlignmentFile(output.bam, 'wb', in_bam) as out_bam:
                    for record in in_bam.fetch():
                        record.query_sequence = in_fa.fetch(record.query_name)
                        out_bam.write(record)

        # Index
        shell("""samtools index {output.bam}""")

# gt_contig_mask_fasta
#
# Hard-mask contig sequences.
rule gt_contig_mask_fasta:
    input:
        fasta='temp/contigs/hardmask/contigs_nomask.fasta',
        fai='temp/contigs/hardmask/contigs_nomask.fasta.fai',
        bed='contigs/hardmask/regions.bed'
    output:
        fasta='contigs/contigs.fasta',
        fai='contigs/contigs.fasta.fai'
    shell:
        """bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output.fasta}; """
        """samtools faidx {output.fasta}"""


# gt_contig_hardmask_bed
#
# Get a BED of regions to be masked.
rule gt_contig_hardmask_bed:
    input:
        bed='temp/contigs/hardmask/regions.bed',
        sizes='contigs/contigs.sizes'
    output:
        bed='contigs/hardmask/regions.bed'
    shell:
        """awk -vOFS="\\t" '{{print $1, "0", $2}}' {input.sizes} | """
        """bedtools subtract -a stdin -b {input.bed} """
        """>{output.bed}"""

# gt_contig_hardmask_flank
#
# Add flank to variant locations.
rule gt_contig_hardmask_flank:
    input:
        bed='temp/contigs/hardmask/regions_noslop.bed',
        sizes='contigs/contigs.sizes'
    output:
        bed=temp('temp/contigs/hardmask/regions.bed')
    params:
        flank=HARDMASK_FLANK
    shell:
        """bedtools slop -i {input.bed} -g {input.sizes} -b {params.flank} """
        """>{output.bed}"""

# gt_contig_hardmask_regions
#
# Get contig regions using SV breakpoints.
rule gt_contig_hardmask_regions:
    input:
        bed='sv_calls/sv_calls.bed'
    output:
        bed=temp('temp/contigs/hardmask/regions_noslop.bed')
    run:

        # Read table
        df = pd.read_table(input.bed, header=0, usecols=('CONTIG', 'CONTIG_START', 'CONTIG_END'))
        df = df.loc[:, ('CONTIG', 'CONTIG_START', 'CONTIG_END')]
        df.columns = ('#CHROM', 'POS', 'END')

        # Sort
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False)
