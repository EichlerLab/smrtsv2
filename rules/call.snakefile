import os

localrules: call_variants

SV_TYPES = ("insertion", "deletion")
INDEL_TYPES = ("insertion", "deletion")
LOCAL_ASSEMBLY_ALIGNMENTS = config.get("local_assembly_alignments", "local_assembly_alignments.bam")
VARIANTS = config.get("variants", "variants.bed")
MIN_CONTIG_LENGTH = 40000


###################
### Definitions ###
###################

def _get_call_comparison_action(wildcards):
    if wildcards.sv_type == "insertion":
        return "window"
    else:
        return "intersect"

def _get_repeat_species(wildcards):
    if "species" in config:
        return config["species"]
    else:
        return "Homo sapiens"


#############
### Rules ###
#############


#
# All variant calls
#

rule call_variant_vcf:
    input:
        'call/sv_calls.vcf',
        'call/indel_calls.vcf',
        'call/inversions.vcf'
    output:
        VARIANTS
    shell:
        """grep "^##" {input[0]} | sed '/INFO/d' > {output}; """
        """grep -h "^##" {input} | grep INFO | sort | uniq >> {output}; """
        """grep -h "^#CHROM" {input[0]} >> {output}; """
        """sed '/^#/d' {input} | sort -k 1,1 -k 2,2n >> {output}"""


#
# SNVs
#

rule call_merge_snvs_calls:
    input:
        'call/snvs.bed'
    output:
        'call/merged_snvs.bed'
    shell:
        """cut -f 1-5 {input} | """
        """sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5 | """
        """groupBy -i stdin -g 1,2,3,4,5 -c 4 -o count | """
        """sort -k 1,1 -k 2,2n -k 3,3n -k 6,6rn | """
        """groupBy -i stdin -g 1,2,3 -c 4,5,6 -o first,first,first | """
        """awk '$6 > 1' """
        """> {output}"""

rule call_find_snvs_alignments:
    input:
        reference=config['reference'],
        alignments=LOCAL_ASSEMBLY_ALIGNMENTS
    output:
        'call/snvs.bed'
    params:
        min_contig_length=str(MIN_CONTIG_LENGTH)
    shell:
        """samtools view {input.alignments} | """
        """{SNAKEMAKE_DIR}/scripts/PrintGaps.py {input.reference} /dev/stdin """
            """--minLength 0 --maxLength 0 --minContigLength {params.min_contig_length} """
            """--outFile /dev/null --snv {output}"""


#
# Small insertion/deletion (indel) calls
#

rule call_convert_indel_bed_to_vcf:
    input:
        'call/indel_calls.bed',
        config['reference']
    output:
        'call/indel_calls.vcf'
    params:
        sample=config.get('sample', 'UnnamedSample')
    shell:
        """{SNAKEMAKE_DIR}/scripts/variants_bed_to_vcf.py {input} {output} {params.sample} indel"""

rule call_call_indels:
    input:
        expand('call/indel_calls/{indel_type}.tab', indel_type=INDEL_TYPES)
    output:
        'call/indel_calls.bed'
    shell:
        """sort -k 1,1 -k 2,2n {input} > {output}"""

rule call_combine_annotations_for_indel_type:
    input:
        'call/indel_calls/{indel_type}/filtered_gaps.bed',
        'call/indel_calls/{indel_type}/support_by_assembled_contigs.txt',
        'call/indel_calls/{indel_type}/assembled_contigs_coverage.txt',
        'call/indel_calls/{indel_type}/read_coverage.txt'
    output:
        'call/indel_calls/{indel_type}.tab'
    shell:
        """paste {input} | """
        """awk 'OFS="\\t" {{ print $0,"{wildcards.indel_type}" }}' """
        """> {output}"""

rule call_annotate_support_from_assembled_contigs_for_indels:
    input:
        'call/indel_calls/{indel_type}/gaps_2bp_or_more_without_homopolymers.bed'
    output:
        'call/indel_calls/{indel_type}/support_by_assembled_contigs.txt'
    shell:
        """cut -f 6 {input} > {output}"""

rule call_annotate_strs_in_indels:
    input:
        'call/indel_calls/{indel_type}/gaps_2bp_or_more_without_homopolymers.bed',
        'call/strs_in_reference.bed'
    output:
        'indel_calls/{indel_type}/strs.txt'
    shell:
        """cut -f 1-3 {input[0]} | """
        """bedtools intersect -a stdin -b {input[1]} -loj -sorted | """
        """sed 's/\\t\./\\t0/g' | """
        """bedtools groupby -c 7 -o first -full | """
        """cut -f 7 """
        """> {output}"""

rule call_annotate_coverage_of_pacbio_reads_for_indels:
    input:
        'call/indel_calls/{indel_type}/gaps_2bp_or_more_without_homopolymers.bed',
        'call/coverage.bed'
    output:
        'call/indel_calls/{indel_type}/read_coverage.txt'
    shell:
        """cut -f 1-3 {input[0]} | """
        """bedtools intersect -a stdin -b {input[1]} -loj -sorted | """
        """sed 's/\\t\./\\t0/g' | """
        """bedtools groupby -c 7 -o mean -full | """
        """cut -f 8 | """
        """awk '{{ printf("%2.2f\\n", $1) }}' """
        """> {output}"""

rule call_annotate_coverage_of_assembled_contigs_for_indels:
    input:
        'call/indel_calls/{indel_type}/gaps_2bp_or_more_without_homopolymers.bed',
        'call/assembled_contigs.depth.bed'
    output:
        'call/indel_calls/{indel_type}/assembled_contigs_coverage.txt'
    shell:
        """cut -f 1-3 {input[0]} | """
        """bedtools intersect -a stdin -b {input[1]} -loj -sorted | """
        """sed 's/\\t\./\\t0/g' | """
        """bedtools groupby -c 7 -o max -full | """
        """cut -f 8 """
        """> {output}"""

rule call_calculate_coverage_from_assembled_contigs:
    input:
        reference=config['reference'],
        alignments=LOCAL_ASSEMBLY_ALIGNMENTS
    output:
        'call/assembled_contigs.depth.bed'
    shell:
        """bedtools bamtobed -i {input.alignments} | """
        """{SNAKEMAKE_DIR}/scripts/BedIntervalsToDepth.py /dev/stdin {input.reference} --out /dev/stdout | """
        """sort -k 1,1 -k 2,2n > {output}"""

rule call_cut_indel_events_by_columns:
    input:
        'call/indel_calls/{indel_type}/gaps_2bp_or_more_without_homopolymers.bed'
    output:
        'call/indel_calls/{indel_type}/filtered_gaps.bed'
    shell:
        """cut -f 1-5 {input} > {output}"""

rule call_filter_indel_events_by_size:
    input:
        'call/indel_calls/{indel_type}/gaps_without_homopolymers.bed'
    output:
        'call/indel_calls/{indel_type}/gaps_2bp_or_more_without_homopolymers.bed'
    shell:
        """awk '$4 >= 2' {input} | bedtools groupby -c 6 -o max -full > {output}"""

rule call_remove_indel_events_in_homopolymers:
    input:
        'call/indel_calls/{indel_type}/gaps.bed'
    output:
        'call/indel_calls/{indel_type}/gaps_without_homopolymers.bed'
    run:
        command = (
            """awk '$11 == "F"' {input} | """
            """cut -f 1,2,3,5,6,8 | """
            """sort -k 1,1 -k 2,2n -k 4,4n | """
            """{SNAKEMAKE_DIR}/scripts/PrintSNVSupport.py /dev/stdin /dev/stdout"""
        )

        if wildcards.indel_type == "insertion":
            command = "%s | {SNAKEMAKE_DIR}/scripts/BedMod.py /dev/stdin {output} --leftjustify 1" % command
        else:
            command = "%s > {output}" % command

        shell(command)

rule call_split_indels_by_type:
    input:
        'call/indel_calls/gaps.tiled.bed'
    output:
        'call/indel_calls/{indel_type}/gaps.bed'
    shell:
        """grep {wildcards.indel_type} {input} | sort -k 1,1 -k 2,2n > {output}"""

rule call_filter_indel_gaps_by_tiling_path:
    input:
        'call/indel_calls/gaps.bed',
        'call/tiling_contigs.tab'
    output:
        'call/indel_calls/gaps.tiled.bed',
        'call/indel_calls/gaps.tiled.log'
    shell:
        """{SNAKEMAKE_DIR}/scripts/FilterGapsByTilingPath.py {input} """
        """> {output[0]} 2> {output[1]}"""

rule call_find_indel_gaps_in_alignments:
    input:
        reference=config['reference'],
        alignments=LOCAL_ASSEMBLY_ALIGNMENTS
    output:
        'call/indel_calls/gaps.bed'
    params:
        indel_pack_distance='0',
        min_contig_length=str(MIN_CONTIG_LENGTH)
    shell:
        """samtools view {input.alignments} | """
        """{SNAKEMAKE_DIR}/scripts/PrintGaps.py {input.reference} /dev/stdin """
            """--minLength 0 --maxLength 50 --context 6 --removeAdjacentIndels --onTarget """
            """--minContigLength {params.min_contig_length} --condense {params.indel_pack_distance} """
            """--outFile {output}"""


#
# Structural variation (SV) calls
#

rule call_convert_sv_bed_to_vcf:
    input:
        'call/sv_calls/sv_calls_with_repeats.bed',
        config['reference']
    output:
        'call/sv_calls.vcf'
    params:
        sample=config.get('sample', 'UnnamedSample')
    shell:
        """{SNAKEMAKE_DIR}/scripts/variants_bed_to_vcf.py {input} {output} {params.sample} sv"""

rule call_collect_all_summarized_sv_calls:
    input:
        expand('call/sv_calls/repeat_classified_{sv_type}.bed', sv_type=SV_TYPES)
    output:
        'call/sv_calls/sv_calls_with_repeats.bed'
    shell:
        """sort -k 1,1 -k 2,2n -m {input} | uniq > {output}"""

rule call_collect_summarized_sv_calls_within_type:
    input:
        'call/sv_calls/summarized_{sv_type}'
    output:
        'call/sv_calls/repeat_classified_{sv_type}.bed'
    shell:
        """for file in {input}/*.bed; do """
            """repeat_type=`basename ${{file/.bed/}} | """
            """sed 's/\./_/g'`; """
            """awk -v repeat_type=$repeat_type 'OFS="\\t" {{ print $0,repeat_type }}' $file; """
        """done | """
        """sort -k 1,1 -k 2,2n | """
        """uniq > """
        """{output}"""

rule call_summarize_calls_by_repeat_type:
    input:
        'call/sv_calls/all_annotated_with_trf.{sv_type}.bed'
    output:
        'call/sv_calls/summarized_{sv_type}'
    shell:
        """mkdir -p {output}; """
        """awk '$20 > 0.8' {input} > {output}/TRF.bed; """
        """awk '$20 <= 0.8' {input} > sv_calls/not_trf.bed; """
        """{SNAKEMAKE_DIR}/scripts/FixMasked.py 17 < sv_calls/not_trf.bed | """
            """awk '$18 < 0.7' > {output}/NotMasked.bed; """
        """{SNAKEMAKE_DIR}/scripts/FixMasked.py 17 < sv_calls/not_trf.bed | """
            """awk '$18 >= 0.7' > sv_calls/repeat.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py sv_calls/repeat.bed --prefix AluY --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder {output}/1.bed > {output}/AluY.simple.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/1.bed --prefix AluS --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder {output}/2.bed > {output}/AluS.simple.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/2.bed --minSTR 1  --maxNotPrefix 0 --remainder {output}/4.bed > {output}/STR.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/4.bed --prefix L1HS  --maxNotPrefix 0 --remainder {output}/5.bed > {output}/L1HS.simple.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/5.bed --prefix Alu  --minPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder {output}/6.bed > {output}/Alu.Mosaic.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/6.bed --prefix Alu  --minSTR 1 --minPrefix 1 --maxNotPrefix 0 --remainder {output}/7.bed > {output}/Alu.STR.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/7.bed --prefix ALR   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/8.bed > {output}/ALR.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/8.bed --prefix SVA   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/9.bed > {output}/SVA.simple.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/9.bed --prefix HERV   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/10.bed > {output}/HERV.simple.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/10.bed --prefix L1P   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/11.bed > {output}/L1P.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/11.bed --prefix BSR/Beta   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/12.bed > {output}/Beta.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/12.bed --prefix HSAT   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/13.bed > {output}/HSAT.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/13.bed --prefix MER   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/14.bed > {output}/MER.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/14.bed --prefix L1   --minPrefix 1 --maxNotPrefix 0 --remainder {output}/15.bed > {output}/L1.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/15.bed --prefix LTR  --minPrefix 1 --maxNotPrefix 0 --remainder {output}/16.bed > {output}/LTR.bed; """
        """{SNAKEMAKE_DIR}/scripts/PrintUniqueEvents.py {output}/16.bed --max 1 --remainder {output}/17.bed > {output}/Singletons.bed; """
        """mv -f {output}/17.bed {output}/Complex.bed; """
        """rm -f {output}/[0-9]*.bed"""

rule call_merge_sv_calls:
    input:
        expand('call/sv_calls/all_annotated.{sv_type}.bed', sv_type=SV_TYPES)
    output:
        'call/sv_calls.bed'
    shell:
        "sort -k 1,1 -k 2,2n {input} > {output}"

# Filter out non-repetitive SVs or SVs with less than 70% masked sequence.
rule call_filter_annotated_sv_calls:
    input:
        'call/sv_calls/all_annotated.{sv_type}.bed'
    output:
        'call/sv_calls/annotated.{sv_type}.bed'
    shell:
        """"grep -v NONE {input} | """
        """awk '$18 >= 0.70' > {output}"""

rule call_annotate_sv_calls_with_trf:
    input:
        'call/sv_calls/all_annotated.{sv_type}.bed',
        'call/sv_calls/{sv_type}/rm/{sv_type}.fasta.trf'
    output:
        'call/sv_calls/all_annotated_with_trf.{sv_type}.bed'
    shell:
        """{SNAKEMAKE_DIR}/scripts/AnnotateWithTRF.py {input} {output}"""

rule call_annotate_sv_calls_with_repeatmasker:
    input:
        calls='call/sv_calls/calls.{sv_type}.bed',
        repeats='call/sv_calls/{sv_type}/rm/{sv_type}.fasta.out',
        masked_fasta='call/sv_calls/{sv_type}/rm/{sv_type}.fasta.masked'
    output:
        'call/sv_calls/all_annotated.{sv_type}.bed'
    shell:
        """{SNAKEMAKE_DIR}/scripts/AnnotateGapBed.py {input.calls} {output} {input.repeats} {input.masked_fasta}"""

rule call_trf_mask_sv_fasta:
    input:
        'call/sv_calls/{sv_type}/{sv_type}.fasta'
    output:
        'call/sv_calls/{sv_type}/rm/{sv_type}.fasta.trf'
    shell:
        """{SNAKEMAKE_DIR}/bin/trf {input} 2 7 7 80 10 20 500 -m -ngs -h """
        """> {output}"""

rule call_repeatmask_sv_fasta:
    input:
        'call/sv_calls/{sv_type}/{sv_type}.fasta'
    output:
        'call/sv_calls/{sv_type}/rm/{sv_type}.fasta.out',
        'call/sv_calls/{sv_type}/rm/{sv_type}.fasta.masked'
    params:
        threads='8',
        species=_get_repeat_species
    shell:
        """RepeatMasker -species "{params.species}" -dir `dirname {output[0]}` -xsmall -no_is -e wublast -s -pa {params.threads} {input}"""

rule call_create_sv_fasta:
    input:
        'call/sv_calls/calls.{sv_type}.bed'
    output:
        'call/sv_calls/{sv_type}/{sv_type}.fasta'
    shell:
        """{SNAKEMAKE_DIR}/scripts/GapBedToFasta.py {input} {output}""""

# call_identify_calls_by_type
#
# Filter input gaps by SV type and omit gaps containing N bases. Insertion
# coordinates are resized to their single-bp representation prior to
# clustering and returned to their original values after clustering. After
# clustering each call is annotated by the coverage of local assemblies.
rule call_identify_calls_by_type:
    input:
        gaps='call/sv_calls/gaps.bed',
        alignments=LOCAL_ASSEMBLY_ALIGNMENTS
    output:
        'call/sv_calls/calls.{sv_type}.bed'
    params:
        call_comparison_action=_get_call_comparison_action,
        window='20',
        overlap='0.5'
    shell:
        """awk '$4 == "{wildcards.sv_type}" && index($6, "N") == 0' {input.gaps} | """
        """awk 'OFS="\\t" {{ if ("{wildcards.sv_type}" == "insertion") {{ $3=$2 + 1 }} print }}' | """
        """python {SNAKEMAKE_DIR}/scripts/cluster_calls.py --window {params.window} --reciprocal_overlap {params.overlap} /dev/stdin {params.call_comparison_action} | """
        """awk 'OFS="\\t" {{ if ("{wildcards.sv_type}" == "insertion") {{ $3=$2 + $5 }} print }}' | """
        """sort -k 1,1 -k 2,2n | """
        """while read line; do set -- $line; coverage=`samtools view -c {input.alignments} $1:$2-$3`; echo -e "$line\\t$coverage"; done """
        """> {output}"""

rule call_find_calls_by_gaps_in_alignments:
    input:
        reference=config['reference'],
        alignments=LOCAL_ASSEMBLY_ALIGNMENTS
    output:
        'call/sv_calls/gaps.bed'
    params:
        tsd_length='20',
        indel_pack_distance='20'
    shell:
        """samtools view -h {input.alignments} | """
        """{SNAKEMAKE_DIR}/scripts/PrintGaps.py {input.reference} /dev/stdin --qpos --condense {params.indel_pack_distance} --tsd {params.tsd_length} | """
        """sort -k 1,1 -k 2,2n """
        """> {output}"""

rule call_tile_contigs_from_alignments:
    input:
        LOCAL_ASSEMBLY_ALIGNMENTS
    output:
        'call/tiling_contigs.tab'
    shell:
        """samtools view -h {input} | """
        """{SNAKEMAKE_DIR}/scripts/TilingPath.py /dev/stdin """
        """> {output}"""


#
# Inversions
#

rule call_convert_inversion_bed_to_vcf:
    input:
        'call/sv_calls/merged_inversions.bed',
        config['reference']
    output:
        'call/inversions.vcf'
    params:
        sample=config.get('sample', 'UnnamedSample')
    shell:
        """{SNAKEMAKE_DIR}/scripts/variants_bed_to_vcf.py {input} {output} {params.sample} inversion"""

rule call_merge_inversions:
    input:
        inversions='call/sv_calls/inversions.bed',
        alignments=LOCAL_ASSEMBLY_ALIGNMENTS
    output:
        'call/sv_calls/merged_inversions.bed'
    run:
        # Only try to sort and merge inversions if any exist.
        if os.stat(input.inversions).st_size > 0:
            shell(
                """sort -k 1,1 -k 2,2n {input.inversions} | """
                """bedtools merge -i stdin -d 0 -c 4 -o count | """
                """while read line; do """
                    """set -- $line; """
                    """coverage=`samtools view -c {input.alignments} $1:$2-$3`; """
                    """echo -e "$line\\t$coverage"; """
                """done | ""
                """awk 'OFS="\\t" {{ print $1,$2,$3,"inversion",$4,$5 }}' > """
                """{output}"""
            )
        else:
            shell("touch {output}")

rule call_find_inversions:
    input:
        alignments=LOCAL_ASSEMBLY_ALIGNMENTS,
        reference=config['reference']
    output:
        bed='call/sv_calls/inversions.bed'
    params:
        reference_window='5000',
        threads='8'
    shell:
        """samtools view {input.alignments} | """
        """{SNAKEMAKE_DIR}/scripts/mcst/screenInversions """
            """/dev/stdin {input.reference} {output.bed} -w {params.reference_window} -r --noClip -j {params.threads}"""
