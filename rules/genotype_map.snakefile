"""
Called by genotype.snakefile to map sample reads and perform post-alt processing on the alignments. The
output BAM file aligns reads to the primary contigs (reference the SVs were called against) and alternate contigs
(local assemblies the SVs were found on).
"""

import pandas as pd
import socket


###################
### Definitions ###
###################

# Get parameters
SAMPLE = config['sample']                      # Sample name
SAMPLE_ALN = config['sample_aln']              # SAM/BAM/CRAM of aligned input reads
SV_REF = config['sv_ref']                      # Augmented reference with local assemblies as alternate contigs.
SV_REF_ALT = config['sv_ref_alt']              # Reference ALT file (.alt); aligns alternate to primary contigs.
MAP_REGIONS_BED = config['map_regions_bed']    # Filter alignments outside these loci
OUTPUT_ALN = config['output_aln']              # Output SAM/BAM/CRAM file with reads mapped to the primary contigs and alternate contigs.
OUTPUT_ALN_INDEX = config['output_aln_index']  # Output SAM/BAM/CRAM file index.

THREADS = config['threads']                  # Number of alignment threads
MAPQ = config['mapq']                        # Minimum mapping quality of reads once re-aligned to the augmented reference.
SMRTSV_DIR = config['smrtsv_dir']            # Directory where SMRTSV is installed.
POSTALT_PATH = config['postalt_path']        # Full path to bwa-postalt.js

BM_DIR = os.path.join(os.path.dirname(OUTPUT_ALN), 'bm')  # Path to benchmark logs

# Get hostname
HOSTNAME = socket.gethostname()

# Create subdirectories of the temp directories

os.makedirs('collate', exist_ok=True)  # For samtools collate
os.makedirs('sort', exist_ok=True)  # For samtools sort


#############
### Rules ###
#############

# gt_map_align_sample_reads
#
# Map reads to the augmented reference (primary + local-assembly-contigs).
rule gt_map_align_sample_reads:
    output:
        aln=OUTPUT_ALN,
        aln_index=OUTPUT_ALN_INDEX
    benchmark:
        os.path.join(BM_DIR, 'map_sample_reads.txt')
    shell:
        """echo "Running alignment for sample {SAMPLE}"; """
        """echo "Input Reads: {SAMPLE_ALN}"; """
        """echo "Temp: $(pwd)"; """
        """echo "Host: {HOSTNAME}"; """
        """samtools collate {SAMPLE_ALN} -Ou -n 2048 - | """
        """samtools bam2fq - | """
        """seqtk dropse - | """
        """bwa mem -R '@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}' -p -t {THREADS} {SV_REF} - | """
        """k8 {POSTALT_PATH} {SV_REF_ALT} | """
        """grep -Ev '^\s*NaN\s*' | """
        """samtools view -h -q {MAPQ} -L {MAP_REGIONS_BED} | """
        """samtools sort -@ 4 -m1G -T sort/sample_ | """
        """samtools view -T {SV_REF} -o {output.aln}; """
        """echo "Indexing mapped reads..."; """
        """samtools index {output.aln}; """
        """echo "Done mapping sample: {SAMPLE}" """
