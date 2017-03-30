"""
Defines a dictionary of command-line options.
"""

args_dict = dict()

# alignment_parameters
args_dict['alignment_parameters'] = {
    'help': 'BLASR parameters for raw read alignments',
    'default':
        '-bestn 2 '
        '-maxAnchorsPerPosition 100 '
        '-advanceExactMatches 10 '
        '-affineAlign '
        '-affineOpen 100 '
        '-affineExtend 0 '
        '-insertion 5 '
        '-deletion 5 '
        '-extend '
        '-maxExtendDropoff 50'
}

# alignments
args_dict['alignments'] = {
    'help': 'text file with one absolute path to a BLASR alignments file (.bam) per line',
    'default': 'alignments.fofn'
}

# alignments_dir
args_dict['alignments_dir'] = {
    'help': 'absolute path of directory for BLASR alignment files',
    'default': 'alignments'
}

# asm_alignment_parameters
args_dict['asm_alignment_parameters'] = {
    'default':
        '-affineAlign '
        '-affineOpen 8 '
        '-affineExtend 0 '
        '-bestn 1 '
        '-maxMatch 30 '
        '-sdpTupleSize 13',
    'help': 'BLASR parameters to use to align local assemblies',
}

# assembly_alignments
args_dict['assembly_alignments'] = {
    help: 'BAM file with BLASR alignments of local assemblies against the reference'
}

# assembly_log
args_dict['assembly_log'] = {
    'default': 'assembly.log',
    'help': 'name of log file for local assemblies'
}

# assembly_window_size
args_dict['assembly_window_size'] = {
    'type': int,
    'help': 'size of reference window for local assemblies',
    'default': 60000
}

# assembly_window_slide
args_dict['assembly_window_slide'] = {
    'type': int,
    'help': 'size of reference window slide for local assemblies',
    'default': 20000
}

# batches
args_dict['batches'] = {
    'help': 'number of batches to split input reads into such that there will be one BAM output file per batch',
    'type': int,
    'default': 1
}

# candidates
args_dict['candidates'] = {
    'help': 'BED file of candidates detected in read alignments'
}

# exclude
args_dict['exclude'] = {
    'help': 'BED file of regions to exclude from local assembly (e.g., heterochromatic sequences, etc.)'
}

# mapping_quality
args_dict['mapping_quality'] = {
    'type': int,
    'help': 'minimum mapping quality of raw reads to use for local assembly',
    'default': 30
}

# max_candidate_length
args_dict['max_candidate_length'] = {
    'type': int,
    'default': 60000,
    'help': 'maximum length allowed for an SV candidate region'
}

# max_coverage
args_dict['max_coverage'] = {
    'type': int,
    'default': 100,
    'help': 'maximum number of total reads allowed to flag a region as an SV candidate'
}

# max_support
args_dict['max_support'] = {
    'type': int,
    'default': 100,
    'help': 'maximum number of supporting reads allowed to flag a region as an SV candidate'
}

# min_coverage
args_dict['min_coverage'] = {
    'type': int,
    'default': 5,
    'help': 'minimum number of total reads required to flag a region as an SV candidate'
}

# genotyped_variants
args_dict['genotyped_variants'] = {
    'help': 'VCF of SMRT SV variant genotypes for the given sample-level BAMs'
}

# genotyper_config
args_dict['genotyper_config'] = {
    'help':
        'JSON configuration file with SV reference paths, samples to genotype as BAMs, '
        'and their corresponding references'
}

# min_hardstop_support
args_dict['min_hardstop_support'] = {
    'type': int,
    'default': 11,
    'help': 'minimum number of reads with hardstops required to flag a region as an SV candidate'
}

# min_length
args_dict['min_length'] = {
    'type': int,
    'default': 50,
    'help': 'minimum length required for SV candidates'
}

# min_support
args_dict['min_support'] = {
    'type': int,
    'default': 5,
    'help': 'minimum number of supporting reads required to flag a region as an SV candidate'
}

# minutes_to_delay_jobs
args_dict['minutes_to_delay_jobs'] = {
    'type': int,
    'default': 1,
    'help': 'maximum number of minutes to delay local assembly jobs to limit simultaneous I/O on shared storage'
}

# reads
args_dict['reads'] = {
    'help': 'text file with one absolute path to a PacBio reads file (.bax.h5) per line'
}

# rebuild_regions
args_dict['rebuild_regions'] = {
    'action': 'store_true',
    'help': 'rebuild subset of regions to assemble'
}

# reference
args_dict['reference'] = {
    'help': 'FASTA file of reference to index'
}

# threads
args_dict['threads'] = {
    'help': 'number of threads to use for each alignment job',
    'type': int,
    'default': 1
}

# sample
args_dict['sample'] = {
    'default': 'UnnamedSample',
    'help': 'Sample name to use in final variant calls'
}

# species
args_dict['species'] = {
    'default': 'human',
    'help': 'Common or scientific species name to pass to RepeatMasker'
}

# variants
args_dict['variants'] = {
    'help': 'VCF of variants called by local assembly alignments'
}
