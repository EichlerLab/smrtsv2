"""
Structural variant caller for PacBio reads.

See also: https://github.com/EichlerLab/pacbio_variant_caller
"""
import math
import os

#
# Define internal constants.
#
BLASR_BIN = "/net/eichler/vol20/projects/pacbio/nobackups/users/jlhudd/blasr_jlhudd/alignment/bin/blasr"
CWD = os.getcwd()

#
# Load user variables.
#
configfile: "config.json"
TMP_DIR = config["tmp_dir"]
EVENT_TYPES = ("insertion", "deletion")

#
# Include rules.
#

# TODO: fix bug caused by Snakemake not understanding more than one dynamic
# output type per file.
# include: "rules/alignment.rules"
include: "rules/sv_candidates.rules"

if config["assembly"]["assembler"] == "celera":
    include: "rules/local_assembly.celera.rules"
else:
    include: "rules/local_assembly.rules"

#
# Define rules.
#

# Create list of all final outputs.
rule all:
    input:
        "sv_candidate_lengths.pdf",
        "sv_candidate_support.pdf",
        "assembly_candidates.bed",
        "sv_assemblies.txt"
