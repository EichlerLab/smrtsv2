"""
Structural variant caller for PacBio reads.

See also: https://github.com/EichlerLab/pacbio_variant_caller
"""
import math
import os
import tempfile

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

# Always set the environment
LD_LIBRARY_PATH = config.get("ld_path")
PATH = config.get("path")
PERL5LIB = config.get("perl5lib")

shell.prefix("export PATH=\"%s\"; export LD_LIBRARY_PATH=\"%s\"; export PERL5LIB=\"%s\"; " % (PATH, LD_LIBRARY_PATH, PERL5LIB))

#
# Define internal constants.
#
CWD = os.getcwd()

#
# Load user variables.
#

# If the user has a config file in the current working directory, use
# that. Otherwise, use SMRT SV defaults.
if os.path.exists("config.json"):
    configfile: "config.json"
else:
    configfile: "%s/config.template.json" % SNAKEMAKE_DIR

TMP_DIR = config.get("tmp_dir", tempfile.gettempdir())
EVENT_TYPES = ("insertion", "deletion")

#
# Include rules.
#

if not "genotyper_config" in config:
    CHROMOSOME_LENGTHS = config.get("reference_index", "%s.fai" % config["reference"])

    # TODO: fix bug caused by Snakemake not understanding more than one dynamic
    # output type per file.
    include: "rules/prepare_reference.rules"

    # Only include alignment rules if alignments aren't defined already or don't
    # exist yet.
    #if config.get("alignments") is None or not os.path.exists(config.get("alignments")):

    include: "rules/alignment.rules"
    include: "rules/sv_candidates.rules"
    include: "rules/local_assembly.mhap_celera.rules"
    include: "rules/variant_caller.rules"
else:
    include: "rules/genotyper.rules"

#
# Determine which outputs to create.
#

OUTPUTS = []

if config.get("detection"):
    OUTPUTS.extend([
        "sv_candidate_lengths.pdf",
        "sv_candidate_support.pdf",
        "assembly_candidates.bed"
    ])

# if config.get("assembly") and config["assembly"].get("regions_to_assemble"):
#     OUTPUTS.append("sv_assemblies.txt")

# if config.get("gap_extension") and config["gap_extension"].get("regions_to_assemble"):
#     OUTPUTS.append("gap_assemblies.txt")

#
# Define rules.
#

# Create list of all final outputs.
rule all:
    input: OUTPUTS
