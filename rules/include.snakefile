"""
Sets up common constructs needed by Snakefile called by Snakemake.
"""

import os
import tempfile


#
# Set SMRTSV locations
#

RULES_DIR = os.path.dirname(workflow.snakefile)
SMRTSV_DIR = os.path.dirname(RULES_DIR)


# Always set the environment
LD_LIBRARY_PATH = config.get("ld_path")
PATH = config.get("path")
PERL5LIB = config.get("perl5lib")

#shell.prefix("export PATH=\"%s\"; export LD_LIBRARY_PATH=\"%s\"; export PERL5LIB=\"%s\"; " % (PATH, LD_LIBRARY_PATH, PERL5LIB))

shell.prefix('set -euo pipefail')

#
# Define internal constants.
#

WORKING_DIR = os.getcwd()


#
# Load tiered configurations
#

# Set locations
CONFIG_DEFAULT = os.path.join(SMRTSV_DIR, 'config.default.json')
CONFIG_LOCAL = os.path.join(WORKING_DIR, 'config.json')

# Load default
configfile: CONFIG_DEFAULT

# Load local configuration
if os.path.exists(CONFIG_LOCAL):
    configfile: CONFIG_LOCAL

#
# Get temporary directory
#

TEMP_DIR = config.get('tempdir', None)

if TEMP_DIR is None or TEMP_DIR == '':
    TEMP_DIR = tempfile.gettempdir()

os.makedirs(TEMP_DIR, exist_ok=True)
