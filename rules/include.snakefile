"""
Sets up common constructs needed by Snakefile called by Snakemake.
"""

import os
import sys
import tempfile


############
### Init ###
############

#
# Set SMRTSV locations
#

# Must be done before imports so that smrtsvlib can be found

RULES_DIR = os.path.dirname(workflow.snakefile)

if os.path.basename(workflow.snakefile) == 'Snakefile':
    RULES_DIR = os.path.join(RULES_DIR, 'rules')

SMRTSV_DIR = os.path.dirname(RULES_DIR)

WORKING_DIR = os.getcwd()

sys.path.append(SMRTSV_DIR)


#
# Imports with modified system path
#

from smrtsvlib.args import args_dict
from smrtsvlib.args import get_arg



###################
### Definitions ###
###################

#
# Load project configuration
#

CONFIG_LOCAL = os.path.join(WORKING_DIR, 'config.json')

if os.path.exists(CONFIG_LOCAL):
    configfile: CONFIG_LOCAL


#
# Snakemake directives
#

# Explicitly BASH safemode
shell.prefix('set -euo pipefail; ')


#
# Universal functions
#

def get_config_param(param_name, as_is=False):
    """
    Get a configuration parameter or the default if it was not set. Explicitly getting the default from
    `smrtsvlib.args.args_dict` enables rules to be run by `Snakefile` without requiring every parameter be set
    manually.

    :param param_name: Name of the configuration parameter.
    :param as_is: Leave the default type as it is defined; do not translate to a string.
    """

    # Get configured parameter
    if param_name in config:
        return config[param_name]

    # Get default parameter
    val = get_arg(param_name)

    if not as_is:
        val = str(val)

    return val


#
# Set universal constants
#

SVTYPES = ['INS', 'DEL', 'INV']
INSDEL= ['INS', 'DEL']


#
# Set environment
#

# Get PATHs carried over from the caller. This ensures that PATH is not modified and that LD_LIBRARY_PATH
# is not cleared when Snakemake jobs are distributed, which Grid Engine does for security reasons.

LD_LIBRARY_PATH = config.get('ld_path', None)
PATH = config.get('path', None)

if LD_LIBRARY_PATH is not None:
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH

if PATH is not None:
    os.environ['PATH'] = PATH

PROCESS_ENV = os.environ.copy()


#
# Get temporary directory
#

TEMP_DIR = config.get('tempdir', None)

if TEMP_DIR is None or TEMP_DIR == '':
    TEMP_DIR = tempfile.gettempdir()
else:
    TEMP_DIR = os.path.abspath(TEMP_DIR)

if os.path.samefile(TEMP_DIR, '.'):
    # Defaults to local directory if a temp cannot be found. This
    # should not occur on real systems, but don't clutter the working
    # directory if it does.
    TEMP_DIR = os.path.join(TEMP_DIR, 'temp')




#
# Flag include
#

INCLUDE_SNAKEFILE = True
