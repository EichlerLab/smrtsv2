"""
Snakefile for troubleshooting and manually running SMRTSV. Typical usage should be executed through `smrtsv.py`.
Individual rules may be run using this Snakefile and manually setting up the parameters rules will require.
"""


#############
### Rules ###
#############

# Base inclusion
if not 'INCLUDE_SNAKEFILE' in globals():
    include: 'rules/include.snakefile'

# Pipeline components
include: 'rules/reference.snakefile'
include: 'rules/align.snakefile'
include: 'rules/detect.snakefile'

# Plots
include: 'rules/plots/detect.snakefile'
