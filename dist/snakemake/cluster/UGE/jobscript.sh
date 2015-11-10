#!/bin/sh
# properties = {properties}
#$ -S /bin/bash
#ulimit -c 0
#$ -V
#$ -cwd
#$ -e ./log
#$ -o ./log
#$ -q all.q
set -o pipefail

# Trap signals from UGE and exit gracefully.
trap 'echo "Job killed by UGE: too much wallclock time"; exit 2' SIGUSR1
trap 'echo "Job killed by UGE: too much memory"; exit 2' SIGXCPU

{exec_job}
exit 0
