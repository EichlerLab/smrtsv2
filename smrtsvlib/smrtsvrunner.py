#!/usr/bin/env python3

"""
Utilities for `smrtsv.py`.
"""

import os
import subprocess
import sys


# List of relative paths for PATH
INSTALL_PATH = [
    'bin',
    'dist/miniconda/envs/python2/bin',
    'dist/miniconda/envs/python3/bin',
    'dist/miniconda/bin',
    'dist/celera/wgs-8.3rc2/Linux-amd64/bin/',
    'dist/amos-3.1.0/bin',
    'canu/Linux-amd64/bin'
]

# List of relative paths for LD_LIBRARY_PATH
INSTALL_LD_PATH = [
    'dist/hdf5/lib'
]

# Cluster parameters
CLUSTER_SETTINGS = ' -V -cwd -e ./log -o ./log {cluster.params} -w n -S /bin/bash'
CLUSTER_FLAG = ('--drmaa', CLUSTER_SETTINGS, '-w', '120')

# Empty arguments
class EmptyArgs:
    pass


def get_env(install_dir):
    """
    Get environment variables

    :param install_dir: The directory where SMRTSV is installed (contains 'smrtsv.py').

    :return: A hash of all environment variables with modifications made by this function.
    """

    # Setup environment for executing commands
    process_env = os.environ.copy()

    # Prepend PATH
    process_env_path = ':'.join([os.path.join(install_dir, this_path) for this_path in INSTALL_PATH])

    if 'PATH' in process_env:
        process_env['PATH'] = process_env_path + ':' + process_env['PATH']
    else:
        process_env['PATH'] = process_env_path

    # Prepend LD_LIBRARY_PATH
    process_env_ld_path = ':'.join([os.path.join(install_dir, this_path) for this_path in INSTALL_LD_PATH])

    if 'LD_LIBRARY_PATH' in process_env:
        process_env['LD_LIBRARY_PATH'] = process_env_ld_path + ':' + process_env['LD_LIBRARY_PATH']
    else:
        process_env['LD_LIBRARY_PATH'] = process_env_ld_path

    os.environ['LD_LIBRARY_PATH'] = process_env['LD_LIBRARY_PATH']

    # Return environment variables
    return process_env


def run_cmd(args, process_env, stdout=None, stderr=None):
    """
    Run a command with the proper environment set.

    :param args: A tuple of arguments starting with the command name.
    :param process_env: A dictionary of environment variables to be set for the process.
    :param stdout: Specify the output stream. This argument is passed directly to `subprocess.Popen`.
    :param stderr: Specify the error stream. This argument is passed directly to `subprocess.Popen`.

    :return: Return code. Negative numbers is the negation of a POSIX signal that killed the process. -1024 is returned
        if the subprocess module did not give a return code.
    """

    sys.stdout.flush()

    p = subprocess.Popen(args=args, env=process_env, stdout=stdout, stderr=stderr)

    p.wait()

    ret_code = p.returncode

    return ret_code if ret_code is not None else -1024


def run_snake_target(snakefile, args, process_env, smrtsv_dir, cmd, stdout=None, stderr=None):
    """
    Run a snakemake target.

    :param args: Arguments processed from the command line.
    :param cmd: The command to run as a tuple starting with the name of the snakemake target.
    :param stdout: Specify the output stream. This argument is passed directly to `subprocess.Popen`.
    :param stderr: Specify the error stream. This argument is passed directly to `subprocess.Popen`.

    :return: Return code from snakemake.
    """

    cmd = list(cmd)

    if args is None:
        args = EmptyArgs()

    # Use the user-defined cluster config path if one is given. Otherwise, use
    # an empty config that comes with the SMRT-SV distribution.
    if hasattr(args, 'cluster_config') and args.cluster_config is not None:
        cluster_config_path = args.cluster_config
    else:
        cluster_config_path = os.path.join(smrtsv_dir, 'cluster.template.json')

    # Setup snakemake command
    prefix = [
        'snakemake',
        '--snakefile', os.path.join(smrtsv_dir, snakefile),
        '-T',
        '--rerun-incomplete'
    ]

    # Set jobs
    if hasattr(args, 'jobs'):
        prefix.extend(['-j', str(args.jobs)])

    # Set dryrun
    if hasattr(args, 'dryrun') and args.dryrun:
        prefix.append('-n')

    if hasattr(args, 'nt') and args.nt:
        prefix.append('--nt')

    # Set distribute
    if hasattr(args, 'distribute') and args.distribute:
        prefix.extend(['--cluster-config', cluster_config_path])
        prefix.extend(CLUSTER_FLAG)

        if not hasattr(args, 'jobs'):
            prefix.extend(['-j', '1'])

    # Append command
    prefix.extend(cmd)

    # Append path and ld_path
    if '--config' not in cmd:
        prefix.append('--config')

    prefix.extend([
        'ld_path={}'.format(process_env['LD_LIBRARY_PATH']),
        'path={}'.format(process_env['PATH'])
    ])

    # Report (verbose)
    if hasattr(args, 'verbose') and args.verbose:
        print('Running snakemake command: {}'.format(' '.join(prefix)))

    # Run snakemake command
    return run_cmd(prefix, process_env, stdout=stdout, stderr=stderr)
