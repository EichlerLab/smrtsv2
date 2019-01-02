#!/usr/bin/env python3

"""
Utilities for `smrtsv.py`.
"""

import os
import subprocess
import sys

from smrtsvlib.args import get_arg


# List of relative paths for PATH
INSTALL_PATH = [
    'dep/bin'
]

# List of relative paths for LD_LIBRARY_PATH
INSTALL_LD_PATH = [
    'dep/lib'
]


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

    # Return environment variables
    return process_env


def run_cmd(args, process_env, stdout=None, stderr=None, cwd=None):
    """
    Run a command with the proper environment set.

    :param args: A tuple of arguments starting with the command name.
    :param process_env: A dictionary of environment variables to be set for the process.
    :param stdout: Specify the output stream. This argument is passed directly to `subprocess.Popen`.
    :param stderr: Specify the error stream. This argument is passed directly to `subprocess.Popen`.
    :param cwd: Switch to this directory before executing the command. If `None`, uses current working directory.

    :return: Return code. Negative numbers is the negation of a POSIX signal that killed the process. -1024 is returned
        if the subprocess module did not give a return code.
    """

    sys.stdout.flush()

    p = subprocess.Popen(args=args, env=process_env, stdout=stdout, stderr=stderr, cwd=cwd)

    p.wait()

    ret_code = p.returncode

    return ret_code if ret_code is not None else -1024


def run_snake_target(snakefile, args, process_env, smrtsv_dir, cmd,
                     stdout=None, stderr=None,
                     cwd=None,
                     cmd_log=None,
                     resources=None
                     ):
    """
    Run a snakemake target.

    :param args: Arguments processed from the command line.
    :param cmd: The command to run as a tuple starting with the name of the snakemake target.
    :param stdout: Specify the output stream. This argument is passed directly to `subprocess.Popen`.
    :param stderr: Specify the error stream. This argument is passed directly to `subprocess.Popen`.
    :param cwd: Switch to this directory before executing the command. If `None`, uses current working directory.
    :param cmd_log: Default to this location for cluster logs if `log` is not explicitly set in `args`. This gives each
        SMRT-SV command a
        way to set it's default log file directory and separate cluster logs
    :param resources: Append "--resources" to command.

    :return: Return code from snakemake.
    """

    cmd = list(cmd)

    # Get args and attributes
    if args is None:
        args = EmptyArgs()

    wait_time = get_arg('wait_time', args)
    dry_run = get_arg('dryrun', args)
    cluster_params = get_arg('cluster_params', args)
    cluster_config_path = get_arg('cluster_config', args, default_none=True)
    job_prefix = get_arg('job_prefix', args)

    if cluster_config_path is None:
        cluster_config_path = os.path.join(smrtsv_dir, 'cluster.template.json')

    # Get log directory
    log = get_arg('log', args, default_none=True)

    if log is None:
        if cmd_log is None:
            log = 'log'
        else:
            log = cmd_log

    # Setup snakemake command
    prefix = [
        'snakemake',
        '--snakefile', os.path.join(smrtsv_dir, snakefile),
        '--rerun-incomplete'
    ]

    # Set keep-going flag
    if hasattr(args, 'keep_going') and args.keep_going:
        prefix.append('--keep-going')

    # Set jobs
    if hasattr(args, 'jobs'):
        prefix.extend(['-j', str(args.jobs)])

    # Set dryrun
    if dry_run:
        prefix.append('-n')

    if hasattr(args, 'nt') and args.nt:
        prefix.append('--nt')

    # Set distribute
    if hasattr(args, 'distribute') and args.distribute:

        cluster_params = cluster_params.format(**{'log': log})

        prefix.extend(['--cluster-config', cluster_config_path])
        prefix.extend(('--drmaa', cluster_params, '-w', str(wait_time)))

        if not hasattr(args, 'jobs'):
            prefix.extend(['-j', '1'])

        # Set job name
        if job_prefix is None:
            prefix.extend(['--jobname', '{rulename}.{jobid}'])
        else:
            prefix.extend(['--jobname', '{}{{rulename}}.{{jobid}}'.format(job_prefix)])

        # Make log directory for distributed jobs
        if not dry_run and not os.path.isdir(log):
            os.makedirs(log, exist_ok=True)

    # Append command
    prefix.extend(cmd)

    # Append path, ld_path, and temp
    if '--config' not in cmd:
        prefix.append('--config')

    temp_dir = args.tempdir.strip() if hasattr(args, 'tempdir') and args.tempdir is not None else ''

    prefix.extend([
        'ld_path={}'.format(process_env['LD_LIBRARY_PATH']),
        'path={}'.format(process_env['PATH']),
        'tempdir={}'.format(temp_dir)
    ])

    # Append resources
    if resources is not None:
        prefix.append('--resources')
        prefix.extend(list(resources))

    # Report (verbose)
    if hasattr(args, 'verbose') and args.verbose:
        print('Running snakemake command:\n\t* {}'.format('\n\t* '.join([element if element is not None else 'Nonetype' for element in prefix])))

    # Run snakemake command
    return run_cmd(prefix, process_env, stdout=stdout, stderr=stderr, cwd=cwd)
