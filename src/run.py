#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from shutil import copytree, copyfile
import os, re, json, sys, subprocess



def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def require(cmds, suggestions, path=None):
    """Enforces an executable is in $PATH
    @param cmds list[<str>]:
        List of executable names to check
    @param suggestions list[<str>]:
        Name of module to suggest loading for a given index
        in param cmd.
    @param path list[<str>]]:
        Optional list of PATHs to check [default: $PATH]
    """
    error = False
    for i in range(len(cmds)):
        available = which(cmds[i])
        if not available:
            error = True
            err("""\x1b[6;37;41m\n\tFatal: {} is not in $PATH and is required during runtime!
            └── Solution: please 'module load {}' and run again!\x1b[0m""".format(cmds[i], suggestions[i])
            )

    if error: fatal()

    return 


def exists(testpath):
    """Checks if file exists on the local filesystem.
    @param parser <argparse.ArgumentParser() object>:
        argparse parser object
    @param testpath <str>:
        Name of file/directory to check
    @return does_exist <boolean>:
        True when file/directory exists, False when file/directory does not exist
    """
    does_exist = True
    if not os.path.exists(testpath):
        does_exist = False # File or directory does not exist on the filesystem

    return does_exist


def setup(sub_args, repo_path, output_path, links=[]):
    """Setup the pipeline for execution and creates config file from templates
    @param sub_args <parser.parse_args() object>:
        Parsed arguments for run sub-command
    @param repo_path <str>:
        Path to RNA-seek source code and its templates
    @param output_path <str>:
        Pipeline output path, created if it does not exist
    @param create_nidap_folder_YN <str>:
        yes or no
    @return config <dict>:
         Config dictionary containing metadata to run the pipeline
    """
    # Check for mixed inputs,
    # inputs which are a mixture
    # of FastQ and BAM files 

    hpcget = subprocess.run(
        "scontrol show config", shell=True, capture_output=True, text=True
    )
    hpcname = ""

    if "biowulf" in hpcget.stdout:
        shorthostname = "biowulf"
        print("Thank you for running XAVIER on Biowulf")
    elif "fsitgl" in hpcget.stdout:
        shorthostname = "frce"
        print("Thank you for running XAVIER on FRCE")
    else:
        shorthostname = "biowulf"
        print("%s unknown host. Configuration files for references may not be correct. Defaulting to Biowulf config"%(hpcget))


def runner(mode, outdir, alt_cache, logger, additional_bind_paths = None, 
    threads=2,  jobname='pl:LOGAN', submission_script='runner',
    tmp_dir = '/lscratch/$SLURM_JOBID/', wait = ''):
    """Runs the pipeline via selected executor: local or slurm.
    If 'local' is selected, the pipeline is executed locally on a compute node/instance.
    If 'slurm' is selected, jobs will be submited to the cluster using SLURM job scheduler.
    Support for additional job schedulers (i.e. PBS, SGE, LSF) may be added in the future.
    @param outdir <str>:
        Pipeline output PATH
    @param mode <str>:
        Execution method or mode:
            local runs serially a compute instance without submitting to the cluster.
            slurm will submit jobs to the cluster using the SLURM job scheduler.
    @param additional_bind_paths <str>:
        Additional paths to bind to container filesystem (i.e. input file paths)
    @param alt_cache <str>:
        Alternative singularity cache location
    @param logger <file-handle>:
        An open file handle for writing
    @param threads <str>:
        Number of threads to use for local execution method
    @param masterjob <str>:
        Name of the master job
    @param wait <str>:
        "--wait" or "" ... used only while submitting job via HPC API
    @return masterjob <subprocess.Popen() object>:
    """
    # Add additional singularity bind PATHs
    # to mount the local filesystem to the 
    # containers filesystem, NOTE: these 
    # PATHs must be an absolute PATHs
    outdir = os.path.abspath(outdir)
    # Add any default PATHs to bind to 
    # the container's filesystem, like 
    # tmp directories, /lscratch
    bindpaths = "{},{}".format(outdir, os.path.dirname(tmp_dir.rstrip('/')))
    # Set ENV variable 'SINGULARITY_CACHEDIR' 
    # to output directory
    my_env = {}; my_env.update(os.environ)
    cache = os.path.join(outdir, ".singularity")
    my_env['SINGULARITY_CACHEDIR'] = cache
    if alt_cache:
        # Override the pipeline's default 
        # cache location
        my_env['SINGULARITY_CACHEDIR'] = alt_cache
        cache = alt_cache

    if additional_bind_paths:
        # Add Bind PATHs for rawdata directories
        bindpaths = "{},{}".format(additional_bind_paths,bindpaths)

    if not exists(os.path.join(outdir, 'logfiles')):
        # Create directory for logfiles
        os.makedirs(os.path.join(outdir, 'logfiles'))
    
    # Create .singularity directory for 
    # installations of snakemake without
    # setuid which creates a sandbox in
    # the SINGULARITY_CACHEDIR
    if not exists(cache):
        # Create directory for sandbox 
        # and image layers
        os.makedirs(cache)

    # Run on compute node or instance
    # without submitting jobs to a scheduler
    if mode == 'local':
        # Run pipeline's main process
        # Look into later: it maybe worth 
        # replacing Popen subprocess with a direct
        # snakemake API call: https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html
        masterjob = subprocess.Popen([
                'snakemake', '-pr', '--rerun-incomplete',
                '--use-singularity',
                '--singularity-args', "'-B {}'".format(bindpaths),
                '--cores', str(threads),
                '--configfile=config.json'
            ], cwd = outdir, stderr=subprocess.STDOUT, stdout=logger, env=my_env)

    # Submitting jobs to cluster via SLURM's job scheduler
    elif mode == 'biowulf':
        if wait=='':
            masterjob = subprocess.Popen([
                    str(os.path.join(outdir, 'resources', str(submission_script))), mode,
                    '-j', jobname, '-b', str(bindpaths),
                    '-o', str(outdir), '-c', str(cache), 
                    '-t', "'{}'".format(tmp_dir)
                ], cwd = outdir, stderr=subprocess.STDOUT, stdout=logger, env=my_env)
        else:
            masterjob = subprocess.Popen([
                    str(os.path.join(outdir, 'resources', str(submission_script))), mode,
                    '-j', jobname, '-b', str(bindpaths),
                    '-o', str(outdir), '-c', str(cache), str(wait),
                    '-t', "'{}'".format(tmp_dir)
                ], cwd = outdir, stderr=subprocess.STDOUT, stdout=logger, env=my_env)            

    return masterjob
