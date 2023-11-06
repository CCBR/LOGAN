# LOGAN 🔬 [![Docker Pulls](https://img.shields.io/docker/pulls/nciccbr/ccbr_wes_base)](https://hub.docker.com/r/nciccbr/ccbr_wes_base) [![GitHub issues](https://img.shields.io/github/issues/ccbr/LOGAN?color=brightgreen)](https://github.com/ccbr/LOGAN/issues)  [![GitHub license](https://img.shields.io/github/license/ccbr/LOGAN)](https://github.com/ccbr/LOGAN/blob/master/LICENSE) 

> **_LOGAN-whoLe genOme-sequencinG Analysis pipeliNe_**. This is the home of the LOGAN Pipeline. Its long-term goals: to accurately call germline and somatic variants, to infer SVs and CNVs, and to annotate variants!

## Overview
Welcome to LOGAN! Before getting started, we highly recommend reading through [LOGAN's documentation](https://ccbr.github.io/logan).

The **`LOGAN`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>logan <b>--fastq</b></code>](https://ccbr.github.io/LOGAN/usage/run/): Run the WGS pipeline with your fastq input files.

LOGAN is a comprehensive whole genome-sequencing pipeline following the Broad's set of best practices. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Nextflow<sup>2</sup>](https://nextflow.io/), a flexible and scalable workflow management system, to submit jobs to a cluster or cloud provider.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FASTQ or BAM files and can be run locally on a compute instance, on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM.

Before getting started, we highly recommend reading through the [usage](https://ccbr.github.io/LOGAN/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/CCBR/LOGAN/issues).

Original pipelining and code forked heavily from the CCBR Exome-seek Pipeline [Exome-seek](https://github.com/CCBR/XAVIER)

## Dependencies
**Requires:** `singularity>=3.5`  `nextflow>=22.10.2`

[singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step relies on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Nextflow uses singularity to pull these images onto the local filesystem prior to job execution, and as so, nextflow and singularity are the only two dependencies.

## Installation
Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/ccbr/LOGAN.git
# Change your working directory
cd LOGAN/
```

## Contribute 
This site is a living document, created for and by members like you. LOGAN is maintained by the members of CCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [repository](https://github.com/ccbr/LOGAN/pulls).


## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
