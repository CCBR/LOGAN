# LOGAN 🔬 [![Docker Pulls](https://img.shields.io/docker/pulls/nciccbr/ccbr_wes_base)](https://hub.docker.com/r/nciccbr/ccbr_wes_base) [![GitHub issues](https://img.shields.io/github/issues/ccbr/LOGAN?color=brightgreen)](https://github.com/ccbr/LOGAN/issues)  [![GitHub license](https://img.shields.io/github/license/ccbr/LOGAN)](https://github.com/ccbr/LOGAN/blob/master/LICENSE) 

> **_LOGAN-whoLe genOme-sequencinG Analysis pipeliNe_**. Call germline and somatic variants, CNVs, and SVs and  annotate variants!

## Overview
Welcome to LOGAN! Before getting started, we highly recommend reading through [LOGAN's documentation](https://ccbr.github.io/LOGAN).

LOGAN is a comprehensive whole genome-sequencing pipeline following the Broad's set of best practices. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Nextflow<sup>2</sup>](https://nextflow.io/), a flexible and scalable workflow management system, to submit jobs to a cluster or cloud provider.

Before getting started, we highly recommend reading through the [usage](https://ccbr.github.io/LOGAN/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/ccbr/LOGAN/issues).

Original pipelining and code forked from the CCBR Exome-seek Pipeline [Exome-seek](https://github.com/CCBR/XAVIER) and [OpenOmics](https://github.com/openOmics/genome-seek)

## Dependencies
**Requires:** `singularity>=3.5`  `nextflow>=22.10.2`

[singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step relies on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Nextflow uses singularity to pull these images onto the local filesystem prior to job execution, and as so, nextflow and singularity are the only two dependencies.

## Setup
LOGAN is installed on the Biowulf in the ccbrpipeliner module.
Please clone this repository to your local filesystem using the following command:
```bash
# start an interactive node
sinteractive --mem=2g --cpus-per-task=2 --gres=lscratch:200
# load the ccbrpipeliener module
module load ccbrpipeliner
```

## Usage
LOGAN supports either

### Input Files
LOGAN supports inputs of either 
1) paired end fastq files

`--fastq_input`- A glob can be used to include all FASTQ files. Like `--fastq_input "*R{1,2}.fastq.gz"`. Globbing requires quotes.

2) Pre aligned BAM files with BAI indices 

`--bam_input`- A glob can be used to include all FASTQ files. Like `--bam_input "*.bam"`. Globbing requires quotes.

3) A sheet that indicates the sample name and either FASTQs or BAM file locations

`--fastq_file_input`-  A headerless tab delimited sheet that has the sample name, R1, and R2 file locations

`--bam_file_input` -  A headerless tab delimited sheet that has the sample name, bam, and bam index (bai) file locations

### Operating Modes

#### 1.  Paired Tumor/Normal Mode 

Required for Paired Tumor/Normal Mode

`--sample_sheet` In Paired mode a sample sheet must be provided with the basename of the Tumor and Normal samples. This sheet must be Tab separated with a header for Tumor and Normal.


#### 2.  Tumor only mode

No flags are required

#### Calling Mode

Adding flags determines SNV (germline and/or somatic), SV, and/or CNV calling modes

`--vc`- Enables somatic SNV calling using mutect2, vardict, varscan, octopus, sage, MUSE (TN only), and lofreq (TN only)
  

`--germline`- Enables germline using DV

`--sv`- Enables somatic SV calling using Manta and SVABA

`--cnv`- Enables somatic CNV calling using FREEC, Sequenza, and Purple (hg38 only)



#### Optional Arguments
`--indelrealign` - Enables indel realignment when running alignment steps. May be helpful for certain callers (VarScan, VarDict)

`--callers`- Comma separated argument for callers, the default is to use all available.  
Example: `--callers mutect2,octopus`

`--cnvcallers`- - Comma separated argument for cnvcallers. Adding flag allows only certain callers to run.  
Example: `--cnvcallers purple`


## Running LOGAN
Example of Tumor_Normal calling mode 
```bash
# copy the logan config files to your current directory
logan init
# preview the logan jobs that will run 
logan run --mode local -profile ci_stub --genome hg38 --sample_sheet samplesheet.tsv --outdir out --fastq_input "*R{1,2}.fastq.gz" -preview --vc --sv --cnv
# run a stub/dryrun of the logan jobs 
logan run --mode local -profile ci_stub --genome hg38 --sample_sheet samplesheet.tsv --outdir out --fastq_input "*R{1,2}.fastq.gz" -stub --vc --sv --cnv
# launch a logan run on slurm with the test dataset
logan run --mode slurm -profile biowulf,slurm --genome hg38 --sample_sheet samplesheet.tsv --outdir out --fastq_input "*R{1,2}.fastq.gz" --vc --sv --cnv 
```

Example of Tumor only calling mode 
```bash
# copy the logan config files to your current directory
logan init
# preview the logan jobs that will run 
logan run --mode local -profile ci_stub --genome hg38 --outdir out --fastq_input "*R{1,2}.fastq.gz" --callers octopus,mutect2 -preview --vc --sv --cnv
# run a stub/dryrun of the logan jobs 
logan run --mode local -profile ci_stub --genome hg38 --outdir out --fastq_input "*R{1,2}.fastq.gz" --callers octopus,mutect2 -stub --vc --sv --cnv
# launch a logan run on slurm with the test dataset
logan run --mode slurm -profile biowulf,slurm --genome hg38 --outdir out --fastq_input "*R{1,2}.fastq.gz" --callers octopus,mutect2 --vc --sv --cnv
```

We currently support the hg38, hg19 (in progress), and mm10 genomes. 




## Contribute 
This site is a living document, created for and by members like you. LOGAN is maintained by the members of CCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [repository](https://github.com/ccbr/LOGAN/pulls).


## References
This repo was originally generated from the [CCBR Nextflow Template](https://github.com/CCBR/CCBR_NextflowTemplate).

<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
