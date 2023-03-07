# How to run WGS-Seek

## Guide

* `./run_nf_fq` - Starts a next nextflow run
Supports runs from Fastq and either Tumor-Normal or Tumor-only Sequencing

## Running Nextflow
Multiple options required for running

## Code
`./run_nf_fq  --fqinput "Samples/Sample_R{1,2}.fastq.gz" --output 'B2' --sample_sheet sample.tsv --paired T --profile`

`./run_nf_fq  --fileinput samples.tsv --output 'B2' --sample_sheet sample.tsv --paired`

### Arguments
Input selection can either be  
`--fqinput`
1) A wildcard expansion of Fastq files
 "Samples/Sample_*_R{1,2}.fastq.gz" which finds all Samples in the directory with the head Sample_  
OR  
`--fileinput`
2) A tab separated file with 3 columns Sample Name, Fastq1 Full path, Fastq2 Full Path

`--output` - Output Directory

`--sample_sheet`- Tab separated file for Tumor Normal delination

`--profile` Biowulf or Local Run

`--resume` Resume previous nextflow run

`--submit`- Submit job to Biowulf?

`--paired`- Are Samples paired Tumor-Normal


