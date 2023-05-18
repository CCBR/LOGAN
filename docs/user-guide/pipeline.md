# How to run WGS-Seek

## Guide

* `./wgs-seek` - Starts a next nextflow run
Supports runs from Fastq and either Tumor-Normal or Tumor-only Sequencing

## Running Nextflow
Multiple options required for running

## Code
`./wgs-seek  --fastq "Samples/Sample_R{1,2}.fastq.gz" --output 'B2' --sample_sheet sample.tsv --paired T --profile biowulf`


### Arguments
Input selection can either be  
`--fastq`
1) A wildcard expansion of Fastq files
 "Samples/Sample_*_R{1,2}.fastq.gz" which finds all Samples in the directory with the head Sample_  
OR  
`--filelist`
2a) A tab separated file with 3 columns Sample Name, Fastq1 Full path, Fastq2 Full Path if using fastq files or
2b) A tab separated file with 2 columns Sample Name, BAM file path

`--output` - Output Directory

`--sample_sheet`- Tab separated file for Normal and Tumor delination with a header for "Normal" and "Tumor"

`--profile` Biowulf or Local Run

`--resume` Resume previous nextflow run

`--submit`- Submit job to Biowulf?

`--paired`- Are Samples paired Tumor-Normal


