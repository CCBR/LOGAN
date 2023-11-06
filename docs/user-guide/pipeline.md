# How to run WGS-Seek

## Guide

* `./logan` - Starts a next nextflow run
Supports runs from fastq or bam files and either Tumor-Normal or Tumor-only Sequencing

Multiple options required for running
## Code
`./logan  --fastq "Samples/Sample_R{1,2}.fastq.gz" --output 'WGSRun' --sample_sheet sample.tsv --profile biowulf --submit`


### Arguments
Input selection can either be  
`--fastq`
1) A wildcard expansion of Fastq files
 "Samples/Sample_*_R{1,2}.fastq.gz" which finds all Samples in the directory with the head Sample_  
OR  

`--file_input`  
2a) A tab separated file with 3 columns Sample Name, fastq R1 Full path, fastq R2 Full Path if using fastq files   (no header required)

> 3_LM8_p5_S3 Sample_3_LM8_p5/3_LM8_p5_S3_R1_001.fastq.gz   Sample_3_LM8_p5/3_LM8_p5_S3_R2_001.fastq.gz  
6_Spleen_S6     Sample_6_Spleen/6_Spleen_S6_R1_001.fastq.gz   Sample_6_Spleen/6_Spleen_S6_R2_001.fastq.gz

2b) A tab separated file with 2 columns Sample Name, BAM file path  (no header required)

> 3_LM8_p5_S3     3_LM8_p5_S3.bam  
6_Spleen_S6     6_Spleen_S6.bam

`--genome`  
Currently supported genomes are hg38 and mm10

`--output`   
Output Directory

### Modes
`--vc`  
Perform Variant Calling using Mutect2, Vardict, Varscan, Octopus, Strelka (TN Pairs only)

`--sv`  
Perform SV Calling using manta, SvAba, Delly

`--cnv`  
Perform CNV Calling using hg38-Purple or mm10- FREEC, Sequenza


### Running Tumor-Normal
`--sample_sheet`  
Tab separated file for Normal and Tumor delination with a header for "Normal" and "Tumor"
Required for Tumor Normal Pairings  

> Tumor  Normal  
3_LM8_p5_S3  6_Spleen_S6

`--submit`  
Submit job to Biowulf

`--stub`  
Dry/Stub run 

#### Additional Options
`--profile`   
Biowulf or Local Run (Defaults to Biowulf for HPC)

`--resume` Resume previous nextflow run




