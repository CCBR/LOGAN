# How to run LOGAN 

## Guide

### Preview Run
```bash
git clone https://github.com/CCBR/LOGAN
module load nextflow
# Starts a next nextflow preview run to see the processes that will run
nextflow run LOGAN/main.nf -profile ci_stub -preview`  
```


##Example run 
## Usage

### Input Files
LOGAN supports inputs of either 
1) paired end fastq files

`--fastq_input`- A glob can be used to include all FASTQ files. Like `--fastq_input "*R{1,2}.fastq.gz"` quotes.

2) Pre aligned BAM files with BAI indices 

`--bam_input`- A glob can be used to include all FASTQ files. Like `--bam_input *.bam`

3) A sheet that indicates the sample name and either FASTQs or BAM file locations

`--fastq_file_input`-  A headerless tab delimited sheet that has the sample name, R1, and R2 file locations

Example
```bash
c130863309_TUMOR   /data/nousomedr/c130863309_TUMOR.R1_001.fastq.gz  /data/nousomedr/c130863309_TUMOR.R2_001.fastq.gz
c130889189_PBMC  /data/nousomedr/c130889189_PBMC.R1_001.fastq.gz  /data/nousomedr/c130889189_PBMC.R2_001.fastq.gz
```


`--bam_file_input` -  A headerless Tab delimited sheet that has the sample name, bam, and bam index (bai) file locations

Example
```bash
c130863309_TUMOR   /data/nousomedr/c130863309_TUMOR.bam  /data/nousomedr/c130863309_TUMOR.bam.bai
c130889189_PBMC  /data/nousomedr/c130889189_PBMC.bam  /data/nousomedr/c130889189_PBMC.bam.bai
```

### Genome
`--genome` - A flag to indicate which genome to run for alignment/variant calling/etc. Like `--genome hg38` to run the hg38 genome

`--genome hg19` and `--genome mm10` are also supported 

#### hg38 has options for either  
`--genome hg38` - Based off the GRCh38.d1.vd1.fa which is consistent with TCGA and other GDC processing pipelines  

`--genome hg38_sf` - Based off the Homo_sapiens_assembly38.fasta which is derived from the Broad Institute/NCI Sequencing Facility.  
The biggest difference between the two is that GRCh38.d1.vd1.fa has fewer contigs (especially related to HLA regions), so reads should map to chr6 vs the HLA contig directly


### Operating Modes

#### 1.  Paired Tumor/Normal Mode 

Required for Paired Tumor/Normal Mode

`--sample_sheet` In Paired mode a sample sheet must be provided with the basename of the Tumor and Normal samples. This sheet must be Tab separated with a header for Tumor and Normal.  

Example
```bash
Tumor  Normal
c130863309_TUMOR  c130863309_PBMC
c130889189_TUMOR  c130889189_PBMC
```

#### 2.  Tumor only mode

No addtional flags for sample sheet are required as all samples will be used to call variants

#### Calling Mode

Adding flags determines SNV (germline and/or somatic), SV, and/or CNV calling modes

`--vc`- Enables somatic SNV calling using mutect2, vardict, varscan, octopus, strelka (TN only), MUSE (TN only), and lofreq (TN only)

`--germline`- Enables germline calling using Deepvariant

`--sv`- Enables somatic SV calling using Manta, GRIDSS, and SVABA

`--cnv`- Enables somatic CNV calling using FREEC, Sequenza, ASCAT, CNVKit, and Purple (hg19/hg38 only)



#### Optional Arguments
`--indelrealign` - Enables indel realignment when running alignment steps. May be helpful for certain callers (VarScan, VarDict)

`--callers`- Comma separated argument for selecting only specified callers, the default is to use all available.  
Example: `--callers mutect2,octopus`

`--cnvcallers`- - Comma separated argument for selecting only specified CNV callers. Adding flag allows only certain callers to run.  
Example: `--cnvcallers purple`

`--svcallers`- - Comma separated argument for selecting only specified SV vallers. Adding flag allows only certain callers to run.  
Example: `--svcallers gridss`

## Running LOGAN
Example of Tumor_Normal calling mode 
```bash
# preview the logan jobs that will run 
nextflow run LOGAN/main.nf --mode local -profile ci_stub --genome hg38 --sample_sheet samplesheet.tsv --outdir out --fastq_input "*R{1,2}.fastq.gz" -preview --vc --sv --cnv
# run a stub/dryrun of the logan jobs 
nextflow run LOGAN/main.nf --mode local -profile ci_stub --genome hg38 --sample_sheet samplesheet.tsv --outdir out --fastq_input "*R{1,2}.fastq.gz" -stub --vc --sv --cnv
# launch a logan run on slurm with the test dataset
nextflow run LOGAN/main.nf --mode slurm -profile biowulf,slurm --genome hg38 --sample_sheet samplesheet.tsv --outdir out --fastq_input "*R{1,2}.fastq.gz" --vc --sv --cnv 
```

Example of Tumor only calling mode 
```bash
# preview the logan jobs that will run 
nextflow run LOGAN/main.nf --mode local -profile ci_stub --genome hg38 --outdir out --fastq_input "*R{1,2}.fastq.gz" --callers octopus,mutect2 -preview --vc --sv --cnv
# run a stub/dryrun of the logan jobs 
nextflow run LOGAN/main.nf --mode local -profile ci_stub --genome hg38 --outdir out --fastq_input "*R{1,2}.fastq.gz" --callers octopus,mutect2 -stub --vc --sv --cnv
# launch a logan run on slurm with the test dataset
nextflow run LOGAN/main.nf --mode slurm -profile biowulf,slurm --genome hg38 --outdir out --fastq_input "*R{1,2}.fastq.gz" --callers octopus,mutect2 --vc --sv --cnv
```

