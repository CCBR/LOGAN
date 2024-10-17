#!/usr/bin/env Rscript

#######
#
#R Script for ASCAT
#
######
library(ASCAT)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
tumor_bam=args[1]
normal_bam=args[2]
tumor_name=gsub(".bam","",basename(tumor_bam))
normal_name=gsub(".bam","",basename(normal_bam))
genome="hg38"


##DETERMINE SEX
system(sprintf('/usr/local/apps/alleleCount/4.2.1/bin/alleleCounter -l /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/ASCAT/G1000_loci/gender_chr.loci -b %s -o %s_temp_gender.out',
normal_bam,normal_name))
s=read.table(sprintf("%s_temp_gender.out",normal_name))
gender=ifelse(sum(s$V7)>5,"XY","XX")

ascat.prepareHTS(
  tumourseqfile = tumor_bam,
  normalseqfile = normal_bam,
  tumourname = tumor_name,
  normalname = normal_name,
  allelecounter_exe = "/usr/local/apps/alleleCount/4.2.1/bin/alleleCounter",
  alleles.prefix = "/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/ASCAT/G1000_alleles/G1000_alleles_hg38_chr",
  loci.prefix = "/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/ASCAT/G1000_loci/G1000_loci_hg38_chr",
  gender = gender,
  genomeVersion = genome,
  nthreads = 8,
  tumourLogR_file = sprintf("%s_LogR.txt",tumor_name),
  tumourBAF_file = sprintf("%s_BAF.txt",tumor_name),
  normalLogR_file = sprintf("%s_LogR.txt",normal_name),
  normalBAF_file = sprintf("%s_BAF.txt",normal_name))

ascat.bc = ascat.loadData(Tumor_LogR_file = sprintf("%s_LogR.txt",tumor_name), 
    Tumor_BAF_file = sprintf("%s_BAF.txt",tumor_name), 
    Germline_LogR_file = sprintf("%s_LogR.txt",normal_name), Germline_BAF_file = sprintf("%s_BAF.txt",normal_name), 
    gender = gender, genomeVersion = genome)

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/ASCAT/GC_G1000/GC_G1000_hg38.txt", 
  replictimingfile = "/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/ASCAT/RT_G1000/RT_G1000_hg38.txt")
  ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = sprintf('%s_vs_%s_ascat.Rdata',tumor_name,normal_name))


#####