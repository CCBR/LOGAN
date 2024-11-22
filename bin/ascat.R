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
tumor_name=args[2]
normal_bam=args[3]
normal_name=args[4]
genome=args[5]
bed=args[6]
#chroms=scan(text=args[4],sep=",",quiet=T)
cpus=as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
cpus=ifelse(is.na(cpus),2,cpus)

##DETERMINE SEX
system(sprintf('alleleCounter -l /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/%s/ASCAT/gender_chr.loci -b %s -c chrX -o %s_temp_gender.out',
  genome,normal_bam,normal_name))
s=read.table(sprintf("%s_temp_gender.out",normal_name))
gender=ifelse(sum(s$V7)>5,"XY","XX")
print(gender)

ascat.prepareHTS(
  tumourseqfile = tumor_bam,
  normalseqfile = normal_bam,
  tumourname = tumor_name,
  normalname = normal_name,
  allelecounter_exe = "alleleCounter",
  alleles.prefix = sprintf("/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/%s/ASCAT/G1000_alleles/G1000_alleles_%s_chr",genome,genome),
  loci.prefix = sprintf("/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/%s/ASCAT/G1000_loci/G1000_loci_%s_chr",genome,genome),
  gender = gender,
  genomeVersion = genome,
  nthreads = cpus,
  tumourLogR_file = sprintf("%s_LogR.txt",tumor_name),
  tumourBAF_file = sprintf("%s_BAF.txt",tumor_name),
  normalLogR_file = sprintf("%s_LogR.txt",normal_name),
  normalBAF_file = sprintf("%s_BAF.txt",normal_name),
  BED_file=bed)

ascat.bc = ascat.loadData(Tumor_LogR_file = sprintf("%s_LogR.txt",tumor_name), 
    Tumor_BAF_file = sprintf("%s_BAF.txt",tumor_name), 
    Germline_LogR_file = sprintf("%s_LogR.txt",normal_name), Germline_BAF_file = sprintf("%s_BAF.txt",normal_name), 
    gender = gender, genomeVersion = genome)

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, 
  GCcontentfile = sprintf("/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/%s/ASCAT/GC_G1000/GC_G1000_%s.txt",genome,genome), 
  replictimingfile = sprintf("/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/%s/ASCAT/RT_G1000/RT_G1000_%s.txt",genome,genome))
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
write.table(QC,sprintf("%s.qc.txt",paste0(tumor_name,"_vs_",normal_name)))
save(ascat.bc, ascat.output, QC, file = sprintf('%s_vs_%s_ascat.Rdata',tumor_name,normal_name))
