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

genome="hg38"

ascat.prepareHTS(
  tumourseqfile = tumor_bam,
  normalseqfile = normal_bam,
  tumourname = tumor_name,
  normalname = normal_name,
  allelecounter_exe = "/PATH/TO/allelecounter",
  alleles.prefix = "/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/ASCAT/G1000_alleles",
  loci.prefix = "/data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/ASCAT_G1000_loci",
  nthreads = 10
  gender = "XX",
  genomeVersion = genome,
  nthreads = 8,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt")

ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", 
    Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = 'XX', genomeVersion = "hg19")

ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_file.txt", replictimingfile = "RT_file.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')


#####