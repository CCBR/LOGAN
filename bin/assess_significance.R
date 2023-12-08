#!/usr/bin/env Rscript

library(rtracklayer)

args <- commandArgs()

dataTable <-read.table(args[5], header=TRUE);
ratio<-data.frame(dataTable)

dataTable <- read.table(args[4], header=FALSE)
cnvs<- data.frame(dataTable)

ratio$Ratio[which(ratio$Ratio==-1)]=NA

cnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))  
ratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio)

overlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)
normals <- setdiff(ratio.bed,cnvs.bed)
normals <- subsetByOverlaps(ratio.bed,normals)

#mu <- mean(score(normals),na.rm=TRUE)
#sigma<- sd(score(normals),na.rm=TRUE)

#hist(score(normals),n=500,xlim=c(0,2))
#hist(log(score(normals)),n=500,xlim=c(-1,1))

#shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])
#qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))
#qqline(score(normals)[which(!is.na(score(normals)))], col = 2)

#shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])
#qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))
#qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)

numberOfCol=length(cnvs)
wscore=c()
kscore=c()
for (i in c(1:length(cnvs[,1]))) {
values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))
resultw <- class(try(wilcox.test(values,score(normals)), silent = TRUE))
ifelse(resultw == "try-error", wscore <- c(wscore, "NA"), wscore <- c(wscore, wilcox.test(values,score(normals))$p.value))
resultks <- class(try(ks.test(values,score(normals)), silent = TRUE))
ifelse(resultks == "try-error",kscore <- c(kscore, "NA"),kscore <- c(kscore, ks.test(values,score(normals))$p.value))
}
cnvs = cbind(cnvs, as.numeric(wscore), as.numeric(kscore))

if (numberOfCol==7) {
  names(cnvs)=c("chr","start","end","copy number","status","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
if (numberOfCol==9) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
if (numberOfCol==11) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","somatic/germline","precentageOfGermline","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
write.table(cnvs, file=paste(args[4],".p.value.txt",sep=""),sep="\t",quote=F,row.names=F)


