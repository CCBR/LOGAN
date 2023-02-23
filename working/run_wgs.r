##Create PON
s="4_CL0157_N1D_E
5_GUHA_MC_0292_001
6_PERU_MC_0013_001
10_GUHA_MC_0310_0040
3_GUHA_MC_0650_1
1_BR7_Normal_DNA
2_BR10_Normal_DNA
4_GUHA_MC_0566_001
9_NCI0422_N1D_E
21_GUHA_MC_0497_002
18_GUHA_MC_0440_007
5_AU_19_117_M
6_AN_16_35_H
1_GUHA_MC_0714_001_N1D"

s1=unlist(strsplit(s,"\n"))

alldir=Sys.glob("../SCLC_BrainMets/BMSCLC*/*")
#basename(alldir) %in% s1

norm=sapply(s1,function(x)grep(x,basename(alldir)))
normdir=alldir[norm]
finame=Sys.glob(paste0(normdir,"/*R1_001.fastq.gz"))
finame=gsub("_R1_001.fastq.gz","",basename(finame))
#Sys.glob("../SCLC_BrainMets/BMSCLC01/Sample_4_CL0157_N1D_E/*.gz")
#lapply(paste0("ls ",normdir,"/*gz"),system)

s2=data.frame(Tumor=finame,Normal="")
write.table(s2,"norm_sampsheet.tsv",quote=F,row.names=F,col.names=T)

lapply(paste0("ls ",normdir,"/*gz"),system)
paste0("../SCLC_BrainMets/*/{",paste(basename(normdir),collapse=","),"}/*_R{1,2}_001.fastq.gz")
head -n3 norm_sampsheet.tsv >test.tsv
module load nextflow singularity

nextflow run /data/nousomedr/wgs/wgs-seek/wgs-seek_tumoronly.nf  -profile biowulf \
--input "../SCLC_BrainMets/*/{Sample_4_CL0157_N1D_E,Sample_5_GUHA_MC_0292_001,Sample_6_PERU_MC_0013_001,Sample_10_GUHA_MC_0310_0040,Sample_3_GUHA_MC_0650_1,Sample_1_BR7_Normal_DNA,Sample_2_BR10_Normal_DNA,Sample_4_GUHA_MC_0566_001,Sample_9_NCI0422_N1D_E,Sample_21_GUHA_MC_0497_002,Sample_18_GUHA_MC_0440_007,Sample_5_AU_19_117_M,Sample_6_AN_16_35_H,Sample_1_GUHA_MC_0714_001_N1D}/*_R{1,2}_001.fastq.gz" \
--sample_sheet norm_sampsheet.tsv --output . -resume





'../SCLC_BrainMets/BMSCLC01/Sample_4_CL0157_N1D_E/*gz' '../SCLC_BrainMets/BMSCLC02/Sample_5_GUHA_MC_0292_001/*gz' ,'../SCLC_BrainMets/BMSCLC04/Sample_6_PERU_MC_0013_001/*gz',\
'../SCLC_BrainMets/BMSCLC05/Sample_10_GUHA_MC_0310_0040/*gz','../SCLC_BrainMets/BMSCLC06/Sample_3_GUHA_MC_0650_1/*gz','../SCLC_BrainMets/BMSCLC08/Sample_1_BR7_Normal_DNA/*gz',\
'../SCLC_BrainMets/BMSCLC09/Sample_2_BR10_Normal_DNA/*gz','../SCLC_BrainMets/BMSCLC32/Sample_4_GUHA_MC_0566_001/*gz','../SCLC_BrainMets/BMSCLC33/Sample_9_NCI0422_N1D_E/*gz',\
'../SCLC_BrainMets/BMSCLC34/Sample_21_GUHA_MC_0497_002/*gz','../SCLC_BrainMets/BMSCLC35/Sample_18_GUHA_MC_0440_007/*gz','../SCLC_BrainMets/BMSCLC36/Sample_5_AU_19_117_M/*gz',\
'../SCLC_BrainMets/BMSCLC41/Sample_6_AN_16_35_H/*gz','../SCLC_BrainMets/BMSCLC55/Sample_1_GUHA_MC_0714_001_N1D/*gz