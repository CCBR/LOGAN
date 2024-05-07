SPLIT_BED=file(params.splitbed)
SPLIT_REGIONS=params.split_regions


// Split Bed Step to create the path 
process splitinterval {
    //Keep Process Local
    executor="local"
    cpus= '2'
    memory=2.GB

    input:
    path(BED_IN)

    output:
    path('bedout/*.bed')

    script:

    """
    mkdir -p bedout
    python $SPLIT_BED -infile ${BED_IN} -num ${SPLIT_REGIONS} -out 'bedout/bed'
    """
}

/*
Code to convert beds to interval list
awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' genome.fa.fai
bedtools subtract -a GRCh38.primary_assembly.genome.bed -b ../hg38.blacklist.bed > GRCh38.primary_assembly.genome.interval.bed

gatk BedToIntervalList -I GRCh38.primary_assembly.genome.interval.bed -O \ 
GRCh38.primary_assembly.genome.interval_list -SD GRCh38.primary_assembly.genome.dict

#hg19
awk -F '\t' '{printf("%s\t0\t%s\n",$1,$2);}' /data/CCBR_Pipeliner/db/PipeDB/lib/hg19.with_extra.fa.fai >hg19_all.bed
bedtools subtract -a hg19_all.bed -b hg19-blacklist.v2.bed > hg19_noblacklist.bed
bedtools sort -i hg19_noblacklist.bed -chrThenSizeD  >hg19_noblacklistsort.bed
awk '/^chr[0-9,X,Y,M]*\t/ {printf("%s\t%s\t%s\n",$1,$2,$3);}' hg19_noblacklistsort.bed  > hg19_noblacklistsort_vc.bed
*/
