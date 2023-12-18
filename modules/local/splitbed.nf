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
*/
