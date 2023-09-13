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
