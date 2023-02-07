// Split Bed Step to create the path 
process splitinterval {
    
    input:
    path(BED_IN)

    output:
    path('bedout/*.bed')

    script:
    
    SPLIT_BED=file(params.splitbed)

    """
    mkdir bedout
    python ${SPLIT_BED} -infile ${BED_IN} -num 32 -out 'bedout/bed'
    """
}
