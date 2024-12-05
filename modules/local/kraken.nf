BACDB=file(params.genomes[params.genome].KRAKENBACDB)

process kraken {
    /*
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interactive krona report
    */
    container = "${params.containers.loganqc}"
    label 'process_high'

    input:
        tuple val(samplename),
        path(fqs)

    output:
        tuple val(samplename),
        //path("${samplename}.trimmed.kraken_bacteria.out.txt"),
        path("${samplename}.trimmed.kraken_bacteria.taxa.txt"),
        path("${samplename}.trimmed.kraken_bacteria.krona.html")


    script:
    """
    #Setups temporary directory for
    #intermediate files with built-in
    #mechanism for deletion on exit


    # Copy kraken2 db to local node storage to reduce filesystem strain
    cp -rv $BACDB .
    kdb_base=\$(basename $BACDB)

    kraken2 --db $BACDB \
        --threads 16 --report ${samplename}.trimmed.kraken_bacteria.taxa.txt \
        --output - \
        --gzip-compressed \
        --paired ${fqs[0]} ${fqs[1]}
    # Generate Krona Report
    cut -f2,3 ${samplename}.trimmed.kraken_bacteria.taxa.txt | \
        ktImportTaxonomy - -o ${samplename}.trimmed.kraken_bacteria.krona.html
    """

    stub:
    """
    touch  ${samplename}.trimmed.kraken_bacteria.taxa.txt ${samplename}.trimmed.kraken_bacteria.krona.html
    """

}