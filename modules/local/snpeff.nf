SNPEFF_GENOME = params.genomes[params.genome].snpeff_genome
SNPEFF_CONFIG = file(params.genomes[params.genome].snpeff_config)
SNPEFF_BUNDLE = file(params.genomes[params.genome].snpeff_bundle)

process snpeff {
    /*
    Data processing and quality-control step to annotate variants, predict its
    functional effects, and collect various summary statistics about variants and
    their annotations. Please see SnpEff's documentation for more information:
    https://pcingola.github.io/SnpEff/
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    */
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(samplename), path("${samplename}.gvcf.gz"), path("${samplename}.gvcf.gz.tbi")
    output:
        tuple path("${samplename}.germline.snpeff.ann.vcf"),
        path("${samplename}.germline.snpeff.ann.csv"),
        path("${samplename}.germline.snpeff.ann.html")

    script:
    """
        java -Xmx12g -jar \$SNPEFF_JAR \
        -v -canon -c $SNPEFF_CONFIG \
        -csvstats ${samplename}.germline.snpeff.ann.csv \
        -stats ${samplename}.germline.snpeff.ann.html \
        $SNPEFF_GENOME \
        ${samplename}.gvcf.gz > ${samplename}.germline.snpeff.ann.vcf
    """

    stub:

    """
    touch ${samplename}.germline.snpeff.ann.vcf
    touch ${samplename}.germline.snpeff.ann.csv
    touch ${samplename}.germline.snpeff.ann.html
    """
}
