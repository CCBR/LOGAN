GENOMEREF=file(params.genomes[params.genome].genome)
DBSNP=file(params.genomes[params.genome].dbsnp) //dbsnp_138.hg38.vcf.gz"


process gatk_varianteval {
    /*
    Quality-control step to calculate various quality control metrics from a
    variant callset. These metrics include the number of raw or filtered SNP
    counts; ratio of transition mutations to transversions; concordance of a
    particular sample's calls to a genotyping chip; number of s per sample.
    Please see GATK's documentation for more information:
    https://gatk.broadinstitute.org/hc/en-us/articles/360040507171-VariantEval
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    */
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(samplename), path("${samplename}.gvcf.gz") ,path("${samplename}.gvcf.gz.tbi")
    output:
        path("${samplename}.germline.eval.grp")
    script:
    """
    gatk --java-options '-Xmx12g -XX:ParallelGCThreads=16' VariantEval \
        -R $GENOMEREF \
        -O ${samplename}.germline.eval.grp \
        --dbsnp $DBSNP \
        --eval ${samplename}.gvcf.gz
    """

    stub:

    """
    touch ${samplename}.germline.eval.grp
    """

}

process collectvariantcallmetrics {
    /*
    Quality-control step to collect summary metrics about snps and indels
    called in a multisample VCF file. Please see the Broad's documentation
    for more information about each field in the generated log file:
    https://broadinstitute.github.io/picard/picard-metric-definitions.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a collection of metrics relating to snps and indels
    */
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple path(germlinevcf),path(germlinetbi)

    output:
        tuple path("raw_variants.variant_calling_detail_metrics"),
        path("raw_variants.variant_calling_summary_metrics")


    script:
    """
    java -Xmx24g -jar \${PICARDJARPATH}/picard.jar \
        CollectVariantCallingMetrics \
        INPUT=${germlinevcf} \
        OUTPUT= "raw_variants" \
        DBSNP=$DBSNP Validation_Stringency=SILENT
    """

    stub:
    """
    touch raw_variants.variant_calling_detail_metrics raw_variants.variant_calling_summary_metrics
    """

}
