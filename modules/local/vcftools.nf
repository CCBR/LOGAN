
process vcftools {
    /*
    Quality-control step to calculates a measure of heterozygosity on
    a per-individual basis. The inbreeding coefficient, F, is estimated
    for each individual using a method of moments. Please see VCFtools
    documentation for more information:
    https://vcftools.github.io/man_latest.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a measure of heterozygosity
    */
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple path(germlinevcf),path(germlinetbi)

    output:
       path("variants_raw_variants.het")


    script:
    """
    vcftools --gzvcf ${germlinevcf} --het --out variants_raw_variants
    """

    stub:
    """
    touch variants_raw_variants.het
    """
}