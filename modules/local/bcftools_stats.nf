
process bcftools_stats {
    /*
    Quality-control step to collect summary statistics from bcftools stats.
    When bcftools stats is run with one VCF file then stats by non-reference
    allele frequency, depth distribution, stats by quality and per-sample
    counts, singleton statsistics are calculated. Please see bcftools'
    documentation for more information:
    http://samtools.github.io/bcftools/bcftools.html#stats
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Text file containing a collection of summary statistics
    */
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(samplename),  path("${samplename}.gvcf.gz"),path("${samplename}.gvcf.gz.tbi")
    output:
        path("${samplename}.germline.bcftools_stats.txt")

    script:
    """
    bcftools stats ${samplename}.gvcf.gz > ${samplename}.germline.bcftools_stats.txt
    """

    stub:
    """
    touch ${samplename}.germline.bcftools_stats.txt
    """

}