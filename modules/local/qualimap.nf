
process qualimap_bamqc {
    /*
    Quality-control step to assess various post-alignment metrics
    and a secondary method to calculate insert size. Please see
    QualiMap's website for more information about BAM QC:
    http://qualimap.conesalab.org/
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Report containing post-aligment quality-control metrics
    */
    container = "${params.containers.loganqc}"
    label 'process_medium'

    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        tuple path("${samplename}_genome_results.txt"), path("${samplename}_qualimapReport.html")

    script:
    """
    unset DISPLAY
    qualimap bamqc -bam ${bam} \
        --java-mem-size=112G \
        -c -ip \
        -outdir ${samplename} \
        -outformat HTML \
        -nt 8 \
        --gd HUMAN \
        --skip-duplicated \
        -nw 500 \
        -p NON-STRAND-SPECIFIC
    mv ${samplename}/genome_results.txt ${samplename}_genome_results.txt
    mv ${samplename}/qualimapReport.html ${samplename}_qualimapReport.html
    """

    stub:
    """
    touch ${samplename}_genome_results.txt ${samplename}_qualimapReport.html
    """
}