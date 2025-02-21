
process samtools_flagstats {
    /*
    Quality-control step to assess alignment quality. Flagstat provides
    counts for each of 13 categories based primarily on bit flags in the
    FLAG field. Information on the meaning of the flags is given in the
    SAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Text file containing alignment statistics
    */
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        path("${samplename}.samtools_flagstat.txt")

    script:
    """
    samtools flagstat ${bam} > ${samplename}.samtools_flagstat.txt
    """

    stub:
    """
    touch ${samplename}.samtools_flagstat.txt
    """
}
