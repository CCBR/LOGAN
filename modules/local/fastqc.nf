
process fastqc {
    """
    Quality-control step to assess sequencing quality of each sample.
    FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        FastQC report and zip file containing sequencing quality information
    """
    container = "${params.containers.loganqc}"
    label 'process_medium'

    input:
        tuple val(samplename), path(bam), path(bai)
    output:
        tuple val(samplename), path("${samplename}_fastqc.html"), path("${samplename}_fastqc.zip")

    script:

    """
    mkdir -p fastqc
    fastqc -t $task.cpus \
        -f bam \
        -o fastqc \
        ${bam}
    mv fastqc/${bam.simpleName}_fastqc.html ${samplename}_fastqc.html
    mv fastqc/${bam.simpleName}_fastqc.zip ${samplename}_fastqc.zip
    """

    stub:
    """
    touch ${samplename}_fastqc.html ${samplename}_fastqc.zip
    """
}
