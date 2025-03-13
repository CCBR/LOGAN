FASTQ_SCREEN_CONF=file(params.fastq_screen_conf)

process fastq_screen {
    //Uses Trimmed Files
    container = "${params.containers.loganqc}"
    label 'process_medium'

    input:
    tuple val(samplename),
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz")

    output:
        tuple path("${samplename}.R1.trimmed_screen.html"),
        path("${samplename}.R1.trimmed_screen.png"),
        path("${samplename}.R1.trimmed_screen.txt"),
        path("${samplename}.R2.trimmed_screen.html"),
        path("${samplename}.R2.trimmed_screen.png"),
        path("${samplename}.R2.trimmed_screen.txt")

    script:
        FASTQ_SCREEN_CONF=file(params.fastq_screen_conf)

        """
    fastq_screen --conf $FASTQ_SCREEN_CONF \
        --outdir . \
        --threads 8 \
        --subset 1000000 \
        --aligner bowtie2 \
        --force \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz

        """

    stub:
    """
    touch ${samplename}.R1.trimmed_screen.html ${samplename}.R1.trimmed_screen.png
    touch ${samplename}.R1.trimmed_screen.txt ${samplename}.R2.trimmed_screen.html
    touch ${samplename}.R2.trimmed_screen.png ${samplename}.R2.trimmed_screen.txt
    """
}