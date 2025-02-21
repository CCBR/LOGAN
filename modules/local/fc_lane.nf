process fc_lane {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(samplename), path(fqs)

    output:
        tuple val(samplename),
        path("${samplename}.fastq.info.txt")

    script:
    GET_FLOWCELL_LANES=file(params.get_flowcell_lanes)

    """
    python $GET_FLOWCELL_LANES \
    ${fqs[0]}  \
    ${samplename} > ${samplename}.fastq.info.txt
    """

    stub:
    """
    touch ${samplename}.fastq.info.txt
    """
}
