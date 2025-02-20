
process mosdepth {
    /*
    Quality-control step to assess depth
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        `{prefix}.mosdepth.global.dist.txt`
        `{prefix}.mosdepth.summary.txt`
        `{prefix}.mosdepth.region.dist.txt` (if --by is specified)
        `{prefix}.per-base.bed.gz|per-base.d4` (unless -n/--no-per-base is specified)
        `{prefix}.regions.bed.gz` (if --by is specified)
        `{prefix}.quantized.bed.gz` (if --quantize is specified)
        `{prefix}.thresholds.bed.gz` (if --thresholds is specified)
    */
    container = "${params.containers.loganqc}"
    label 'process_medium'

    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        tuple path("${samplename}.mosdepth.region.dist.txt"),
        path("${samplename}.mosdepth.summary.txt"),
        path("${samplename}.regions.bed.gz"),
        path("${samplename}.regions.bed.gz.csi")


    script:
    """
    mosdepth -n --fast-mode --by 500  ${samplename} ${bam} -t $task.cpus
    """

    stub:
    """
    touch "${samplename}.mosdepth.region.dist.txt"
    touch "${samplename}.mosdepth.summary.txt"
    touch "${samplename}.regions.bed.gz"
    touch "${samplename}.regions.bed.gz.csi"
    """
}


process mosdepth_exome {
    container = "${params.containers.loganqc}"
    label 'process_medium'

    input:
        tuple val(samplename), path(bam), path(bai), path(bed)

    output:
        tuple path("${samplename}.mosdepth.global.dist.txt"),
        path("${samplename}.mosdepth.region.dist.txt"),
        path("${samplename}.mosdepth.summary.txt"),
        path("${samplename}.per-base.bed.gz"),
        path("${samplename}.per-base.bed.gz.csi"),
        path("${samplename}.regions.bed.gz"),
        path("${samplename}.regions.bed.gz.csi")

    script:
    """
    mosdepth --by ${bed} $samplename $bam -t $task.cpus
    """

    stub:
    """
    touch "${samplename}.mosdepth.global.dist.txt"
    touch "${samplename}.mosdepth.region.dist.txt"
    touch "${samplename}.mosdepth.summary.txt"
    touch "${samplename}.per-base.bed.gz"
    touch "${samplename}.per-base.bed.gz.csi"
    touch "${samplename}.regions.bed.gz"
    touch "${samplename}.regions.bed.gz.csi"
    """
}
