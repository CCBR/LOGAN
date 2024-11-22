
BWAGENOME=file(params.genomes[params.genome].bwagenome)
INDELREF=file(params.genomes[params.genome].INDELREF)


process svaba_somatic {
    container = "${params.containers.logan}"
    label 'process_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname),
        path("${tumor.simpleName}.bps.txt.gz"),
        path("${tumor.simpleName}.contigs.bam"),
        path("${tumor.simpleName}.discordant.txt.gz"),
        path("${tumor.simpleName}.alignments.txt.gz"),
        path("${tumor.simpleName}.svaba.germline.indel.vcf"),
        path("${tumor.simpleName}.svaba.germline.sv.vcf"),
        path("${tumor.simpleName}.svaba.somatic.indel.vcf"),
        path("${tumor.simpleName}.svaba.somatic.sv.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.germline.indel.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.germline.sv.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.somatic.indel.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.somatic.sv.vcf"),
        path("${tumor.simpleName}.log")


    script:
    """
    svaba run -t ${tumor} -n ${normal} -p $task.cpus -D $INDELREF -a ${tumor.simpleName} -G $BWAGENOME
    """

    stub:

    """
    touch "${tumor.simpleName}.bps.txt.gz"
    touch "${tumor.simpleName}.contigs.bam"
    touch "${tumor.simpleName}.discordant.txt.gz"
    touch "${tumor.simpleName}.alignments.txt.gz"
    touch "${tumor.simpleName}.svaba.germline.indel.vcf"
    touch "${tumor.simpleName}.svaba.germline.sv.vcf"
    touch "${tumor.simpleName}.svaba.somatic.indel.vcf"
    touch "${tumor.simpleName}.svaba.somatic.sv.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.germline.indel.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.germline.sv.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.somatic.indel.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.somatic.sv.vcf"
    touch "${tumor.simpleName}.log"

    """
}



process svaba_tonly {
    container = "${params.containers.logan}"
    label 'process_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        tuple val(tumorname),
        path("${tumor.simpleName}.bps.txt.gz"),
        path("${tumor.simpleName}.contigs.bam"),
        path("${tumor.simpleName}.discordant.txt.gz"),
        path("${tumor.simpleName}.alignments.txt.gz"),
        path("${tumor.simpleName}.svaba.indel.vcf"),
        path("${tumor.simpleName}.svaba.sv.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.indel.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.sv.vcf"),
        path("${tumor.simpleName}.log")


    script:
    """
    svaba run -t ${tumor} -p $task.cpus -D $INDELREF -a ${tumor.simpleName} -G $BWAGENOME
    """

    stub:

    """
    touch "${tumor.simpleName}.bps.txt.gz"
    touch "${tumor.simpleName}.contigs.bam"
    touch "${tumor.simpleName}.discordant.txt.gz"
    touch "${tumor.simpleName}.alignments.txt.gz"
    touch "${tumor.simpleName}.svaba.indel.vcf"
    touch "${tumor.simpleName}.svaba.sv.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.indel.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.sv.vcf"
    touch "${tumor.simpleName}.log"

    """
}
