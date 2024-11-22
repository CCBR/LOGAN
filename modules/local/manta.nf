GENOMEREF=file(params.genomes[params.genome].genome)

process manta_somatic {
    container = "${params.containers.logan}"
    label 'process_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname),
        path("${tumor.simpleName}.diplodSV.vcf.gz"),
        path("${tumor.simpleName}.somaticSV.vcf.gz"),
        path("${tumor.simpleName}.candidateSV.vcf.gz"),
        path("${tumor.simpleName}.candidateSmallIndels.vcf.gz")

    script:
    """
    mkdir -p wd

    configManta.py \
        --normalBam=${normal} \
        --tumorBam=${tumor} \
        --referenceFasta=$GENOMEREF \
        --runDir=wd

    wd/runWorkflow.py -m local -j $task.cpus

    mv wd/results/variants/diploidSV.vcf.gz ${tumor.simpleName}.diplodSV.vcf.gz
    mv wd/results/variants/somaticSV.vcf.gz ${tumor.simpleName}.somaticSV.vcf.gz
    mv wd/results/variants/candidateSV.vcf.gz ${tumor.simpleName}.candidateSV.vcf.gz
    mv wd/results/variants/candidateSmallIndels.vcf.gz ${tumor.simpleName}.candidateSmallIndels.vcf.gz

    """

    stub:

    """
    touch ${tumor.simpleName}.diplodSV.vcf.gz
    touch ${tumor.simpleName}.somaticSV.vcf.gz
    touch ${tumor.simpleName}.candidateSV.vcf.gz
    touch ${tumor.simpleName}.candidateSmallIndels.vcf.gz
    """
}





process manta_tonly {
    container = "${params.containers.logan}"
    label 'process_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        tuple val(tumorname),
        path("${tumor.simpleName}.candidateSV.vcf.gz"),
        path("${tumor.simpleName}.candidateSmallIndels.vcf.gz"),
        path("${tumor.simpleName}.tumorSV.vcf.gz")


    script:
    """
    mkdir -p wd

    configManta.py \
        --tumorBam=${tumor} \
        --referenceFasta=$GENOMEREF \
        --runDir=wd

    wd/runWorkflow.py -m local -j $task.cpus

    mv wd/results/variants/candidateSV.vcf.gz ${tumor.simpleName}.candidateSV.vcf.gz
    mv wd/results/variants/candidateSmallIndels.vcf.gz ${tumor.simpleName}.candidateSmallIndels.vcf.gz
    mv wd/results/variants/tumorSV.vcf.gz ${tumor.simpleName}.tumorSV.vcf.gz

    """

    stub:

    """
    touch ${tumor.simpleName}.candidateSV.vcf.gz
    touch ${tumor.simpleName}.candidateSmallIndels.vcf.gz
    touch ${tumor.simpleName}.tumorSV.vcf.gz

    """
}
