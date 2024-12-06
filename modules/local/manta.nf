GENOMEREF=file(params.genomes[params.genome].genome)

process manta_somatic {
    container = "${params.containers.logan}"
    label 'process_high'

    input:
        tuple val(tumorname), path(tumorbam), path(tumorbai), 
        val(normalname), path(normalbam), path(normalbai)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}.diplodSV.vcf.gz"), path("${tumorname}_vs_${normalname}.diplodSV.vcf.gz.tbi"),
        path("${tumorname}_vs_${normalname}.somaticSV.vcf.gz"), path("${tumorname}_vs_${normalname}.somaticSV.vcf.gz.tbi"),
        path("${tumorname}_vs_${normalname}.candidateSV.vcf.gz"), path("${tumorname}_vs_${normalname}.candidateSV.vcf.gz.tbi"),
        path("${tumorname}_vs_${normalname}.candidateSmallIndels.vcf.gz"), path("${tumorname}_vs_${normalname}.candidateSmallIndels.vcf.gz.tbi")

    script:
    """
    mkdir -p wd

    configManta.py \
        --normalBam=${normalbam} \
        --tumorBam=${tumorbam} \
        --referenceFasta=$GENOMEREF \
        --runDir=wd

    wd/runWorkflow.py -m local -j $task.cpus

    mv wd/results/variants/diploidSV.vcf.gz ${tumorname}_vs_${normalname}.diplodSV.vcf.gz
    mv wd/results/variants/somaticSV.vcf.gz ${tumorname}_vs_${normalname}.somaticSV.vcf.gz
    mv wd/results/variants/candidateSV.vcf.gz ${tumorname}_vs_${normalname}.candidateSV.vcf.gz
    mv wd/results/variants/candidateSmallIndels.vcf.gz ${tumorname}_vs_${normalname}.candidateSmallIndels.vcf.gz

    bcftools index -t ${tumorname}_vs_${normalname}.diplodSV.vcf.gz
    bcftools index -t ${tumorname}_vs_${normalname}.somaticSV.vcf.gz
    bcftools index -t ${tumorname}_vs_${normalname}.candidateSV.vcf.gz
    bcftools index -t ${tumorname}_vs_${normalname}.candidateSmallIndels.vcf.gz
    """

    stub:

    """
    touch ${tumorname}_vs_${normalname}.diplodSV.vcf.gz ${tumorname}_vs_${normalname}.diplodSV.vcf.gz.tbi
    touch ${tumorname}_vs_${normalname}.somaticSV.vcf.gz ${tumorname}_vs_${normalname}.somaticSV.vcf.gz.tbi
    touch ${tumorname}_vs_${normalname}.candidateSV.vcf.gz ${tumorname}_vs_${normalname}.candidateSV.vcf.gz.tbi
    touch ${tumorname}_vs_${normalname}.candidateSmallIndels.vcf.gz ${tumorname}_vs_${normalname}.candidateSmallIndels.vcf.gz.tbi
    """
}





process manta_tonly {
    container = "${params.containers.logan}"
    label 'process_high'

    input:
        tuple val(tumorname), path(tumorbam), path(tumorbai)

    output:
        tuple val(tumorname),
        path("${tumorname}.candidateSV.vcf.gz"), path("${tumorname}.candidateSV.vcf.gz.tbi"),
        path("${tumorname}.candidateSmallIndels.vcf.gz"), path("${tumorname}.candidateSmallIndels.vcf.gz.tbi"),
        path("${tumorname}.tumorSV.vcf.gz"), path("${tumorname}.tumorSV.vcf.gz.tbi")


    script:
    """
    mkdir -p wd

    configManta.py \
        --tumorBam=${tumorbam} \
        --referenceFasta=$GENOMEREF \
        --runDir=wd

    wd/runWorkflow.py -m local -j $task.cpus

    mv wd/results/variants/candidateSV.vcf.gz ${tumorname}.candidateSV.vcf.gz
    mv wd/results/variants/candidateSmallIndels.vcf.gz ${tumorname}.candidateSmallIndels.vcf.gz
    mv wd/results/variants/tumorSV.vcf.gz ${tumorname}.tumorSV.vcf.gz

    bcftools index -t ${tumorname}.candidateSV.vcf.gz
    bcftools index -t ${tumorname}.candidateSmallIndels.vcf.gz
    bcftools index -t ${tumorname}.tumorSV.vcf.gz

    """

    stub:

    """
    touch ${tumorname}.candidateSV.vcf.gz ${tumorname}.candidateSV.vcf.gz.tbi
    touch ${tumorname}.candidateSmallIndels.vcf.gz ${tumorname}.candidateSmallIndels.vcf.gz.tbi
    touch ${tumorname}.tumorSV.vcf.gz ${tumorname}.tumorSV.vcf.gz.tbi

    """
}
