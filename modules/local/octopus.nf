/* TODO do not set global variables like this.
    These files should be passed as inputs to the processes that use them
*/
//References
GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEFAI=file(params.genomes[params.genome].genomefai)

//Octopus
SOMATIC_FOREST=params.genomes[params.genome].octopus_sforest
GERMLINE_FOREST=params.genomes[params.genome].octopus_gforest



process octopus_tn {
    container "${params.containers.octopus}"
    label 'process_somaticcaller_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val("${tumorname}_vs_${normalname}"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf.gz")

    script:
    """
    octopus \\
        -R ${GENOMEREF} -I ${normal} ${tumor} \\
        --normal-sample ${normalname} \\
        -C cancer \\
        --annotations AF AC AD DP SB \\
        -t ${bed} \\
        --threads ${task.cpus} \\
        ${GERMLINE_FOREST} \\
        ${SOMATIC_FOREST} \\
        -B ${task.memory.toGiga()}Gb \\
        -o ${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf.gz
    """

    stub:
    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf.gz"
    """

}



process bcftools_index_octopus {
    tag { vcf.getBaseName() }
    container "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(tumor), path(vcf)

    output:
        tuple val(tumor), path(vcf), path("${vcf}.tbi")

    script:
    """
    bcftools index -t ${vcf}
    """

    stub:
    """
    touch ${vcf}.tbi
    """

}



process octopus_convertvcf {
    container "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(tumor), val(normal),
        val(oct), path(vcf), path(vcfindex)

    output:
        tuple val(tumor), val(normal), path("${tumor}.octopus.norm.vcf.gz"),
        path("${tumor}.octopus.norm.vcf.gz.tbi")


    script:
    """
    zcat ${vcf}  | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > ${tumor}_temp.octopus.norm.vcf
    bgzip ${tumor}_temp.octopus.norm.vcf
    mv ${tumor}_temp.octopus.norm.vcf.gz ${tumor}.octopus.norm.vcf.gz
    bcftools index -t ${tumor}.octopus.norm.vcf.gz -f
    """

    stub:
    """
    touch ${tumor}.octopus.norm.vcf.gz ${tumor}.octopus.norm.vcf.gz.tbi
    """
}




process octopus_tonly {
    container "${params.containers.octopus}"
    label 'process_somaticcaller_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)

    output:
        tuple val(tumorname),
        path("${tumorname}_${bed.simpleName}.tonly.octopus.vcf.gz")

    script:
    """
    octopus -R ${GENOMEREF} \\
        -C cancer \\
        -I ${tumor} \\
        --annotations AF AC AD DP SB \\
        -B 92Gb \\
        -t ${bed} \\
        --threads ${task.cpus} \\
        ${SOMATIC_FOREST} \\
        -o ${tumorname}_${bed.simpleName}.tonly.octopus.vcf.gz
    """

    stub:
    """
    touch ${tumorname}_${bed.simpleName}.tonly.octopus.vcf.gz
    """
}




process octopus_convertvcf_tonly {
    container "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(tumor), val(oct), path(vcf), path(vcfindex)

    output:
        tuple val(tumor), path("${tumor}.octopus_tonly.norm.vcf.gz"),
        path("${tumor}.octopus_tonly.norm.vcf.gz.tbi")


    script:
    """
    zcat ${vcf}  | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > ${tumor}_temp.octopus_tonly.norm.vcf
    bgzip ${tumor}_temp.octopus_tonly.norm.vcf
    mv ${tumor}_temp.octopus_tonly.norm.vcf.gz ${tumor}.octopus_tonly.norm.vcf.gz
    bcftools index -t ${tumor}.octopus_tonly.norm.vcf.gz -f
    """

    stub:
    """
    touch ${tumor}.octopus_tonly.norm.vcf.gz ${tumor}.octopus_tonly.norm.vcf.gz.tbi
    """
}
