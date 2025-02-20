GENOMEREF=file(params.genomes[params.genome].genome)
ANNOTSVGENOME=params.genomes[params.genome].annotsvgenome


process annotsv_tn {
     //AnnotSV for Manta/Svaba/GRIDSS works with either vcf.gz or .vcf files
     //Requires bedtools,bcftools
    errorStrategy = 'ignore'
    container = "${params.containers.annotcnvsv}"

    input:
        tuple val(tumorname), path(somaticvcf), val(sv)

    output:
        tuple val(tumorname),
        path("${sv}/${tumorname}.tsv"),
        path("${sv}/${tumorname}.unannotated.tsv")


    script:
    """
    mkdir ${sv}

    AnnotSV -SVinputFile ${somaticvcf} \
    -genomeBuild $ANNOTSVGENOME \
    -SVinputInfo 1 -outputFile ${tumorname} \
    -outputDir ${sv}

    """

    stub:
    """
    mkdir ${sv}

    touch "${sv}/${tumorname}.tsv"
    touch "${sv}/${tumorname}.unannotated.tsv"
    """
}




process gunzip {
    label 'process_single'

    input:
        tuple val(tumorname),
        path(vcf), val(sv)

    output:
        tuple val(tumorname),
        path("${tumorname}.tumorSV_${sv}.vcf"), val(sv)

    script:
    """
    gunzip -f ${vcf} > ${tumorname}.tumorSV_${sv}.vcf
    """

    stub:

    """
    touch ${tumorname}.tumorSV_${sv}.vcf
    """

}


process survivor_sv {
    container = "${params.containers.annotcnvsv}"

    input:
        tuple val(tumorname),
        path(vcfs), val(svs)

    output:
        tuple val(tumorname),
        path("${tumorname}_merged.vcf"),
        val("survivor")


    script:
    strin = vcfs.join("\\n")

    """
    echo -e '$strin' > filelistin
    SURVIVOR merge filelistin 1000 2 1 1 1 30 ${tumorname}_merged.vcf
    """

    stub:
    strin = vcfs.join("\\n")
    """
    echo -e '$strin' > filelistin
    touch "${tumorname}_merged.vcf"
    """
}


process annotsv_tonly {
     //AnnotSV for Manta/Svaba works with either vcf.gz or .vcf files
     //Requires bedtools,bcftools
    errorStrategy = 'ignore'

    container = "${params.containers.annotcnvsv}"

    input:
        tuple val(tumorname), path(somaticvcf), val(sv)

    output:
        tuple val(tumorname),
        path("${sv}/${tumorname}.tsv"),
        path("${sv}/${tumorname}.unannotated.tsv")


    script:
    """
    mkdir ${sv}

    AnnotSV -SVinputFile ${somaticvcf} \
    -genomeBuild $ANNOTSVGENOME \
    -SVinputInfo 1 -outputFile ${tumorname} \
    -outputDir ${sv}

    """

    stub:
    """
    mkdir ${sv}

    touch "${sv}/${tumorname}.tsv"
    touch "${sv}/${tumorname}.unannotated.tsv"
    """
}
