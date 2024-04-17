GENOMEREF=file(params.genomes[params.genome].genome)
ANNOTSVGENOME=file(params.genomes[params.genome].annotsvgenome)
BWAGENOME=file(params.genomes[params.genome].bwagenome)
INDELREF=file(params.genomes[params.genome].INDELREF)



process svaba_somatic {
    container = "${params.containers.logan}"
    label 'process_highcpu'

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



process manta_somatic {
    container = "${params.containers.logan}"
    label 'process_highcpu'

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

    wd/runWorkflow.py -m local -j 10 -g 10

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


process annotsv_tn {
     //AnnotSV for Manta/Svaba works with either vcf.gz or .vcf files
     //Requires bedtools,bcftools
    container = "${params.containers.annotcnvsv}"
    process
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


process manta_tonly {
    container = "${params.containers.logan}"
    label 'process_highcpu'

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

    wd/runWorkflow.py -m local -j 10 -g 10

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



process svaba_tonly {
    container = "${params.containers.logan}"
    label 'process_highcpu'

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


process gunzip {

    input:
        tuple val(tumorname),
        path(vcf), val(sv)

    output:
        tuple val(tumorname),
        path("${tumorname}.tumorSV.vcf"), val(sv)

    script:
    """
    gunzip -f ${vcf} > ${tumorname}.tumorSV.vcf
    """

    stub:

    """
    touch ${tumorname}.tumorSV.vcf
    """

}


process survivor_sv {
    container = "${params.containers.annotcnvsv}"

    input:
        tuple val(tumorname),
        path(vcfs),val(svs)

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
