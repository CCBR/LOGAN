GENOMEREF=file(params.genomes[params.genome].genome)
ANNOTSVGENOME=params.genomes[params.genome].annotsvgenome
BWAGENOME=file(params.genomes[params.genome].bwagenome)
INDELREF=file(params.genomes[params.genome].INDELREF)
BLACKLIST=file(params.genomes[params.genome].GRIDSSBLACKLIST)

if (params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
    GENOMEVER = params.genomes[params.genome].GENOMEVER
    PONSGL = file(params.genomes[params.genome].PONSGL)
    PONSV = file(params.genomes[params.genome].PONSV)
    SVHOTSPOT = file(params.genomes[params.genome].SVHOTSPOT)
    REPEATMASK = file(params.genomes[params.genome].REPEATMASK)
}

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


process gridss_somatic {
    label 'process_high'
    container = "${params.containers.sv}"

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), 
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}.vcf.gz"),
        path("${tumorname}_vs_${normalname}.vcf.gz.tbi"),
        path("${tumorname}_vs_${normalname}.vcf.gz.assembly.bam"),
        path("${tumorname}.gripss.vcf.gz"),
        path("${tumorname}.gripss.vcf.gz.tbi"),
        path("${tumorname}.gripss.filtered.vcf.gz"),
        path("${tumorname}.gripss.filtered.vcf.gz.tbi")

    script:
    """
    gridss --jar /opt2/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    -r $BWAGENOME \
    -l ${normalname},${tumorname} \
    -o ${tumorname}_vs_${normalname}.vcf.gz -b $BLACKLIST \
    --picardoptions VALIDATION_STRINGENCY=LENIENT \
    --jvmheap 90g \
    --otherjvmheap 64g \
    -t $task.cpus \
    ${normal} ${tumor} 

    mkdir -p ${tumorname}_vs_${normalname}

    java -jar /opt2/hmftools/gripss.jar \
        -sample ${tumorname} \
        -reference ${normalname} \
        -ref_genome_version $GENOMEVER \
        -ref_genome $GENOMEREF  \
        -pon_sgl_file $PONSGL \
        -pon_sv_file $PONSV \
        -known_hotspot_file $SVHOTSPOT \
        -repeat_mask_file $REPEATMASK \
        -vcf ${tumorname}_vs_${normalname}.vcf.gz \
        -output_dir ${tumorname}_vs_${normalname}

    mv ${tumorname}_vs_${normalname}/* .
    """


    stub:

    """
    touch "${tumorname}_vs_${normalname}.vcf.gz"
    touch "${tumorname}_vs_${normalname}.vcf.gz.tbi"
    touch "${tumorname}_vs_${normalname}.vcf.gz.assembly.bam"
    touch "${tumorname}.gripss.vcf.gz"
    touch "${tumorname}.gripss.vcf.gz.tbi"
    touch "${tumorname}.gripss.filtered.vcf.gz"
    touch "${tumorname}.gripss.filtered.vcf.gz.tbi"
    """
}



process gridss_tonly {
    label 'process_high'
    container = "${params.containers.sv}"

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        tuple val(tumorname),
        path("${tumorname}.vcf.gz"),
        path("${tumorname}.vcf.gz.tbi"),
        path("${tumorname}.vcf.gz.assembly.bam"),
        path("${tumorname}.gripss.vcf.gz"),
        path("${tumorname}.gripss.vcf.gz.tbi"),
        path("${tumorname}.gripss.filtered.vcf.gz"),
        path("${tumorname}.gripss.filtered.vcf.gz.tbi")

    script:
    """
    gridss --jar /opt2/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    -r $BWAGENOME \
    -l ${tumorname} \
    -o ${tumorname}.vcf.gz -b $BLACKLIST \
    --picardoptions VALIDATION_STRINGENCY=LENIENT \
    ${tumor} -t $task.cpus

    mkdir -p ${tumorname}

    java -jar /opt2/hmftools/gripss.jar \
        -sample ${tumorname} \
        -ref_genome_version $GENOMEVER \
        -ref_genome $GENOMEREF  \
        -pon_sgl_file $PONSGL \
        -pon_sv_file $PONSV \
        -known_hotspot_file $SVHOTSPOT \
        -repeat_mask_file $REPEATMASK \
        -vcf ${tumorname}.vcf.gz \
        -output_dir ${tumorname}

    mv ${tumorname}/* .
    """


    stub:

    """
    touch "${tumorname}.vcf.gz"
    touch "${tumorname}.vcf.gz.tbi"
    touch "${tumorname}.vcf.gz.assembly.bam"
    touch "${tumorname}.gripss.vcf.gz"
    touch "${tumorname}.gripss.vcf.gz.tbi"
    touch "${tumorname}.gripss.filtered.vcf.gz"
    touch "${tumorname}.gripss.filtered.vcf.gz.tbi"
    """
}


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
