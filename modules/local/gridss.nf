BWAGENOME=file(params.genomes[params.genome].bwagenome)
BLACKLIST=file(params.genomes[params.genome].GRIDSSBLACKLIST)
GENOMEREF=file(params.genomes[params.genome].genome)

if (params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
    GENOMEVER = params.genomes[params.genome].GENOMEVER
    PONSGL = file(params.genomes[params.genome].PONSGL)
    PONSV = file(params.genomes[params.genome].PONSV)
    SVHOTSPOT = file(params.genomes[params.genome].SVHOTSPOT)
    REPEATMASK = file(params.genomes[params.genome].REPEATMASK)
}



process gridss_somatic {
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
