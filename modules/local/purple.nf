//PURPLE
GENOMEREF=file(params.genomes[params.genome].genome)
HMFGENOMEREF = file(params.genomes[params.genome].HMFGENOME)
GENOMEVER = params.genomes[params.genome].GENOMEVER
GCPROFILE = file(params.genomes[params.genome].GCPROFILE)
GERMLINEHET = file(params.genomes[params.genome].GERMLINEHET)
DIPLODREG = file(params.genomes[params.genome].DIPLODREG)
ENSEMBLCACHE = params.genomes[params.genome].ENSEMBLCACHE
DRIVERS = file(params.genomes[params.genome].DRIVERS)
SOMATICHOTSPOTS = file(params.genomes[params.genome].SOMATICHOTSPOTS)
GERMLINEHOTSPOTS = file(params.genomes[params.genome].GERMLINEHOTSPOTS)
//if (params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
//}



process amber_tonly {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'


    input:
        tuple val(tumorname), path(tumor), path(tumorbai)


    output:
        tuple val(tumorname), path("${tumorname}_amber")
      
    script:

    """

    java -Xmx32G -cp /opt2/hmftools/amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -tumor ${tumorname} -tumor_bam ${tumor} \
    -output_dir ${tumorname}_amber \
    -threads $task.cpus \
    -ref_genome_version $GENOMEVER \
    -loci $GERMLINEHET

    """

    stub:

    """
    mkdir ${tumorname}_amber
    touch ${tumorname}_amber/${tumorname}.amber.baf.tsv.gz ${tumorname}_amber/${tumorname}.amber.baf.pcf ${tumorname}_amber/${tumorname}.amber.qc
    """
}

process amber_tn {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val("${tumorname}_vs_${normalname}"),
        val(tumorname), val(normalname), path("${tumorname}_vs_${normalname}_amber")
      
    script:

    """

    java -Xmx32G -cp /opt2/hmftools/amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -tumor ${tumorname} -tumor_bam ${tumor} \
    -reference ${normalname} -reference_bam ${normal} \
    -output_dir ${tumorname}_vs_${normalname}_amber \
    -threads $task.cpus \
    -ref_genome_version $GENOMEVER \
    -loci $GERMLINEHET

    """

    stub:

    """
    mkdir ${tumorname}_vs_${normalname}_amber
    touch ${tumorname}_vs_${normalname}_amber/${tumorname}.amber.baf.tsv.gz ${tumorname}_vs_${normalname}_amber/${tumorname}.amber.baf.pcf ${tumorname}_vs_${normalname}_amber/${tumorname}.amber.qc
    """
}

process cobalt_tonly {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        tuple val(tumorname), path("${tumorname}_cobalt")

    script:

    """

    java -jar -Xmx8G /opt2/hmftools/cobalt.jar \
    -tumor ${tumorname} -tumor_bam ${tumor} \
    -output_dir ${tumorname}_cobalt \
    -threads $task.cpus \
    -tumor_only_diploid_bed $DIPLODREG \
    -gc_profile $GCPROFILE

    """

    stub:

    """
    mkdir ${tumorname}_cobalt
    touch ${tumorname}_cobalt/${tumorname}.cobalt.ratio.tsv.gz ${tumorname}_cobalt/${tumorname}.cobalt.ratio.pcf ${tumorname}_cobalt/${tumorname}.cobalt.gc.median.tsv
    """
}

process cobalt_tn {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val("${tumorname}_vs_${normalname}"),
        val(tumorname), val(normalname), path("${tumorname}_vs_${normalname}_cobalt")

    script:

    """
    java -jar -Xmx8G /opt2/hmftools/cobalt.jar \
    -tumor ${tumorname} -tumor_bam ${tumor} \
    -reference ${normalname} -reference_bam ${normal} \
    -output_dir ${tumorname}_vs_${normalname}_cobalt \
    -threads $task.cpus \
    -gc_profile $GCPROFILE

    """

    stub:

    """
    mkdir ${tumorname}_vs_${normalname}_cobalt
    touch ${tumorname}_vs_${normalname}_cobalt/${tumorname}.cobalt.ratio.tsv.gz ${tumorname}_vs_${normalname}_cobalt/${tumorname}.cobalt.ratio.pcf ${tumorname}_vs_${normalname}_cobalt/${tumorname}.cobalt.gc.median.tsv
    """
}


process purple {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
        tuple val(id), val(tumorname), val(normalname),
        path(amberin), path(cobaltin),
        path(somaticvcf), path(somaticvcfindex)

    output:
        tuple val(id), path("${id}")

    script:

    """
    java -jar /opt2/hmftools/purple.jar \
    -tumor ${tumorname} \
    -reference ${normalname} \
    -amber ${amberin} \
    -cobalt ${cobaltin} \
    -gc_profile $GCPROFILE \
    -ref_genome_version $GENOMEVER \
    -ref_genome $GENOMEREF \
    $ENSEMBLCACHE \
    -somatic_vcf ${somaticvcf} \
    -driver_gene_panel $DRIVERS \
    -somatic_hotspots $SOMATICHOTSPOTS \
    -threads $task.cpus \
    -output_dir ${id}
    """

    stub:

    """
    mkdir ${id}
    touch ${id}/${id}.purple.cnv.somatic.tsv ${id}/${id}.purple.cnv.gene.tsv ${id}/${id}.driver.catalog.somatic.tsv
    """

}


process purple_novc {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
        tuple val(id), val(tumorname), val(normalname),
        path(amberin), path(cobaltin)

    output:
        tuple val(id), val(tumorname), val(normalname), 
            path("${id}")

    script:

    """
    java -jar /opt2/hmftools/purple.jar \
    -tumor ${tumorname} \
    -reference ${normalname} \
    -amber ${amberin} \
    -cobalt ${cobaltin} \
    -gc_profile $GCPROFILE \
    -ref_genome_version $GENOMEVER \
    -ref_genome $HMFGENOMEREF \
    $ENSEMBLCACHE \
    -threads $task.cpus \
    -output_dir ${id}

    """

    stub:

    """
    mkdir ${id}
    touch ${id}/${id}.purple.cnv.somatic.tsv ${id}/${id}.purple.cnv.gene.tsv ${id}/${id}.driver.catalog.somatic.tsv
    """

}


process purple_tonly {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), 
        path(amberin), path(cobaltin),
        path(somaticvcf), path(somaticvcfindex)

    output:
        tuple val(tumorname), path("${tumorname}")

    script:

    """
    java -jar /opt2/hmftools/purple.jar \
    -tumor ${tumorname} \
    -amber ${amberin} \
    -cobalt ${cobaltin} \
    -gc_profile $GCPROFILE \
    -ref_genome_version $GENOMEVER \
    -ref_genome $GENOMEREF \
    $ENSEMBLCACHE \
    -somatic_vcf ${somaticvcf} \
    -driver_gene_panel $DRIVERS \
    -somatic_hotspots $HOTSPOTS \
    -threads $task.cpus \
    -output_dir ${tumorname}
    """

    stub:

    """
    mkdir ${tumorname}
    touch ${tumorname}/${tumorname}.purple.cnv.somatic.tsv ${tumorname}/${tumorname}.purple.cnv.gene.tsv ${tumorname}/${tumorname}.driver.catalog.somatic.tsv
    """

}


process purple_tonly_novc {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'
    
    input:
        tuple val(tumorname), val(normalname),
        path(cobaltin), path(amberin)

    output:
        tuple val(tumorname), path("${tumorname}")

    script:

    """
    java -jar /opt2/hmftools/purple.jar \
    -tumor ${tumorname} \
    -amber ${amberin} \
    -cobalt ${cobaltin} \
    -gc_profile $GCPROFILE \
    -ref_genome_version $GENOMEVER \
    -ref_genome $GENOMEREF \
    $ENSEMBLCACHE \
    -threads $task.cpus \
    -output_dir ${tumorname}
    """

    stub:

    """
    mkdir ${tumorname}
    touch ${tumorname}/${tumorname}.purple.cnv.somatic.tsv ${tumorname}/${tumorname}.purple.cnv.gene.tsv ${tumorname}/${tumorname}.driver.catalog.somatic.tsv
    """

}

