GENOMEREF = file(params.genomes[params.genome].genome)
SEQUENZAGC = file(params.genomes[params.genome].SEQUENZAGC)
SEQUENZA_SCRIPT = params.script_sequenza

if (params.genome=="mm10"){
    FREECLENGTHS = params.genomes[params.genome].FREEC.FREECLENGTHS
    FREECCHROMS = params.genomes[params.genome].FREEC.FREECCHROMS
    FREECPILEUP = params.genomes[params.genome].FREEC.FREECPILEUP
    FREECSNPS = params.genomes[params.genome].FREEC.FREECSNPS
    FREECTARGETS = params.genomes[params.genome].intervals
    FREECSCRIPT = params.script_freec
    FREECPAIR_SCRIPT = params.script_freecpaired
    FREECSIGNIFICANCE = params.freec_significance
    FREECPLOT = params.freec_plot
}

if (params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
    HMFGENOMEREF = file(params.genomes[params.genome].HMFGENOME)
    GENOMEVER = params.genomes[params.genome].GENOMEVER
    GCPROFILE = file(params.genomes[params.genome].GCPROFILE)
    GERMLINEHET = file(params.genomes[params.genome].GERMLINEHET)
    DIPLODREG = file(params.genomes[params.genome].DIPLODREG)
    ENSEMBLCACHE = params.genomes[params.genome].ENSEMBLCACHE
    DRIVERS = file(params.genomes[params.genome].DRIVERS)
    SOMATICHOTSPOTS = file(params.genomes[params.genome].SOMATICHOTSPOTS)
    GERMLINEHOTSPOTS = file(params.genomes[params.genome].GERMLINEHOTSPOTS)
}

//mm10 Paired-Sequenza, FREEC-tumor only
process seqz_sequenza_bychr {
    container = "${params.containers.logan}"
    label 'process_long'

    input:
        tuple val(pairid), val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), val(chr)

    output:
        tuple val(pairid), path("${tumorname}_${normalname}_${chr}.seqz.gz")

    script:
    """
        sequenza-utils bam2seqz \
        -gc ${SEQUENZAGC} \
        -F $GENOMEREF \
        -C ${chr} \
        -n ${normal} \
        -t ${tumor} | gzip > "${tumorname}_${normalname}_${chr}.seqz.gz"

    """

    stub:
    """
    touch "${tumorname}_${normalname}_${chr}.seqz.gz"
    """
}

process pileup_sequenza {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(pairid), val(name), 
        path(bam), path(bai), path(bed)

    output:
        tuple val(pairid), path("${name}_${bed}.mpileup.gz"), path("${name}_${bed}.mpileup.gz.tbi") 

    script:
    //Q20 is default in sequenza
    """
        samtools mpileup -f $GENOMEREF -R ${bed} -Q 20 ${bam} |gzip > ${name}_${bed}.mpileup.gz
        tabix -s1 -b2 -e2 ${name}_${bed}.mpileup.gz
    """

    stub:
    """
    touch "${name}_${bed}.mpileup.gz"
    touch "${name}_${bed}.mpileup.gz.tbi"
    """
}

process seqz_sequenza_reg {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(pairid), val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(pairid), path("${tumorname}_${normalname}_${chr}.seqz.gz")

    script:
    """
        sequenza-utils bam2seqz \
        -gc ${SEQUENZAGC} \
        -p \
        -F $GENOMEREF \
        -n ${normal} \
        -t ${tumor} | gzip > "${tumorname}_${normalname}_${bed}.seqz.gz"

    """

    stub:
    """
    touch "${tumorname}_${normalname}_${chr}.seqz.gz"
    """
}

process seqz_sequenza {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(pairid), val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(pairid), path("${tumorname}_${normalname}_${chr}.seqz.gz")

    script:
    """
        sequenza-utils bam2seqz \
        -gc ${SEQUENZAGC} \
        -p \
        -F $GENOMEREF \
        -n ${normal} \
        -t ${tumor} | gzip > "${tumorname}_${normalname}_${bed}.seqz.gz"

    """

    stub:
    """
    touch "${tumorname}_${normalname}_${chr}.seqz.gz"
    """
}




process sequenza {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(pairid), path(seqz)

    output:
        tuple val(pairid),
        path("${pairid}_alternative_solutions.txt"),
        path("${pairid}_alternative_fit.pdf"),
        path("${pairid}_model_fit.pdf"),
        path("${pairid}_confints_CP.txt"),
        path("${pairid}_CN_bars.pdf"),
        path("${pairid}_genome_view.pdf"),
        path("${pairid}_chromosome_view.pdf"),
        path("${pairid}_mutations.txt"),
        path("${pairid}_segments.txt"),
        path("${pairid}_CP_contours.pdf"),
        path("${pairid}_sequenza_cp_table.RData"),
        path("${pairid}_chromosome_depths.pdf"),
        path("${pairid}_gc_plots.pdf"),
        path("${pairid}_sequenza_extract.RData")


    shell:
    '''

    zcat !{seqz} | awk '{if (NR==1) {print $0} else {if ($1!="chromosome"){print $0}}}' |\
    sequenza-utils seqz_binning \
        -w 100 \
        -s - > !{pairid}.bin100.seqz

    Rscript !{SEQUENZA_SCRIPT} \
        !{pairid}.bin100.seqz \
        . \
        !{pairid} \
        !{task.cpus}

    '''

    stub:

    """
    touch "${pairid}_alternative_solutions.txt"
    touch "${pairid}_alternative_fit.pdf"
    touch "${pairid}_model_fit.pdf"
    touch "${pairid}_confints_CP.txt"
    touch "${pairid}_CN_bars.pdf"
    touch "${pairid}_genome_view.pdf"
    touch "${pairid}_chromosome_view.pdf"
    touch "${pairid}_mutations.txt"
    touch "${pairid}_segments.txt"
    touch "${pairid}_CP_contours.pdf"
    touch "${pairid}_sequenza_cp_table.RData"
    touch "${pairid}_chromosome_depths.pdf"
    touch "${pairid}_gc_plots.pdf"
    touch "${pairid}_sequenza_extract.RData"

    """

}


process freec_paired {
    container = "${params.containers.logan}"
    label 'process_long'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_CNVs.p.value.txt"),
        path("${tumorname}_vs_${normalname}_ratio.txt"),
        path("${tumorname}_vs_${normalname}_BAF.txt"),
        path("${tumorname}_vs_${normalname}_ratio.txt.log2.png"),
        path("${tumorname}_vs_${normalname}_ratio.txt.png")

    shell:
    """

    perl $FREECPAIR_SCRIPT \
        . \
        $FREECLENGTHS \
        $FREECCHROMS \
        ${tumor} \
        ${normal} \
        $FREECPILEUP \
        $GENOMEREF \
        $FREECSNPS \
        $FREECTARGETS

    freec -conf freec_genome_config.txt

    cat $FREECSIGNIFICANCE | \
        R --slave \
        --args ${tumor}_CNVs \
        ${tumor}_ratio.txt

    cat $FREECPLOT | \
        R --slave \
        --args 2 \
        ${tumor}_ratio.txt \
        ${tumor}_BAF.txt

    mv ${tumor}_CNVs.p.value.txt ${tumorname}_vs_${normalname}_CNVs.p.value.txt
    mv ${tumor}_ratio.txt ${tumorname}_vs_${normalname}_ratio.txt
    mv ${tumor}_BAF.txt ${tumorname}_vs_${normalname}_BAF.txt
    mv ${tumor}_BAF.txt.png ${tumorname}_vs_${normalname}_BAF.txt.png
    mv ${tumor}_ratio.txt.log2.png ${tumorname}_vs_${normalname}_ratio.txt.log2.png
    mv ${tumor}_ratio.txt.png ${tumorname}_vs_${normalname}_ratio.txt.png

    """

    stub:
    """
    touch ${tumorname}_vs_${normalname}_CNVs.p.value.txt
    touch ${tumorname}_vs_${normalname}_ratio.txt
    touch ${tumorname}_vs_${normalname}_BAF.txt
    touch ${tumorname}_vs_${normalname}_BAF.txt.png
    touch ${tumorname}_vs_${normalname}_ratio.txt.log2.png
    touch ${tumorname}_vs_${normalname}_ratio.txt.png

    """
}


process freec {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        tuple val(tumorname),
        path("${tumorname}_CNVs.p.value.txt"),
        path("${tumorname}_ratio.txt"),
        path("${tumorname}_BAF.txt"),
        path("${tumorname}_ratio.txt.log2.png"),
        path("${tumorname}_ratio.txt.png")


    shell: """

    perl $FREECSCRIPT \
        . \
        $FREECLENGTHS \
        $FREECCHROMS \
        ${tumor} \
        $FREECPILEUP \
        $GENOMEREF \
        $FREECSNPS \
        $FREECTARGETS

    freec -conf freec_genome_config.txt

    cat $FREECSIGNIFICANCE | \
        R --slave \
        --args ${tumor}_CNVs \
        ${tumor}_ratio.txt

    cat $FREECPLOT | \
        R --slave \
        --args 2 \
        ${tumor}_ratio.txt \
        ${tumor}_BAF.txt

    mv ${tumor}_CNVs.p.value.txt ${tumorname}_CNVs.p.value.txt
    mv ${tumor}_ratio.txt ${tumorname}_ratio.txt
    mv ${tumor}_BAF.txt ${tumorname}_BAF.txt
    mv ${tumor}_BAF.txt.png ${tumorname}_BAF.txt.png
    mv ${tumor}_ratio.txt.log2.png ${tumorname}_ratio.txt.log2.png
    mv ${tumor}_ratio.txt.png ${tumorname}_ratio.txt.png

    """

    stub:
    """
    touch ${tumorname}_CNVs.p.value.txt
    touch ${tumorname}_ratio.txt
    touch ${tumorname}_BAF.txt
    touch ${tumorname}_BAF.txt.png
    touch ${tumorname}_ratio.txt.log2.png
    touch ${tumorname}_ratio.txt.png

    """
}


process amber_tonly {
    container = "${params.containers.logan}"

    label 'process_medium'

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

    input:
        tuple val(id), val(tumorname), val(normalname),
        path(amberin), path(cobaltin)

    output:
        tuple val(id), val(tumorname), val(normalname), path("${id}")

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
    mkdir ${tumorname}
    touch ${tumorname}/${tumorname}.purple.cnv.somatic.tsv ${tumorname}/${tumorname}.purple.cnv.gene.tsv ${tumorname}/${tumorname}.driver.catalog.somatic.tsv
    """

}


process purple_tonly {
    container = "${params.containers.logan}"
    label 'process_medium'

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

/*
process ascat_tn {
    module=["java/12.0.1","R/3.6.3"]

    input:
        tuple val(samplename), path(cobaltin), path(amberin), path("${samplename}.tonly.final.mut2.vcf.gz")

    output:
        tuple val(samplename), path("${samplename}")

    script:

    """
    Rscript ${ascatR}
    """

    stub:

    """
    touch ${prefix}.after_correction.gc_rt.test.tumour.germline.png
    touch ${prefix}.after_correction.gc_rt.test.tumour.tumour.png
    touch ${prefix}.before_correction.test.tumour.germline.png
    touch ${prefix}.before_correction.test.tumour.tumour.png
    touch ${prefix}.cnvs.txt
    touch ${prefix}.metrics.txt
    touch ${prefix}.normal_alleleFrequencies_chr21.txt
    touch ${prefix}.normal_alleleFrequencies_chr22.txt
    touch ${prefix}.purityploidy.txt
    touch ${prefix}.segments.txt
    touch ${prefix}.tumour.ASPCF.png
    touch ${prefix}.tumour.sunrise.png
    touch ${prefix}.tumour_alleleFrequencies_chr21.txt
    touch ${prefix}.tumour_alleleFrequencies_chr22.txt
    touch ${prefix}.tumour_normalBAF.txt
    touch ${prefix}.tumour_normalLogR.txt
    touch ${prefix}.tumour_tumourBAF.txt
    touch ${prefix}.tumour_tumourLogR.txt
        """

}

*/
