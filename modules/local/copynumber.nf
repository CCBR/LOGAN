GENOMEREF=file(params.genomes[params.genome].genome)
SEQUENZAGC=file(params.genomes[params.genome].SEQUENZAGC)
SEQUENZA_SCRIPT=params.script_sequenza

if (params.genome=="mm10"){
FREECLENGTHS=file(params.genomes[params.genome].FREEC.FREECLENGTHS)
FREECCHROMS=file(params.genomes[params.genome].FREEC.FREECCHROMS)
FREECPILEUP=file(params.genomes[params.genome].FREEC.FREECPILEUP)
FREECSNPS = file(params.genomes[params.genome].FREEC.FREECSNPS)
FREECTARGETS=file(params.genomes[params.genome].intervals)
FREECSCRIPT = params.script_freec
FREECPAIR_SCRIPT = params.script_freecpaired
FREECSIGNIFICANCE = params.freec_significance
FREECPLOT = params.freec_plot
}

GERMLINEHET="/data/SCLC-BRAINMETS/cn/copy_number/GermlineHetPon.38.vcf.gz"
GCPROFILE='/data/SCLC-BRAINMETS/cn/copy_number/GC_profile.1000bp.38.cnp'
DIPLODREG='/data/SCLC-BRAINMETS/cn/copy_number/DiploidRegions.38.bed.gz'
ENSEMBLCACHE='/data/SCLC-BRAINMETS/cn/common/ensembl_data'
DRIVERS='/data/SCLC-BRAINMETS/cn/common/DriverGenePanel.38.tsv'
HOTSPOTS='/data/SCLC-BRAINMETS/cn/variants/KnownHotspots.somatic.38.vcf.gz'

//DBSNP_INDEL=file(params.genomes[params.genome].KNOWNINDELS) 
//ascatR=


//mm10 Paired-Sequenza, FREEC-tumor only 
process seqz_sequenza_bychr {
    label 'process_low'

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




process sequenza {
    label 'process_highcpu'

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

    //samtools mpileup ${tumor} -f $GENOMEREF -Q 20 |gzip > ${tumorname}.mpileup.gz
    //samtools mpileup ${normal} -f $GENOMEREF -Q 20 |gzip > ${normalname}.mpileup.gz
    //sequenza-utils seqz_binning --seqz --window 50 -o ${sample}_bin50.seqz.gz

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
    label 'process_highcpu'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai)

    shell: """

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
        ${tumorname}_vs_${normalname}.bam_ratio.txt \
        ${tumorname}_vs_${normalname}.bam_BAF.txt

    """      

    stub:
    """
    touch ${tumorname}_vs_${normalname}.bam_CNVs.p.value.txt  
    touch ${tumorname}_vs_${normalname}.bam_ratio.txt 
    touch ${tumorname}_vs_${normalname}.bam_BAF.txt
    touch ${tumorname}_vs_${normalname}.bam_ratio.txt.log2.png
    touch ${tumorname}_vs_${normalname}.bam_ratio.txt.png

    """
}


process freec {
    label 'process_mid'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

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

    """      

    stub:
    """
    touch ${tumor}_CNVs.p.value.txt  
    touch ${tumor}_ratio.txt 
    touch ${tumor}_BAF.txt 
    touch ${tumor}_ratio.txt.log2.png
    touch ${tumor}_ratio.txt.png

    """
}


process amber_tonly {
    label 'process_mid'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)
       

    output:
        tuple val(tumorname), path("${tumorname}_amber")
        //path("${samplename}.amber.baf.tsv.gz"),
        //path("${samplename}.amber.baf.pcf"),
        //path("${samplename}.amber.qc")
        //path("${samplename}.amber.contamination.vcf.gz") Contamination maybe only with tumor

    script:

    """

    java -Xmx32G -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -tumor ${tumorname} -tumor_bam ${tumor} \
    -output_dir ${tumorname}_amber \
    -threads $task.cpus \
    -ref_genome_version V38 \
    -loci $GERMLINEHET

    """

    stub:

    """
    mkdir ${tumorname}_amber
    touch ${tumorname}_amber/${tumorname}.amber.baf.tsv.gz ${tumorname}_amber/${tumorname}.amber.baf.pcf ${tumorname}_amber/${tumorname}.amber.qc 
    """
}

process amber_tn {
    label 'process_mid'
    
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname), path("${tumorname}_vs_${normalname}_amber")
        //path("${samplename}.amber.baf.tsv.gz"),
        //path("${samplename}.amber.baf.pcf"),
        //path("${samplename}.amber.qc")
        //path("${samplename}.amber.contamination.vcf.gz") Contamination maybe only with tumor

    script:

    """

    java -Xmx32G -cp amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -tumor ${tumorname} -tumor_bam ${tumor} \
    -reference ${normalname} -reference_bam ${normal} \
    -output_dir ${tumorname}_vs_${normalname}_amber \
    -threads $task.cpus \
    -ref_genome_version V38 \
    -loci $GERMLINEHET

    """

    stub:

    """
    mkdir ${tumorname}_vs_${normalname}_amber
    touch ${tumorname}_vs_${normalname}_amber/${tumorname}.amber.baf.tsv.gz ${tumorname}_vs_${normalname}_amber/${tumorname}.amber.baf.pcf ${tumorname}_vs_${normalname}_amber/${tumorname}.amber.qc 
    """
}

process cobalt_tonly {
    label "process_mid"

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        tuple val(tumorname), path("${tumorname}_cobalt")
        //path("${samplename}/${samplename}.cobalt.ratio.tsv.gz"), 
        //path("${samplename}/${samplename}.cobalt.ratio.pcf"),
        //path("${samplename}/${samplename}.cobalt.gc.median.tsv")

    script:

    """

    java -jar -Xmx8G cobalt.jar \
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
    label "process_mid"

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname), path("${tumorname}_vs_${normalname}_cobalt")
        //path("${samplename}/${samplename}.cobalt.ratio.tsv.gz"), 
        //path("${samplename}/${samplename}.cobalt.ratio.pcf"),
        //path("${samplename}/${samplename}.cobalt.gc.median.tsv")

    script:

    """

    java -jar -Xmx8G cobalt.jar \
    -tumor ${tumorname} -tumor_bam ${tumorname} \
    -reference ${normalname} -reference_bam ${normal} \
    -output_dir ${tumorname}_vs_${normalname}_cobalt \
    -threads $task.cpus \
    -tumor_only_diploid_bed $DIPLODREG \
    -gc_profile $GCPROFILE

    """

    stub:

    """
    mkdir ${tumorname}_vs_${normalname}_cobalt
    touch ${tumorname}_vs_${normalname}_cobalt/${tumorname}.cobalt.ratio.tsv.gz ${tumorname}_vs_${normalname}_cobalt/${tumorname}.cobalt.ratio.pcf ${tumorname}_vs_${normalname}_cobalt/${tumorname}.cobalt.gc.median.tsv
    """
}


process purple {
    label 'process_mid'
    publishDir("${outdir}/cnv/purple", mode: 'copy')

    input:
        tuple val(tumorname),
        path(cobaltin), 
        path(amberin),
        path(somaticvcf),
        path(somaticvcfindex)

    output:
        tuple val(tumorname), path("${tumorname}")

    script:

    """
    java -jar purple.jar \
    -tumor ${tumorname} \
    -amber ${amberin} \
    -cobalt ${cobaltin} \
    -gc_profile $GCPROFILE \
    -ref_genome_version 38 \
    -ref_genome $GENOME \
    -ensembl_data_dir $ENSEMBLCACHE \
    -somatic_vcf ${somaticvcf} \
    -driver_gene_panel $DRIVERS \
    -somatic_hotspots $HOTSPOTS \
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

    publishDir("${outdir}/purple", mode: 'copy')

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

