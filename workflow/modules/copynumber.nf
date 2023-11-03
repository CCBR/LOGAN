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
FREECSIGNIFICANCE = params.freec_significance
FREECPLOT = params.freec_plot
}

cobalt="/data/SCLC-BRAINMETS/cn/cobalt_v1.14.jar"
amber="/data/SCLC-BRAINMETS/cn/amber-3.9.jar"
purple="/data/SCLC-BRAINMETS/cn/purple_v3.8.2.jar"
GERMLINEHET="/data/SCLC-BRAINMETS/cn/copy_number/GermlineHetPon.38.vcf.gz"
GCPROFILE='/data/SCLC-BRAINMETS/cn/copy_number/GC_profile.1000bp.38.cnp'
DIPLODREG='/data/SCLC-BRAINMETS/cn/copy_number/DiploidRegions.38.bed.gz'
ENSEMBLCACHE='/data/SCLC-BRAINMETS/cn/common/ensembl_data'
DRIVERS='/data/SCLC-BRAINMETS/cn/common/DriverGenePanel.38.tsv'
HOTSPOTS='/data/SCLC-BRAINMETS/cn/variants/KnownHotspots.somatic.38.vcf.gz'



//DBSNP_INDEL=file(params.genomes[params.genome].KNOWNINDELS) 
//ascatR=

outdir=file(params.output)

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
    publishDir("${outdir}/cnv/sequenza", mode: 'copy')

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


process freec {
    label 'process_highcpu'
    publishDir("${outdir}/cnv/freec", mode: 'copy')

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    shell: """

    perl $FREECSCRIPT \
        . \
        $FREECLENGTHS \
        $FREECCHROMS \
        ${tumor} \
        $FREECPILE \
        $GENOMEREF \
        $FREECSNPS \
        $FREECTARGETS

    freec -conf freec_genome_config.txt

    cat $FREECSIGNIFICANCE | \
        R --slave \
        --args ${tumorname}.bam_CNVs \
        ${tumorname}.bam_ratio.txt

    cat $FREECPLOT | \
        R --slave \
        --args 2 \
        ${tumorname}.bam_ratio.txt \
        ${tumorname}.bam_BAF.txt

    """      

    stub:
    """
    touch ${tumorname}.bam_CNVs.p.value.txt  
    touch ${tumorname}.bam_ratio.txt 
    touch ${tumorname}.bam_BAF.txt 
    """
}



process amber {
    module=["java/12.0.1","R/3.6.3"]
    publishDir("${outdir}/amber", mode: 'copy')
    input:
        tuple val(samplename), path("${samplename}.bqsr.bam"), path("${samplename}.bqsr.bam.bai")

    output:
        tuple val(samplename), path("${samplename}_amber")
        //path("${samplename}.amber.baf.tsv.gz"),
        //path("${samplename}.amber.baf.pcf"),
        //path("${samplename}.amber.qc")
        //path("${samplename}.amber.contamination.vcf.gz") Contamination maybe only with tumor

    script:

    """

    java -Xmx32G -cp $amber com.hartwig.hmftools.amber.AmberApplication \
    -tumor ${samplename} -tumor_bam ${samplename}.bqsr.bam \
    -output_dir ${samplename}_amber \
    -threads 16 \
    -ref_genome_version V38 \
    -loci $GERMLINEHET

    """

    stub:

    """
    mkdir ${samplename}_amber
    touch ${samplename}_amber/${samplename}.amber.baf.tsv.gz ${samplename}_amber/${samplename}.amber.baf.pcf ${samplename}_amber/${samplename}.amber.qc 
    """
}

process cobalt {
    module=["java/12.0.1","R/3.6.3"]
    publishDir("${outdir}/cobalt", mode: 'copy')

    input:
        tuple val(samplename), path("${samplename}.bqsr.bam"), path("${samplename}.bqsr.bam.bai")

    output:
        tuple val(samplename), path("${samplename}_cobalt")
        //path("${samplename}/${samplename}.cobalt.ratio.tsv.gz"), 
        //path("${samplename}/${samplename}.cobalt.ratio.pcf"),
        //path("${samplename}/${samplename}.cobalt.gc.median.tsv")

    script:

    """

    java -jar -Xmx8G $cobalt \
    -tumor ${samplename} -tumor_bam ${samplename}.bqsr.bam \
    -output_dir ${samplename}_cobalt \
    -threads 16 \
    -tumor_only_diploid_bed $DIPLODREG \
    -gc_profile $GCPROFILE

    """

    stub:

    """
    mkdir ${samplename}_cobalt
    touch ${samplename}_cobalt/${samplename}.cobalt.ratio.tsv.gz ${samplename}_cobalt/${samplename}.cobalt.ratio.pcf ${samplename}_cobalt/${samplename}.cobalt.gc.median.tsv
    """
}


process purple {
    module=["java/12.0.1","R/3.6.3"]

    publishDir("${outdir}/purple", mode: 'copy')

    input:
        tuple val(samplename), path(cobaltin), path(amberin), path("${samplename}.tonly.final.mut2.vcf.gz")

    output:
        tuple val(samplename), path("${samplename}")

    script:

    """

    java -jar $purple \
    -tumor ${samplename} \
    -amber ${amberin} \
    -cobalt ${cobaltin} \
    -gc_profile $GCPROFILE \
    -ref_genome_version 38 \
    -ref_genome $GENOME \
    -ensembl_data_dir $ENSEMBLCACHE \
    -somatic_vcf ${samplename}.tonly.final.mut2.vcf.gz \
    -driver_gene_panel $DRIVERS \
    -somatic_hotspots $HOTSPOTS \
    -output_dir ${samplename}

    """

    stub:

    """
    mkdir ${samplename}
    touch ${samplename}/${samplename}.purple.cnv.somatic.tsv ${samplename}/${samplename}.purple.cnv.gene.tsv ${samplename}/${samplename}.driver.catalog.somatic.tsv
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

