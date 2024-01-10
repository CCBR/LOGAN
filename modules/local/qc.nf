///References to assign
GENOMEREF=file(params.genomes[params.genome].genome)
DBSNP=file(params.genomes[params.genome].dbsnp) //dbsnp_138.hg38.vcf.gz"
FASTQ_SCREEN_CONF=file(params.fastq_screen_conf)
BACDB=file(params.genomes[params.genome].KRAKENBACDB)
SNPEFF_GENOME = params.genomes[params.genome].snpeff_genome
SNPEFF_CONFIG = file(params.genomes[params.genome].snpeff_config)
SNPEFF_BUNDLE = file(params.genomes[params.genome].snpeff_bundle)

//SOMALIER
SITES_VCF= file(params.genomes[params.genome].sites_vcf)
ANCESTRY_DB=file(params.genomes[params.genome].somalier_ancestrydb)
SCRIPT_PATH_GENDER = file(params.script_genderPrediction)
SCRIPT_PATH_SAMPLES = file(params.script_combineSamples)
SCRIPT_PATH_PCA = file(params.script_ancestry)


//OUTPUT DIRECTORY
process fc_lane {
    container = "${params.containers.logan}"

    label 'process_low'

    input:
        tuple val(samplename), path(fqs)

    output:
        tuple val(samplename),
        path("${samplename}.fastq.info.txt")

    script:
    GET_FLOWCELL_LANES=file(params.get_flowcell_lanes)

    """
    python $GET_FLOWCELL_LANES \
    ${fqs[0]}  \
    ${samplename} > ${samplename}.fastq.info.txt
    """

    stub:
    """
    touch ${samplename}.fastq.info.txt
    """
}


process fastq_screen {
    //Uses Trimmed Files


    input:
    tuple val(samplename),
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz"),
        path("${samplename}.fastp.json"),
        path("${samplename}.fastp.html")

    output:
        tuple path("${samplename}.R1.trimmed_screen.html"),
        path("${samplename}.R1.trimmed_screen.png"),
        path("${samplename}.R1.trimmed_screen.txt"),
        path("${samplename}.R2.trimmed_screen.html"),
        path("${samplename}.R2.trimmed_screen.png"),
        path("${samplename}.R2.trimmed_screen.txt")

    script:
        FASTQ_SCREEN_CONF=file(params.fastq_screen_conf)

        """
    fastq_screen --conf $FASTQ_SCREEN_CONF \
        --outdir . \
        --threads 8 \
        --subset 1000000 \
        --aligner bowtie2 \
        --force \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz

        """

    stub:
    """
    touch ${samplename}.R1.trimmed_screen.html ${samplename}.R1.trimmed_screen.png
    touch ${samplename}.R1.trimmed_screen.txt ${samplename}.R2.trimmed_screen.html
    touch ${samplename}.R2.trimmed_screen.png ${samplename}.R2.trimmed_screen.txt
    """
}

process kraken {
    /*
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interactive krona report
    */

    input:
        tuple val(samplename),
        path(fqs)

    output:
        tuple val(samplename),
        //path("${samplename}.trimmed.kraken_bacteria.out.txt"),
        path("${samplename}.trimmed.kraken_bacteria.taxa.txt"),
        path("${samplename}.trimmed.kraken_bacteria.krona.html")


    script:
    """
    #Setups temporary directory for
    #intermediate files with built-in
    #mechanism for deletion on exit


    # Copy kraken2 db to local node storage to reduce filesystem strain
    cp -rv $BACDB .
    kdb_base=\$(basename $BACDB)

    kraken2 --db $BACDB \
        --threads 16 --report ${samplename}.trimmed.kraken_bacteria.taxa.txt \
        --output - \
        --gzip-compressed \
        --paired ${fqs[0]} ${fqs[1]}
    # Generate Krona Report
    cut -f2,3 ${samplename}.trimmed.kraken_bacteria.taxa.txt} | \
        ktImportTaxonomy - -o ${samplename}.trimmed.kraken_bacteria.krona.html
    """

    stub:
    """
    touch  ${samplename}.trimmed.kraken_bacteria.taxa.txt ${samplename}.trimmed.kraken_bacteria.krona.html
    """

}

process fastqc {
    """
    Quality-control step to assess sequencing quality of each sample.
    FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        FastQC report and zip file containing sequencing quality information
    """


    input:
        tuple val(samplename), path("${samplename}.bqsr.bam"), path("${samplename}.bqsr.bai")
    output:
        tuple val(samplename), path("${samplename}_fastqc.html"), path("${samplename}_fastqc.zip")

    //message: "Running FastQC with {threads} threads on '{input}' input file"
    //threads: 8
    //module=['fastqc/0.11.9']

    script:
    """
    mkdir -p fastqc
    fastqc -t 8 \
        -f bam \
        -o fastqc \
        ${samplename}.bqsr.bam
    mv fastqc/${samplename}.bqsr_fastqc.html ${samplename}_fastqc.html
    mv fastqc/${samplename}.bqsr_fastqc.zip ${samplename}_fastqc.zip
    """

    stub:
    """
    touch  ${samplename}_fastqc.html ${samplename}_fastqc.zip
    """
}

process qualimap_bamqc {
    /*
    Quality-control step to assess various post-alignment metrics
    and a secondary method to calculate insert size. Please see
    QualiMap's website for more information about BAM QC:
    http://qualimap.conesalab.org/
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Report containing post-aligment quality-control metrics
    */

    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        tuple path("${samplename}_genome_results.txt"), path("${samplename}_qualimapReport.html")

    script:
    """
    unset DISPLAY
    qualimap bamqc -bam ${bam} \
        --java-mem-size=112G \
        -c -ip \
        -outdir ${samplename} \
        -outformat HTML \
        -nt 8 \
        --gd HUMAN \
        --skip-duplicated \
        -nw 500 \
        -p NON-STRAND-SPECIFIC
    mv ${samplename}/genome_results.txt ${samplename}_genome_results.txt
    mv ${samplename}/qualimapReport.html ${samplename}_qualimapReport.html
    """

    stub:
    """
    touch ${samplename}_genome_results.txt ${samplename}_qualimapReport.html
    """
}

process samtools_flagstats {
    /*
    Quality-control step to assess alignment quality. Flagstat provides
    counts for each of 13 categories based primarily on bit flags in the
    FLAG field. Information on the meaning of the flags is given in the
    SAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Text file containing alignment statistics
    */
    label 'process_medium'

    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        path("${samplename}.samtools_flagstat.txt")

    script:
    """
    samtools flagstat ${bam} > ${samplename}.samtools_flagstat.txt
    """

    stub:
    """
    touch ${samplename}.samtools_flagstat.txt
    """
}


process mosdepth {
    /*
    Quality-control step to assess depth
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        `{prefix}.mosdepth.global.dist.txt`
        `{prefix}.mosdepth.summary.txt`
        `{prefix}.mosdepth.region.dist.txt` (if --by is specified)
        `{prefix}.per-base.bed.gz|per-base.d4` (unless -n/--no-per-base is specified)
        `{prefix}.regions.bed.gz` (if --by is specified)
        `{prefix}.quantized.bed.gz` (if --quantize is specified)
        `{prefix}.thresholds.bed.gz` (if --thresholds is specified)
    */
    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        path("${samplename}.mosdepth.region.dist.txt"),
        path("${samplename}.mosdepth.summary.txt"),
        path("${samplename}.regions.bed.gz"),
        path("${samplename}.regions.bed.gz.csi")


    script:
    """
    mosdepth -n --fast-mode --by 500  ${samplename} ${bam} -t $task.cpus
    """

    stub:
    """
    touch "${samplename}.mosdepth.region.dist.txt"
    touch "${samplename}.mosdepth.summary.txt"
    touch "${samplename}.regions.bed.gz"
    touch "${samplename}.regions.bed.gz.csi"
    """
}

process vcftools {
    /*
    Quality-control step to calculates a measure of heterozygosity on
    a per-individual basis. The inbreeding coefficient, F, is estimated
    for each individual using a method of moments. Please see VCFtools
    documentation for more information:
    https://vcftools.github.io/man_latest.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a measure of heterozygosity
    */
    label 'process_medium'


    input:
        tuple path(germlinevcf),path(germlinetbi)
    output:
       path("variants_raw_variants.het")


    script:
    """
    vcftools --gzvcf ${germlinevcf} --het --out variants_raw_variants
    """

    stub:
    """
    touch variants_raw_variants.het
    """
}

process collectvariantcallmetrics {
    /*
    Quality-control step to collect summary metrics about snps and indels
    called in a multisample VCF file. Please see the Broad's documentation
    for more information about each field in the generated log file:
    https://broadinstitute.github.io/picard/picard-metric-definitions.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a collection of metrics relating to snps and indels
    */
    input:
        tuple path(germlinevcf),path(germlinetbi)

    output:
        tuple path("raw_variants.variant_calling_detail_metrics"),
        path("raw_variants.variant_calling_summary_metrics")


    script:
    """
    java -Xmx24g -jar \${PICARDJARPATH}/picard.jar \
        CollectVariantCallingMetrics \
        INPUT=${germlinevcf} \
        OUTPUT= "raw_variants" \
        DBSNP=$DBSNP Validation_Stringency=SILENT
    """

    stub:
    """
    touch raw_variants.variant_calling_detail_metrics raw_variants.variant_calling_summary_metrics
    """

}


process bcftools_stats {
    /*
    Quality-control step to collect summary statistics from bcftools stats.
    When bcftools stats is run with one VCF file then stats by non-reference
    allele frequency, depth distribution, stats by quality and per-sample
    counts, singleton statsistics are calculated. Please see bcftools'
    documentation for more information:
    http://samtools.github.io/bcftools/bcftools.html#stats
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Text file containing a collection of summary statistics
    */

    label 'process_medium'

    input:
        tuple val(samplename),  path("${samplename}.gvcf.gz"),path("${samplename}.gvcf.gz.tbi")
    output:
        path("${samplename}.germline.bcftools_stats.txt")

    script:
    """
    bcftools stats ${samplename}.gvcf.gz > ${samplename}.germline.bcftools_stats.txt
    """

    stub:
    """
    touch ${samplename}.germline.bcftools_stats.txt
    """

}

process gatk_varianteval {
    /*
    Quality-control step to calculate various quality control metrics from a
    variant callset. These metrics include the number of raw or filtered SNP
    counts; ratio of transition mutations to transversions; concordance of a
    particular sample's calls to a genotyping chip; number of s per sample.
    Please see GATK's documentation for more information:
    https://gatk.broadinstitute.org/hc/en-us/articles/360040507171-VariantEval
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    */
    label 'process_medium'

    input:
        tuple val(samplename), path("${samplename}.gvcf.gz") ,path("${samplename}.gvcf.gz.tbi")
    output:
        path("${samplename}.germline.eval.grp")
    //params:
     //   rname    = "vareval",
      //  genome   = config['references']['GENOME'],
       // dbsnp    = config['references']['DBSNP'],
      //  ver_gatk = config['tools']['gatk4']['version']
    //message: "Running GATK4 VariantEval on '{input.vcf}' input file"
    //container: config['images']['wes_base']
    //threads: 16
    script:
    """
    gatk --java-options '-Xmx12g -XX:ParallelGCThreads=16' VariantEval \
        -R $GENOMEREF \
        -O ${samplename}.germline.eval.grp \
        --dbsnp $DBSNP \
        --eval ${samplename}.gvcf.gz
    """

    stub:

    """
    touch ${samplename}.germline.eval.grp
    """

}

process snpeff {
    /*
    Data processing and quality-control step to annotate variants, predict its
    functional effects, and collect various summary statistics about variants and
    their annotations. Please see SnpEff's documentation for more information:
    https://pcingola.github.io/SnpEff/
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    */
    label 'process_medium'

    input:
        tuple val(samplename), path("${samplename}.gvcf.gz"), path("${samplename}.gvcf.gz.tbi")
    output:
        tuple path("${samplename}.germline.snpeff.ann.vcf"),
        path("${samplename}.germline.snpeff.ann.csv"),
        path("${samplename}.germline.snpeff.ann.html")

    script:
    """
        java -Xmx12g -jar \$SNPEFF_JAR \
        -v -canon -c $SNPEFF_CONFIG \
        -csvstats ${samplename}.germline.snpeff.ann.csv \
        -stats ${samplename}.germline.snpeff.ann.html \
        $SNPEFF_GENOME \
        ${samplename}.gvcf.gz > ${samplename}.germline.snpeff.ann.vcf
    """

    stub:

    """
    touch ${samplename}.germline.snpeff.ann.vcf
    touch ${samplename}.germline.snpeff.ann.csv
    touch ${samplename}.germline.snpeff.ann.html
    """
}

process somalier_extract {
    /*
    To estimate ancestry, Somalier first extracts known sites from mapped reads
    @Input:
        Mapped and pre-processed BAM file
    @Output:
        Exracted sites in (binary) somalier format
    */
    label 'process_low'

    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")
    output:
        path("output/${samplename}.somalier")
    //params:
    //    sites_vcf = config['references']['SOMALIER']['SITES_VCF'],
    //    genomeFasta = config['references']['GENOME'],
    //    rname = 'somalier_extract'
    //container: config['images']['wes_base']
    script:
    """
    mkdir -p output
    somalier extract \
        -d output \
        --sites $SITES_VCF \
        -f $GENOMEREF \
        ${samplename}.bam
    """

    stub:
    """
    mkdir -p output
    touch output/${samplename}.somalier
    """
}

process somalier_analysis_human {
    /*
    To estimate relatedness, Somalier uses extracted site information to
    compare across all samples. This step also runs the ancestry estimation
    function in Somalier.
    @Input:
        Exracted sites in (binary) somalier format for ALL samples in the cohort
    @Output:
        Separate tab-separated value (TSV) files with relatedness and ancestry outputs

    */
    label 'process_low'


    input:
        path(somalierin)

    output:
        tuple path("relatedness.pairs.tsv"), path("relatedness.samples.tsv"),
        path("ancestry.somalier-ancestry.tsv"), path("predicted.genders.tsv"),
        path("predicted.pairs.tsv"),
        path("sampleAncestryPCAPlot.html"),
        path("predictedPairsAncestry.pdf")

    script:
    """
    echo "Estimating relatedness"
    somalier relate \
        -o "relatedness" \
        $somalierin

    echo "Estimating ancestry"
    somalier ancestry \
        -o "ancestry" \
        --labels $ANCESTRY_DB/ancestry-labels-1kg.tsv \
        $ANCESTRY_DB/*.somalier ++ \
        $somalierin

    Rscript $SCRIPT_PATH_GENDER \
        relatedness.samples.tsv \
        predicted.genders.tsv

    Rscript $SCRIPT_PATH_SAMPLES \
        relatedness.pairs.tsv \
        predicted.pairs.tsv

    Rscript $SCRIPT_PATH_PCA \
        ancestry.somalier-ancestry.tsv \
        predicted.pairs.tsv \
        sampleAncestryPCAPlot.html \
        predictedPairsAncestry.pdf
    """

    stub:

    """
    touch relatedness.pairs.tsv
    touch relatedness.samples.tsv
    touch ancestry.somalier-ancestry.tsv predicted.genders.tsv
    touch predicted.pairs.tsv sampleAncestryPCAPlot.html
    touch predictedPairsAncestry.pdf
    """
}

process somalier_analysis_mouse {
    /*
    To estimate relatedness, Somalier uses extracted site information to
    compare across all samples. This step also runs the ancestry estimation
    function in Somalier.
    @Input:
        Exracted sites in (binary) somalier format for ALL samples in the cohort
    @Output:
        Separate tab-separated value (TSV) files with relatedness and ancestry outputs

    */
    label 'process_low'

    input:
        path(somalierin)

    output:
        tuple path("relatedness.pairs.tsv"),
        path("relatedness.samples.tsv"),
        path("predicted.genders.tsv"),
        path("predicted.pairs.tsv")

    script:
    """
    echo "Estimating relatedness"
    somalier relate \
        -o "relatedness" \
        $somalierin

    Rscript $SCRIPT_PATH_GENDER \
        relatedness.samples.tsv \
        predicted.genders.tsv

    Rscript $SCRIPT_PATH_SAMPLES \
        relatedness.pairs.tsv \
        predicted.pairs.tsv

    """

    stub:

    """
    touch relatedness.pairs.tsv
    touch relatedness.samples.tsv
    touch predicted.genders.tsv
    touch predicted.pairs.tsv

    """
}

process multiqc {

    """
    Reporting step to aggregate sample summary statistics and quality-control
    information across all samples. This will be one of the last steps of the
    pipeline. The inputs listed here are to ensure that this step runs last.
    During runtime, MultiQC will recursively crawl through the working directory
    and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """

    input:
        path(allqcin)

    output:
        path("MultiQC_Report.html")

    script:

    """
    multiqc . \
    -f --interactive \
    -n "MultiQC_Report.html" \
    """

    stub:

    """
    touch MultiQC_Report.html
    """
}
