#!/usr/bin/env nextflow

/*
QC Processes for WGS Data
1-FC LAne
2-FQ Screen
3-Kraken
4-fastqc_Bam ##If starig from BAM
5-RERUN Ttargets bed
6-Qualimap QC
7-samtools flagstats
8-vcftools
9-collect vc metrics
10-bcftools stats
11-gatk-varianteval
12snpeff
13-somalier extract
14-multiqc~
*/


process fc_lane {
    
    input:
        tuple val(samplename),
        path(fqs)

    output: 
        tuple val(samplename),
        path("${samplename}.fc_lane")    
    
    script:
    GET_FLOWCELL_LANES=file(params.get_flowcell_lanes)

    """
    python ${GET_FLOWCELL_LANES} \
    ${fqs[0]}  \
    ${samplename} > ${samplename}.fc_lane
    """
}

//Uses Trimmed Files
process fastq_screen {

    input:
        tuple val(samplename), 
        path(fqs)

    output:
        tuple val(samplename),
        path("${samplename}_fastq_screen")

    script: 
    
    FASTQ_SCREEN_CONF=file(params.fastq_screen_conf)

    """
    fastq_screen --conf ${FASTQ_SCREEN_CONF} \
        --outdir ${samplename}_fastq_screen \
        --threads 8 \
        --subset 1000000 \
        --aligner bowtie2 \
        --force \
        ${fqs[0]} ${fqs[1]}
    """

}

///kraken
rule kraken:
    """
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interative krona report
    """
    input:
        tuple val(samplename), 
        path(fqs)

    output:
        tuple val(samplename),
        path("${samplename}_fastq_screen"),
        path("${samplename}_fastq_screen")
        path("${samplename}_fastq_screen")
        out  = os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.out.txt"),
        taxa = os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.taxa.txt"),
        html = os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.krona.html"),
    params:
        rname  ='kraken',
        outdir = os.path.join(output_qcdir, "kraken"),
        bacdb  = config['references']['KRAKENBACDB'],
        localdisk = config['input_params']['tmpdisk']
    envmodules:
        'kraken/2.1.2', 
        'kronatools/2.8'
    container: config['images']['kraken']
    threads: 24
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    #tmp=$(mktemp -d -p "{params.localdisk}")
    #trap 'rm -rf "${{tmp}}"' EXIT
    # Copy kraken2 db to local node storage to reduce filesystem strain
    #cp -rv {params.bacdb} ${{tmp}}
    #kdb_base=$(basename {params.bacdb})
    
    kraken2 --db ${{tmp}}/${{kdb_base}} \\
        --threads {threads} --report {output.taxa} \\
        --output {output.out} \\
        --gzip-compressed \\
        --paired {input.fq1} {input.fq2}
    # Generate Krona Report
    cut -f2,3 {output.out} | \\
        ktImportTaxonomy - -o {output.html}
    """


rule qualimap_bamqc:
    """
    Quality-control step to assess various post-alignment metrics 
    and a secondary method to calculate insert size. Please see
    QualiMap's website for more information about BAM QC:
    http://qualimap.conesalab.org/
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Report containing post-aligment quality-control metrics
    """
    input:
        bam  = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
        bed=os.path.join(output_qcdir, "exome_targets.bed"),
    output: 
        txt  = os.path.join(output_qcdir,"{samples}","genome_results.txt"),
        html = os.path.join(output_qcdir,"{samples}","qualimapReport.html")
    params:
        outdir = os.path.join(output_qcdir, "{samples}"),
        rname  = "qualibam"
    message: "Running QualiMap BAM QC with {threads} threads on '{input}' input file"
    envmodules: 'qualimap/2.2.1'
    container: config['images']['qualimap']
    threads: 8
    shell: """
    unset DISPLAY
    qualimap bamqc -bam {input.bam} \\
        --java-mem-size=48G \\
        -c -ip \\
        -gff {input.bed} \\
        -outdir {params.outdir} \\
        -outformat HTML \\
        -nt {threads} \\
        --skip-duplicated \\
        -nw 500 \\
        -p NON-STRAND-SPECIFIC
    """


rule samtools_flagstats:
    """
    Quality-control step to assess alignment quality. Flagstat provides 
    counts for each of 13 categories based primarily on bit flags in the 
    FLAG field. Information on the meaning of the flags is given in the 
    SAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Text file containing alignment statistics
    """
    input:
        bam  = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
    output:
        txt  = os.path.join(output_qcdir,"{samples}.samtools_flagstat.txt")
    params: 
        rname = "samtools_flagstats"
    message: "Running SAMtools flagstat on '{input}' input file"
    envmodules: 'samtools/1.12'
    container: config['images']['wes_base']
    shell: """
    samtools flagstat {input.bam} > {output.txt}
    """


rule vcftools:
    """
    Quality-control step to calculates a measure of heterozygosity on 
    a per-individual basis. The inbreeding coefficient, F, is estimated
    for each individual using a method of moments. Please see VCFtools
    documentation for more information: 
    https://vcftools.github.io/man_latest.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a measure of heterozygosity
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
    output: 
        het = os.path.join(output_qcdir,"raw_variants.het"),
    params: 
        prefix = os.path.join(output_qcdir,"raw_variants"),
        rname  = "vcftools",
    message: "Running VCFtools on '{input.vcf}' input file"
    envmodules: 'vcftools/0.1.16'
    container: config['images']['wes_base']
    shell: """
    vcftools --gzvcf {input.vcf} --het --out {params.prefix}
    """


rule collectvariantcallmetrics:
    """
    Quality-control step to collect summary metrics about snps and indels
    called in a multisample VCF file. Please see the Broad's documentation
    for more information about each field in the generated log file:
    https://broadinstitute.github.io/picard/picard-metric-definitions.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a collection of metrics relating to snps and indels 
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
    output: 
        metrics = os.path.join(output_qcdir,"raw_variants.variant_calling_detail_metrics"),
    params: 
        dbsnp=config['references']['DBSNP'],
        prefix = os.path.join(output_qcdir,"raw_variants"),
        rname="varcallmetrics",
    message: "Running Picard CollectVariantCallingMetrics on '{input.vcf}' input file"
    envmodules: 'picard/2.20.8'
    container: config['images']['picard']
    shell: """
    java -Xmx24g -jar ${{PICARDJARPATH}}/picard.jar \\
        CollectVariantCallingMetrics \\
        INPUT={input.vcf} \\
        OUTPUT={params.prefix} \\
        DBSNP={params.dbsnp} Validation_Stringency=SILENT
    """


rule bcftools_stats:
    """
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
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz"),
    output: 
        txt = os.path.join(output_qcdir,"{samples}.germline.bcftools_stats.txt"),
    params: 
        rname="bcfstats",
    message: "Running BCFtools on '{input.vcf}' input file"
    envmodules: 'bcftools/1.9'
    container: config['images']['wes_base']
    shell: """
    bcftools stats {input.vcf} > {output.txt}
    """


rule gatk_varianteval:
    """
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
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz"), 
    output: 
        grp = os.path.join(output_qcdir,"{samples}.germline.eval.grp"),
    params:
        rname    = "vareval",
        genome   = config['references']['GENOME'],
        dbsnp    = config['references']['DBSNP'],
        ver_gatk = config['tools']['gatk4']['version']
    message: "Running GATK4 VariantEval on '{input.vcf}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    threads: 16
    shell: """
    gatk --java-options '-Xmx12g -XX:ParallelGCThreads={threads}' VariantEval \\
        -R {params.genome} \\
        -O {output.grp} \\
        --dbsnp {params.dbsnp} \\
        --eval {input.vcf} 
    """


rule snpeff:
    """
    Data processing and quality-control step to annotate variants, predict its
    functional effects, and collect various summary statistics about variants and
    their annotations. Please see SnpEff's documentation for more information: 
    https://pcingola.github.io/SnpEff/
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    """
    input:  
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz")
    output: 
        vcf  = os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.vcf"),
        csv  = os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.csv"),
        html = os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.html"),
    params: 
        rname  = "snpeff",
        genome = config['references']['SNPEFF_GENOME'],
        config = config['references']['SNPEFF_CONFIG'],
        bundle = config['references']['SNPEFF_BUNDLE'],
    envmodules: 'snpEff/4.3t'
    container: config['images']['wes_base']
    shell: """
    java -Xmx12g -jar $SNPEFF_JAR \\
        -v -canon -c {params.config} \\
        -csvstats {output.csv} \\
        -stats {output.html} \\
        {params.genome} \\
        {input.vcf} > {output.vcf}
    """

rule somalier_extract:
    """
    To estimate ancestry, Somalier first extracts known sites from mapped reads
    @Input:
        Mapped and pre-processed BAM file
    @Output:
        Exracted sites in (binary) somalier format
    """
    input:
        bam = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
        bai = os.path.join(output_bamdir,"final_bams","{samples}.bai"),
    output: 
        somalierOut = os.path.join(output_germline_base,"somalier","{samples}.somalier")
    params:
        sites_vcf = config['references']['SOMALIER']['SITES_VCF'],
        genomeFasta = config['references']['GENOME'],
        rname = 'somalier_extract'
    container: config['images']['wes_base']
    shell: """ 
    echo "Extracting sites to estimate ancestry"
    somalier extract \\
        -d "$(dirname {output.somalierOut})" \\
        --sites {params.sites_vcf} \\
        -f {params.genomeFasta} \\
        {input.bam}
    """


rule somalier_analysis:
    """
    To estimate relatedness, Somalier uses extracted site information to
    compare across all samples. This step also runs the ancestry estimation
    function in Somalier.
    @Input:
        Exracted sites in (binary) somalier format for ALL samples in the cohort
    @Output:
        Separate tab-separated value (TSV) files with relatedness and ancestry outputs
    """
    input:
        somalier=expand(os.path.join(output_germline_base,"somalier","{samples}.somalier"), samples=samples),
    output:
        relatedness = os.path.join(output_germline_base,"somalier","relatedness.pairs.tsv"),
        relatednessSamples = os.path.join(output_germline_base,"somalier","relatedness.samples.tsv"),
        ancestry = os.path.join(output_germline_base,"somalier","ancestry.somalier-ancestry.tsv"),
        finalFileGender = os.path.join(output_germline_base,"predicted.genders.tsv"),
        finalFilePairs = os.path.join(output_germline_base,"predicted.pairs.tsv"),
        ancestoryPlot = os.path.join(output_germline_base,"sampleAncestryPCAPlot.html"),
        pairAncestoryHist = os.path.join(output_germline_base,"predictedPairsAncestry.pdf"),
    params:
        ancestry_db = config['references']['SOMALIER']['ANCESTRY_DB'],
        sites_vcf = config['references']['SOMALIER']['SITES_VCF'],
        genomeFasta = config['references']['GENOME'],
        script_path_gender = config['scripts']['genderPrediction'],
        script_path_samples = config['scripts']['combineSamples'],
        script_path_pca = config['scripts']['ancestry'],
        rname = 'somalier_analysis'
    container: config['images']['wes_base']
    shell: """ 
    echo "Estimating relatedness"
    somalier relate \\
        -o "$(dirname {output.relatedness})/relatedness" \\
        {input.somalier}
    
    echo "Estimating ancestry"
    somalier ancestry \\
        -o "$(dirname {output.relatedness})/ancestry" \\
        --labels {params.ancestry_db}/ancestry-labels-1kg.tsv \\
        {params.ancestry_db}/*.somalier ++ \\
        {input.somalier}
    Rscript {params.script_path_gender} \\
        {output.relatednessSamples} \\
        {output.finalFileGender}    
    
    Rscript {params.script_path_samples} \\
        {output.relatedness} \\
        {output.finalFilePairs}
    
    Rscript {params.script_path_pca} \\
        {output.ancestry} \\
        {output.finalFilePairs} \\
        {output.ancestoryPlot} \\
        {output.pairAncestoryHist}
    """



rule multiqc:
    """
    Reporting step to aggregate sample summary statistics and quality-control
    information across all samples. This will be one of the last steps of the 
    pipeline. The inputs listed here are to ensure that this step runs last. 
    During runtime, MultiQC will recurively crawl through the working directory
    and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """
    input:  
        expand(os.path.join(output_fqdir,"{samples}.fastq.info.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"FQscreen","{samples}.R2.trimmed_screen.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.krona.html"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}_fastqc.zip"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}","genome_results.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.samtools_flagstat.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.bcftools_stats.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.eval.grp"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.html"), samples=samples),
        os.path.join(output_qcdir,"raw_variants.het"), 
        os.path.join(output_qcdir,"raw_variants.variant_calling_detail_metrics"),
        os.path.join(output_germline_base,"somalier","ancestry.somalier-ancestry.tsv"),
    output: 
        report  = os.path.join(output_qcdir,"finalQC","MultiQC_Report.html"),
    params: 
        rname  = "multiqc",
        workdir = os.path.join(BASEDIR)
    envmodules: 'multiqc/1.11'
    container: config['images']['multiqc']
    shell: """
    multiqc --ignore '*/.singularity/*' \\
        --ignore '*/*/*/*/*/*/*/*/pyflow.data/*' \\
        --ignore 'slurmfiles/' \\
        -f --interactive \\
        -n {output.report} \\
        {params.workdir}
    """