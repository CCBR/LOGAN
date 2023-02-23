
get_flowcell_lanes = os.path.join("workflow", "scripts", "get_flowcell_lanes.py"),

process fastq_screen {
    tag { name }
    module=['fastq_screen/0.15.2','bowtie/2-2.4.5']

    input:
        tuple val(samplename), 
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz")

    output:
        tuple val(samplename),
        path("${samplename}_fastq_screen")


    script: 
    """
    fastq_screen --conf ${FASTQ_SCREEN_CONF} \
        --outdir ${samplename}_fastq_screen \
        --threads 8 \
        --subset 1000000 \
        --aligner bowtie2 \
        --force \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz 
    """

}

/*
process kraken {
//Kraken and Krona Rerpot
    tag { name }
    module=['kraken/2.1.2','kronatools/2.8']

    input:
        tuple val(samplename), 
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz")

    output:
        tuple val(samplename),
        path("${samplename}.trimmed.kraken_bacteria.out.txt"),
        path("${samplename}.trimmed.kraken_bacteria.taxa.txt"),
        path("${samplename}.trimmed.kraken_bacteria.krona.html")


    script:
    """
    kraken2 --db \${tmp}/\${kdb_base} \
        --threads {threads} --report {output.taxa} \
        --output {output.out} \
        --gzip-compressed \
        --paired {input.fq1} {input.fq2}

    # Generate Krona Report
    cut -f2,3 {output.out} | \
        ktImportTaxonomy - -o {output.html}
    
    """
}
*/

process qualimap_bam{
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

modules: 'qualimap/2.2.1'
//container: config['images']['qualimap']
threads: 8
script:
 """
 
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

}
rule qualimap_bamqc:

   
    params:
        outdir = os.path.join(output_qcdir, "{samples}"),
        rname  = "qualibam"
    message: "Running QualiMap BAM QC with {threads} threads on '{input}' input file"
  

  process multiqc{

    input:

    output:
        path("QC/MultiQC_Report.html")

    script:
        multiqc --ignore '*/.singularity/*' \\
        --ignore 'slurmfiles/' \\
        --exclude peddy \\
        -f --interactive \\
        -n {output.report} \\
        {params.workdir}
    
  }