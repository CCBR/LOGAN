#!/usr/bin/env nextflow

/*
QC Processes for WGS Data
1- FC LAne
2- FQ Screen
3-Kraken
4-fastqc_Bam ##If starig from BAM
5-RERUN Ttargets bed
6-QUalimap QC
7-samtools flgastats
8-vcftools
9-collect vc metrics
10-bcftoosl stats
11-gatk-varianteval
12snpeff
13-somalier extract
14-multiqc~~
*/



process fc_lane {
    tag { name }
    
    module 'python/3.8'

    input:

    output: 
        tuple val(samplename),
        path("${samplename}.fc_lane")    
    script:

"""
    python {params.get_flowcell_lanes} \
    {input.r1} \
    $ > {output.txt}
"""
}


process fastq_screen {
    tag { name }
    module ['fastq_screen/0.15.2','bowtie/2-2.4.5']
    fastq_screen_conf=params.fastq_screen_conf

    input:
        tuple val(samplename), 
        baseDir
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz")

    output:
        tuple val(samplename),
        path("${samplename}_fastq_screen")

//fastq_screen uses directory output
    script: 
    """
    fastq_screen --conf resources/fastq_screen.conf \
        --outdir ${samplename}_fastq_screen \
        --threads 8 \
        --subset 1000000 \
        --aligner bowtie2 \
        --force \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz 
    """

}


process kraken{
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
        path("${samplename}.trimmed.kraken_bacteria.krona.html"),


    script:
    """
    kraken2 --db ${{tmp}}/${{kdb_base}} \
        --threads {threads} --report {output.taxa} \
        --output {output.out} \
        --gzip-compressed \
        --paired {input.fq1} {input.fq2}

    # Generate Krona Report
    cut -f2,3 {output.out} | \
        ktImportTaxonomy - -o {output.html}
    
    """
}
