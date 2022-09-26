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

script:
//Use Python Container here
"""
python {params.get_flowcell_lanes} \\
        {input.r1} \\
        {wildcards.samples} > {output.txt}
"""
}


process fastq_screen {

script: 
    """
    fastq_screen --conf {params.fastq_screen_config} \
        --outdir {params.outdir} \
        --threads {threads} \
        --subset 1000000 \
        --aligner bowtie2 \
        --force \
        {input.fq1} {input.fq2}
    """

}


process kraken{


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
