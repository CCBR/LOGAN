#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
C A L L I N G S  -  N F    v 2.1 
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
results  : $params.results
"""


// ##MUTECT2 part only using the BQSR s

GENOME = file(params.genome)
GOLD1 = file(params.gold_indels1)
GOLD2 = file(params.gold_indels2)
KNOWNS= file(params.known)
GERMLINE = file{params.germline} // /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON = file{params.pom} ///data/CCBR_Pipeliner/Exome-seek/hg38/PON/hg38.noCOSMIC_ClinVar.pon.vcf.gz

inputFile = file(params.infile)
inputFileHeader = params.infile_header

// Read Pairs.tsv file for Tumor/normal pairings
sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> 
                                    // Optionally allow user to specify a bam location
                                    // if (bam_folder != "") {
                                    row.bam = "${bam_folder}/${row.bam}"
                                    // }           
                                    [row.strain,
                                     file("${row.bam}", checkIfExists: true),
                                     file("${row.bam}.bai", checkIfExists: true)] }

//Channel.from(inputFile)
       //.splitCsv(sep: '\t', header: inputFileHeader)
   //       .into { readPairsFastQC; readPairsFastqToSam }
//#

//### Start with BQSR to split the file sinto the correct intervals
//###
process bqsr {
    
     tag { "${strain}:${region}" }

    input:
    val(chr)
    path(vcf)
    path(vcf_index)

    output:
    tuple path("${prefix}.${chr}.bam"), path("${prefix}.${chr}.vcf.gz.tbi")

    
    script:

    """
    gatk --java-options '-Xmx{params.memory}g' BaseRecalibrator \
    --input ${inputBAM} \
    --reference ${GENOME} \
    ${params.knowns} \
    --output {output.recal} \
    ${intervals}

    gatk --java-options '-Xmx48g' ApplyBQSR \
        --reference ${GENOME} \
        --input ${inputBAM} \
        --bqsr-recal-file {output.re} \
        --output {output.bam} \
        --use-jdk-inflater \
        --use-jdk-deflater

    """
}


process mutect2 {
    
     tag { "${strain}:${region}" }

    input:
    val(chr)
    path(vcf)
    path(vcf_index)

    output:
    tuple path("${prefix}.${chr}.bam"), path("${prefix}.${chr}.vcf.gz.tbi")

    
    script:

    """
    java -Xmx8G -jar \$GATK Mutect2 \
    --reference params.GENOME \
    --intervals ${splitbed} \
    --input ${tumor_recal.bam} \
    --input ${normal_recal.bam} \
    --normal-sample ${normal_name} \
    --tumor-sample ${tumor_name} \
    --native-pair-hmm-threads 1 \
    --germline-resource /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals /data/CCBR_Pipeliner/Exome-seek/hg38/PON/hg38.noCOSMIC_ClinVar.pon.vcf.gz \
    --output ${prefix}_${intervalbed}.vcf.gz
    
    java -Xmx8G -jar \$GATK FilterMutectCalls \
    --variant ${prefix}_${intervalbed}.vcf.gz \
    --reference /a97664a825fa4983ab480a796f46136d/Homo_sapiens_assembly38.fasta \
    --output ${prefix}_${intervalbed_final}.bam
    """


}


proces merge_mutect2 {
    
    input:
    val(chr)
    path(vcf)
    path(vcf_index)

    output:
    tuple path("${prefix}.${intervalbed}.vcf.tz"), path("${prefix}.${intervalbed_final}.vcf.gz")

    
    script:
    singularity exec --cleanenv  --bind /data/CCBR_Pipeliner/Exome-seek/hg38/genome:/b741c14d026c477f8c567de54bf500ae --bind /data/desaipa/ccbr_analysis/WGS/pipeliner/variant:/7c37dbfa13a44b8ab86ae259450a181d --bind /data/CCBR_Pipeliner/Exome-seek/hg38/ docker://lethalfang/somaticseq:3.7.3 \
    concat.py --bgzip-output -infiles \

}
