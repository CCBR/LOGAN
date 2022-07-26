#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )


//params.outputdir = "${workflow.projectDir}"
//params.inputdir = "${workflow.inputDir}"
GENOME=file(params.genome)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"// file(params.gold_indels1) //
SHAPEITINDEL=file(params.shapeitindel) //params.shapeitindel =  "/data/OpenOmics/references/genome-seek/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz" //file(params.gold_indels2) //
DBSNP=file(params.dbsnp) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //= '/data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz' // /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 


baminput= Channel.fromPath(params.inputdir, checkIfExists: true, type: 'file')



workflow {
    applybqsr(baminput)
}



process applybqsr {
    /*
    Base quality recalibration for all samples to 
    */    
    input:
        path baminput 

    output:
        path("${baminput.simpleName}.bqsr.bam")

    script:
    """

    gatk --java-options '-Xmx32g' ApplyBQSR \
        --reference ${GENOME} \
        --input ${baminput} \
        --bqsr-recal-file  /data/CCBR/projects/ccbr1160/nf/testnf/recal/${baminput.simpleName}.recal_data.grp \
        --output ${baminput.simpleName}.bqsr.bam \
        --use-jdk-inflater \
        --use-jdk-deflater

    """
}



