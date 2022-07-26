#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )


GENOME=file(params.genome)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"// file(params.gold_indels1) //
SHAPEITINDEL=file(params.shapeitindel) //params.shapeitindel =  "/data/OpenOmics/references/genome-seek/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz" //file(params.gold_indels2) //
DBSNP=file(params.dbsnp) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //= '/data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz' // /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 


log.info """\
        W G S  C C B R   P I P E L I N E
         ===================================
         genome: ${GENOME.baseName}
         """
         .stripIndent()

baminput= Channel.fromPath(params.inputdir, checkIfExists: true, type: 'file')
intervalbed = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')


//Turn Interval into name and then path 
//allsamples=combine(Channel.from("${params.sample_sheet}")).view()
//intervalbed.view()
//allsamples.view()
//sample_sheet.view()
//inbam = Channel.fromPath(params.inputdir).map {
  //  tuple( it.name.split('.realign.bam')[0], it )
//}


bambyinterval=baminput.combine(intervalbed).view()
bamb1=bambyinterval.groupTuple().view()

workflow {
    bqsr(bambyinterval) 
    
}


process bqsr {
    /*
    Base quality recalibration for all samples to 
    */    
    input:
        tuple path(bam),  path(bed)

    output:
          tuple path(bam), path(bed),emit: bqsrby
          //path("${bam.simpleName}_${bed.simpleName}.recal_data.grp")
        //path("${bam.simpleName}_${bed.simpleName}.recal_data.grp"), emit: bqsrby

    script:
    """
    gatk --java-options '-Xmx32g' BaseRecalibrator \
    --input ${bam} \
    --reference ${GENOME} \
    --known-sites ${MILLSINDEL} --known-sites ${SHAPEITINDEL} \
    --output ${bam.simpleName}_${bed.simpleName}.recal_data.grp \
    --intervals ${bed}
    """
}

