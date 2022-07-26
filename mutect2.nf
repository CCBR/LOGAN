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
                .map { row -> tuple(
                        row.simpleName,row
                       )}
intervalbed = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor,
                       row.Normal
                       )
                                  }




//bambyinterval=bamwithsample.combine(intervalbed).view()

workflow {
    //Create index and then 
    
    samtoolsindex(baminput)
    bamwithsample=samtoolsindex.out.join(sample_sheet).map{it.swap(3,0)}.join(samtoolsindex.out).map{it.swap(3,0)}
    bambyinterval=bamwithsample.combine(intervalbed)
    mutect2(bambyinterval)
    mutect2_tonly(bambyinterval)
}

process samtoolsindex{
  input:
    tuple val(bamname), path(bam)
    
    output:
    tuple val(bamname), path(bam), path("${bam}.bai")

    
    script:
    """
    samtools index -@ 4 ${bam} ${bam}.bai
    """

}


process mutect2 {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        path("${tumor.simpleName}.${bed.simpleName}.mut2.vcf.gz")

    
    script:

    """
    gatk Mutect2 \
    --reference ${GENOME} \
    --intervals ${bed} \
    --input ${tumor} \
    --input ${normal} \
    --normal-sample ${normal.simpleName} \
    --tumor-sample ${tumor.simpleName} \
    --native-pair-hmm-threads 1 \
    --germline-resource ${GNOMAD} \
    --panel-of-normals ${PON} \
    --output ${tumor.simpleName}.${bed.simpleName}.mut2.vcf.gz
    """
}



process mutect2_tonly {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        path("${tumor.simpleName}.${bed.simpleName}.tonly.mut2.vcf.gz")

    
    script:

    """
    gatk Mutect2 \
    --reference ${GENOME} \
    --intervals ${bed} \
    --input ${tumor} \
    --tumor-sample ${tumor.simpleName} \
    --native-pair-hmm-threads 1 \
    --germline-resource ${GNOMAD} \
    --panel-of-normals ${PON} \
    --output ${tumor.simpleName}.${bed.simpleName}.tonly.mut2.vcf.gz
    """
}
