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



//baminput= Channel.fromPath(params.inputdir, checkIfExists: true, type: 'file')
                //.map { row -> tuple(
                        //row.simpleName,row
                       //)}

fastqinput=Channel.fromFilePairs(params.fastqs,checkIfExists: true,type:'file')
intervalbed = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor,
                       row.Normal
                       )
                                  }


workflow {
    fastp(fastqinput) 
}

process fastp{
  input:
    tuple val(samplename), path(fqs)
    
    output:
    tuple val(samplename)
    path({$fqs[0].simpleName}.R1.trimmed.fastq.gz)
    path({$fqs[0].simpleName}.R2.trimmed.fastq.gz)
    path({$fqs[0].simpleName}.fastp.json)
    path({$fqs[0].simpleName}.fastp.html)

    script:
    """
    fastp -w 4 \
        --detect_adapter_for_pe \
        --in1 ${fqs[0]} \
        --in2 ${fqs[1]} \
        --out1 ${fqs[0].simpleName}.R1.trimmed.fastq.gz \
        --out2 ${fqs[0].simpleName}.R2.trimmed.fastq.gz  \
        --json ${fqs[0].simpleName}.fastp.json \
        --html ${fqs[0].simpleName}.fastp.html
    """

}



process bwamem2{
    input:
        tuple val(samplename)
        path({$outfq.simpleName}.R1.trimmed.fastq.gz)
        path({$outfq.simpleName}.R2.trimmed.fastq.gz)
        path({$outfq.simpleName}.fastp.json)
        path({$outfq.simpleName}.fastp.html)

    output:
        path(outbam)

    script:
    adfadfj;ajdfkl;jakldsf
    skl;fkl;
}

process realign {
    input:
    path(bam)

    output:
    path(outbam)

    script: 
    
    """
    java -Xmx{params.memory}g -jar ${{GATK_JAR}} -T RealignerTargetCreator \\
        -I {input.bam} \\
        -R {params.genome} \\
        -o {output.intervals} \\
        {params.knowns}
    
    java -Xmx{params.memory}g -jar ${{GATK_JAR}} -T IndelRealigner \\
        -R {params.genome} \\
        -I {input.bam} \\
        {params.knowns} \\
        --use_jdk_inflater \\
        --use_jdk_deflater \\
        -targetIntervals {output.intervals} \\
        -o {output.bam}
    """
    
}