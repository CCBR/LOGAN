
process fastp{
    tag { name }
    module 'fastp/0.23.2'

    input:
    tuple val(samplename), path(fqs)
    
    output:
    tuple val(samplename),
    path("${samplename}.R1.trimmed.fastq.gz"),
    path("${samplename}.R2.trimmed.fastq.gz"),
    path("${samplename}.fastp.json"),
    path("${samplename}.fastp.html")

    script:
    """
    fastp -w 4 \
        --detect_adapter_for_pe \
        --in1 ${fqs[0]} \
        --in2 ${fqs[1]} \
        --out1 ${samplename}.R1.trimmed.fastq.gz \
        --out2 ${samplename}.R2.trimmed.fastq.gz  \
        --json ${samplename}.fastp.json \
        --html ${samplename}.fastp.html
    """

}


process bwamem2{
    tag { name }
    module=['bwa-mem2/2.2.1','samblaster/0.1.26','samtools/1.15.1']
    input:
        tuple val(samplename), 
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz"),
        path("${samplename}.fastp.json"),
        path("${samplename}.fastp.html")
        
    output:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")

    script:
    //BWAmem2/samblaster/samtools sort for marking duplicates;
    """

     bwa-mem2 mem -M \
        -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:illumina\\tLB:${samplename}\\tPU:${samplename}\\tCN:hgsc\\tDS:wes' \
        -t 12 \
        ${GENOME} \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -@12 -m 4G - -o ${samplename}.bam

    samtools index -@ 8 ${samplename}.bam ${samplename}.bai

    """
}

process indelrealign {
    tag { name }
    module=['GATK/3.8-1']
    
    input:
    tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")

    output:
    tuple val(samplename), path("${samplename}.ir.bam")

    script: 
    
    """
    /usr/bin/java -Xmx32g -jar \${GATK_JAR} -T RealignerTargetCreator \
        -I  ${samplename}.bam \
        -R ${GENOME} \
        -o ${samplename}.intervals \
        -known ${MILLSINDEL} -known ${SHAPEITINDEL} 
    
    /usr/bin/java -Xmx32g -jar \${GATK_JAR} -T IndelRealigner \
        -R ${GENOME} \
        -I ${samplename}.bam \
        -known ${MILLSINDEL} -known ${SHAPEITINDEL} \
        --use_jdk_inflater \
        --use_jdk_deflater \
        -targetIntervals ${samplename}.intervals \
        -o  ${samplename}.ir.bam
    """
    
}



process bqsr {
    /*
    Base quality recalibration for all samples to 
    */    
    input:
        tuple val(samplename), path("${samplename}.ir.bam"), path(bed)

    output:
        tuple val(samplename),path("${samplename}_${bed.simpleName}.recal_data.grp"),emit: bqsrby
        //path("${bam.simpleName}_${bed.simpleName}.recal_data.grp"), emit: bqsrby

    script:
    """
    gatk --java-options '-Xmx32g' BaseRecalibrator \
    --input ${samplename}.ir.bam \
    --reference ${GENOME} \
    --known-sites ${MILLSINDEL} --known-sites ${SHAPEITINDEL} \
    --output ${samplename}_${bed.simpleName}.recal_data.grp \
    --intervals ${bed}

    """
}

process gatherbqsr {

    input: 
        tuple val(samplename), path(recalgroups)
    output:
        tuple val(samplename), path("${samplename}.recal_data.grp")
    script:
    
    strin = recalgroups.join(" --input ")

    """
    gatk --java-options '-Xmx32g' GatherBQSRReports \
    --input ${strin} \
    --output ${samplename}.recal_data.grp

    """
}


process applybqsr {
    /*
    Base quality recalibration for all samples to 
    */    
    input:
        tuple val(samplename),path("${samplename}.ir.bam"), path("${samplename}.recal_data.grp")

    output:
        tuple val(samplename),path("${samplename}.bqsr.bam")

    script:
    """

    gatk --java-options '-Xmx32g' ApplyBQSR \
        --reference ${GENOME} \
        --input ${samplename}.ir.bam \
        --bqsr-recal-file  ${samplename}.recal_data.grp \
        --output ${samplename}.bqsr.bam \
        --use-jdk-inflater \
        --use-jdk-deflater

    """
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