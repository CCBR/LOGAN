GENOME=file(params.genome)
GENOMEDICT=file(params.genomedict)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
SHAPEITINDEL=file(params.shapeitindel) //ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
KGP=file(params.kgp) //1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.dbsnp) //dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 
outdir=file(params.output)


process fastp{
    tag { name }

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

    stub:
    """
    touch ${samplename}.R1.trimmed.fastq.gz
    touch ${samplename}.R2.trimmed.fastq.gz
    touch ${samplename}.fastp.json
    touch ${samplename}.fastp.html

    """
}


process bwamem2 {
    tag { name }
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
        -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:illumina\\tLB:${samplename}\\tPU:${samplename}\\tCN:hgsc\\tDS:wgs' \
        -t 16 \
        ${GENOME} \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -@12 -m 4G - -o ${samplename}.bam

    samtools index -@ 8 ${samplename}.bam ${samplename}.bai

    """

    stub:
    """
    touch ${samplename}.bam ${samplename}.bai
    """
}

process indelrealign {
    /*
    Briefly, RealignerTargetCreator runs faster with increasing -nt threads, 
    while IndelRealigner shows diminishing returns for increasing scatter
    */
    tag { name }
    
    input:
    tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")

    output:
    tuple val(samplename), path("${samplename}.ir.bam")

    script: 
    
    """
    /usr/bin/java -Xmx32g -jar \${GATK_JAR} -T RealignerTargetCreator \
        -I ${samplename}.bam \
        -R ${GENOME} \
        -o ${samplename}.intervals \
        -nt 16 \
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
    

    stub:
    """
    touch ${samplename}.ir.bam 
    """

}



process bqsr {
    /*
    Base quality recalibration for all samples 
    */    
    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai"), path(bed)

    output:
        tuple val(samplename),path("${samplename}_${bed.simpleName}.recal_data.grp"),emit: bqsrby
        //path("${bam.simpleName}_${bed.simpleName}.recal_data.grp"), emit: bqsrby

    script:
    """
    gatk --java-options '-Xmx32g' BaseRecalibrator \
    --input ${samplename}.bam \
    --reference ${GENOME} \
    --known-sites ${MILLSINDEL} --known-sites ${SHAPEITINDEL} \
    --output ${samplename}_${bed.simpleName}.recal_data.grp \
    --intervals ${bed}
    """

    stub:
    """
    touch ${samplename}_${bed.simpleName}.recal_data.grp 
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

    stub:
    """
    touch ${samplename}.recal_data.grp
    """
}


process applybqsr {
    /*
    Base quality recalibration for all samples to 
    */   
    publishDir(path: "${outdir}/bams/BQSR", mode: 'copy') 

    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai"), path("${samplename}.recal_data.grp")

    output:
        tuple val(samplename), path("${samplename}.bqsr.bam"),  path("${samplename}.bqsr.bai")

    script:
    """

    gatk --java-options '-Xmx32g' ApplyBQSR \
        --reference ${GENOME} \
        --input ${samplename}.bam \
        --bqsr-recal-file ${samplename}.recal_data.grp \
        --output ${samplename}.bqsr.bam \
        --use-jdk-inflater \
        --use-jdk-deflater

    """

    stub:
    """
    touch ${samplename}.bqsr.bam ${samplename}.bqsr.bai
    """

}



process samtoolsindex {
    publishDir(path: "${outdir}/bams/BQSR", mode: 'copy') 
    
    input:
    tuple val(bamname), path(bam)
    
    output:
    tuple val(bamname), path(bam), path("${bam}.bai")

    script:
    """
    samtools index -@ 4 ${bam} ${bam}.bai
    """

    stub:
    """
    touch ${bam}.bai
    """

}

//Save to CRAM for output and publish
process bamtocram_tonly{
    
    input: 
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        path("${sample}.cram")

    script:
    """
        samtools view -@ 4 -C -T $GENOME -o ${sample}.cram {$tumor}.bam
    """
}