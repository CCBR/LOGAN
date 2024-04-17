GENOMEREF = file(params.genomes[params.genome].genome)
KNOWNRECAL = params.genomes[params.genome].KNOWNRECAL
KNOWNINDELS = params.genomes[params.genome].KNOWNINDELS

process fastp {
    container = "${params.containers.logan}"
    label 'process_medium'
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
    fastp -w $task.cpus \
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
    container = "${params.containers.logan}"
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
    """

     bwa-mem2 mem -M \
        -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:illumina\\tLB:${samplename}\\tPU:${samplename}\\tCN:hgsc\\tDS:wgs' \
        -t $task.cpus \
        ${GENOMEREF} \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -@ $task.cpus -m 4G - -o ${samplename}.bam

    samtools index -@ $task.cpus ${samplename}.bam ${samplename}.bai

    """

    stub:
    """
    touch ${samplename}.bam ${samplename}.bai
    """
}



process indelrealign {
    container "${params.containers.logan}"
    label 'process_long'

    input:
    tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")

    output:
    tuple val(samplename), path("${samplename}.ir.bam"), path("${samplename}.ir.bai")

    script:

    """
    /usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx16g -jar \$GATK_JAR -T RealignerTargetCreator \
        -I ${samplename}.bam \
        -R ${GENOMEREF} \
        -o ${samplename}.intervals \
        -nt $task.cpus \
        ${KNOWNINDELS}

    /usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx16g -jar \$GATK_JAR -T IndelRealigner \
        -R ${GENOMEREF} \
        -I ${samplename}.bam \
        ${KNOWNINDELS} \
        -targetIntervals ${samplename}.intervals \
        -o ${samplename}.ir.bam
    """

    stub:
    """
    touch ${samplename}.ir.bam ${samplename}.ir.bai
    """

}


process bqsr_ir {
    /*
    Base quality recalibration for all samples
    */
    container = "${params.containers.logan}"
    label 'process_low'
    input:
        tuple val(samplename), path("${samplename}.ir.bam"), path("${samplename}.ir.bai"), path(bed)

    output:
        tuple val(samplename), path("${samplename}_${bed.simpleName}.recal_data.grp")

    script:
    """
    gatk --java-options '-Xmx16g' BaseRecalibrator \
    --input ${samplename}.ir.bam \
    --reference ${GENOMEREF} \
    ${KNOWNRECAL} \
    --output ${samplename}_${bed.simpleName}.recal_data.grp \
    --intervals ${bed}
    """

    stub:
    """
    touch ${samplename}_${bed.simpleName}.recal_data.grp
    """
}

process bqsr {
    /*
    Base quality recalibration for all samples
    */
    container = "${params.containers.logan}"
    label 'process_low'
    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai"), path(bed)

    output:
        tuple val(samplename),path("${samplename}_${bed.simpleName}.recal_data.grp"), emit: bqsrby

    script:
    """
    gatk --java-options '-Xmx16g' BaseRecalibrator \
    --input ${samplename}.bam \
    --reference ${GENOMEREF} \
    ${KNOWNRECAL} \
    --output ${samplename}_${bed.simpleName}.recal_data.grp \
    --intervals ${bed}
    """

    stub:
    """
    touch ${samplename}_${bed.simpleName}.recal_data.grp
    """
}

process gatherbqsr {
    container = "${params.containers.logan}"
    label 'process_low'
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
    container = "${params.containers.logan}"
    label 'process_long'

    input:
        tuple val(samplename), path(bam), path(bai), path("${samplename}.recal_data.grp")

    output:
        tuple val(samplename), path("${samplename}.bqsr.bam"),  path("${samplename}.bqsr.bai")

    script:

    """
    gatk --java-options '-Xmx32g' ApplyBQSR \
        --reference ${GENOMEREF} \
        --input ${bam} \
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
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
    tuple val(bamname), path(bam)

    output:
    tuple val(bamname), path(bam), path("${bam}.bai")

    script:
    """
    samtools index -@ $task.cpus ${bam} ${bam}.bai
    """

    stub:
    """
    touch ${bam}.bai
    """

}

//Save to CRAM for output
process bamtocram_tonly {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        path("${tumorname}.cram"), path("${tumorname}.cram.crai")


    script:
    """
        samtools view -@ $task.cpus -C -T $GENOMEREF -o ${sample}.cram $tumor
        samtools index ${tumorname}.cram -@ $task.cpus
    """
    
    stub:
    """
    touch ${tumorname}.cram ${tumorname}.cram.crai
    """
}




