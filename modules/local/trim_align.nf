GENOMEREF = file(params.genomes[params.genome].genome)
KNOWNRECAL = params.genomes[params.genome].KNOWNRECAL
KNOWNINDELS = params.genomes[params.genome].KNOWNINDELS



process indelrealign {
    container "${params.containers.logan}"
    label 'process_long'

    input:
    tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bam.bai")

    output:
    tuple val(samplename), path("${samplename}.ir.bam"), path("${samplename}.ir.bam.bai")

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
    touch ${samplename}.ir.bam ${samplename}.ir.bam.bai
    """

}


process bqsr_ir {
    /*
    Base quality recalibration for all samples
    */
    errorStrategy 'ignore'
    container = "${params.containers.logan}"
    label 'process_low'
    input:
        tuple val(samplename), path("${samplename}.ir.bam"), path("${samplename}.ir.bam.bai"), path(bed)

    output:
        tuple val(samplename), path("${samplename}_${bed.simpleName}.recal_data.grp")

    script:
    """
    gatk --java-options '-Xmx10g' BaseRecalibrator \
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
    errorStrategy 'ignore'
    container = "${params.containers.logan}"
    label 'process_low'
    input:
        tuple val(samplename), path(bam), path(bai), path(bed)

    output:
        tuple val(samplename), path("${samplename}_${bed.simpleName}.recal_data.grp"), optional: true

    script:
    """
    gatk --java-options '-Xmx10g' BaseRecalibrator \
    --input ${bam} \
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
    gatk --java-options '-Xmx10g' GatherBQSRReports \
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
    label 'process_highmem'

    input:
        tuple val(samplename), path(bam), path(bai), path("${samplename}.recal_data.grp")

    output:
        tuple val(samplename), path("${samplename}.bqsr.bam"), path("${samplename}.bqsr.bam.bai")

    script:

    """
    gatk --java-options '-Xmx20g' ApplyBQSR \
        --reference ${GENOMEREF} \
        --input ${bam} \
        --bqsr-recal-file ${samplename}.recal_data.grp \
        --output ${samplename}.bqsr.bam \
        --use-jdk-inflater \
        --use-jdk-deflater
    mv ${samplename}.bqsr.bai ${samplename}.bqsr.bam.bai
    """

    stub:

    """
    touch ${samplename}.bqsr.bam ${samplename}.bqsr.bam.bai
    """

}


process samtoolsindex {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
    tuple val(bamname), path(bam)

    output:
    tuple val(bamname), path(bam), path("${bam.simpleName}.bam.bai")

    script:
    """
    samtools index -@ $task.cpus ${bam} ${bam.simpleName}.bam.bai
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
        tuple val(id), path(bam), path(bai)

    output:
        tuple val(id), path("${id}.cram"), path("${id}.cram.crai")


    script:
    """
        samtools view -@ $task.cpus -C -T $GENOMEREF -o ${id}.cram $tumor
        samtools index ${id}.cram -@ $task.cpus
    """
    
    stub:
    """
    touch ${id}.cram ${id}.cram.crai
    """
}


process samtools2fq {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(id), path(bam), path(bai)

    output:
        tuple val(id), path("${id}.R1.fastq"), path("${id}.R2.fastq")

    script:
    """
    samtools fastq -@ $task.cpus \
        -1 ${id}.R1.fastq -2 ${id}.R2.fastq -0 /dev/null -s /dev/null \
        -n $bam
    """
    
    stub:
    """
    touch ${id}.R1.fastq ${id}.R2.fastq 
    """
}

