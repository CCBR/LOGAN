GENOMEREF = file(params.genomes[params.genome].genome)

process bwamem2 {
    container = "${params.containers.logan}"
    tag { name }
    errorStrategy { task.exitStatus in [137,140,143] ? 'retry' : 'terminate' }
    maxRetries 2

    memory { 
        if (task.attempt == 2) return '180 GB'
        else if (task.attempt == 3) return '200 GB'
    }
    

    input:
        tuple val(samplename),
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz")

    output:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bam.bai")

    script:
    sub_cpus = "$task.cpus".toInteger()/2
    def sub_mem= (task.memory.toGiga().toInteger()/task.cpus.toInteger())-2

    """
    # Check for AVX512 then AVX2 then SSE41
    #if lscpu | grep -q avx512; then
    #    BWA_BINARY="bwa-mem2.avx512bw"
    if lscpu | grep -q avx2; then
        BWA_BINARY="bwa-mem2.avx2"
    elif lscpu | grep -q sse4_1; then
        BWA_BINARY="bwa-mem2.sse41"
    else
        BWA_BINARY="bwa-mem2"
    fi
    
    mkdir -p tmp
    \$BWA_BINARY mem -M \
        -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:illumina\\tLB:${samplename}\\tPU:${samplename}\\tCN:hgsc\\tDS:wgs' \
        -t $task.cpus \
        ${GENOMEREF} \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz |\
    samblaster -M | \
    samtools sort -T tmp/ -@ $sub_cpus -m ${sub_mem}G - --write-index -o ${samplename}.bam##idx##${samplename}.bam.bai
    """

    stub:
    """
    touch ${samplename}.bam ${samplename}.bam.bai
    """
}


process BWAMEM2_SPLIT {
    container = "${params.containers.logan}"
    tag { name }
    errorStrategy { task.exitStatus in [137,140,143] ? 'retry' : 'terminate' }
    maxRetries 2

    memory { 
        if (task.attempt == 2) return '48 GB'
        else if (task.attempt == 3) return '64 GB'
    }
    

    input:
        tuple val(samplename),
        path(reads), val(chunk)

    output:
        tuple val(samplename), 
		path("${samplename}_${chunk}.bam")

    script:

    """
    # Check for AVX512 then AVX2 then SSE41
    #if lscpu | grep -q avx512; then
    #    BWA_BINARY="bwa-mem2.avx512bw"
    if lscpu | grep -q avx2; then
        BWA_BINARY="bwa-mem2.avx2"
    elif lscpu | grep -q sse4_1; then
        BWA_BINARY="bwa-mem2.sse41"
    else
        BWA_BINARY="bwa-mem2"
    fi

    \$BWA_BINARY mem -M \
        -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:illumina\\tLB:${samplename}\\tPU:${samplename}\\tCN:hgsc\\tDS:wgs' \
        -t $task.cpus \
        ${GENOMEREF} \
        ${reads[0]} ${reads[1]} | \
	samtools view - -b -o ${samplename}_${chunk}.bam

    """

    stub:
    """
    touch ${samplename}_${chunk}.bam 
    """
}




process COMBINE_ALIGNMENTS {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
    tuple val(samplename), path(bam)

    output:
    tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bam.bai"), emit: bams
    path("${samplename}.metrics")

    script:
    allbams = bam.join(" ")
	sub_cpus = "$task.cpus".toInteger()/2

    """
    mkdir tmp
    samtools cat -@ $task.cpus $allbams -o ${samplename}_tmp.sam
    samblaster -M ${samplename}.metrics ${samplename}_tmp.sam  |\
        samtools sort -T tmp/ -@ $sub_cpus -m 10G - --write-index -o ${samplename}.bam##idx##${samplename}.bam.bai

    """

    stub:
    """
    touch "${samplename}.bam" "${samplename}.bam.bai"
    """

}