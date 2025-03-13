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

process fastp_split {
    container = "${params.containers.logan}"
    label 'process_medium'
    tag { name }

    input:
    tuple val(samplename), path(fqs)

    output:
    tuple val(samplename), path("*_R{1,2}.trimmed.fastq.gz"),
    path("${samplename}.fastp.json"),
    path("${samplename}.fastp.html")
    
        

    script:
    """
    fastp -w $task.cpus \
        --detect_adapter_for_pe \
		-S $params.split_fastq \
        --in1 ${fqs[0]} \
        --in2 ${fqs[1]} \
        --out1 ${samplename}_R1.trimmed.fastq.gz \
        --out2 ${samplename}_R2.trimmed.fastq.gz  \
        --json ${samplename}.fastp.json \
        --html ${samplename}.fastp.html
    """

    stub:
    """
    touch 0001.${samplename}_R1.trimmed.fastq.gz
    touch 0001.${samplename}_R2.trimmed.fastq.gz
    touch ${samplename}.fastp.json
    touch ${samplename}.fastp.html

    """
}
