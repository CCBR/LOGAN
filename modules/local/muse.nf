//References
GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEFAI=file(params.genomes[params.genome].genomefai)
GENOMEDICT=file(params.genomes[params.genome].genomedict)
GERMLINE_RESOURCE=file(params.genomes[params.genome].germline_resource)

//DBSNP for LOFREQ/MUSE
DBSNP=file(params.genomes[params.genome].dbsnp)


process muse_tn {
    container "${params.containers.logan}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai),
        val(mode)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}.vcf.gz")

    script:

    """
    MuSE call -f $GENOMEREF -O ${tumorname}_vs_${normalname} -n $task.cpus $tumor $normal
    MuSE sump -I ${tumorname}_vs_${normalname}.MuSE.txt \
        -O ${tumorname}_vs_${normalname}.vcf -n $task.cpus -D $DBSNP $mode

    bcftools view ${tumorname}_vs_${normalname}.vcf -Oz -o ${tumorname}_vs_${normalname}_temp.vcf.gz

    printf "NORMAL\t${normalname}\nTUMOR\t${tumorname}\n" > sampname

    bcftools reheader -s sampname ${tumorname}_vs_${normalname}_temp.vcf.gz \
        | bcftools view -Oz -o ${tumorname}_vs_${normalname}.vcf.gz

    """

    stub:

    """
    touch "${tumorname}_vs_${normalname}.vcf.gz"
    """

}


process muse_tn_arch {
    container "${params.containers.logan}"
    label 'process_somaticcaller'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)


    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}.vcf.gz")

    script:

    """
    MuSE call -f $GENOMEREF -O ${tumorname}_vs_${normalname} -n $task.cpus $tumor $normal
    MuSE sump -I ${tumorname}_vs_${normalname}.MuSE.txt \
        -O ${tumorname}_vs_${normalname}.vcf -n $task.cpus -D $DBSNP -G

    bcftools view ${tumorname}_vs_${normalname}.vcf -Oz -o ${tumorname}_vs_${normalname}_temp.vcf.gz

    printf "NORMAL\t${normalname}\nTUMOR\t${tumorname}\n" > sampname

    bcftools reheader -s sampname ${tumorname}_vs_${normalname}_temp.vcf.gz \
        | bcftools view -Oz -o ${tumorname}_vs_${normalname}.vcf.gz

    """

    stub:

    """
    touch "${tumorname}_vs_${normalname}.vcf.gz"
    """

}
