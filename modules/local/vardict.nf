//References
GENOMEREF=file(params.genomes[params.genome].genome)


process vardict_tn {
    container "${params.containers.logan}"
    label 'process_somaticcaller_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.vardict.vcf.gz")
    //bcbio notes of vardict filtering var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M” and
    //filtered with “((AF*DP < 6) && ((MQ < 55.0 && NM > 1.0) || (MQ < 60.0 && NM > 2.0) || (DP < 10) || (QUAL < 45)))”
    script:

    """
    bedtools makewindows -b ${bed} -w 50150 -s 50000 > temp_${bed}

    VarDict -G $GENOMEREF \
        -f 0.01 \
        --nosv \
        -b "${tumor}|${normal}" --fisher \
        -t -Q 20 -c 1 -S 2 -E 3 \
        --th $task.cpus temp_${bed} \
        | var2vcf_paired.pl \
            -N "${tumor}|${normal}" \
            -Q 20 \
            -d 10 \
            -v 6 \
            -S \
            -f 0.05 >  ${tumorname}_vs_${normalname}_${bed.simpleName}.vardict.vcf

    printf "${normal.Name}\t${normalname}\n${tumor.Name}\t${tumorname}\n" > sampname

    bcftools reheader -s sampname ${tumorname}_vs_${normalname}_${bed.simpleName}.vardict.vcf \
        | bcftools view -Oz -o ${tumorname}_vs_${normalname}_${bed.simpleName}.vardict.vcf.gz


    """

    stub:

    """
    touch ${tumorname}_vs_${normalname}_${bed.simpleName}.vardict.vcf.gz
    """
}



process vardict_tonly {
    container "${params.containers.logan}"
    label 'process_somaticcaller_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)

    output:
        tuple val(tumorname),
        path("${tumorname}_${bed.simpleName}.tonly.vardict.vcf.gz")

    script:

    """
    bedtools makewindows -b ${bed} -w 50150 -s 50000 > temp_${bed}

    VarDict -G $GENOMEREF \
        -f 0.01 \
        -x 500 \
        --nosv \
        -b ${tumor} --fisher \
        -t -Q 20 -c 1 -S 2 -E 3 --th ${task.cpus} \
        temp_${bed} | var2vcf_valid.pl \
            -N ${tumor} \
            -Q 20 \
            -d 10 \
            -v 6 \
            -S \
            -E \
            -f 0.05 >  ${tumorname}_${bed.simpleName}_temp.tonly.vardict.vcf

    printf "${tumor.Name}\t${tumorname}\n" > sampname

    bcftools reheader -s sampname ${tumorname}_${bed.simpleName}_temp.tonly.vardict.vcf \
        | bcftools view -Oz -o ${tumorname}_${bed.simpleName}.tonly.vardict.vcf.gz

    """

    stub:

    """
    touch ${tumorname}_${bed.simpleName}.tonly.vardict.vcf.gz

    """

}