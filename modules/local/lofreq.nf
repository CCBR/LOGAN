GENOMEREF=file(params.genomes[params.genome].genome)

//DBSNP for LOFREQ/MUSE
DBSNP=file(params.genomes[params.genome].dbsnp)
//HelperScripts
LOFREQ_CONVERT=params.lofreq_convert

process lofreq_tn {
    container "${params.containers.lofreq}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)


    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.snvs.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.indels.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz.tbi")

    script:

    """
    lofreq somatic -f $GENOMEREF -n ${normal} -t ${tumor} \
        -d $DBSNP \
        --threads $task.cpus \
        -l ${bed} \
        --call-indels \
        -o ${tumorname}_vs_${normalname}_${bed.simpleName}_

    bcftools concat ${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz \
        ${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz --threads $task.cpus -Oz -o \
        ${tumorname}_vs_${normalname}_${bed.simpleName}_temp_lofreq.vcf.gz

    $LOFREQ_CONVERT -i ${tumorname}_vs_${normalname}_${bed.simpleName}_temp_lofreq.vcf.gz -g 1/0 \
        -n ${tumorname} -o ${tumorname}_vs_${normalname}_${bed.simpleName}_temp1_lofreq.vcf.gz

    bcftools view -h ${tumorname}_vs_${normalname}_${bed.simpleName}_temp1_lofreq.vcf.gz >temphead

    sed 's/^##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">/##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref\\/fwd, ref\\/rev, var\\/fwd, var\\/rev">/' temphead > temphead1
    bcftools reheader ${tumorname}_vs_${normalname}_${bed.simpleName}_temp1_lofreq.vcf.gz -h temphead1 |\
        bcftools view -Oz -o ${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz

    bcftools index -t ${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz

    """

    stub:

    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.snvs.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz" "${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz.tbi"

    """
}


