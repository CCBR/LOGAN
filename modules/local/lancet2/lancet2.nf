GENOMEREF=file(params.genomes[params.genome].genome)

process lancet2_tn {
    container "${params.containers.lancet}"
    label 'process_somaticcaller'
    errorStrategy 'ignore'

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
	Lancet2 pipeline \
	--normal ${normal} \
	--tumor ${tumor} \
	--reference $GENOMEREF \
	--num-threads $task.cpus -R ${bed} \
	--out-vcfgz ${tumorname}_vs_${normalname}_${bed.simpleName}_temp.vcf.gz

	python3 score_variants.py \
    ${tumorname}_vs_${normalname}_${bed.simpleName}_temp.vcf.gz somatic_ebm.lancet_6ef7ba445a.v1.pkl > ${tumorname}_vs_${normalname}_${bed.simpleName}_scored.vcf
	
	bcftools view ${tumorname}_vs_${normalname}_${bed.simpleName}_scored.vcf -Oz -o ${tumorname}_vs_${normalname}_${bed.simpleName}_lancet.vcf.gz
    bcftools index -t ${tumorname}_vs_${normalname}_${bed.simpleName}_lancet.vcf.gz

    """

    stub:

    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz" "${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz.tbi"

    """
}


