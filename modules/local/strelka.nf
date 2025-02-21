//References
GENOMEREF=file(params.genomes[params.genome].genome)

//HelperScripts
STRELKA_CONVERT=params.strelka_convert


process strelka_tn {
    container "${params.containers.logan}"
    label 'process_high'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.snvs.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.snvs.vcf.gz.tbi"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.indels.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.indels.vcf.gz.tbi")

    script:

    """
    mkdir -p wd

    bgzip ${bed}
    tabix ${bed}.gz

    configureStrelkaSomaticWorkflow.py \
        --ref=$GENOMEREF \
        --tumor=${tumor} \
        --normal=${normal} \
        --runDir=wd \
        --callRegions ${bed}.gz
    ./wd/runWorkflow.py -m local -j $task.cpus
    mv wd/results/variants/somatic.snvs.vcf.gz  ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic_temp.snvs.vcf.gz
    mv wd/results/variants/somatic.indels.vcf.gz  ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic_temp.indels.vcf.gz

    printf %s "NORMAL\t${normalname}\nTUMOR\t${tumorname}\n" >sampname

    bcftools reheader -s sampname ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic_temp.snvs.vcf.gz \
        | bcftools view -Oz -o ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.snvs.vcf.gz
    bcftools reheader -s sampname ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic_temp.indels.vcf.gz \
        | bcftools view -Oz -o ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.indels.vcf.gz

    bcftools index -t ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.snvs.vcf.gz
    bcftools index -t ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.indels.vcf.gz

    """

    stub:

    """
    touch ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.snvs.vcf.gz  ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.snvs.vcf.gz.tbi
    touch ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.indels.vcf.gz ${tumorname}_vs_${normalname}_${bed.simpleName}.somatic.indels.vcf.gz.tbi

    """

}

process convert_strelka {
    //Add GT/AD column to Strelka Variants
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumor), val(normal), val(vc),
        path(strelkavcf), path(strelkaindex)

    output:
        tuple val(tumor), val(normal), val("strelka"),
        path("${tumor}_vs_${normal}.filtered.strelka-fixed.vcf.gz"), 
        path("${tumor}_vs_${normal}.filtered.strelka-fixed.vcf.gz.tbi")


    script:

    """
    python $STRELKA_CONVERT ${strelkavcf} ${tumor}_vs_${normal}.filtered.strelka-fixed.vcf.gz
    bcftools index -t ${tumor}_vs_${normal}.filtered.strelka-fixed.vcf.gz
    """

    stub:

    """
    touch ${tumor}_vs_${normal}.filtered.strelka-fixed.vcf.gz ${tumor}_vs_${normal}.filtered.strelka-fixed.vcf.gz.tbi
    """

}

process combineVariants_strelka {
    //Concat all somatic snvs/indels across all files, strelka separates snv/indels
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(sample),
        path(strelkasnvs), path(snvindex),
        path(strelkaindels), path(indelindex)

    output:
        tuple val(sample),
        path("${sample}.strelka.vcf.gz"), path("${sample}.strelka.vcf.gz.tbi"),
        path("${sample}.filtered.strelka.vcf.gz"), path("${sample}.filtered.strelka.vcf.gz.tbi")


    script:

    vcfin = strelkasnvs.join(" ")
    indelsin = strelkaindels.join(" ")
    samplist=sample.split('_vs_')
    if(samplist.size()>1){
        samporder = samplist.join(",")
    }else{
        samporder = sample
    }
    """
    bcftools concat $vcfin $indelsin --threads $task.cpus -Oz -o ${sample}.temp.strelka.vcf.gz -a
    bcftools norm ${sample}.temp.strelka.vcf.gz -m- --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' > ${sample}.temp1.strelka.vcf.gz

    bcftools sort ${sample}.temp1.strelka.vcf.gz |bcftools view - -s $samporder -Oz -o ${sample}.strelka.vcf.gz

    bcftools view ${sample}.strelka.vcf.gz --threads $task.cpus -f PASS -Oz -o ${sample}.filtered.strelka.vcf.gz

    bcftools index ${sample}.strelka.vcf.gz -t
    bcftools index ${sample}.filtered.strelka.vcf.gz -t
    """

    stub:

    """
    touch ${sample}.strelka.vcf.gz ${sample}.strelka.vcf.gz.tbi
    touch ${sample}.filtered.strelka.vcf.gz ${sample}.filtered.strelka.vcf.gz.tbi

    """

}
