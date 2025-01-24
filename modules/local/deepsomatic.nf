GENOMEREF=file(params.genomes[params.genome].genome)

if(params.exome && params.ffpe) {
    DS_MODEL = "/opt/models/deepsomatic/ffpe_wes"
}else if(params.ffpe){
    DS_MODEL = "/opt/models/deepsomatic/ffpe_wgs"
}else if(params.exome){
    DS_MODEL = "/opt/models/deepsomatic/wes"
}else{
    DS_MODEL = "/opt/models/deepsomatic/wgs"
}

process deepsomatic_tn_step1 {
    container = "${params.containers.deepsomatic}"
    label 'process_somaticcaller'

    input:
        tuple val(tname), path(tbam), path(tbai), 
        val(nname), path(nbam), path(nbai),
        path(bed)

    output:
        tuple val(tname), val(nname),
        path("${tname}_vs_${nname}.tfrecord_${bed}.gz"),
        path("${tname}_vs_${nname}.tfrecord_${bed}.gz.example_info.json"),
        path(bed)

    script:
    """
    make_examples_somatic \
    --mode calling \
    --ref $GENOMEREF \
    --regions ${bed} \
    --checkpoint $DS_MODEL \
    --population_vcfs "/opt/models/deepsomatic/pons/AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz" \
    --vsc_max_fraction_indels_for_non_target_sample "0.5" --vsc_max_fraction_snps_for_non_target_sample "0.5" --vsc_min_fraction_indels "0.05" --vsc_min_fraction_snps "0.029" \
    --reads_tumor ${tbam} \
    --reads_normal ${nbam} \
    --sample_name_tumor ${tname} \
    --sample_name_normal ${nname} \
    --examples ${tname}_vs_${nname}.tfrecord_${bed}.gz \
    """

    stub:
    """
    touch ${tname}_vs_${nname}.tfrecord_${bed}.gz
    touch ${tname}_vs_${nname}.tfrecord_${bed}.gz.example_info.json
    """

}


process deepsomatic_tonly_step1 {
    container = "${params.containers.deepsomatic}"
    label 'process_somaticcaller'

    input:
        tuple val(tname), path(tbam), path(tbai), 
        path(bed)

    output:
        tuple val(tname),
        path("${tname}.tfrecord_${bed}.gz"),
        path("${tname}.tfrecord_${bed}.gz.example_info.json"),
        path(bed)

    script:
    """
    make_examples_somatic \
    --mode calling \
    --ref $GENOMEREF \
    --regions ${bed} \
    --checkpoint /opt/models/deepsomatic/wgs_tumor_only \
    --population_vcfs "/opt/models/deepsomatic/pons/AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz" \
    --vsc_max_fraction_indels_for_non_target_sample "0.5" --vsc_max_fraction_snps_for_non_target_sample "0.5" --vsc_min_fraction_indels "0.07" --vsc_min_fraction_snps "0.05" \
    --reads_tumor ${tbam} \
    --sample_name_tumor ${tname} \
    --examples ${tname}.tfrecord_${bed}.gz
    """

    stub:
    """
    touch ${tname}.tfrecord_${bed}.gz
    touch ${tname}.tfrecord_${bed}.gz.example_info.json
    """

}


//Step 2 can run in CPU or GPU mode for now use only CPUs
process deepsomatic_step2 {
    container = "${params.containers.deepsomatic}"
    label 'process_somaticcaller'
    errorStrategy { task.exitStatus == 1 ? 'ignore' : 'terminate' }

    input:
        tuple val(samplename), path(tfrecords), path(json), path(bed)

    output:
        tuple val(samplename),
        path(tfrecords),
        path("outds/*"), path(bed)

    script:
    sub_cpus = "$task.cpus".toInteger() - 1

    """
    mkdir -p outds/
    call_variants \
    --examples $tfrecords \
    --outfile outds/${samplename}_call_variants_output.tfrecord.gz \
    --checkpoint $DS_MODEL \
    --num_readers $sub_cpus
    """

    stub:
    """
    mkdir -p outds
    touch "outds/${samplename}_call_variants_output.tfrecord.gz"
    """
}

process deepsomatic_tonly_step2 {
    container = "${params.containers.deepsomatic}"
    label 'process_somaticcaller'
    errorStrategy { task.exitStatus == 1 ? 'ignore' : 'terminate' }

    input:
        tuple val(samplename), path(tfrecords), path(json), path(bed)

    output:
        tuple val(samplename),
        path(tfrecords),
        path("outds/*"), path(bed)

    script:
    sub_cpus = "$task.cpus".toInteger() - 1

    """
    mkdir -p outds/
    call_variants \
    --examples $tfrecords \
    --outfile outds/${samplename}_call_variants_output.tfrecord.gz \
    --checkpoint /opt/models/deepsomatic/wgs_tumor_only \
    --num_readers $sub_cpus
    """

    stub:
    """
    mkdir -p outds
    touch "outds/${samplename}_call_variants_output.tfrecord.gz"
    """
}


//Step 3 DV
process deepsomatic_step3 {
    container = "${params.containers.deepsomatic}"
    label 'process_somaticcaller'

    input:
        tuple val(samplename), path(tfrecords),
        path("outds/*"), path(bed)

    output:
        tuple val(samplename), path("${samplename}_${bed}.vcf.gz"), path("${samplename}_${bed}.vcf.gz.tbi")       


    script:
    """
    postprocess_variants \
        --ref $GENOMEREF \
        -j $task.cpus \
        --process_somatic=true --pon_filtering "/opt/models/deepsomatic/pons/PON_dbsnp138_gnomad_ILMN1000g_pon.vcf.gz" \
        --infile outds/${samplename}_call_variants_output.tfrecord.gz \
        --outfile ${samplename}_${bed}.vcf.gz
    """

    stub:
    """
    touch ${samplename}_${bed}.vcf.gz ${samplename}_${bed}.vcf.gz.tbi
    """

}


process bcfconcat {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(samplename), path(vcf), path(index), val(type)

    output:
        tuple val(samplename), path("${samplename}.${type}.gz"), path("${samplename}.${type}.gz.tbi")

    script:
        vcfin=vcf.join(" ")

    """
    bcftools concat $vcfin --write-index -Oz -o ${samplename}.${type}.gz##idx##${samplename}.${type}.gz.tbi
    """

    stub:
    """
    touch ${samplename}.${type}.gz ${samplename}.${type}.gz.tbi
    """

}

