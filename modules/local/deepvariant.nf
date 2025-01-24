GENOMEREF=file(params.genomes[params.genome].genome)
MODEL="/opt/models/wgs/"


//Processes
//Deep Variant
process deepvariant_step1 {
    container = "${params.containers.deepvariant}"
    label 'process_somaticcaller'

    input:
        tuple val(samplename), path(bam), path(bai), path(bed)

    output:
        tuple val(samplename), path("${samplename}.tfrecord_${bed}.gz"),
        path("${samplename}.tfrecord_${bed}.gz.example_info.json"),
        path("${samplename}.gvcf.tfrecord_${bed}.gz"), path(bed)

    script:
    """
    make_examples \
    --mode calling \
    --ref $GENOMEREF \
    --regions ${bed} \
    --reads ${bam} \
    --channels insert_size \
    --examples ${samplename}.tfrecord_${bed}.gz \
    --gvcf ${samplename}.gvcf.tfrecord_${bed}.gz
    """

    stub:
    """
    touch ${samplename}.tfrecord_${bed}.gz
    touch ${samplename}.tfrecord_${bed}.gz.example_info.json
    touch ${samplename}.gvcf.tfrecord_${bed}.gz
    """

}

//Step 2 requires GPU
process deepvariant_step2 {
    container = "${params.containers.deepvariant}"
    //clusterOptions '--gres=lscratch:100,gpu:p100:1  --partition=gpu'
    label 'process_somaticcaller'

    input:
        tuple val(samplename), path(tfrecords), path(json), path(tfgvcf), path(bed)

    output:
        tuple val(samplename), path(tfrecords),
        path(tfgvcf), path("outdv/*"), path(bed)

    script:
    sub_cpus = "$task.cpus".toInteger() - 1

    """
    mkdir -p outdv/
    call_variants \
    --examples $tfrecords \
    --outfile outdv/${samplename}_call_variants_output.tfrecord.gz \
    --checkpoint $MODEL \
    --writer_threads $sub_cpus
    """

    stub:
    """
    mkdir -p outdv
    touch "outdv/${samplename}_call_variants_output.tfrecord.gz"
    """
}


//Step 3 DV
process deepvariant_step3 {
    container = "${params.containers.deepvariant}"
    label 'process_somaticcaller'
    
    input:
        tuple val(samplename), path(tfrecords),
        path(tfgvcf), path("outdv/*"), path(bed)

    output:
        tuple val(samplename), path("${samplename}_${bed}.vcf.gz"), path("${samplename}_${bed}.vcf.gz.tbi"),
        path("${samplename}_${bed}.gvcf.gz"), path("${samplename}_${bed}.gvcf.gz.tbi")


    script:
    """
    postprocess_variants \
        --ref $GENOMEREF \
        --infile outdv/${samplename}_call_variants_output.tfrecord.gz \
        --outfile ${samplename}_${bed}.vcf.gz \
        --gvcf_outfile ${samplename}_${bed}.gvcf.gz \
        --nonvariant_site_tfrecord_path .
    """

    stub:
    """
    touch ${samplename}_${bed}.vcf.gz ${samplename}_${bed}.vcf.gz.tbi
    touch ${samplename}_${bed}.gvcf.gz ${samplename}_${bed}.gvcf.gz.tbi

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


process glnexus {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        path(gvcfs)

    output:
        tuple path("germline.v.bcf"),
        path("germline.norm.vcf.gz"), path("germline.norm.vcf.gz.tbi")

    script:

    """
    glnexus_cli --config DeepVariant_unfiltered \
        *.gvcf.gz --threads 8 > germline.v.bcf

    bcftools norm \
        -m - \
        -Oz \
        --threads 8 \
        -f $GENOMEREF \
        -o germline.norm.vcf.gz \
        germline.v.bcf

    bcftools index \
        -f -t \
        --threads 8 \
        germline.norm.vcf.gz
    """

    stub:
    """
        touch germline.v.bcf
        touch germline.norm.vcf.gz
        touch germline.norm.vcf.gz.tbi
    """
}




//Combined DeepVariant
process deepvariant_combined {
    module = ['deepvariant/1.6.0']

    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        tuple val(samplename), path("${samplename}.gvcf.gz"), path("${samplename}.gvcf.gz.tbi"),
        path("${samplename}.vcf.gz"), path("${samplename}.vcf.gz.tbi")


    script:
    """
    run_deepvariant \
        --model_type=WGS \
        --ref=$GENOMEREF \
        --reads=${bam} \
        --output_gvcf=${samplename}.gvcf.gz \
        --output_vcf=${samplename}.vcf.gz \
        --num_shards=16
    """


    stub:
    """
    touch ${samplename}.vcf.gz ${samplename}.vcf.gz.tbi
    touch ${samplename}.gvcf.gz  ${samplename}.gvcf.gz.tbi
    """


}