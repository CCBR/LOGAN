GENOMEREF=file(params.genomes[params.genome].genome)
MODEL="/opt/models/wgs/model.ckpt"


//Processes
//Deep Variant
process deepvariant_step1 {

    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai"), path(bed)

    output:
        tuple val(samplename), path("outputshard/${samplename}.tfrecord_${bed}.gz"),
        path("gvcf/${samplename}.gvcf.tfrecord_${bed}.gz")

    script:
    """
    mkdir -p outputshard
    mkdir -p gvcf
    make_examples \
    --mode calling \
    --ref $GENOMEREF \
    --regions ${bed} \
    --reads ${samplename}.bam \
    --channels insert_size \
    --examples outputshard/${samplename}.tfrecord_${bed}.gz \
    --gvcf gvcf/${samplename}.gvcf.tfrecord_${bed}.gz
    """

    stub:
    """
    mkdir -p outputshard
    mkdir -p gvcf
    touch outputshard/${samplename}.tfrecord_${bed}.gz
    touch gvcf/${samplename}.gvcf.tfrecord_${bed}.gz
    """

}

//Step 2 requires GPU
process deepvariant_step2 {

    input:
        tuple val(samplename), path(tfrecords), path(tfgvcf)

    output:
        tuple val(samplename), path(tfrecords),
        path("${samplename}_call_variants_output.tfrecord.gz"), path(tfgvcf)

    script:

    """
    call_variants \
    --examples "${samplename}.tfrecord_*.gz" \
    --outfile ${samplename}_call_variants_output.tfrecord.gz \
    --checkpoint $MODEL \
    --num_readers 16
    """

    stub:
    """
    touch ${samplename}_call_variants_output.tfrecord.gz
    """

}


//Step 3 DV
process deepvariant_step3 {

    input:
        tuple val(samplename), path(tfrecords), path("${samplename}_call_variants_output.tfrecord.gz"),
        path(tfgvcf)

    output:
        tuple val(samplename), path("${samplename}.vcf.gz"), path("${samplename}.vcf.gz.tbi"),
        path("${samplename}.gvcf.gz"), path("${samplename}.gvcf.gz.tbi")


    script:
    """
   postprocess_variants \
    --ref $GENOMEREF \
    --infile ${samplename}_call_variants_output.tfrecord.gz \
    --outfile ${samplename}.vcf.gz \
    --gvcf_outfile ${samplename}.gvcf.gz \
    --nonvariant_site_tfrecord_path .
    """

    stub:
    """
    touch ${samplename}.vcf.gz ${samplename}.vcf.gz.tbi
    touch ${samplename}.gvcf.gz ${samplename}.gvcf.gz.tbi

    """

}

//Combined DeepVariant
process deepvariant_combined {

    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")

    output:
        tuple val(samplename), path("${samplename}.gvcf.gz"), path("${samplename}.gvcf.gz.tbi"),
        path("${samplename}.vcf.gz"), path("${samplename}.vcf.gz.tbi")


    script:
    """
    run_deepvariant \
        --model_type=WGS \
        --ref=$GENOMEREF \
        --reads=${samplename}.bam \
        --output_gvcf= ${samplename}.gvcf.gz \
        --output_vcf=${samplename}.vcf.gz \
        --num_shards=16
    """


    stub:
    """
    touch ${samplename}.vcf.gz ${samplename}.vcf.gz.tbi
    touch ${samplename}.gvcf.gz  ${samplename}.gvcf.gz.tbi

    """


}

process glnexus {

    input:
        path(gvcfs)

    output:
        tuple path("germline.v.bcf"),
        path("germline.norm.vcf.gz"),path("germline.norm.vcf.gz.tbi")

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
