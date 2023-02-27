//References
GENOME=file(params.genome)
MODEL="/opt/models/wgs/model.ckpt"
intervalbedin = file(params.intervals)




//Output Directory
outdir=file(params.output)

//Processes
//Deep Variant
process deepvariant_step1 {
    module=['deepvariant/1.4.0']
    
    //publishDir("${outdir}/deepvariant", mode: 'copy')

    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai"), path(bed)
    
    output:
        tuple val(samplename), path("outputshard/${samplename}.tfrecord_${bed}.gz"), 
        path("gvcf/${samplename}.gvcf.tfrecord_${bed}.gz")

    script:     
    """
    mkdir outputshard gvcf
    make_examples \
    --mode calling \
    --ref $GENOME \
    --regions ${bed} \
    --reads ${samplename}.bam \
    --channels insert_size \
    --examples outputshard/${samplename}.tfrecord_${bed}.gz \
    --gvcf gvcf/${samplename}.gvcf.tfrecord_${bed}.gz 
    """

}

//Step 2 requires GPU
process deepvariant_step2 {
    
    module=['deepvariant/1.4.0']
    
    //publishDir("${outdir}/deepvariant", mode: 'copy')
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
}


//Step 3 DV
process deepvariant_step3 {
    scratch '/lscratch/$SLURM_JOB_ID/dv'
    publishDir("${outdir}/deepvariant", mode: 'copy')

    module=['deepvariant/1.4.0']
    
    input:
        tuple val(samplename), path(tfrecords), path("${samplename}_call_variants_output.tfrecord.gz"),
        path(tfgvcf)
    
    output:
        tuple val(samplename), path("${samplename}.vcf.gz"), path("${samplename}.vcf.gz.tbi"),
        path("${samplename}.gvcf.gz"), path("${samplename}.gvcf.gz.tbi")


    script: 
    """
   postprocess_variants \
    --ref $GENOME \
    --infile ${samplename}_call_variants_output.tfrecord.gz \
    --outfile ${samplename}.vcf.gz \
    --gvcf_outfile ${samplename}.gvcf.gz \
    --nonvariant_site_tfrecord_path .
    """
}

//Combined DeepVariant
process deepvariant_combined {
    module=['deepvariant/1.4.0']
    scratch '/lscratch/$SLURM_JOB_ID/dv'

    publishDir("${outdir}/germline_deepvariant", mode: 'copy')

    input:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai"), path(bed)
    
    output:
        tuple val(samplename), path("${samplename}.gvcf.gz"), path("${samplename}.gvcf.gz.tbi"),
        path("${samplename}.vcf.gz"), path("${samplename}.vcf.gz.tbi")


    script:     
    """
    run_deepvariant \
        --model_type=WGS \
        --ref=$GENOME \
        --reads=${samplename}.bam \
        --output_gvcf= ${samplename}.gvcf.gz \
        --output_vcf=${samplename}.vcf.gz \
        --num_shards=16 
    """

}

process glnexus {
    //scratch '/lscratch/$SLURM_JOB_ID/dv'
    publishDir("${outdir}/deepvariant", mode: 'copy')

    module=['glnexus','bcftools']
    
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
        -f $GENOME \
        -o germline.norm.vcf.gz \
        germline.v.bcf

    bcftools index \
        -f -t \
        --threads 8 \
        germline.norm.vcf.gz
    
    """
}





