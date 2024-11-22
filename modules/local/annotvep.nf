//References
GENOMEREF=file(params.genomes[params.genome].genome)

//VEP
VEPCACHEDIR=file(params.genomes[params.genome].vepcache)
VEPSPECIES=params.genomes[params.genome].vepspecies
VEPBUILD=params.genomes[params.genome].vepbuild

process annotvep_tn {
    label 'process_medium'
    container "${params.containers.vcf2maf}"

    input:
        tuple val(tumorsample), val(normalsample),
        val(vc), path(tumorvcf), path(vcfindex)

    output:
        path("paired/${vc}/${tumorsample}_vs_${normalsample}.maf")

    shell:

    '''
    VCF_SAMPLE_IDS=($(bcftools query -l !{tumorvcf}))
    TID_IDX=0
    NID_IDX=""
    VCF_NID=""
    NORM_VCF_ID_ARG=""
    NSAMPLES=${#VCF_SAMPLE_IDS[@]}
    if [ $NSAMPLES -gt 1 ]; then
        # Assign tumor, normal IDs
        # Look through column names and
        # see if they match provided IDs
        for (( i = 0; i < $NSAMPLES; i++ )); do
            echo "${VCF_SAMPLE_IDS[$i]}"
            if [ "${VCF_SAMPLE_IDS[$i]}" == !{tumorsample} ]; then
                TID_IDX=$i
            fi

            if [ "${VCF_SAMPLE_IDS[$i]}" == !{normalsample} ]; then
                NID_IDX=$i
            fi
        done

        if [ ! -z $NID_IDX ]; then
            VCF_NID=${VCF_SAMPLE_IDS[$NID_IDX]}
            NORM_VCF_ID_ARG="--vcf-normal-id $VCF_NID"
        fi
    fi
    VCF_TID=${VCF_SAMPLE_IDS[$TID_IDX]}

    zcat !{tumorvcf} > !{tumorvcf.baseName}

    mkdir -p paired/!{vc}

    vcf2maf.pl \
    --vep-forks !{task.cpus} --input-vcf !{tumorvcf.baseName} \
    --output-maf paired/!{vc}/!{tumorsample}_vs_!{normalsample}.maf \
    --tumor-id !{tumorsample} \
    --normal-id !{normalsample} \
    --vep-path /opt/vep/src/ensembl-vep \
    --vep-data !{VEPCACHEDIR} \
    --ncbi-build !{VEPBUILD} --species !{VEPSPECIES} --ref-fasta !{GENOMEREF} \
    --retain-info "set" \
    --vep-overwrite

    '''

    stub:
    """
    mkdir -p paired/${vc}
    touch paired/${vc}/${tumorsample}_vs_${normalsample}.maf
    """
}


process annotvep_tonly {
    container "${params.containers.vcf2maf}"
    label 'process_medium'

    input:
        tuple val(tumorsample),
        val(vc), path(tumorvcf),
        path(vcfindex)


    output:
        path("tumor_only/${vc}/${tumorsample}.tonly.maf")

    shell:

    '''
    VCF_SAMPLE_IDS=($(bcftools query -l !{tumorvcf}))
    TID_IDX=0
    NID_IDX=""
    VCF_NID=""
    NORM_VCF_ID_ARG=""
    NSAMPLES=${#VCF_SAMPLE_IDS[@]}
    if [ $NSAMPLES -gt 1 ]; then
        # Assign tumor, normal IDs
        # Look through column names and
        # see if they match provided IDs
        for (( i = 0; i < $NSAMPLES; i++ )); do
            echo "${VCF_SAMPLE_IDS[$i]}"
            if [ "${VCF_SAMPLE_IDS[$i]}" == !{tumorsample} ]; then
                TID_IDX=$i
            fi

        done

        if [ ! -z $NID_IDX ]; then
            VCF_NID=${VCF_SAMPLE_IDS[$NID_IDX]}
            NORM_VCF_ID_ARG="--vcf-normal-id $VCF_NID"
        fi
    fi
    VCF_TID=${VCF_SAMPLE_IDS[$TID_IDX]}

    zcat !{tumorvcf} > !{tumorvcf.baseName}

    mkdir -p tumor_only/!{vc}

    vcf2maf.pl \
    --vep-forks !{task.cpus} --input-vcf !{tumorvcf.baseName} \
    --output-maf tumor_only/!{vc}/!{tumorsample}.tonly.maf \
    --tumor-id !{tumorsample} \
    --vep-path /opt/vep/src/ensembl-vep \
    --vep-data !{VEPCACHEDIR} \
    --ncbi-build !{VEPBUILD} --species !{VEPSPECIES} --ref-fasta !{GENOMEREF} \
    --retain-info "set" \
    --vep-overwrite


    '''

    stub:
    """
    mkdir -p tumor_only/${vc}
    touch tumor_only/${vc}/${tumorsample}.tonly.maf
    """
}
