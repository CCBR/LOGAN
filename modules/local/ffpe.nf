
process sobdetect_pass1 {
    container = "${params.containers.ffpe}"
    label 'process_medium'

    input:
    tuple val(sample), path(vcf), path(bam), val(vc)

    output:
    tuple val(sample),
        path("${vc}/pass1/${sample}.pass1.sobdetect.vcf"),
        path("${vc}/pass1/${sample}.info")

    script:
    """
    mkdir -p ${vc}/pass1
    java -jar $SOB_JAR \
        --input-type VCF \
        --input-variants ${vcf} \
        --input-bam ${bam} \
        --output-variants ${sample}.pass1.sobdetect.vcf \
        --only-passed false 

    bcftools query \
        -f '%INFO/numF1R2Alt\t%INFO/numF2R1Alt\t%INFO/numF1R2Ref\t%INFO/numF2R1Ref\t%INFO/numF1R2Other\t%INFO/numF2R1Other\t%INFO/SOB\n' \
        ${sample}.sobdetect.vcf \
        | awk '{if (\$1 != "."){tum_alt=\$1+\$2; tum_depth=\$1+\$2+\$3+\$4+\$5+\$6; if (tum_depth==0){tum_af=1} else {tum_af=tum_alt/tum_depth }; print tum_alt,tum_depth,tum_af,\$7}}' \
        > ${sample}.info
    
    mv ${sample}.pass1.sobdetect.vcf ${vc}/pass1
    mv ${sample}.info ${vc}/pass1
 
    """

    stub:
    """
    mkdir -p ${vc}/pass1
    touch ${vc}/pass1/${sample}.pass1.sobdetect.vcf
    touch ${vc}/pass1/${sample}.info
    """
}



// Cohort parameter calculation
process sobdetect_cohort_params {
    container = "${params.containers.ffpe}"
    label 'process_medium'

    input:
        path info_files

    output:
        tuple path("all_samples.info"), path("cohort_params.txt")

    script:
    allinfos = info_file.join(" ")

    """
    echo -e "#TUMOR.alt\tTUMOR.depth\tTUMOR.AF\tSOB\tFS\tSOR\tTLOD\tReadPosRankSum" > all_samples.info
    cat ${allinfos} >> all_samples.info


    grep -v '^#' all_samples.info \
        | awk '{ total1 += \$1; ss1 += \$1^2; total2 += \$2; ss2 += \$2^2; total3 += \$3; ss3 += \$3^2; total4 += \$4; ss4 += \$4^2 } END { print total1/NR,total2/NR,total3/NR,total4/NR; print sqrt(ss1/NR-(total1/NR)^2),sqrt(ss2/NR-(total2/NR)^2),sqrt(ss3/NR-(total3/NR)^3),sqrt(ss4/NR-(total4/NR)^2) }' > cohort_params.txt
    """
    
    stub:
    """
    touch all_samples.info cohort_params.txt
    """
}   

process sobdetect_pass2 {
    container = "${params.containers.ffpe}"
    label 'process_medium'

    input:
    tuple val(sample), path(vcf), path(bam), val(vc), path(sample_info), path(params_file)
        
    output:
    tuple val(sample), 
        path("${vc}/pass2/${sample}.pass2.sobdetect.vcf"),
        path("${vc}/pass2/${sample}.info"),
        path("${vc}/pass2/${sample}_${vc}.artifact_filtered.vcf.gz"),
        path("${vc}/pass2/${sample}_${vc}.artifact_filtered.vcf.gz.tbi")

    script:
    """
    mkdir -p ${vc}/pass2
    java -jar $SOB_JAR \
        --input-type VCF \
        --input-variants ${vcf} \
        --input-bam ${bam} \
        --output-variants pass2_output.vcf \
        --only-passed true \
        --standardization-parameters ${params_file}

    bcftools query \
        -f '%INFO/numF1R2Alt\t%INFO/numF2R1Alt\t%INFO/numF1R2Ref\t%INFO/numF2R1Ref\t%INFO/numF1R2Other\t%INFO/numF2R1Other\t%INFO/SOB\n' \
        pass2_output.vcf \
        | awk '{if (\$1 != "."){tum_alt=\$1+\$2; tum_depth=\$1+\$2+\$3+\$4+\$5+\$6; if (tum_depth==0){tum_af=1} else {tum_af=tum_alt/tum_depth }; print tum_alt,tum_depth,tum_af,\$7}}' \
        > F{sample}.info

    # Artifact filtering
    bcftools filter \
        -e 'INFO/pArtifact < 0.05' \
        -Oz \
        -o ${sample}.artifact_filtered.vcf.gz ${sample}.sobdetect.vcf

    bcftools index -f -t ${sample}.artifact_filtered.vcf.gz

    mv ${sample}.pass2.sobdetect.vcf ${vc}/pass2
    mv ${sample}.info ${vc}/pass2
    mv ${sample}.artifact_filtered.vcf.gz ${vc}/pass2
    """

    stub:
    """
    mkdir -p ${vc}/pass2
    touch ${vc}/pass2/${sample}.pass2.sobdetect.vcf
    touch ${vc}/pass2/${sample}.info
    touch ${vc}/pass2/${sample}_${vc}.artifact_filtered.vcf.gz
    touch ${vc}/pass2/${sample}_${vc}.artifact_filtered.vcf.gz.tbi

    """
}

// Metrics calculation
process sobdetect_metrics {
    container = "${params.containers.ffpe}"
    label 'process_medium'

    input:
        path (pass1_vcfs)
        path (pass2_vcfs)

    output:
    tuple path("variant_count_table.txt"), 
        path("all_metrics.txt")

    script:
    """
    echo -e "#ID\tDefaultParam\tCohortParam\tTotalVariants" > variant_count_table.txt
    echo -e "#SAMPLE_ID\tParam\tCHROM\tPOS\tnumF1R2Alt\tnumF2R1Alt\tnumF1R2Ref\tnumF2R1Ref\tnumF1R2Other\tnumF2R1Other\tSOB\tpArtifact\tFS\tSOR\tTLOD\tReadPosRankSum" > all_metrics.txt

    P1FILES=(\$(echo ${pass1_vcfs}))
    P2FILES=(\$(echo ${pass2_vcfs}))
    for (( i=0; i<\${#P1FILES[@]}; i++ )); do
        MYID=\$(basename -s ".sobdetect.vcf" \${P1FILES[\$i]})
        
        total_count=\$(grep -v ^# \${P1FILES[\$i]} | wc -l) || total_count=0
        count_1p=\$(bcftools query -f '%INFO/pArtifact\n' \${P1FILES[\$i]} | awk '{if (\$1 != "." && \$1 < 0.05){print}}' | wc -l)
        count_2p=\$(bcftools query -f '%INFO/pArtifact\n' \${P2FILES[\$i]} | awk '{if (\$1 != "." && \$1 < 0.05){print}}' | wc -l)

        echo -e "\$MYID\t\$count_1p\t\$count_2p\t\$total_count" >> variant_count_table.txt

        bcftools query -f '%CHROM\t%POS\t%INFO/numF1R2Alt\t%INFO/numF2R1Alt\t%INFO/numF1R2Ref\t%INFO/numF2R1Ref\t%INFO/numF1R2Other\t%INFO/numF2R1Other\t%INFO/SOB\t%INFO/pArtifact\n' \${P1FILES[\$i]} | awk -v id=\$MYID 'BEGIN{OFS="\t"}{print id,"PASS_1",\$0}' >> all_metrics.txt
        bcftools query -f '%CHROM\t%POS\t%INFO/numF1R2Alt\t%INFO/numF2R1Alt\t%INFO/numF1R2Ref\t%INFO/numF2R1Ref\t%INFO/numF1R2Other\t%INFO/numF2R1Other\t%INFO/SOB\t%INFO/pArtifact\n' \${P2FILES[\$i]} | awk -v id=\$MYID 'BEGIN{OFS="\t"}{print id,"PASS_2",\$0}' >> all_metrics.txt
    done
    """

    stub:
    """
    touch variant_count_table.txt all_metrics.txt
    """

}

