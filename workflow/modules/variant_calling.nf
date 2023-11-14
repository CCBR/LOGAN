GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEDICT=file(params.genomes[params.genome].genomedict)
KGPGERMLINE=params.genomes[params.genome].kgp 
DBSNP=file(params.genomes[params.genome].dbsnp) 
GNOMADGERMLINE=params.genomes[params.genome].gnomad 
PON=file(params.genomes[params.genome].pon) 
VEPCACHEDIR=file(params.genomes[params.genome].vepcache)
VEPSPECIES=params.genomes[params.genome].vepspecies
VEPBUILD=params.genomes[params.genome].vepbuild
SOMATIC_FOREST=params.genomes[params.genome].octopus_sforest
GERMLINE_FOREST=params.genomes[params.genome].octopus_gforest

//Output
outdir=file(params.output)


process mutect2 {
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz.stats")

    
    script:
    """
    gatk Mutect2 \
    --reference $GENOMEREF \
    --intervals ${bed} \
    --input ${tumor} \
    --input ${normal} \
    --normal-sample ${normal.simpleName} \
    --tumor-sample ${tumor.simpleName} \
    $GNOMADGERMLINE \
    --panel-of-normals ${PON} \
    --output ${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz \
    --f1r2-tar-gz ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates
    """

    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz
    touch ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz
    touch ${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz.stats
    """
}


process pileup_paired_t {
    label 'process_highmem'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table")

    script:
    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${tumor} \
        -V $KGPGERMLINE \
        -L ${bed} \
        -O ${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table 

    """

    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table
    """

}


process pileup_paired_n {
    label 'process_highmem'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.normal.pileup.table")

    script:
    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${normal} \
        -V $KGPGERMLINE \
        -L ${bed} \
        -O ${tumor.simpleName}_${bed.simpleName}.normal.pileup.table 

    """
    
    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.normal.pileup.table
    """

}


process contamination_paired {
    label 'process_highmem'

    publishDir(path: "${outdir}/vcfs/mutect2", mode: 'copy')

    input:
        tuple val(tumorname),
        path(tumor_pileups),
        path(normal_pileups)
    
    output:
        tuple val(tumorname),
        path("${tumorname}_allpileups.table"),
        path("${tumorname}_normal.allpileups.table"),
        path("${tumorname}.contamination.table"),
        path("${tumorname}_normal.contamination.table")

    script:
    //Gather all the Pileup summaries first for Tumor and Also for NORMAL and then run!
    alltumor = tumor_pileups.join(" -I ")
    allnormal = normal_pileups.join(" -I ")


    """
    gatk GatherPileupSummaries \
    --sequence-dictionary $GENOMEDICT \
    -I ${alltumor} -O ${tumorname}_allpileups.table
    
    gatk GatherPileupSummaries \
    --sequence-dictionary $GENOMEDICT \
    -I ${allnormal} -O ${tumorname}_normal.allpileups.table

    gatk CalculateContamination \
        -I ${tumorname}_allpileups.table \
        --matched-normal ${tumorname}_normal.allpileups.table \
        -O ${tumorname}.contamination.table

    gatk CalculateContamination \
        -I ${tumorname}_normal.allpileups.table \
        -O ${tumorname}_normal.contamination.table

    """

    stub:
    """
    touch ${tumorname}_allpileups.table
    touch ${tumorname}_normal.allpileups.table
    touch ${tumorname}.contamination.table
    touch ${tumorname}_normal.contamination.table
    """

    
}


process learnreadorientationmodel {
    label 'process_highmem'

    publishDir(path: "${outdir}/vcfs/mutect2", mode: 'copy')

    input:
        tuple val(sample), path(f1r2)
      
    output:
        tuple val(sample), path("${sample}.read-orientation-model.tar.gz")

    script: 
    f1r2in = f1r2.join(" --input ")

    """
    gatk LearnReadOrientationModel \
        --output ${sample}.read-orientation-model.tar.gz \
        --input ${f1r2in}
    """

    stub:
    """
    touch ${sample}.read-orientation-model.tar.gz
    """
}


process mergemut2stats {
    label 'process_low'

    publishDir(path: "${outdir}/vcfs/mutect2", mode: 'copy')

    input:
        tuple val(sample), path(stats)
      
    output:
        tuple val(sample), path("${sample}.final.stats")

    script: 
    statsin = stats.join(" --stats ")

    """
    gatk MergeMutectStats \
        --stats ${statsin} \
        -O ${sample}.final.stats
    """

    stub:
    """
    touch ${sample}.final.stats
    """

}


process mutect2filter {
    label 'process_mid'
        
    publishDir(path: "${outdir}/vcfs/mutect2", mode: 'copy')

    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups), path(normal_pileups),path(tumorcontamination),path(normalcontamination)
    output:
        tuple val(sample), path("${sample}.mut2.marked.vcf.gz"), 
        path("${sample}.mut2.norm.vcf.gz"), 
        path("${sample}.mut2.marked.vcf.gz.filteringStats.tsv")

    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")


    """
    gatk GatherVcfs -I ${mut2in} -O ${sample}.concat.vcf.gz 
    gatk IndexFeatureFile -I ${sample}.concat.vcf.gz 
    gatk FilterMutectCalls \
        -R $GENOMEREF \
        -V ${sample}.concat.vcf.gz \
        --ob-priors ${obs} \
        --contamination-table ${tumorcontamination} \
        --stats ${stats} \
        -O ${sample}.mut2.marked.vcf.gz


    gatk SelectVariants \
        -R $GENOMEREF \
        --variant ${sample}.mut2.marked.vcf.gz \
        --exclude-filtered \
        --output ${sample}.mut2.final.vcf.gz
    
    bcftools sort ${sample}.mut2.final.vcf.gz |\
    bcftools norm --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M]\\y/,"N",\$4); OFS = "\\t"; print}}' |\
        sed '/^\$/d' > ${sample}.mut2.norm.vcf |\
    bcftools view - -Oz -o  ${sample}.mut2.norm.vcf.gz
    """

    stub:
    """
    touch ${sample}.mut2.marked.vcf.gz
    touch ${sample}.mut2.norm.vcf.gz
    touch ${sample}.mut2.marked.vcf.gz.filteringStats.tsv
    """


}


process strelka_tn {
    label 'process_highcpu'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz")
    
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
    mv wd/results/variants/somatic.snvs.vcf.gz  ${tumor.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz
    mv wd/results/variants/somatic.indels.vcf.gz  ${tumor.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz

    """

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz
    touch ${tumor.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz

    """


}


process vardict_tn {
    label 'process_highcpu'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.vardict.vcf")
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
            -f 0.05 >  ${tumor.simpleName}_${bed.simpleName}.vardict.vcf

    """

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.vardict.vcf

    """


}


process varscan_tn {
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), 
        val(normalname), path(normal), path(normalbai), path(bed),
        path(tumorpileup), path(normalpileup), path(tumor_con_table), path(normal_con_table)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.varscan.vcf")
    
    shell:

    '''
    tumor_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !{tumor_con_table} | cut -f2 ))" | bc -l)
    normal_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !{normal_con_table} | cut -f2 ))" | bc -l)
    varscan_opts="--strand-filter 1 --min-var-freq 0.01 --min-avg-qual 30 --somatic-p-value 0.05 --output-vcf 1 --normal-purity $normal_purity --tumor-purity $tumor_purity"
    varscan somatic < samtools mpileup -d 10000 -q 15 -Q 15 -f !GENOME -l !{bed.simpleName} !{normal} !{tumor} !{tumor.simpleName}_{bed.simpleName}.vardict.vcf $varscan_opts --mpileup 1 
    '''

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.varscan.vcf
    
    """

}

process octopus_tn {
    //label 'process_highcpu' Using separate docker for octopus

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), 
        val(normalname), path(normal), path(normalbai), path(bed)
    

    output:
        tuple val(tumorname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf")
    
    script:

    """
    octopus -R $GENOMEREF -I ${normal} ${tumor} --normal-sample ${normalname} \
    -C cancer \
    --annotations AC AD DP -t ${bed} \
    --threads $task.cpus \
    $GERMLINE_FOREST \
    $SOMATIC_FOREST \
    -o ${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf 
    """

    stub:
    
    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf"
    """

} 

process lofreq_tn {
    label 'process_somaticcaller' 
    module=["lofreq/2.1.5","bcftools/1.17"]

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), 
        val(normalname), path(normal), path(normalbai), path(bed)
    

    output:
        tuple val(tumorname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.snvs.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.indels.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz")
    
    script:

    """
    lofreq -f $GENOMEREF -n ${normal} -t ${tumor} \
        -d $DBSNP \
        --threads $task.cpus \
        -l ${bed} \
        --call-indels \
        -o ${tumorname}_vs_${normalname}_${bed.simpleName}
    
    bcftools concat ${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz \
        ${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz" --threads $task.cpus -Oz -o \
        ${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz"


    """

    stub:
    
    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.snvs.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz"
    
    """
} 



process muse_tn {
    label 'process_somaticcaller' 
    module=["muse/2.0.1"]

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), 
        val(normalname), path(normal), path(normalbai)
    

    output:
        tuple val(tumorname),
        path("${tumorname}_vs_${normalname}.vcf.gz")
    
    script:

    """
    MuSE call -f $GENOMEREF -O ${tumorname}_vs_${normalname} -n $task.cpus $tumor $normal
    MuSE sump -I ${tumorname}_vs_${normalname}.MuSE.txt \
        -O ${tumorname}_vs_${normalname} -n $task.cpus -D $DBSNP -G
        
    """

    stub:
    
    """
    touch "${tumorname}_vs_${normalname}.vcf.gz"
    """

} 


process combineVariants {
    label 'process_highmem'
    publishDir(path: "${outdir}/vcfs/", mode: 'copy')

    input:
        tuple val(sample), path(inputvcf), val(vc)
    
    output:
        tuple val(sample), 
        path("${vc}/${sample}.${vc}.marked.vcf.gz"), path("${vc}/${sample}.${vc}.norm.vcf.gz")
    
    script:
    vcfin = inputvcf.join(" -I ")
    
    """
    mkdir ${vc}
    gatk --java-options "-Xmx48g" MergeVcfs \
        -O ${sample}.${vc}.temp.vcf.gz \
        -D $GENOMEDICT \
        -I $vcfin
    bcftools sort ${sample}.${vc}.temp.vcf.gz -Oz -o ${sample}.${vc}.marked.vcf.gz
    bcftools norm ${sample}.${vc}.marked.vcf.gz --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M]\\y/,"N",\$4); OFS = "\\t"; print}}' |\
        sed '/^\$/d' > ${sample}.${vc}.temp.vcf

    bcftools view ${sample}.${vc}.temp.vcf -f PASS -Oz -o ${vc}/${sample}.${vc}.norm.vcf.gz

    mv ${sample}.${vc}.marked.vcf.gz ${vc}
    """

    stub:

    """
    mkdir ${vc}
    touch ${vc}/${sample}.${vc}.marked.vcf.gz
    touch ${vc}/${sample}.${vc}.norm.vcf.gz
    
    """

}



process bcftools_index_octopus {
    label 'process_low'

    input:
        tuple val(sample),
        path(vcf)

    output:
        tuple val(sample), 
        path(vcf), 
        path("${vcf}.tbi")
    
    script:    
    """
    bcftools index -t ${vcf}
    """

    stub:
    """
    touch ${vcf}
    touch ${vcf}.tbi
    """

}

process combineVariants_octopus {
    label 'process_highmem'
    publishDir(path: "${outdir}/vcfs/", mode: 'copy')

    input:
        tuple val(sample), path(vcfs), path(vcfsindex), val(vc)
    
    output:
        tuple val(sample), 
        path("${vc}/${sample}.${vc}.marked.vcf.gz"), path("${vc}/${sample}.${vc}.norm.vcf.gz")
    
    script:
    vcfin = vcfs.join(" ")
    
    """
    mkdir ${vc}
    bcftools concat $vcfin -a -Oz -o ${sample}.${vc}.temp.vcf.gz
    bcftools sort ${sample}.${vc}.temp.vcf.gz -Oz -o ${sample}.${vc}.marked.vcf.gz
    bcftools norm ${sample}.${vc}.marked.vcf.gz --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M]\\y/,"N",\$4); OFS = "\\t"; print}}' |\
        sed '/^\$/d' > ${sample}.${vc}.temp.vcf

    bcftools view ${sample}.${vc}.temp.vcf -f PASS -Oz -o ${vc}/${sample}.${vc}.norm.vcf.gz

    mv ${sample}.${vc}.marked.vcf.gz ${vc}
    """

    stub:

    """
    mkdir ${vc}
    touch ${vc}/${sample}.${vc}.marked.vcf.gz
    touch ${vc}/${sample}.${vc}.norm.vcf.gz
    
    """

}






process combineVariants_strelka {
    //Concat all somatic snvs/indels across all files, strelka separates snv/indels
    label 'process_mid'
    publishDir(path: "${outdir}/vcfs/strelka", mode: 'copy')

    input:
        tuple val(sample), path(strelkasnvs), path(strelkaindels)
    
    output:
        tuple val(sample), path("${sample}.strelka.vcf.gz"),path("${sample}.filtered.strelka.vcf.gz")
    
    
    script:
    
    vcfin = strelkasnvs.join(" ")
    indelsin = strelkaindels.join(" ")


    """
    bcftools concat $vcfin $indelsin --threads $task.cpus -Oz -o ${sample}.temp.strelka.vcf.gz
    bcftools sort ${sample}.temp.strelka.vcf.gz -Oz -o ${sample}.strelka.vcf.gz 

    bcftools view ${sample}.strelka.vcf.gz --threads $task.cpus -f PASS -Oz -o ${sample}.filtered.strelka.vcf.gz

    """

    stub:

    """
    touch ${sample}.strelka.vcf.gz
    touch ${sample}.filtered.strelka.vcf.gz
    
    """

}


process annotvep_tn {    
    publishDir(path: "${outdir}/mafs/", mode: 'copy')

    input:
        tuple val(tumorsample), val(normalsample), 
        val(vc), path(tumorvcf) 

    output:
        path("paired/${vc}/${tumorsample}.maf")

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
    --output-maf paired/!{vc}/!{tumorsample}.maf \
    --tumor-id !{tumorsample} \
    --normal-id !{normalsample} \
    --vep-path /opt/vep/src/ensembl-vep \
    --vep-data !{VEPCACHEDIR} \
    --ncbi-build !{VEPBUILD} --species !{VEPSPECIES} --ref-fasta !{GENOMEREF} \
    --vep-overwrite


    '''

    stub:
    """
    mkdir -p paired/${vc}
    touch paired/${vc}/${tumorsample}.maf
    """
}




process combinemafs_tn {
    label 'process_low'
    publishDir(path: "${outdir}/mafs/paired", mode: 'copy')

    input: 
        path(allmafs)

    output:
        path("final_tn.maf")

    shell:
    mafin= allmafs.join(" ")

    """
    echo "Combining MAFs..."
    head -2 ${allmafs[0]} > final_tn.maf
    awk 'FNR>2 {{print}}' ${mafin}  >> final_tn.maf
    """

    stub:
    """
    touch final_tn.maf
    """
}



/*
process combineVariants_allcallers {

    publishDir(path: "${outdir}/vcfs/", mode: 'copy')

    input:
        tuple val(sample), path(inputvcf), val(vc)
    
    output:
        tuple val(sample), 
        path("${vc}/${sample}.${vc}.marked.vcf.gz"), path("${vc}/${sample}.${vc}.norm.vcf.gz")

}
*/