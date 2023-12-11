GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEFAI=file(params.genomes[params.genome].genomefai)
GENOMEDICT=file(params.genomes[params.genome].genomedict)
KGPGERMLINE=params.genomes[params.genome].kgp //1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.genomes[params.genome].dbsnp) //dbsnp_138.hg38.vcf.gz"
GNOMADGERMLINE=params.genomes[params.genome].gnomad //somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.genomes[params.genome].pon) 
VEPCACHEDIR=file(params.genomes[params.genome].vepcache)
VEPSPECIES=params.genomes[params.genome].vepspecies
VEPBUILD=params.genomes[params.genome].vepbuild
SOMATIC_FOREST=params.genomes[params.genome].octopus_sforest
GERMLINE_FOREST=params.genomes[params.genome].octopus_gforest

//Output
outdir=file(params.output)



process pileup_paired_tonly {
    label 'process_highmem'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)
    
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


process contamination_tumoronly {
    label 'process_highmem'
    publishDir(path: "${outdir}/vcfs/mutect2/", mode: 'copy')

    input:
        tuple val(tumorname),
        path(tumor_pileups)

    output:
        tuple val(tumorname),
        path("${tumorname}_allpileups.table"),
        path("${tumorname}.contamination.table")

    script:
    //Gather all the Pileup summaries first for Tumor and Also for NORMAL and then run!
    alltumor = tumor_pileups.join(" -I ")


    """
    gatk GatherPileupSummaries \
    --sequence-dictionary $GENOMEDICT \
    -I ${alltumor} -O ${tumorname}_allpileups.table
    
    gatk CalculateContamination \
        -I ${tumorname}_allpileups.table \
        -O ${tumorname}.contamination.table

    """

    stub:
    """
    touch ${tumorname}_allpileups.table
    touch ${tumorname}.contamination.table
    """

}



process learnreadorientationmodel_tonly {
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





process mergemut2stats_tonly {
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



process mutect2_t_tonly {
    label 'process_somaticcaller'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz.stats")
    
    script:

    """
    gatk Mutect2 \
    --reference $GENOMEREF \
    --intervals ${bed} \
    --input ${tumor} \
    --tumor-sample ${tumor.simpleName} \
    $GNOMADGERMLINE \
    --panel-of-normals $PON \
    --output ${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz \
    --f1r2-tar-gz ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates    
    """

    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz
    touch ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz
    touch ${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz.stats
    """


}



process mutect2filter_tonly {
    label 'process_mid'
    publishDir(path: "${outdir}/vcfs/mutect2_tonly", mode: 'copy')

    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups),path(tumorcontamination)
    output:
        tuple val(sample), 
        path("${sample}.tonly.mut2.marked.vcf.gz"),path("${sample}.tonly.mut2.marked.vcf.gz.tbi"), 
        path("${sample}.tonly.mut2.norm.vcf.gz"),path("${sample}.tonly.mut2.norm.vcf.gz.tbi"), 
        path("${sample}.tonly.mut2.marked.vcf.gz.filteringStats.tsv")

    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")


    """
    gatk GatherVcfs -I ${mut2in} -O ${sample}.tonly.concat.vcf.gz 
    gatk IndexFeatureFile -I ${sample}.tonly.concat.vcf.gz 
    gatk FilterMutectCalls \
        -R $GENOMEREF \
        -V ${sample}.tonly.concat.vcf.gz \
        --ob-priors ${obs} \
        --contamination-table ${tumorcontamination} \
        --stats ${stats} \
        -O ${sample}.tonly.mut2.marked.vcf.gz

    gatk SelectVariants \
        -R $GENOMEREF \
        --variant ${sample}.tonly.mut2.marked.vcf.gz \
        --exclude-filtered \
        --output ${sample}.tonly.mut2.final.vcf.gz

    bcftools sort ${sample}.tonly.mut2.final.vcf.gz |\
    bcftools norm --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' |\
    bcftools view - -Oz -o  ${sample}.tonly.mut2.norm.vcf.gz
    bcftools index -t ${sample}.tonly.mut2.norm.vcf.gz

    """

    stub:
    """
    touch ${sample}.tonly.mut2.marked.vcf.gz ${sample}.tonly.mut2.marked.vcf.gz.tbi
    touch ${sample}.tonly.mut2.norm.vcf.gz ${sample}.tonly.mut2.norm.vcf.gz.tbi
    touch ${sample}.tonly.mut2.marked.vcf.gz.filteringStats.tsv
    """
}


process varscan_tonly {
    label 'process_somaticcaller'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), 
        path(bed),
        path(tumorpileup),  path(tumor_con_table)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.tonly.varscan.vcf.gz")
    
    shell:

    '''
    varscan_opts="--strand-filter 0 --min-var-freq 0.01 --output-vcf 1 --variants 1"
    pileup_cmd="samtools mpileup -d 100000 -q 15 -Q 15 -f !{GENOMEREF} !{tumor}"
    varscan_cmd="varscan mpileup2cns <($pileup_cmd) $varscan_opts"

    eval "$varscan_cmd > !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf"

    printf "TUMOR\t!{tumorname}\n" > sampname 
    
    bcftools reheader -s sampname !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf \
        | bcftools view -Oz -o !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf.gz

    '''

    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.tonly.varscan.vcf.gz
    """

}


process vardict_tonly {
    label 'process_highcpu'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.tonly.vardict.vcf.gz")
    
    script:

    """
    bedtools makewindows -b ${bed} -w 50150 -s 50000 > temp_${bed}

    VarDict -G $GENOMEREF \
        -f 0.01 \
        -x 500 \
        --nosv \
        -b ${tumor} --fisher \
        -t -Q 20 -c 1 -S 2 -E 3 --th $task.cpus \
        temp_${bed} | var2vcf_valid.pl \
            -N ${tumor} \
            -Q 20 \
            -d 10 \
            -v 6 \
            -S \
            -E \
            -f 0.05 >  ${tumor.simpleName}_${bed.simpleName}.tonly.vardict.vcf

    printf "${tumor.Name}\t${tumorname}\n" > sampname 
    
    bcftools reheader -s sampname ${tumor.simpleName}_${bed.simpleName}.tonly.vardict.vcf \
        | bcftools view -Oz -o ${tumor.simpleName}_${bed.simpleName}.tonly.vardict.vcf.gz

    """

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.tonly.vardict.vcf.gz

    """

}


process octopus_tonly {
    //label 'process_highcpu'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumorname}_${bed.simpleName}.tonly.octopus.vcf.gz")
    
    script:

    """
    octopus -R $GENOMEREF -C cancer -I ${tumor} \
    --annotations AC AD DP \
    --target-working-memory 64Gb \
    -t ${bed} \
    $SOMATIC_FOREST \
    -o ${tumorname}_${bed.simpleName}.tonly.octopus.vcf.gz --threads $task.cpus


    """

    stub:
    
    """
    touch ${tumorname}_${bed.simpleName}.tonly.octopus.vcf.gz
    """
}



process somaticcombine_tonly {
    label 'process_mid'
    publishDir(path: "${outdir}/vcfs/combined_tonly", mode: 'copy')

    input: 
        tuple val(tumorsample), 
        val(callers),
        path(vcfs), path(vcfindex)

    output:
        tuple val(tumorsample),
        path("${tumorsample}_combined_tonly.vcf.gz"),
        path("${tumorsample}_combined_tonly.vcf.gz.tbi")

    script:
        vcfin1=[callers, vcfs].transpose().collect { a, b -> a + " " + b }
        vcfin2="-V:" + vcfin1.join(" -V:")

    """
    java -jar \$DISCVRSeq_JAR MergeVcfsAndGenotypes \
        -R $GENOMEREF \
        --genotypeMergeOption PRIORITIZE \
        --priority_list mutect2_tonly,octopus_tonly,vardict_tonly,varscan_tonly \
        --filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED \
        -O ${tumorsample}_combined_tonly.vcf.gz \
        $vcfin2
    """

    stub:
    """
    touch ${tumorsample}_combined_tonly.vcf.gz ${tumorsample}_combined_tonly.vcf.gz.tbi
    """

}

process annotvep_tonly {
    publishDir("${outdir}/mafs", mode: "copy")

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
    --vep-overwrite


    '''

    stub:
    """
    mkdir -p tumor_only/${vc}
    touch tumor_only/${vc}/${tumorsample}.tonly.maf
    """
}

process combinemafs_tonly {
    label 'process_low'
    publishDir(path: "${outdir}/mafs/tumor_only", mode: 'copy')

    input: 
        path(allmafs)

    output:
        path("final_tonly.maf")

    shell:
    mafin= allmafs.join(" ")
    
    """
    echo "Combining MAFs..."
    head -2 ${allmafs[0]} > final_tonly.maf
    awk 'FNR>2 {{print}}' ${mafin}  >> final_tonly.maf
    """

    stub:
    """
    touch final_tonly.maf
    """
}



