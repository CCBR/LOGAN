GENOME=file(params.genome)
GENOMEDICT=file(params.genomedict)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"// file(params.gold_indels1) //
SHAPEITINDEL=file(params.shapeitindel) //params.shapeitindel =  "/data/OpenOmics/references/genome-seek/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz" //file(params.gold_indels2) //
KGP=file(params.kgp) ///data/CCBR_Pipeliner/Exome-seek/hg38/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.dbsnp) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //= '/data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz' // /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 

//Output
outdir=file(params.output)


process mutect2 {
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
    --reference ${GENOME} \
    --intervals ${bed} \
    --input ${tumor} \
    --input ${normal} \
    --normal-sample ${normal.simpleName} \
    --tumor-sample ${tumor.simpleName} \
    --germline-resource ${GNOMAD} \
    --panel-of-normals ${PON} \
    --output ${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz \
    --f1r2-tar-gz ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates
    """
}


process pileup_paired_t {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table")

    script:

    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${tumor} \
        -V ${KGP} \
        -L ${bed} \
        -O ${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table 

    """
}


process pileup_paired_n {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.normal.pileup.table")

    script:

    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${normal} \
        -V ${KGP} \
        -L ${bed} \
        -O ${tumor.simpleName}_${bed.simpleName}.normal.pileup.table 

    """
}


process contamination_paired {
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
    --sequence-dictionary ${GENOMEDICT} \
    -I ${alltumor} -O ${tumorname}_allpileups.table
    
    gatk GatherPileupSummaries \
    --sequence-dictionary ${GENOMEDICT} \
    -I ${allnormal} -O ${tumorname}_normal.allpileups.table

    gatk CalculateContamination \
        -I ${tumorname}_allpileups.table \
        --matched-normal ${tumorname}_normal.allpileups.table \
        -O ${tumorname}.contamination.table
    gatk CalculateContamination \
        -I ${tumorname}_normal.allpileups.table \
        -O ${tumorname}_normal.contamination.table

    """
}

process pileup_paired_tonly {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table")

    script:

    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${tumor} \
        -V ${KGP} \
        -L ${bed} \
        -O ${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table 

    """
}


process contamination_tumoronly {
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
    --sequence-dictionary ${GENOMEDICT} \
    -I ${alltumor} -O ${tumorname}_allpileups.table
    
    gatk CalculateContamination \
        -I ${tumorname}_allpileups.table \
        -O ${tumorname}.contamination.table

    """
}


process learnreadorientationmodel {
    input:
    tuple val(sample),
    path(f1r2)
      
    output:
    tuple val(sample),
    path("${sample}.read-orientation-model.tar.gz")

    script: 
    f1r2in = f1r2.join(" --input ")

    """
    gatk LearnReadOrientationModel \
        --output ${sample}.read-orientation-model.tar.gz \
        --input ${f1r2in}
    """

}


process learnreadorientationmodel_tonly {
    input:
    tuple val(sample),
    path(f1r2)
      
    output:
    tuple val(sample),
    path("${sample}.read-orientation-model.tar.gz")

    script: 
    f1r2in = f1r2.join(" --input ")

    """
    gatk LearnReadOrientationModel \
        --output ${sample}.read-orientation-model.tar.gz \
        --input ${f1r2in}
    """

}



process mergemut2stats {
    input:
    tuple val(sample),
    path(stats)
      
    output:
    tuple val(sample),
    path("${sample}.final.stats")

    script: 
    statsin = stats.join(" --stats ")

    """
    gatk MergeMutectStats \
        --stats ${statsin} \
        -O ${sample}.final.stats
    """

}



process mergemut2stats_tonly {
    input:
    tuple val(sample),
    path(stats)
      
    output:
    tuple val(sample),
    path("${sample}.final.stats")

    script: 
    statsin = stats.join(" --stats ")

    """
    gatk MergeMutectStats \
        --stats ${statsin} \
        -O ${sample}.final.stats
    """

}



process mutect2filter {
    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups), path(normal_pileups),path(tumorcontamination),path(normalcontamination)
    output:
        tuple val(sample), path("${sample}.marked.vcf.gz"),path("${sample}.final.mut2.vcf.gz"),path("${sample}.marked.vcf.gz.filteringStats.tsv")
    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")


    """
    gatk GatherVcfs -I ${mut2in} -O ${sample}.concat.vcf.gz 
    gatk IndexFeatureFile -I ${sample}.concat.vcf.gz 
    gatk FilterMutectCalls \
        -R ${GENOME} \
        -V ${sample}.concat.vcf.gz \
        --ob-priors ${obs} \
        --contamination-table ${tumorcontamination} \
        --stats ${stats} \
        -O ${sample}.marked.vcf.gz


    gatk SelectVariants \
        -R ${GENOME} \
        --variant ${sample}.marked.vcf.gz \
        --exclude-filtered \
        --output ${sample}.final.mut2.vcf.gz
    """
}



process mutect2_t_tonly {
    
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
    --reference ${GENOME} \
    --intervals ${bed} \
    --input ${tumor} \
    --tumor-sample ${tumor.simpleName} \
    --germline-resource ${GNOMAD} \
    --panel-of-normals ${PON} \
    --output ${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz \
    --f1r2-tar-gz ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates    
    """
}


process mutect2filter_tonly {
    publishDir(path: "${results_dir}/vcfs/mutect", mode: 'copy')

    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups),path(tumorcontamination)
    output:
        tuple val(sample), path("${sample}.tonly.marked.vcf.gz"),path("${sample}.tonly.final.mut2.vcf.gz"),path("${sample}.tonly.marked.vcf.gz.filteringStats.tsv")
    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")


    """
    gatk GatherVcfs -I ${mut2in} -O ${sample}.concat.vcf.gz 
    gatk IndexFeatureFile -I ${sample}.concat.vcf.gz 
    gatk FilterMutectCalls \
        -R ${GENOME} \
        -V ${sample}.concat.vcf.gz \
        --ob-priors ${obs} \
        --contamination-table ${tumorcontamination} \
        --stats ${stats} \
        -O ${sample}.tonly.marked.vcf.gz


    gatk SelectVariants \
        -R ${GENOME} \
        --variant ${sample}.tonly.marked.vcf.gz \
        --exclude-filtered \
        --output ${sample}.tonly.final.mut2.vcf.gz
    """
}




process annotvep_tn {
    module=['vcf2maf/1.6.21','VEP/102']
    
    publishDir(path: "${results_dir}/mafs/", mode: 'copy')

    input:
        tuple val(tumorsample), 
        path("${tumorsample}.marked.vcf.gz"), 
        path("${tumorsample}.final.mut2.vcf.gz"), 
        path("${tumorsample}.marked.vcf.gz.filteringStats.tsv"), 
        val(normalsample)

    output:
        path("${tumorsample}.maf")

    script:

    """
    
    zcat ${tumorsample}.final.mut2.vcf.gz > ${tumorsample}.final.mut2.vcf

    vcf2maf.pl \
    --vep-forks 16 --input-vcf ${tumorsample}.final.mut2.vcf \
    --output-maf ${tumorsample}.maf \
    --tumor-id ${tumorsample} \
    --normal-id ${normalsample} \
    --vep-path \${VEP_HOME}/bin \
    --vep-data \${VEP_CACHEDIR} \
    --ncbi-build GRCh38 --species homo_sapiens --ref-fasta ${GENOME}

    """
}


process annotvep_tonly {
    module=['vcf2maf/1.6.21','VEP/102']

    publishDir("${results_dir}/mafs/", mode: "copy")

    input:
        tuple val(tumorsample), 
        path("${tumorsample}.tonly.marked.vcf.gz"),
        path("${tumorsample}.tonly.final.mut2.vcf.gz"),
        path("${tumorsample}.tonly.marked.vcf.gz.filteringStats.tsv")

    output:
        path("${tumorsample}.tonly.maf")

    script:

    """
    
    zcat ${tumorsample}.tonly.final.mut2.vcf.gz  > ${tumorsample}.tonly.final.mut2.vcf

    vcf2maf.pl \
    --vep-forks 16 --input-vcf ${tumorsample}.tonly.final.mut2.vcf \
    --output-maf ${tumorsample}.tonly.maf \
    --tumor-id ${tumorsample} \
    --vep-path \${VEP_HOME}/bin \
    --vep-data \${VEP_CACHEDIR} \
    --ncbi-build GRCh38 --species homo_sapiens --ref-fasta ${GENOME}

    """
}