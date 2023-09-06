GENOME=file(params.genome)
GENOMEDICT=file(params.genomedict)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
SHAPEITINDEL=file(params.shapeitindel) //ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
KGP=file(params.kgp) //1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.dbsnp) //dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //somatic-hg38-af-only-gnomad.hg38.vcf.gz
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

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz
    touch ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz
    touch ${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz.stats
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

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table

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
    
    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.normal.pileup.table
    """

}


process contamination_paired {

    publishDir(path: "${outdir}/vcfs/mutect2", mode: 'copy')

    //OUTPUT THE CONTAMINATION TABLE
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

    stub:
    """
    touch ${tumorname}_allpileups.table
    touch ${tumorname}_normal.allpileups.table
    touch ${tumorname}.contamination.table
    touch ${tumorname}_normal.contamination.table
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

    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.tumor.pileup.table

    """

}


process contamination_tumoronly {
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
    --sequence-dictionary ${GENOMEDICT} \
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


process learnreadorientationmodel {
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


process learnreadorientationmodel_tonly {
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



process mergemut2stats_tonly {
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
        
    publishDir(path: "${outdir}/vcfs/mutect2", mode: 'copy')

    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups), path(normal_pileups),path(tumorcontamination),path(normalcontamination)
    output:
        tuple val(sample), path("${sample}.concat.vcf.gz"),path("${sample}.marked.vcf.gz"),path("${sample}.final.mut2.vcf.gz"),path("${sample}.marked.vcf.gz.filteringStats.tsv")

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

    stub:
    """
    touch ${sample}.concat.vcf.gz
    touch ${sample}.marked.vcf.gz
    touch ${sample}.final.mut2.vcf.gz
    touch ${sample}.marked.vcf.gz.filteringStats.tsv
    """


}



process strelka_tn {
    
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz")
    
    script:

    """
    mkdir -p wd

    configureStrelkaSomaticWorkflow.py \
        --ref=${GENOME} \
        --tumor=${tumor} \
        --normal=${normal} \
        --runDir=wd \
        --callRegions ${bed}
    ./wd/runWorkflow.py -m local -j 
    mv wd/results/variants/somatic.snvs.vcf.gz  ${tumor.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz
    mv wd/results/variants/somatic.indels.vcf.gz  ${tumor.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz

    """

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz
    touch ${tumor.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz

    """


}

process combineVariants_strelka {
    //Concat all somatic snvs/indesl across all beds

    publishDir(path: "${outdir}/vcfs/strelka", mode: 'copy')

    input:
        tuple val(sample), path(strelkasnvs), path(strelkaindels)
    
    output:
        tuple val(sample), path("${sample}.final.strelka.vcf.gz")
    
    script:
    
    vcfin = strelkavcfs.join(" ")
    indelsin = strelkaindels.join(" ")


    """
    bcftools concat $vcfin $indelsin -Oz -o ${sample}.final.strelka.vcf.gz
    bcftools view ${sample}.final.strelka.vcf.gz -f PASS -Oz -o ${sample}.unsorted.filtered.strelka.vcf.gz
    bcftools sort ${sample}.unsorted.filtered.strelka.vcf.gz -@ 16

    """

    stub:

    """
    touch ${sample}.final.strelka.vcf.gz
    touch ${sample}.filtered.strelka.vcf.gz
    
    """


}


process vardict {
    
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.vardict.vcf"),
    
    script:

    """
    VarDict -G ${GENOME} \
        -f 0.05 \
=        --nosv \
        -b ${tumor}|${normal} \
        -t -Q 20 -c 1 -S 2 -E 3 
        ${bed} \
        | teststrandbias.R \
        | var2vcf_paired.pl \
            -N ${tumor}|${normal} \
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



process vardict_tonly {
    
    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.vardict.vcf"),
    
    script:

    """
    VarDict -G ${GENOME} \
        -f 0.05 \
        -x 500 \
        --nosv \
        -b ${tumor} \
        -t -Q 20 -c 1 -S 2 -E 3 
        ${bed} \
        | teststrandbias.R \
        | var2vcf_valid.pl \
            -N ${tumor} \
            -Q 20 \
            -d 10 \
            -v 6 \
            -S \
            -E \
            -f 0.05 >  ${tumor.simpleName}_${bed.simpleName}.vardict.vcf

    """

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.vardict.vcf

    """


}



process combineVariants_vardict {
    publishDir(path: "${outdir}/vcfs/vardict", mode: 'copy')

    input:
        tuple val(sample), path(vardictout)
    
    output:
        tuple val(sample), path("${sample}.final.vardict.vcf.gz")
    
    script:
    
    vcfin = vardictout.join(" ")

    //Concat all somatic snvs/indesl across all beds
    """

    bcftools concat $vcfin -Oz -o ${sample}.final.strelka.vcf.gz
    bcftools view ${sample}.final.vardict.vcf.gz -f PASS -Oz -o ${sample}.unsorted.filtered.vardict.vcf.gz
    bcftools sort ${sample}.unsorted.filtered.strelka.vcf.gz -@ 16 -Oz -o ${sample}.filtered.vardict.vcf.gz

    """

    stub:

    """
    touch ${sample}.final.vardict.vcf.gz
    touch ${sample}.filtered.vardict.vcf.gz
    
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

    stub:
    """
    touch ${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz
    touch ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz
    touch ${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz.stats
    """


}



process mutect2filter_tonly {
    publishDir(path: "${outdir}/vcfs/mutect2", mode: 'copy')

    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups),path(tumorcontamination)
    output:
        tuple val(sample), path("${sample}.tonly.concat.vcf.gz"), path("${sample}.tonly.marked.vcf.gz"),path("${sample}.tonly.final.mut2.vcf.gz"),path("${sample}.tonly.marked.vcf.gz.filteringStats.tsv")
    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")


    """
    gatk GatherVcfs -I ${mut2in} -O ${sample}.tonly.concat.vcf.gz 
    gatk IndexFeatureFile -I ${sample}.tonly.concat.vcf.gz 
    gatk FilterMutectCalls \
        -R ${GENOME} \
        -V ${sample}.tonly.concat.vcf.gz \
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

    stub:
    """
    touch ${sample}.tonly.concat.vcf.gz
    touch ${sample}.tonly.marked.vcf.gz
    touch ${sample}.tonly.final.mut2.vcf.gz
    touch ${sample}.tonly.marked.vcf.gz.filteringStats.tsv
    """
}



process varscan_tn {
    module=['samtools/1.9','VarScan/2.4.6']

    input:
        tuple val(tumorname),
        path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai), path(bed)
        path(allpileup.table), path(normalallpileups.table), path(tumor.contamination.table), path(normal.contamination.table)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.vardict.vcf"),
    
    shell:

    """
    tumor_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !normal.contamination.table | cut -f2 ))" | bc -l)
    normal_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !normal.contamination.table | cut -f2 ))" | bc -l)
    varscan_opts="--strand-filter 1 --min-var-freq 0.01 --min-avg-qual 30 --somatic-p-value 0.05 --output-vcf 1 --normal-purity \$normal_purity --tumor-purity \$tumor_purity"
    varscan somatic < samtools mpileup -d 10000 -q 15 -Q 15 -f !GENOME -l !{bed.simpleName} !normal !tumor !{tumor.simpleName}_{bed.simpleName}.vardict.vcf $varscan_opts --mpileup 1 
    """

    stub:
    
    """
    touch ${tumor.simpleName}_${bed.simpleName}.varscan.vcf
    
    """

}

process combineVariants_varscan {

    publishDir(path: "${outdir}/vcfs/varscan", mode: 'copy')

    input:
        tuple val(sample), path(varscanout)
    
    output:
        tuple val(sample), path("${sample}.final.varscan.vcf.gz")
    
    script:
    
    vcfin = varscan.join(" ")
    //Concat all somatic snvs/indesl across all beds
    """
    bcftools concat $vcfin -Oz -o ${sample}.final.strelka.vcf.gz
    bcftools view ${sample}.final.varscan.vcf.gz -f PASS -Oz -o ${sample}.unsorted.filtered.varscan.vcf.gz
    bcftools sort ${sample}.unsorted.filtered.varscan.vcf.gz -@ 16 -Oz -o ${sample}.filtered.varscan.vcf.gz

    """

    stub:

    """
    touch ${sample}.final.vardict.vcf.gz
    touch ${sample}.filtered.vardict.vcf.gz
    
    """


}


process annotvep_tn {
    module=['vcf2maf/1.6.21','VEP/102']
    
    publishDir(path: "${outdir}/mafs/", mode: 'copy')

    input:
        tuple val(tumorsample), 
        path("${sample}.concat.vcf.gz"),
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

    stub:
    """
    touch ${tumorsample}.maf
    """
}


process annotvep_tonly {
    module=['vcf2maf/1.6.21','VEP/102']

    publishDir("${outdir}/mafs", mode: "copy")

    input:
        tuple val(tumorsample), 
        path("${sample}.tonly.concat.vcf.gz"),
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

    stub:
    """
    touch ${tumorsample}.tonly.maf
    """
}