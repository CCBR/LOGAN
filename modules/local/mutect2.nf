//References
GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEDICT=file(params.genomes[params.genome].genomedict)
GERMLINE_RESOURCE=file(params.genomes[params.genome].germline_resource)
GNOMADGERMLINE=params.genomes[params.genome].gnomad
//PON Mutect2
PON=file(params.genomes[params.genome].PON)
TONLYPON=file(params.genomes[params.genome].tonly_PON)



process mutect2 {
    container "${params.containers.logan}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai),
        path(bed)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.mut2.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.f1r2.tar.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.mut2.vcf.gz.stats")


    script:
    """
    gatk Mutect2 \
    --reference $GENOMEREF \
    --intervals ${bed} \
    --input ${tumor} \
    --input ${normal} \
    --normal-sample ${normalname} \
    --tumor-sample ${tumorname} \
    $GNOMADGERMLINE \
    --panel-of-normals ${PON} \
    --output ${tumorname}_vs_${normalname}_${bed.simpleName}.mut2.vcf.gz \
    --f1r2-tar-gz ${tumorname}_vs_${normalname}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates
    """

    stub:
    """
    touch ${tumorname}_vs_${normalname}_${bed.simpleName}.mut2.vcf.gz
    touch ${tumorname}_vs_${normalname}_${bed.simpleName}.f1r2.tar.gz
    touch ${tumorname}_vs_${normalname}_${bed.simpleName}.mut2.vcf.gz.stats
    """
}

process pileup_paired {
    container "${params.containers.logan}"
    label 'process_highmem'

    input:
        tuple val(tumorname),
        path(bam), path(bai),
        path(bed), val(pilename)

    output:
        tuple val(tumorname),
        path("${tumorname}_${bed.simpleName}.${pilename}.table")

    script:
    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${bam} \
        -V $GERMLINE_RESOURCE \
        -L ${bed} \
        -O ${tumorname}_${bed.simpleName}.${pilename}.table

    """

    stub:
    """
    touch ${tumorname}_${bed.simpleName}.${pilename}.table
    """

}


process pileup_paired_t {
    container "${params.containers.logan}"
    label 'process_highmem'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_${bed.simpleName}.tpileup.table")

    script:
    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${tumor} \
        -V $GERMLINE_RESOURCE \
        -L ${bed} \
        -O ${tumorname}_${bed.simpleName}.tpileup.table

    """

    stub:
    """
    touch ${tumorname}_${bed.simpleName}.tpileup.table
    """

}

process pileup_paired_n {
    container "${params.containers.logan}"
    label 'process_highmem'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(tumorname),
        val(normalname),
        path("${normalname}_${bed.simpleName}.npileup.table")

    script:
    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${normal} \
        -V $GERMLINE_RESOURCE \
        -L ${bed} \
        -O ${normalname}_${bed.simpleName}.npileup.table

    """

    stub:
    """
    touch ${normalname}_${bed.simpleName}.npileup.table
    """
}


process contamination_paired {
    container "${params.containers.logan}"
    label 'process_highmem'

    input:
        tuple val(tumorname), val(normalname),
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
    container "${params.containers.logan}"
    label 'process_highmem'

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
    container "${params.containers.logan}"
    label 'process_low'

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
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumor), val(normal),path(mutvcfs), path(stats), path(obs),
        path(pileups), path(normal_pileups), path(tumorcontamination), path(normalcontamination)

    output:
        tuple val("${tumor}_vs_${normal}"),
        path("${tumor}_vs_${normal}.mut2.marked.vcf.gz"), path("${tumor}_vs_${normal}.mut2.marked.vcf.gz.tbi"),
        path("${tumor}_vs_${normal}.mut2.norm.vcf.gz"), path("${tumor}_vs_${normal}.mut2.norm.vcf.gz.tbi"),
        path("${tumor}_vs_${normal}.mut2.marked.vcf.gz.filteringStats.tsv")

    script:
    mut2in = mutvcfs.join(" -I ")

    """
    gatk SortVcf -I ${mut2in} -O ${tumor}_vs_${normal}.concat.vcf.gz --CREATE_INDEX
    gatk FilterMutectCalls \
        -R $GENOMEREF \
        -V ${tumor}_vs_${normal}.concat.vcf.gz \
        --ob-priors ${obs} \
        --contamination-table ${tumorcontamination} \
        --stats ${stats} \
        -O ${tumor}_vs_${normal}.mut2.marked.vcf.gz
    gatk SelectVariants \
        -R $GENOMEREF \
        --variant ${tumor}_vs_${normal}.mut2.marked.vcf.gz \
        --exclude-filtered \
        --output ${tumor}_vs_${normal}.mut2.final.vcf.gz

    bcftools sort ${tumor}_vs_${normal}.mut2.final.vcf.gz |\
    bcftools norm --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' | bcftools view - -Oz -o ${tumor}_vs_${normal}.mut2.norm.vcf.gz
    bcftools index -t ${tumor}_vs_${normal}.mut2.norm.vcf.gz
    """

    stub:
    """
    touch ${tumor}_vs_${normal}.mut2.marked.vcf.gz ${tumor}_vs_${normal}.mut2.marked.vcf.gz.tbi
    touch ${tumor}_vs_${normal}.mut2.norm.vcf.gz ${tumor}_vs_${normal}.mut2.norm.vcf.gz.tbi
    touch ${tumor}_vs_${normal}.mut2.marked.vcf.gz.filteringStats.tsv
    """


}

process pileup_paired_tonly {
    container "${params.containers.logan}"

    label 'process_highmem'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)

    output:
        tuple val(tumorname),
        path("${tumorname}_${bed.simpleName}.tpileup.table")

    script:

    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${tumor} \
        -V $GERMLINE_RESOURCE \
        -L ${bed} \
        -O ${tumorname}_${bed.simpleName}.tpileup.table

    """

    stub:
    """
    touch ${tumorname}_${bed.simpleName}.tpileup.table

    """

}


process contamination_tumoronly {
    container "${params.containers.logan}"

    label 'process_highmem'

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
    container "${params.containers.logan}"

    label 'process_highmem'

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
    container "${params.containers.logan}"
    label 'process_low'

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
    container "${params.containers.logan}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), path(bed)

    output:
        tuple val(tumorname),
        path("${tumorname}_${bed.simpleName}.tonly.mut2.vcf.gz"),
        path("${tumorname}_${bed.simpleName}.f1r2.tar.gz"),
        path("${tumorname}_${bed.simpleName}.tonly.mut2.vcf.gz.stats")

    script:

    """
    gatk Mutect2 \
    --reference $GENOMEREF \
    --intervals ${bed} \
    --input ${tumor} \
    --tumor-sample ${tumorname} \
    $GNOMADGERMLINE \
    --panel-of-normals $TONLYPON \
    --output ${tumorname}_${bed.simpleName}.tonly.mut2.vcf.gz \
    --f1r2-tar-gz ${tumorname}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates
    """

    stub:
    """
    touch ${tumorname}_${bed.simpleName}.tonly.mut2.vcf.gz
    touch ${tumorname}_${bed.simpleName}.f1r2.tar.gz
    touch ${tumorname}_${bed.simpleName}.tonly.mut2.vcf.gz.stats
    """


}



process mutect2filter_tonly {
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups), path(tumorcontamination)
    output:
        tuple val(sample),
        path("${sample}.tonly.mut2.marked.vcf.gz"),path("${sample}.tonly.mut2.marked.vcf.gz.tbi"),
        path("${sample}.tonly.mut2.norm.vcf.gz"),path("${sample}.tonly.mut2.norm.vcf.gz.tbi"),
        path("${sample}.tonly.mut2.marked.vcf.gz.filteringStats.tsv")

    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")

    """
    gatk SortVcf -I ${mut2in} -O ${sample}.tonly.concat.vcf.gz --CREATE_INDEX
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
    bcftools norm --threads ${task.cpus} --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
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



