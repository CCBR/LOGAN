GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEFAI=file(params.genomes[params.genome].genomefai)
GENOMEDICT=file(params.genomes[params.genome].genomedict)
KGPGERMLINE=params.genomes[params.genome].kgp
DBSNP=file(params.genomes[params.genome].dbsnp)
GNOMADGERMLINE=params.genomes[params.genome].gnomad
PON=file(params.genomes[params.genome].pon)
VEPCACHEDIR=file(params.genomes[params.genome].vepcache)
VEPSPECIES=params.genomes[params.genome].vepspecies
VEPBUILD=params.genomes[params.genome].vepbuild
LOFREQ_CONVERT=params.lofreq_convert
//Octopus
SOMATIC_FOREST=params.genomes[params.genome].octopus_sforest
GERMLINE_FOREST=params.genomes[params.genome].octopus_gforest
//HMFTOOLS
HOTSPOTS=params.genomes[params.genome].HOTSPOTS
PANELBED=params.genomes[params.genome].PANELBED
HCBED=params.genomes[params.genome].HCBED
ENSEMBLCACHE=params.genomes[params.genome].ENSEMBLCACHE
GENOMEVER=params.genomes[params.genome].GENOMEVER

process mutect2 {
    container "${params.containers.logan}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai),
        path(bed)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.mut2.vcf.gz"),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.f1r2.tar.gz"),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.mut2.vcf.gz.stats")


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
    --output ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.mut2.vcf.gz \
    --f1r2-tar-gz ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates
    """

    stub:
    """
    touch ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.mut2.vcf.gz
    touch ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.f1r2.tar.gz
    touch ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.mut2.vcf.gz.stats
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
    container "${params.containers.logan}"
    label 'process_highmem'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(tumorname),
        val(normalname),
        path("${normal.simpleName}_${bed.simpleName}.normal.pileup.table")

    script:
    """
    gatk --java-options -Xmx48g GetPileupSummaries \
        -I ${normal} \
        -V $KGPGERMLINE \
        -L ${bed} \
        -O ${normalname}_${bed.simpleName}.normal.pileup.table

    """

    stub:
    """
    touch ${normalname}_${bed.simpleName}.normal.pileup.table
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

process octopus_convertvcf {
    container "${params.containers.logan}"
    label 'process_low'
    
    input:
        tuple val(tumor), val(normal), 
        val(oct), path(vcf), path(vcfindex)

    output:
        tuple val(tumor), val(normal), path("${tumor}.octopus.norm.vcf.gz"), 
        path("${tumor}.octopus.norm.vcf.gz.tbi")


    script:
    """
    zcat ${vcf}  | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > ${tumor}_temp.octopus.norm.vcf
    bgzip ${tumor}_temp.octopus.norm.vcf
    mv ${tumor}_temp.octopus.norm.vcf.gz ${tumor}.octopus.norm.vcf.gz
    bcftools index -t ${tumor}.octopus.norm.vcf.gz -f
    """

    stub:
    """
    touch ${tumor}.octopus.norm.vcf.gz ${tumor}.octopus.norm.vcf.gz.tbi
    """
}

process mutect2filter {
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumor), val(normal),path(mutvcfs), path(stats), path(obs),
        path(pileups), path(normal_pileups),path(tumorcontamination),path(normalcontamination)

    output:
        tuple val("${tumor}_vs_${normal}"),
        path("${tumor}_vs_${normal}.mut2.marked.vcf.gz"),
        path("${tumor}_vs_${normal}.mut2.marked.vcf.gz.tbi"),
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


process strelka_tn {
    container "${params.containers.logan}"
    label 'process_highcpu'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz"),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz.tbi"),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz"),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz.tbi")

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
    mv wd/results/variants/somatic.snvs.vcf.gz  ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic_temp.snvs.vcf.gz
    mv wd/results/variants/somatic.indels.vcf.gz  ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic_temp.indels.vcf.gz

    printf "NORMAL\t${normalname}\nTUMOR\t${tumorname}\n" >sampname

    bcftools reheader -s sampname ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic_temp.snvs.vcf.gz \
        | bcftools view -Oz -o ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz
    bcftools reheader -s sampname ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic_temp.indels.vcf.gz \
        | bcftools view -Oz -o ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz

    bcftools index -t ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz
    bcftools index -t ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz

    """

    stub:

    """
    touch ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz  ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.snvs.vcf.gz.tbi
    touch ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.somatic.indels.vcf.gz.tbi

    """

}


process vardict_tn {
    container "${params.containers.logan}"
    label 'process_somaticcaller_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai), val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.vardict.vcf.gz")
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
            -f 0.05 >  ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.vardict.vcf

    printf "${normal.Name}\t${normalname}\n${tumor.Name}\t${tumorname}\n" > sampname

    bcftools reheader -s sampname ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.vardict.vcf \
        | bcftools view -Oz -o ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.vardict.vcf.gz


    """

    stub:

    """
    touch ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.vardict.vcf.gz

    """


}


process varscan_tn {
    container "${params.containers.logan}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed),
        path(tumorpileup), path(normalpileup),
        path(tumor_con_table), path(normal_con_table)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.varscan.vcf.gz")

    shell:
    '''
    tumor_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !{tumor_con_table} | cut -f2 ))" | bc -l)
    normal_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !{normal_con_table} | cut -f2 ))" | bc -l)
    dual_pileup="samtools mpileup -d 10000 -q 15 -Q 15 -f !{GENOMEREF} -l !{bed} !{normal} !{tumor}"
    varscan_opts="--strand-filter 1 --min-var-freq 0.01 --min-avg-qual 30 --somatic-p-value 0.05 --output-vcf 1 --normal-purity $normal_purity --tumor-purity $tumor_purity"
    varscan_cmd="varscan somatic <($dual_pileup) !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.vcf $varscan_opts --mpileup 1"
    eval "$varscan_cmd"

    awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.vcf.indel \
        | sed '/^$/d' | bcftools view - -Oz -o !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.indel_temp.vcf.gz
    awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.vcf.snp \
        | sed '/^$/d' | bcftools view - -Oz -o !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.snp_temp.vcf.gz

    gatk SortVcf -I !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.snp_temp.vcf.gz \
    -I !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.indel_temp.vcf.gz \
    -R !{GENOMEREF} -SD !{GENOMEDICT} \
    -O !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}_temp.varscan.vcf

    printf "NORMAL\t!{normalname}\nTUMOR\t!{tumorname}\n" > sampname

    bcftools reheader -s sampname !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}_temp.varscan.vcf \
       | bcftools view -Oz -o !{tumor.simpleName}_vs_!{normal.simpleName}_!{bed.simpleName}.varscan.vcf.gz

    '''

    stub:
    """
    touch ${tumor.simpleName}_vs_${normal.simpleName}_${bed.simpleName}.varscan.vcf.gz
    """

}


process octopus_tn {
    container "${params.containers.octopus}"
    label 'process_somaticcaller_high'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val("${tumorname}_vs_${normalname}"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf.gz")

    script:
    """
    octopus -R $GENOMEREF -I ${normal} ${tumor} --normal-sample ${normalname} \
    -C cancer \
    --annotations AC AD DP -t ${bed} \
    --threads $task.cpus \
    $GERMLINE_FOREST \
    $SOMATIC_FOREST \
    -B 92Gb \
    -o ${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf.gz
    """

    stub:
    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}.octopus.vcf.gz"
    """

}


process sage_tn {
    container "${params.containers.hmftools}"
    label 'process_somaticcaller'

     input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

 output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_sage.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_sage.vcf.gz.tbi")

script:
    """
    java -Xms4G -Xmx32G -cp sage.jar com.hartwig.hmftools.sage.SageApplication \
    -tumor $tumorname -tumor_bam $tumorbam \
    -reference $normalname -reference_bam $normalbam \
    -threads $task.cpus \
    -ref_genome_version $GENOMEVER \
    -ref_genome $GENOMEREF \
    $HOTSPOTS $PANELBED $HCBED $ENSEMBLCACHE \
    -output_vcf ${tumorname}_vs_${normalname}_${bed.simpleName}.sage.vcf.gz

    """

    stub:

    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}.sage.vcf.gz" "${tumorname}_vs_${normalname}_${bed.simpleName}.sage.vcf.gz.tbi"
    """
}


process lofreq_tn {
    container "${params.containers.logan}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)


    output:

        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.snvs.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.indels.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz"),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz.tbi")

    script:

    """
    lofreq somatic -f $GENOMEREF -n ${normal} -t ${tumor} \
        -d $DBSNP \
        --threads $task.cpus \
        -l ${bed} \
        --call-indels \
        -o ${tumorname}_vs_${normalname}_${bed.simpleName}_

    bcftools concat ${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz \
        ${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz --threads $task.cpus -Oz -o \
        ${tumorname}_vs_${normalname}_${bed.simpleName}_temp_lofreq.vcf.gz

    $LOFREQ_CONVERT -i ${tumorname}_vs_${normalname}_${bed.simpleName}_temp_lofreq.vcf.gz -g 1/0 \
        -n ${tumorname} -o ${tumorname}_vs_${normalname}_${bed.simpleName}_temp1_lofreq.vcf.gz

    bcftools view -h ${tumorname}_vs_${normalname}_${bed.simpleName}_temp1_lofreq.vcf.gz >temphead

    sed 's/^##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">/##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref\\/fwd, ref\\/rev, var\\/fwd, var\\/rev">/' temphead > temphead1
    bcftools reheader ${tumorname}_vs_${normalname}_${bed.simpleName}_temp1_lofreq.vcf.gz -h temphead1 |\
        bcftools view -Oz -o ${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz

    bcftools index -t ${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz

    """

    stub:

    """
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.snvs.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.snvs.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_somatic_final_minus-dbsnp.indels.vcf.gz"
    touch "${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz" "${tumorname}_vs_${normalname}_${bed.simpleName}_lofreq.vcf.gz.tbi"

    """
}



process muse_tn {
    container "${params.containers.logan}"
    label 'process_somaticcaller'
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)


    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}.vcf.gz")

    script:

    """
    MuSE call -f $GENOMEREF -O ${tumorname}_vs_${normalname} -n $task.cpus $tumor $normal
    MuSE sump -I ${tumorname}_vs_${normalname}.MuSE.txt \
        -O ${tumorname}_vs_${normalname}.vcf -n $task.cpus -D $DBSNP -G

    bcftools view ${tumorname}_vs_${normalname}.vcf -Oz -o ${tumorname}_vs_${normalname}_temp.vcf.gz

    printf "NORMAL\t${normalname}\nTUMOR\t${tumorname}\n" > sampname

    bcftools reheader -s sampname ${tumorname}_vs_${normalname}_temp.vcf.gz \
        | bcftools view -Oz -o ${tumorname}_vs_${normalname}.vcf.gz

    """

    stub:

    """
    touch "${tumorname}_vs_${normalname}.vcf.gz"
    """

}


process combineVariants {
    container "${params.containers.logan}"
    label 'process_highmem'

    input:
        tuple val(sample), path(inputvcf), val(vc)

    output:
        tuple val(sample),
        path("${vc}/${sample}.${vc}.marked.vcf.gz"),
        path("${vc}/${sample}.${vc}.marked.vcf.gz.tbi"),
        path("${vc}/${sample}.${vc}.norm.vcf.gz"),
        path("${vc}/${sample}.${vc}.norm.vcf.gz.tbi")

    script:
    vcfin = inputvcf.join(" -I ")

    """
    mkdir ${vc}
    gatk --java-options "-Xmx48g" SortVcf \
        -O ${sample}.${vc}.marked.vcf.gz \
        -SD $GENOMEDICT \
        -I $vcfin
    bcftools norm ${sample}.${vc}.marked.vcf.gz -m- --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' > ${sample}.${vc}.temp.vcf

    bcftools view ${sample}.${vc}.temp.vcf -f PASS -Oz -o ${vc}/${sample}.${vc}.norm.vcf.gz

    mv ${sample}.${vc}.marked.vcf.gz ${vc}
    mv ${sample}.${vc}.marked.vcf.gz.tbi ${vc}

    bcftools index ${vc}/${sample}.${vc}.norm.vcf.gz -t
    """

    stub:

    """
    mkdir ${vc}
    touch ${vc}/${sample}.${vc}.marked.vcf.gz
    touch ${vc}/${sample}.${vc}.norm.vcf.gz
    touch ${vc}/${sample}.${vc}.marked.vcf.gz.tbi
    touch ${vc}/${sample}.${vc}.norm.vcf.gz.tbi
    """

}




process combineVariants_alternative {
    container "${params.containers.logan}"
    label 'process_highmem'

    input:
        tuple val(sample), path(vcfs), path(vcfsindex), val(vc)

    output:
        tuple val(sample),
        path("${vc}/${sample}.${vc}.marked.vcf.gz"),
        path("${vc}/${sample}.${vc}.marked.vcf.gz.tbi"),
        path("${vc}/${sample}.${vc}.norm.vcf.gz"),
        path("${vc}/${sample}.${vc}.norm.vcf.gz.tbi")

    script:
    vcfin = vcfs.join(" ")

    """
    mkdir ${vc}
    bcftools concat $vcfin -a -Oz -o ${sample}.${vc}.temp1.vcf.gz
    bcftools reheader -f $GENOMEFAI ${sample}.${vc}.temp1.vcf.gz -o ${sample}.${vc}.temp.vcf
    bcftools sort ${sample}.${vc}.temp.vcf -Oz -o ${sample}.${vc}.marked.vcf.gz
    bcftools norm ${sample}.${vc}.marked.vcf.gz -m- --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' > ${sample}.${vc}.temp.vcf

    bcftools view ${sample}.${vc}.temp.vcf -f PASS -Oz -o ${vc}/${sample}.${vc}.norm.vcf.gz

    mv ${sample}.${vc}.marked.vcf.gz ${vc}

    bcftools index ${vc}/${sample}.${vc}.marked.vcf.gz -t
    bcftools index ${vc}/${sample}.${vc}.norm.vcf.gz -t
    """

    stub:

    """
    mkdir ${vc}
    touch ${vc}/${sample}.${vc}.marked.vcf.gz
    touch ${vc}/${sample}.${vc}.norm.vcf.gz
    touch ${vc}/${sample}.${vc}.marked.vcf.gz.tbi
    touch ${vc}/${sample}.${vc}.norm.vcf.gz.tbi

    """

}


process bcftools_index_octopus {
    container "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(tumor),
        path(vcf)

    output:
        tuple val(tumor),
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


process combineVariants_strelka {
    //Concat all somatic snvs/indels across all files, strelka separates snv/indels
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(sample),
        path(strelkasnvs), path(snvindex),
        path(strelkaindels), path(indelindex)

    output:
        tuple val(sample),
        path("${sample}.strelka.vcf.gz"), path("${sample}.strelka.vcf.gz.tbi"),
        path("${sample}.filtered.strelka.vcf.gz"), path("${sample}.filtered.strelka.vcf.gz.tbi")


    script:

    vcfin = strelkasnvs.join(" ")
    indelsin = strelkaindels.join(" ")


    """
    bcftools concat $vcfin $indelsin --threads $task.cpus -Oz -o ${sample}.temp.strelka.vcf.gz -a
    bcftools norm ${sample}.temp.strelka.vcf.gz -m- --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' > ${sample}.temp1.strelka.vcf.gz

    bcftools sort ${sample}.temp1.strelka.vcf.gz -Oz -o ${sample}.strelka.vcf.gz

    bcftools view ${sample}.strelka.vcf.gz --threads $task.cpus -f PASS -Oz -o ${sample}.filtered.strelka.vcf.gz

    bcftools index ${sample}.strelka.vcf.gz -t
    bcftools index ${sample}.filtered.strelka.vcf.gz -t
    """

    stub:

    """
    touch ${sample}.strelka.vcf.gz ${sample}.strelka.vcf.gz.tbi
    touch ${sample}.filtered.strelka.vcf.gz ${sample}.filtered.strelka.vcf.gz.tbi

    """

}


process somaticcombine {
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumorsample), val(normal),
        val(caller),
        path(vcfs), path(vcfindex)

    output:
        tuple val(tumorsample), val(normal),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz"),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz.tbi")

    script:
        vcfin1=[caller, vcfs].transpose().collect { a, b -> a + " " + b }
        vcfin2="-V:" + vcfin1.join(" -V:")

        callerin=caller.join(",")
    """
    /usr/lib/jvm/java-8-openjdk-amd64/bin/java -jar \$GATK_JAR -T CombineVariants  \
        -R $GENOMEREF \
        --genotypemergeoption PRIORITIZE \
        --rod_priority_list $callerin \
        --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
        -o ${tumorsample}_vs_${normal}_combined.vcf.gz \
        $vcfin2
        
    """

    stub:
    vcfin1=[caller, vcfs].transpose().collect { a, b -> a + " " + b }
    vcfin2="-V:" + vcfin1.join(" -V:")

    callerin=caller.join(",")

    """
    touch ${tumorsample}_vs_${normal}_combined.vcf.gz
    touch ${tumorsample}_vs_${normal}_combined.vcf.gz.tbi
    """

}



/*DISCVR
process somaticcombine {
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumorsample), val(normal),
        val(callers),
        path(vcfs), path(vcfindex)

    output:
        tuple val(tumorsample), val(normal),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz"),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz.tbi")

    script:
        vcfin1=[callers, vcfs].transpose().collect { a, b -> a + " " + b }
        vcfin2="-V:" + vcfin1.join(" -V:")

    """
    java -jar \$DISCVRSeq_JAR MergeVcfsAndGenotypes \
        -R $GENOMEREF \
        --genotypeMergeOption PRIORITIZE \
        --priority_list mutect2,strelka,octopus,muse,lofreq,vardict,varscan \
        --filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED \
        -O ${tumorsample}_vs_${normal}_combined.vcf.gz \
        $vcfin2
    """

    stub:
    vcfin1=[callers, vcfs].transpose().collect { a, b -> a + " " + b }
    vcfin2="-V:" + vcfin1.join(" -V:")

    """
    touch ${tumorsample}_vs_${normal}_combined.vcf.gz
    touch ${tumorsample}_vs_${normal}_combined.vcf.gz.tbi
    """

}
*/

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


process combinemafs_tn {
    container "${params.containers.logan}"
    label 'process_low'

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
