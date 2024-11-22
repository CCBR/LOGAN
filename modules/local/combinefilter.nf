//References
GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEFAI=file(params.genomes[params.genome].genomefai)
GENOMEDICT=file(params.genomes[params.genome].genomedict)


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
    //Create Tumor Normal here
    samplist=sample.split('_vs_')
    if(samplist.size()>1){
        samporder = samplist.join(",")
    }else{
        samporder = sample
    }

    """
    mkdir ${vc}
    gatk --java-options "-Xmx48g" SortVcf \
        -O ${sample}.${vc}.markedtemp.vcf.gz \
        -SD $GENOMEDICT \
        -I $vcfin
    
    bcftools view ${sample}.${vc}.markedtemp.vcf.gz -s $samporder -Oz -o ${sample}.${vc}.marked.vcf.gz 
    bcftools index -t ${sample}.${vc}.marked.vcf.gz 

    bcftools norm ${sample}.${vc}.marked.vcf.gz -m- --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' > ${sample}.${vc}.temp.vcf

    bcftools view ${sample}.${vc}.temp.vcf -f PASS -s $samporder -Oz -o ${vc}/${sample}.${vc}.norm.vcf.gz
    bcftools index ${vc}/${sample}.${vc}.norm.vcf.gz -t

    mv ${sample}.${vc}.marked.vcf.gz ${vc}
    mv ${sample}.${vc}.marked.vcf.gz.tbi ${vc}

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
    samplist=sample.split('_vs_')
    if (vc.contains("lofreq") | vc.contains('deepsomatic')) {
        samporder = samplist[0]
    }else if(samplist.size()>1){
        samporder = samplist.join(",")
    }else{
        samporder = sample
    }

    if (vc.contains("octopus")) {
        """
        mkdir ${vc}
        bcftools concat $vcfin -a -Oz -o ${sample}.${vc}.temp1.vcf.gz
        bcftools reheader -f $GENOMEFAI ${sample}.${vc}.temp1.vcf.gz -o ${sample}.${vc}.temp.vcf
        bcftools sort ${sample}.${vc}.temp.vcf | bcftools view - -i "INFO/SOMATIC==1" -s $samporder -Oz -o ${sample}.${vc}.marked.vcf.gz
        bcftools norm ${sample}.${vc}.marked.vcf.gz -m- --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' > ${sample}.${vc}.temp.vcf

        bcftools view ${sample}.${vc}.temp.vcf -f PASS -Oz -o ${vc}/${sample}.${vc}.norm.vcf.gz
        mv ${sample}.${vc}.marked.vcf.gz ${vc}

        bcftools index ${vc}/${sample}.${vc}.marked.vcf.gz -t
        bcftools index ${vc}/${sample}.${vc}.norm.vcf.gz -t
        """
    
    }else{
        """
        mkdir ${vc}
        bcftools concat $vcfin -a -Oz -o ${sample}.${vc}.temp1.vcf.gz
        bcftools reheader -f $GENOMEFAI ${sample}.${vc}.temp1.vcf.gz -o ${sample}.${vc}.temp.vcf
        bcftools sort ${sample}.${vc}.temp.vcf  | bcftools view - -s $samporder -Oz -o ${sample}.${vc}.marked.vcf.gz
        bcftools norm ${sample}.${vc}.marked.vcf.gz -m- --threads $task.cpus --check-ref s -f $GENOMEREF -O v |\
        awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' |\
        sed '/^\$/d' > ${sample}.${vc}.temp.vcf

        bcftools view ${sample}.${vc}.temp.vcf -f PASS -Oz -o ${vc}/${sample}.${vc}.norm.vcf.gz
        mv ${sample}.${vc}.marked.vcf.gz ${vc}

        bcftools index ${vc}/${sample}.${vc}.marked.vcf.gz -t
        bcftools index ${vc}/${sample}.${vc}.norm.vcf.gz -t
        """
    }
   
    stub:

    """
    mkdir ${vc}
    touch ${vc}/${sample}.${vc}.marked.vcf.gz
    touch ${vc}/${sample}.${vc}.norm.vcf.gz
    touch ${vc}/${sample}.${vc}.marked.vcf.gz.tbi
    touch ${vc}/${sample}.${vc}.norm.vcf.gz.tbi

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



process combinemafs_tonly {
    container "${params.containers.logan}"
    label 'process_low'

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




process somaticcombine_tonly {
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumorsample),
        val(caller),
        path(vcfs), path(vcfindex)

    output:
        tuple val(tumorsample),
        path("${tumorsample}_combined_tonly.vcf.gz"),
        path("${tumorsample}_combined_tonly.vcf.gz.tbi")

    script:
        vcfin1=[caller, vcfs].transpose().collect { a, b -> a + " " + b }
        vcfin2="-V:" + vcfin1.join(" -V:")

        callerin=caller.join(",")//.replaceAll("_tonly","")

    """
    /usr/lib/jvm/java-8-openjdk-amd64/bin/java -jar \$GATK_JAR -T CombineVariants  \
        -R $GENOMEREF \
        --genotypemergeoption PRIORITIZE \
        --rod_priority_list $callerin \
        --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
        -o ${tumorsample}_combined_tonly.vcf.gz \
        $vcfin2
    """

    stub:

        vcfin1=[caller, vcfs].transpose().collect { a, b -> a + " " + b }
        vcfin2="-V:" + vcfin1.join(" -V:")
        callerin=caller.join(",")//.replaceAll("_tonly","")
    
    """
    touch ${tumorsample}_combined_tonly.vcf.gz ${tumorsample}_combined_tonly.vcf.gz.tbi
    """ 

}

