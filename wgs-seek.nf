#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

GENOME=file(params.genome)
GENOMEDICT=file(params.genomedict)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"// file(params.gold_indels1) //
SHAPEITINDEL=file(params.shapeitindel) //params.shapeitindel =  "/data/OpenOmics/references/genome-seek/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz" //file(params.gold_indels2) //
KGP=file(params.kgp) ///data/CCBR_Pipeliner/Exome-seek/hg38/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.dbsnp) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //= '/data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz' // /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 


results_dir=params.output


fastqinput=Channel.fromFilePairs(params.fastqs,checkIfExists: true)
intervalbed = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor,
                       row.Normal
                       )
                                  }
//Final Workflow
workflow {
    fastqinput.view()
    fastp(fastqinput)
    bwamem2(fastp.out)
    indelrealign(bwamem2.out)
    indelbambyinterval=indelrealign.out.combine(intervalbed)
    bqsr(indelbambyinterval)
    bqsrs=bqsr.out.groupTuple()
        .map { samplename,beds -> tuple( samplename, 
        beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )
        }
    gatherbqsr(bqsrs)

    tobqsr=indelrealign.out.combine(gatherbqsr.out,by:0)
    applybqsr(tobqsr) 
    samtoolsindex(applybqsr.out)
    
    bamwithsample=samtoolsindex.out.join(sample_sheet).map{it.swap(3,0)}.join(samtoolsindex.out).map{it.swap(3,0)}

    bambyinterval=bamwithsample.combine(intervalbed)

    //Paired Mutect2    
    mutect2(bambyinterval)
    pileup_paired_t(bambyinterval)
    pileup_paired_n(bambyinterval)
    
    pileup_paired_tout=pileup_paired_t.out.groupTuple()
    .map{samplename,pileups-> tuple( samplename,
    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tumor.pileup.table/)[0][1].toInteger() } ,
    )}
    pileup_paired_nout=pileup_paired_n.out.groupTuple()
    .map{samplename,pileups-> tuple( samplename,
    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).normal.pileup.table/)[0][1].toInteger() } ,
    )}


    pileup_paired_all=pileup_paired_tout.join(pileup_paired_nout)
    contamination_paired(pileup_paired_all)

    mut2out_lor=mutect2.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    f1r2.toSorted{ it -> (it.name =~ /${samplename}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } 
    )}

    learnreadorientationmodel(mut2out_lor)

    mut2out_mstats=mutect2.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    stats.toSorted{ it -> (it.name =~ /${samplename}_(.*?).mut2.vcf.gz.stats/)[0][1].toInteger() } 
    )}

    mergemut2stats(mut2out_mstats)

    allmut2tn=mutect2.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    vcfs.toSorted{ it -> (it.name =~ /${samplename}_(.*?).mut2.vcf.gz/)[0][1].toInteger() } 
    )}
    
    mut2tn_filter=allmut2tn
    .join(mergemut2stats.out)
    .join(learnreadorientationmodel.out)
    .join(contamination_paired.out)
    mutect2filter(mut2tn_filter)



    //Tumor Only Calling

    mutect2_tonly(bambyinterval)    
    
    //LOR     
    mut2tout_lor=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    f1r2.toSorted{ it -> (it.name =~ /${samplename}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } 
    )}
    learnreadorientationmodel_tonly(mut2tout_lor)


    //Stats
    mut2tonly_mstats=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    stats.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tonly.mut2.vcf.gz.stats/)[0][1].toInteger() } 
    )}
    mergemut2stats_tonly(mut2out_mstats)


    //Contamination
    contamination_tumoronly(pileup_paired_tout)


    
    //Final TUMOR ONLY FILTER
    allmut2tonly=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    vcfs.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tonly.mut2.vcf.gz/)[0][1].toInteger() } 
    )}
    
    mut2tonly_filter=allmut2tonly
    .join(mergemut2stats_tonly.out)
    .join(learnreadorientationmodel_tonly.out)
    .join(contamination_tumoronly.out)

    mutect2filter_tonly(mut2tonly_filter)


        
    
    //#To implement
        //CNMOPs from the BAM BQSRs
        //##VCF2MAF TO
    tn_vepin=mutect2filter.out
    .join(sample_sheet)

   annotvep_tn(tn_vepin)
   annotvep_tonly(mutect2filter_tonly.out)

}

process fastp{
    tag { name }
    module 'fastp/0.23.2'

    input:
    tuple val(samplename), path(fqs)
    
    output:
    tuple val(samplename),
    path("${samplename}.R1.trimmed.fastq.gz"),
    path("${samplename}.R2.trimmed.fastq.gz"),
    path("${samplename}.fastp.json"),
    path("${samplename}.fastp.html")

    script:
    """
    fastp -w 4 \
        --detect_adapter_for_pe \
        --in1 ${fqs[0]} \
        --in2 ${fqs[1]} \
        --out1 ${samplename}.R1.trimmed.fastq.gz \
        --out2 ${samplename}.R2.trimmed.fastq.gz  \
        --json ${samplename}.fastp.json \
        --html ${samplename}.fastp.html
    """

}


process bwamem2{
    tag { name }
    module=['bwa-mem2/2.2.1','samblaster/0.1.26','samtools/1.15.1']
    input:
        tuple val(samplename), 
        path("${samplename}.R1.trimmed.fastq.gz"),
        path("${samplename}.R2.trimmed.fastq.gz"),
        path("${samplename}.fastp.json"),
        path("${samplename}.fastp.html")
        
    output:
        tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")

    script:
    //BWAmem2/sambalster/samtools sort for marking duplicates;
    """

     bwa-mem2 mem -M \
        -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:illumina\\tLB:${samplename}\\tPU:${samplename}\\tCN:hgsc\\tDS:wes' \
        -t 12 \
        ${GENOME} \
        ${samplename}.R1.trimmed.fastq.gz ${samplename}.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -@12 -m 4G - -o ${samplename}.bam

    samtools index -@ 8 ${samplename}.bam ${samplename}.bai

    """
}

process indelrealign {
    tag { name }
    module=['GATK/3.8-1']
    
    input:
    tuple val(samplename), path("${samplename}.bam"), path("${samplename}.bai")

    output:
    tuple val(samplename), path("${samplename}.ir.bam")

    script: 
    
    """
    /usr/bin/java -Xmx32g -jar \${GATK_JAR} -T RealignerTargetCreator \
        -I  ${samplename}.bam \
        -R ${GENOME} \
        -o ${samplename}.intervals \
        -known ${MILLSINDEL} -known ${SHAPEITINDEL} 
    
    /usr/bin/java -Xmx32g -jar \${GATK_JAR} -T IndelRealigner \
        -R ${GENOME} \
        -I ${samplename}.bam \
        -known ${MILLSINDEL} -known ${SHAPEITINDEL} \
        --use_jdk_inflater \
        --use_jdk_deflater \
        -targetIntervals ${samplename}.intervals \
        -o  ${samplename}.ir.bam
    """
    
}



process bqsr {
    /*
    Base quality recalibration for all samples to 
    */    
    input:
        tuple val(samplename), path("${samplename}.ir.bam"), path(bed)

    output:
        tuple val(samplename),path("${samplename}_${bed.simpleName}.recal_data.grp"),emit: bqsrby
        //path("${bam.simpleName}_${bed.simpleName}.recal_data.grp"), emit: bqsrby

    script:
    """
    gatk --java-options '-Xmx32g' BaseRecalibrator \
    --input ${samplename}.ir.bam \
    --reference ${GENOME} \
    --known-sites ${MILLSINDEL} --known-sites ${SHAPEITINDEL} \
    --output ${samplename}_${bed.simpleName}.recal_data.grp \
    --intervals ${bed}

    """
}

process gatherbqsr {

    input: 
        tuple val(samplename), path(recalgroups)
    output:
        tuple val(samplename), path("${samplename}.recal_data.grp")
    script:
    
    strin = recalgroups.join(" --input ")

    """
    gatk --java-options '-Xmx32g' GatherBQSRReports \
    --input ${strin} \
    --output ${samplename}.recal_data.grp

    """
}


process applybqsr {
    /*
    Base quality recalibration for all samples to 
    */    
    input:
        tuple val(samplename),path("${samplename}.ir.bam"), path("${samplename}.recal_data.grp")

    output:
        tuple val(samplename),path("${samplename}.bqsr.bam")

    script:
    """

    gatk --java-options '-Xmx32g' ApplyBQSR \
        --reference ${GENOME} \
        --input ${samplename}.ir.bam \
        --bqsr-recal-file  ${samplename}.recal_data.grp \
        --output ${samplename}.bqsr.bam \
        --use-jdk-inflater \
        --use-jdk-deflater

    """
}



process samtoolsindex{
  input:
    tuple val(bamname), path(bam)
    
    output:
    tuple val(bamname), path(bam), path("${bam}.bai")

    script:
    """
    samtools index -@ 4 ${bam} ${bam}.bai
    """

}


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



process mutect2_tonly {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
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
    module=['vcf2maf/1.6.21','VEP/106']
    
    publishDir("${results_dir}/mafs/", mode: "copy")

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
    module=['vcf2maf/1.6.21','VEP/106']

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