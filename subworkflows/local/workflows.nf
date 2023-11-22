//All Worksflows in One Place         
intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')


include {fc_lane; fastq_screen;kraken;qualimap_bamqc;fastqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;
    somalier_extract;somalier_analysis;multiqc} from  '../../modules/local/qc.nf'
include {deepvariant_step1;deepvariant_step2;deepvariant_step3;
    deepvariant_combined;glnexus} from '../../modules/local/germline.nf'
include {fastp; bwamem2; //indelrealign; 
    bqsr; gatherbqsr; applybqsr; samtoolsindex} from '../../modules/local/trim_align.nf'
include {mutect2; mutect2filter; pileup_paired_t; pileup_paired_n; 
    contamination_paired; learnreadorientationmodel;mergemut2stats;
    strelka_tn; combineVariants_strelka; 
    varscan_tn; vardict_tn;
    combineVariants as combineVariants_vardict; combineVariants as combineVariants_varscan; 
    combineVariants as combineVariants_vardict_tonly; combineVariants as combineVariants_varscan_tonly
    annotvep_tn as annotvep_tn_mut2; annotvep_tn as annotvep_tn_strelka; annotvep_tn as annotvep_tn_varscan; annotvep_tn as annotvep_tn_vardict;
    combinemafs_tn} from '../../modules/local/variant_calling.nf'
include {mutect2_t_tonly; mutect2filter_tonly; 
    varscan_tonly; vardict_tonly; 
    contamination_tumoronly;
    learnreadorientationmodel_tonly; 
    mergemut2stats_tonly;
    annotvep_tonly as annotvep_tonly_varscan; annotvep_tonly as annotvep_tonly_vardict; annotvep_tonly as annotvep_tonly_mut2;
    combinemafs_tonly} from '../../modules/local/variant_calling_tonly.nf'
include {svaba_somatic} from '../../modules/local/structural_variant.nf'
include {splitinterval} from "../../modules/local/splitbed.nf"



workflow INPUT_PIPE {
    
    if(params.fastq_input){
        fastqinput=Channel.fromFilePairs(params.fastq_input)
    }else if(params.file_input) {
        fastqinput=Channel.fromPath(params.file_input)
                        .splitCsv(header: false, sep: "\t", strip:true)
                        .map{ sample,fq1,fq2 -> 
                        tuple(sample, tuple(file(fq1),file(fq2)))
                                  }
    }
    
    if(params.sample_sheet){
        sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true).view()
                       .ifEmpty { "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t", strip:true)
                       .map { row -> tuple(
                        row.Tumor,
                        row.Normal
                       )
                                  }
    }else{
        sample_sheet=fastqinput.map{samplename,f1 -> tuple (
             samplename)}
        
    }

    emit:
        fastqinput
        sample_sheet

}

workflow TRIM_ALIGN_PIPE {
    take:
        fastqinput
        sample_sheet
    main: 
    fastp(fastqinput)
    splitinterval(intervalbedin)
    
    bwamem2(fastp.out)

    //indelrealign(bwamem2.out)

    //indelbambyinterval=indelrealign.out.combine(splitinterval.out.flatten())
    bqsrbambyinterval=bwamem2.out.combine(splitinterval.out.flatten())
    bambyinterval=bwamem2.out.combine(splitinterval.out.flatten())
    
        
    bqsr(bqsrbambyinterval)
    bqsrs=bqsr.out.groupTuple()
        .map { samplename,beds -> tuple( samplename, 
        beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )
        }
    gatherbqsr(bqsrs)

    tobqsr=bwamem2.out.combine(gatherbqsr.out,by:0)
    applybqsr(tobqsr) 
    //samtoolsindex(applybqsr.out)
    
    //samtoolsindex.out.view()
    //sample_sheet.view()
    bamwithsample=applybqsr.out.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(applybqsr.out,by:0).map{it.swap(3,0)}

    emit:
        bamwithsample
        bambyinterval
        fastpout=fastp.out
        fastqin=fastqinput
        splitout=splitinterval.out
        //indelbambyinterval
        bqsrbambyinterval
        sample_sheet
        bqsrout=applybqsr.out
}

workflow GERMLINE_PIPE {
    //GERMLINE REQUIRES only BAMBYINTERVAL
    take:
        bambyinterval
    main:
    deepvariant_step1(bambyinterval) 
    deepvariant_1_sorted=deepvariant_step1.out.groupTuple()
        .map { samplename,tfbeds,gvcfbed -> tuple( samplename, 
        tfbeds.toSorted{ it -> (it.name =~ /${samplename}.tfrecord_(.*?).bed.gz/)[0][1].toInteger() } ,
        gvcfbed.toSorted{ it -> (it.name =~ /${samplename}.gvcf.tfrecord_(.*?).bed.gz/)[0][1].toInteger() } )
        }
    deepvariant_step2(deepvariant_1_sorted) | deepvariant_step3 
    glin=deepvariant_step3.out.map{samplename,vcf,vcf_tbi,gvcf,gvcf_tbi -> gvcf}.collect()   

    glnexus(glin)
    emit:
        glnexusout=glnexus.out
        bcfout=deepvariant_step3.out
    
}
    
workflow VARIANTCALL_PIPE {
    take:
    //Input is the BAMby interval
        bamwithsample
        splitout
        sample_sheet
        
    main: 
    bambyinterval=bamwithsample.combine(splitout.flatten())

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
    bambyinterval_t=bambyinterval.map{tumorname,tumor,tumorbai,normalname,normalbam,normalbai,bed ->tuple(tumorname,tumor,tumorbai,bed)}
    mutect2_t_tonly(bambyinterval_t)    
    
    //LOR     
    mut2tout_lor=mutect2_t_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    f1r2.toSorted{ it -> (it.name =~ /${samplename}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } 
    )}
    learnreadorientationmodel_tonly(mut2tout_lor)


    //Stats
    mut2tonly_mstats=mutect2_t_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    stats.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tonly.mut2.vcf.gz.stats/)[0][1].toInteger() } 
    )}
    mergemut2stats_tonly(mut2tonly_mstats)


    //Contamination
    contamination_tumoronly(pileup_paired_tout)
    
    //Final TUMOR ONLY FILTER
    allmut2tonly=mutect2_t_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    vcfs.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tonly.mut2.vcf.gz/)[0][1].toInteger() } 
    )}
    
    mut2tonly_filter=allmut2tonly
    .join(mergemut2stats_tonly.out)
    .join(learnreadorientationmodel_tonly.out)
    .join(contamination_tumoronly.out)

    mutect2filter_tonly(mut2tonly_filter)

    //VCF2maf
    mutect2filter.out
    .join(sample_sheet)
    .map{tumor,markedvcf,finalvcf,stats,normal -> tuple(tumor,normal,"mutect2",finalvcf)} | annotvep_tn_mut2

    mutect2filter_tonly.out
    .join(sample_sheet)
    .map{tumor,markedvcf,finalvcf,stats,normal -> tuple(tumor,"mutect2",finalvcf)} | annotvep_tonly_mut2

    //Strelka
    strelka_tn(bambyinterval)
    strelkaout=strelka_tn.out.groupTuple()
    .map { samplename,vcfs,indels -> tuple( samplename,
    vcfs.toSorted{ it -> (it.name =~ /${samplename}_(.*?).somatic.snvs.vcf.gz/)[0][1].toInteger() },
    indels.toSorted{ it -> (it.name =~ /${samplename}_(.*?).somatic.indels.vcf.gz/)[0][1].toInteger() }  
    )}
    combineVariants_strelka(strelkaout)
    combineVariants_strelka.out.join(sample_sheet)
    .map{tumor,markedvcf,finalvcf,normal -> tuple(tumor,normal,"strelka",finalvcf)} | annotvep_tn_strelka

    //Vardict
    vardict_comb=vardict_tn(bambyinterval).groupTuple().map{tumor,vcf-> tuple(tumor,vcf,"vardict")} | combineVariants_vardict
    vardict_comb.join(sample_sheet)
     .map{tumor,marked,normvcf,normal ->tuple(tumor,normal,"vardict",normvcf)} | annotvep_tn_vardict

     //VarDict_tonly
    vardict_tonly_comb=bambyinterval.map{tumorname,tumorbam,tumorbai,normname,normbam,normbai,bed ->
        tuple(tumorname,tumorbam,tumorbai,bed)} 
    vardict_tonly(vardict_tonly_comb).groupTuple().map{tumor,vcf-> tuple(tumor,vcf,"vardict_tonly")} |combineVariants_vardict_tonly
    combineVariants_vardict_tonly.out.join(sample_sheet)
    .map{tumor,marked,normvcf,normal ->tuple(tumor,"vardict_tonly",normvcf)} | annotvep_tonly_vardict


    //VarScan
    varscan_in=bambyinterval.join(contamination_paired.out)
    varscan_comb=varscan_tn(varscan_in).groupTuple().map{tumor,vcf-> tuple(tumor,vcf,"varscan")} | combineVariants_varscan
    varscan_comb.join(sample_sheet)
    .map{tumor,marked,normvcf,normal ->tuple(tumor,normal,"varscan",normvcf)} | annotvep_tn_varscan

    //VarScan_tonly
    varscan_tonly_comb=varscan_in.map{tumor,bam,bai,normal,nbam,nbai,bed,tpile,npile,tumorc,normalc ->
    tuple(tumor,bam,bai,bed,tpile,tumorc)} | varscan_tonly 
    varscan_tonly_comb1=varscan_tonly_comb.groupTuple().map{tumor,vcf-> tuple(tumor,vcf,"varscan_tonly")} | combineVariants_varscan_tonly
    
    varscan_tonly_comb1.join(sample_sheet)
    .map{tumor,marked,normvcf,normal ->tuple(tumor,"varscan_tonly",normvcf)} | annotvep_tonly_varscan

    
    //Combine All MAFs
    annotvep_tn_mut2.out.concat(annotvep_tn_strelka.out).concat(annotvep_tn_vardict.out).concat(annotvep_tn_varscan.out) | combinemafs_tn

    annotvep_tonly_mut2.out.concat(annotvep_tonly_vardict.out).concat(annotvep_tonly_varscan.out) | combinemafs_tonly

    //Implement Copy Number
    //Implement PCGR Annotator/CivIC Next
}


workflow SV_PIPE {
    take:
        bamwithsample
        
    main: 
        //Svaba
        svaba_somatic(bamwithsample)    

        //Manta

}


workflow QC_PIPE {
    take:
        fastqin
        fastpout
        applybqsr
        glnexusout //GLnexus germline output
        bcfout //DV germline output

    main:
    //QC Steps
    fc_lane(fastqin)
    fastq_screen(fastpout)
    kraken(fastqin)
    qualimap_bamqc(applybqsr)
    samtools_flagstats(applybqsr)
    fastqc(applybqsr)
    //Cohort VCF
    glout=glnexusout.map{germlinev,germlinenorm,tbi->tuple(germlinenorm,tbi)}
    vcftools(glout)
    collectvariantcallmetrics(glout)
    //Per sample VCFs
    bcfin=bcfout.map{samplename,vcf,vcf_tbi,gvcf,gvcf_tbi -> tuple(samplename,gvcf,gvcf_tbi)}
    bcftools_stats(bcfin)
    gatk_varianteval(bcfin)
    snpeff(bcfin)
    //Somalier
    somalier_extract(applybqsr) 
    som_in=somalier_extract.out.collect()
    somalier_analysis(som_in)

    //Prep for MultiQC input
    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect() 

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    fastqc_out=fastqc.out.map{samplename,html,zip->tuple(html,zip)}.collect()
    samtools_flagstats_out=samtools_flagstats.out.collect()
    bcftools_stats_out= bcftools_stats.out.collect()
    gatk_varianteval_out= gatk_varianteval.out.collect()
    snpeff_out=snpeff.out.collect()//map{vcf,csv,html->vcf,csv,html}.collect()
    vcftools_out=vcftools.out
    collectvariantcallmetrics_out=collectvariantcallmetrics.out//.map{details,summary->details,summary}
    somalier_analysis_out=somalier_analysis.out.collect()

    conall=fclane_out.concat(fqs_out,kraken_out,qualimap_out,samtools_flagstats_out,bcftools_stats_out,
    gatk_varianteval_out,snpeff_out,vcftools_out,collectvariantcallmetrics_out,somalier_analysis_out).flatten().toList()
    multiqc(conall)
}


//Variant Calling from BAM only
workflow INPUT_BAMVC_PIPE {
    
   if(params.sample_sheet){
        sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t", strip:true)
                       .map { row -> tuple(
                        row.Tumor,
                        row.Normal
                       )
                                  }
    } 
    
    //Either BAM Input or File sheet input 
    if(params.bam_input){
        baminputonly=Channel.fromPath(params.bam_input)
           .map{it-> tuple(it.simpleName,it,file("${it}.bai"))}
    }else if(params.file_input) {
        baminputonly=Channel.fromPath(params.file_input)
                        .splitCsv(header: false, sep: "\t", strip:true)
                        .map{ sample,bam -> 
                        tuple(sample, file(bam),file("${bam}.bai"))
                                  }
    }

    
    splitinterval(intervalbedin)
    
    bamwithsample=baminputonly.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(baminputonly,by:0).map{it.swap(3,0)}
    
    emit:
        bamwithsample
        splitout=splitinterval.out
        sample_sheet

}

