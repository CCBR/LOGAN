//All Worksflows in One Place         
intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')


include {fc_lane; fastq_screen;kraken;qualimap_bamqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;
    somalier_extract;somalier_analysis;multiqc} from  './qc.nf'
include {deepvariant_step1;deepvariant_step2;deepvariant_step3;
    deepvariant_combined;glnexus} from './germline.nf'
include {fastp; bwamem2; indelrealign; bqsr; 
    gatherbqsr; applybqsr; samtoolsindex} from './trim_align.nf'
include {mutect2; mutect2_t_tonly; mutect2filter; mutect2filter_tonly; 
    pileup_paired_t; pileup_paired_n; 
    contamination_paired; contamination_tumoronly;
    learnreadorientationmodel; learnreadorientationmodel_tonly; 
    mergemut2stats; mergemut2stats_tonly;
    annotvep_tn; annotvep_tonly} from './variant_calling.nf'
include {splitinterval} from "./splitbed.nf"



workflow INPUT_PIPE {
    
    if(params.fastq_input){
        fastqinput=Channel.fromFilePairs(params.fastq_input).view()
    }else if(params.file_input) {
        fastqinput=Channel.fromPath(params.file_input).view()
                        .splitCsv(header: false, sep: "\t", strip:true)
                        .map{ sample,fq1,fq2 -> 
                        tuple(sample, tuple(file(fq1),file(fq2)))
                                  }
                       .view()
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
             samplename)}.view()
        
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
    indelrealign(bwamem2.out)

//Indel Realignment after BQSR FOR ALL BAMS 
    indelbambyinterval=indelrealign.out.combine(splitinterval.out.flatten())
    bambyinterval=bwamem2.out.combine(splitinterval.out.flatten())
    
        
    bqsr(indelbambyinterval)
    bqsrs=bqsr.out.groupTuple()
        .map { samplename,beds -> tuple( samplename, 
        beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )
        }
    gatherbqsr(bqsrs)

    tobqsr=indelrealign.out.combine(gatherbqsr.out,by:0)
    applybqsr(tobqsr) 
    samtoolsindex(applybqsr.out)
    
    bamwithsample=samtoolsindex.out.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(samtoolsindex.out,by:0).map{it.swap(3,0)}

    emit:
        bamwithsample
        bambyinterval
        fastpout=fastp.out
        fastqin=fastqinput
        splitout=splitinterval.out
        indelbambyinterval
        sample_sheet
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
    
    glin=deepvariant_step3.out.map{samplename,vcf,vcf_tbi,gvcf,gvcf_tbi -> gvcf}.collect()
    glnexus(glin)

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


    //#To implement
    //CNMOPs from the BAM BQSRs
    //##VCF2MAF TO
    tn_vepin=mutect2filter.out
    .join(sample_sheet)

    annotvep_tn(tn_vepin)
    annotvep_tonly(mutect2filter_tonly.out)
    // Immplment PCGR Annotator/CivIC?
}



workflow QC_PIPE {
    
    main:
    //QC Steps
    fc_lane(fastqinput)
    fastq_screen(fastp.out)
    kraken(fastqinput)
    qualimap_bamqc(bwamem2.out)
    samtools_flagstats(bwamem2.out)
    glout=glnexus.out.map{germlinev,germlinenorm,tbi->tuple(germlinenorm,tbi)}
    vcftools(glout)
    collectvariantcallmetrics(glout)
    //bcfin=deepvariant_combined.out.map{samplename,vcf,vcf_tbi,gvcf,gvcf_tbi -> tuple(samplename,gvcf,gvcf_tbi)}
    bcfin=deepvariant_step3.out.map{samplename,vcf,vcf_tbi,gvcf,gvcf_tbi -> tuple(samplename,gvcf,gvcf_tbi)}
    bcftools_stats(bcfin)
    gatk_varianteval(bcfin)
    snpeff(bcfin)
    somalier_extract(bwamem2.out) 
    som_in=somalier_extract.out.collect()
    somalier_analysis(som_in)
    
    
    //FASTQC DOWN THE LINE ADD
    //MULTIQC
    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect() 
    //map{r1_ht,r1_png,r1_txt,r2_ht,r2_png,r2_txt-> r2_txt}.collect()

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
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
        sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true).view()
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
    
    bamwithsample=baminputonly.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(baminputonly,by:0).map{it.swap(3,0)}.view()
    
    emit:
        bamwithsample
        splitout=splitinterval.out
        sample_sheet

}
