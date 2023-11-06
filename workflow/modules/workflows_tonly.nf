//All Worksflows in One Place  
intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')


include {fc_lane; fastq_screen;kraken;qualimap_bamqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;fastqc;
    somalier_extract;somalier_analysis;multiqc} from  './qc.nf'

include {deepvariant_step1;deepvariant_step2;deepvariant_step3;
    deepvariant_combined;glnexus} from './germline.nf'

include {fastp; bwamem2; 
    bqsr; gatherbqsr; applybqsr; samtoolsindex} from './trim_align.nf'
    
include {mutect2; mutect2filter; pileup_paired_t; pileup_paired_n; 
    contamination_paired; learnreadorientationmodel;mergemut2stats;
    combineVariants as combineVariants_vardict; combineVariants as combineVariants_varscan; 
    combineVariants as combineVariants_vardict_tonly; combineVariants as combineVariants_varscan_tonly
    annotvep_tn as annotvep_tn_mut2; annotvep_tn as annotvep_tn_strelka; annotvep_tn as annotvep_tn_varscan; annotvep_tn as annotvep_tn_vardict;
    combinemafs_tn} from './variant_calling.nf'

include {mutect2_t_tonly; mutect2filter_tonly; pileup_paired_tonly; 
    varscan_tonly; vardict_tonly; 
    contamination_tumoronly;
    learnreadorientationmodel_tonly; 
    mergemut2stats_tonly;
    annotvep_tonly as annotvep_tonly_varscan; annotvep_tonly as annotvep_tonly_vardict; annotvep_tonly as annotvep_tonly_mut2;
    combinemafs_tonly} from './variant_calling_tonly.nf'

include {manta_tonly; svaba_tonly; annotsv_tn as annotsv_manta_tonly; annotsv_tn as annotsv_svaba_tonly} from './structural_variant.nf'

include {freec; purple } from './copynumber.nf'

include {splitinterval} from "./splitbed.nf"



workflow INPUT_TONLY {
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
        sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor
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

workflow ALIGN_TONLY {
    take:
        fastqinput
        sample_sheet

    main:
    fastp(fastqinput)
    splitinterval(intervalbedin)
    
    bwamem2(fastp.out)
    //indelrealign(bwamem2.out)

    bqsrbambyinterval=bwamem2.out.combine(splitinterval.out.flatten())

    
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
    //bamwithsample=samtoolsindex.out.join(sample_sheet).map{it.swap(3,0)}.join(samtoolsindex.out).map{it.swap(3,0)}
    bamwithsample=applybqsr.out.join(sample_sheet)
    .map{samplename,tumor,tumorbai -> tuple( samplename,tumor,tumorbai)
        }
    bambyinterval=bamwithsample.combine(splitinterval.out.flatten())

     emit:
        bamwithsample
        bambyinterval
        fastpout=fastp.out
        fastqin=fastqinput
        splitout=splitinterval.out
        bqsrbambyinterval
        sample_sheet
        bqsrout=applybqsr.out

}

workflow VC_TONLY {
    take:
    //Input is the BAMby interval
        bamwithsample
        splitout
        sample_sheet

    main:
    
    bambyinterval=bamwithsample.combine(splitout.flatten())
    pileup_paired_tonly(bambyinterval)
    pileup_paired_tout=pileup_paired_tonly.out.groupTuple()
    .map{samplename,pileups-> tuple( samplename,
    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tumor.pileup.table/)[0][1].toInteger() } ,
    )}

    mutect2_t_tonly(bambyinterval)    
    
    
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
    
    
    //##VCF2MAF TO

    mutect2filter_tonly.out
    .join(sample_sheet)
    .map{tumor,markedvcf,finalvcf,stats -> tuple(tumor,"mutect2",finalvcf)} | annotvep_tonly_mut2

    //VarDict_tonly
    vardict_tonly(bambyinterval).groupTuple().map{tumor,vcf-> tuple(tumor,vcf,"vardict_tonly")} | combineVariants_vardict_tonly
    combineVariants_vardict_tonly.out.join(sample_sheet)
    .map{tumor,marked,normvcf ->tuple(tumor,"vardict_tonly",normvcf)} | annotvep_tonly_vardict

    //VarScan_tonly
    varscan_in=bambyinterval.join(contamination_tumoronly.out)
    varscan_tonly_comb=varscan_tonly(varscan_in).groupTuple().map{tumor,vcf-> tuple(tumor,vcf,"varscan")} | combineVariants_varscan_tonly
    
    varscan_tonly_comb.join(sample_sheet)
    .map{tumor,marked,normvcf ->tuple(tumor,"varscan_tonly",normvcf)} | annotvep_tonly_varscan

    //Combine All Final
    //annotvep_tonly_mut2.out.concat(annotvep_tonly_vardict.out).concat(annotvep_tonly_varscan.out) | combinemafs_tonly

}


workflow SV_TONLY {
    take:
        bamwithsample
        
    main: 
        //Svaba
        svaba_out=svaba_tonly(bamwithsample)
        .map{ tumor,bps,contigs,discord,alignments,so_indel,so_sv,unfil_so_indel,unfil_sv,log ->
            tuple(tumor,so_sv,"svaba")} 
        annotsv_svaba_tonly(svaba_out).ifEmpty("Empty SV input--No SV annotated")

        //Manta
        manta_out=manta_tonly(bamwithsample)
            .map{tumor,gsv,so_sv,unfil_sv,unfil_indel,tumorSV -> 
            tuple(tumor,so_sv,"manta")} 
        annotsv_manta_tonly(manta_out).ifEmpty("Empty SV input--No SV annotated")

}

workflow CNV_TONLY {
    take:
        bamwithsample
        
    main: 
        //mm10 use sequenza only, hg38 use purple
        if(params.genome=="hg38"){
            purple(bamwithsample)
        }
        if(params.genome=="mm10"){
            freec(bamwithsample)
        } 
       
}

workflow QC_TONLY {
    take:
        fastqin
        fastpout
        bqsrout

    main:
    //QC Steps For Tumor-Only-No Germline Variant QC
    fc_lane(fastqin)
    fastq_screen(fastpout)
    kraken(fastqin)

    //BQSR BAMs 
    fastqc(bqsrout)
    samtools_flagstats(bqsrout)
    qualimap_bamqc(bqsrout)

    somalier_extract(bqsrout) 
    som_in=somalier_extract.out.collect()
    somalier_analysis(som_in)
    
    //Prep for MultiQC input
    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect() 

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    fastqc_out=fastqc.out.map{samplename,html,zip->tuple(html,zip)}.collect()

    samtools_flagstats_out=samtools_flagstats.out.collect()
    somalier_analysis_out=somalier_analysis.out.collect()

    conall=fclane_out.concat(fqs_out,kraken_out,qualimap_out,fastqc_out,
    samtools_flagstats_out,
    somalier_analysis_out).flatten().toList()
    
    multiqc(conall)
}






//Variant Calling from BAM only
workflow INPUT_TONLY_BAM {
    main:
  
    //Either BAM Input or File sheet input 
    if(params.bam_input){
        baminputonly=Channel.fromPath(params.bam_input)
           .map{it-> tuple(it.simpleName,it,file("${it}.bai"))}

        sample_sheet=baminputonly.map{samplename,bam,bai -> tuple (
             samplename)}.view()
    }else if(params.file_input) {
        baminputonly=Channel.fromPath(params.file_input)
                        .splitCsv(header: false, sep: "\t", strip:true)
                        .map{ sample,bam -> 
                        tuple(sample, file(bam),file("${bam}.bai"))
                                  }
        sample_sheet=baminputonly.map{samplename,bam,bai -> tuple (
             samplename)}.view()
        
    }

    splitinterval(intervalbedin)
    
    bamwithsample=baminputonly

    emit:
        bamwithsample
        splitout=splitinterval.out
        sample_sheet
    

}

