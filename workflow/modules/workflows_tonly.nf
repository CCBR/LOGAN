//All Worksflows in One Place         
intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')


include {fc_lane; fastq_screen;kraken;qualimap_bamqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;fastqc;
    somalier_extract;somalier_analysis;multiqc} from  './qc.nf'
include {deepvariant_step1;deepvariant_step2;deepvariant_step3;
    deepvariant_combined;glnexus} from './germline.nf'
include {fastp; bwamem2; 
    bqsr; gatherbqsr; applybqsr; samtoolsindex} from './trim_align.nf'
include {mutect2; mutect2_t_tonly; mutect2filter; mutect2filter_tonly; 
    pileup_paired_t; pileup_paired_n; pileup_paired_tonly;
    contamination_paired; contamination_tumoronly;
    learnreadorientationmodel; learnreadorientationmodel_tonly; 
    mergemut2stats; mergemut2stats_tonly;
    annotvep_tn; annotvep_tonly} from './variant_calling.nf'
include {splitinterval} from "./splitbed.nf"



workflow INPUT_TONLY_PIPE {
    fastqinput=Channel.fromFilePairs(params.fastq_input)

 
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
             samplename)}.view()
        
    }
    
    emit:
        fastqinput
        sample_sheet

}

workflow TRIM_ALIGN_TONLY_PIPE {
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
        //indelbambyinterval
        bqsrbambyinterval
        sample_sheet
        bqsrout=applybqsr.out

}

workflow VARIANT_TONLY_PIPE {
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
    
    //#To implement
    //CNMOPs from the BAM BQSRs
    //##VCF2MAF TO

    annotvep_tonly(mutect2filter_tonly.out)

}


workflow QC_TONLY_PIPE {
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
workflow INPUT_TONLY_BAMVC_PIPE {
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
