//All Worksflows in One Place
// TODO split subworkflows out into one per file  
// TODO: this line should be moved to within a subworkflow or the main workflow
intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')



include {fc_lane; fastq_screen;kraken;qualimap_bamqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;fastqc;
    somalier_extract;somalier_analysis_human;somalier_analysis_mouse;
    multiqc} from  '../../modules/local/qc.nf'

include {deepvariant_step1;deepvariant_step2;deepvariant_step3;
    deepvariant_combined;glnexus} from '../../modules/local/germline.nf'

include {fastp; bwamem2; 
    bqsr; gatherbqsr; applybqsr; samtoolsindex} from '../../modules/local/trim_align.nf'
    
include {mutect2; mutect2filter; pileup_paired_t; pileup_paired_n; 
    bcftools_index_octopus;
    contamination_paired; learnreadorientationmodel; mergemut2stats;
    combineVariants as combineVariants_vardict; combineVariants as combineVariants_varscan; 
    combineVariants as combineVariants_vardict_tonly; combineVariants as combineVariants_varscan_tonly;
    combineVariants_alternative ; 
    annotvep_tn as annotvep_tn_mut2; annotvep_tn as annotvep_tn_strelka; annotvep_tn as annotvep_tn_varscan; annotvep_tn as annotvep_tn_vardict;
    combinemafs_tn} from '../../modules/local/variant_calling.nf'

include {mutect2_t_tonly; mutect2filter_tonly; pileup_paired_tonly; 
    varscan_tonly; vardict_tonly; 
    octopus_tonly; 
    contamination_tumoronly;
    learnreadorientationmodel_tonly; 
    mergemut2stats_tonly; octopus_convertvcf_tonly;
    annotvep_tonly as annotvep_tonly_varscan; annotvep_tonly as annotvep_tonly_vardict; 
    annotvep_tonly as annotvep_tonly_mut2; annotvep_tonly as annotvep_tonly_octopus;
    annotvep_tonly as annotvep_tonly_combined;
    combinemafs_tonly; somaticcombine_tonly} from '../../modules/local/variant_calling_tonly.nf'

include {manta_tonly; svaba_tonly; survivor_sv; gunzip;
annotsv_tonly as annotsv_manta_tonly; annotsv_tonly as annotsv_svaba_tonly;
annotsv_tonly as annotsv_survivor_tonly} from '../../modules/local/structural_variant.nf'

include {freec; amber_tonly; cobalt_tonly; purple  } from '../../modules/local/copynumber.nf'

include {splitinterval} from '../../modules/local/splitbed.nf'


workflow INPUT_TONLY {
    if(params.fastq_input){
        fastqinput=Channel.fromFilePairs(params.fastq_input)

    }else if(params.fastq_file_input) {
        fastqinput=Channel.fromPath(params.fastq_file_input)
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
    //indelrealign(bwamem2.out) Consider indelreaglinement using ABRA?

    bqsrbambyinterval=bwamem2.out.combine(splitinterval.out.flatten())

    
    bqsr(bqsrbambyinterval)
    bqsrs=bqsr.out.groupTuple()
        .map { samplename,beds -> tuple( samplename, 
        beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )
        }
    gatherbqsr(bqsrs)

    tobqsr=bwamem2.out.combine(gatherbqsr.out,by:0)
    applybqsr(tobqsr) 
    
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
    
    mutect2_t_tonly.out.groupTuple()
        | multiMap { tumor,vcfs,f1r2,stats -> 
        mut2tout_lor: tuple(tumor,
                f1r2.toSorted{ it -> (it.name =~ /${tumor}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } )
        mut2tonly_mstats:  tuple( tumor,
                stats.toSorted{ it -> (it.name =~ /${tumor}_(.*?).tonly.mut2.vcf.gz.stats/)[0][1].toInteger() })
        allmut2tonly: tuple(tumor,
                vcfs.toSorted{ it -> (it.name =~ /${tumor}_(.*?).tonly.mut2.vcf.gz/)[0][1].toInteger() } )
        } 
    | set{mut2tonlyout}

  

    learnreadorientationmodel_tonly(mut2tonlyout.mut2tout_lor)
    mergemut2stats_tonly(mut2tonlyout.mut2tonly_mstats)
    contamination_tumoronly(pileup_paired_tout)

    
    mut2tonly_filter=mut2tonlyout.allmut2tonly
        | join(mergemut2stats_tonly.out)
        | join(learnreadorientationmodel_tonly.out)
        | join(contamination_tumoronly.out) 

    mutect2_tonly_in=mutect2filter_tonly(mut2tonly_filter) 
        | join(sample_sheet)
        | map{tumor,markedvcf,markedindex,finalvcf,finalindex,stats -> tuple(tumor,"mutect2_tonly",finalvcf,finalindex)} 
    annotvep_tonly_mut2(mutect2_tonly_in)


    //VarDict
    vardict_in_tonly=vardict_tonly(bambyinterval) | groupTuple()
        | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.vardict.vcf/)[0][1].toInteger()},"vardict_tonly")}
        | combineVariants_vardict_tonly
        | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex ->tuple(tumor,"vardict_tonly",normvcf,normindex)}
    annotvep_tonly_vardict(vardict_in_tonly)

    //VarScan_tonly
    varscan_in_tonly=bambyinterval.combine(contamination_tumoronly.out,by: 0)
        | varscan_tonly | groupTuple() 
        | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.varscan.vcf/)[0][1].toInteger()},"varscan_tonly")}
        | combineVariants_varscan_tonly 
        | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex ->tuple(tumor,"varscan_tonly",normvcf,normindex)} 
    annotvep_tonly_varscan(varscan_in_tonly)

    //Octopus_tonly
    octopus_in_tonly=bambyinterval | octopus_tonly | bcftools_index_octopus
        | groupTuple()
        | map{tumor,vcf,vcfindex -> tuple(tumor,vcf.toSorted{it -> it.name}
                ,vcfindex, "octopus_tonly")} 
        | combineVariants_alternative | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex ->tuple(tumor,"octopus_tonly",normvcf,normindex)} 
    annotvep_tonly_octopus(octopus_in_tonly)
    octopus_in_tonly_sc=octopus_in_tonly | octopus_convertvcf_tonly 
        | map{tumor,normvcf,normindex ->tuple(tumor,"octopus_tonly",normvcf,normindex)} 

    //Combined Variants and Annotated
    mutect2_tonly_in | concat(octopus_in_tonly_sc)
        | concat(vardict_in_tonly) | concat(varscan_in_tonly)
        | groupTuple()
        | somaticcombine_tonly 
        | map{tumor,vcf,index ->tuple(tumor,"combined_tonly",vcf,index)} 
        | annotvep_tonly_combined



    emit:
        somaticcall_input=octopus_in_tonly


}


workflow SV_TONLY {
    take:
        bamwithsample
        
    main: 
        //Svaba
        svaba_out=svaba_tonly(bamwithsample)
        .map{ tumor,bps,contigs,discord,alignments,so_indel,so_sv,unfil_so_indel,unfil_sv,log ->
            tuple(tumor,so_sv,"svaba_tonly")} 
        annotsv_svaba_tonly(svaba_out).ifEmpty("Empty SV input--No SV annotated")

        //Manta
        manta_out=manta_tonly(bamwithsample)
            .map{tumor, sv, indel, tumorsv -> 
            tuple(tumor,tumorsv,"manta_tonly")} 
        annotsv_manta_tonly(manta_out).ifEmpty("Empty SV input--No SV annotated")

        //Delly-WIP

        //Survivor
        gunzip(manta_out).concat(svaba_out).groupTuple()
       | survivor_sv | annotsv_survivor_tonly | ifEmpty("Empty SV input--No SV annotated")
}



workflow CNVmouse_tonly {
    take:
        bamwithsample
        
    main: 
        freec(bamwithsample)
}


workflow CNVhuman_tonly {
    take:
        bamwithsample
        somaticcall_input

    main: 
        //FREEC-Unpaired onlypu
        bamwithsample | freec 
        
        //Purple
        bamwithsample | amber_tonly
        bamwithsample | cobalt_tonly
        purplein=amber_tonly.out.join(cobalt_tonly.out)
        purplein.join(somaticcall_input)| 
        map{t1,amber,cobalt,vc,vcf,index -> tuple(t1,amber,cobalt,vcf,index)}  
            | purple

        
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
    if(params.genome=="hg38"){ 
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome=="mm10"){ 
        somalier_analysis_mouse(som_in)
        somalier_analysis_out=somalier_analysis_mouse.out.collect()
    }
    
    //Prep for MultiQC input
    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect() 

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    fastqc_out=fastqc.out.map{samplename,html,zip->tuple(html,zip)}.collect()

    samtools_flagstats_out=samtools_flagstats.out.collect()



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
        bambai = params.bam_input + ".bai"
        baionly = bambai.replace(".bam", "")
        bamcheck1 = file(bambai)
        bamcheck2 = file(baionly)

        if (bamcheck1.size()>0){
            baminputonly=Channel.fromPath(params.bam_input)
                | map{it-> tuple(it.simpleName,it,file("${it}.bai"))} 
        }else if (bamcheck2.size()>0){
            bai=Channel.from(bamcheck2).map{it -> tuple(it.simpleName,it)}
            baminputonly=Channel.fromPath(params.bam_input)
            | map{it-> tuple(it.simpleName,it)}
            | join(bai)
        }else if (bamcheck1.size==0 && bamcheck2.size==0 ){
            println "Missing BAM Index"
        }

        sample_sheet=baminputonly.map{samplename,bam,bai -> tuple (
             samplename)}

    }else if(params.bam_file_input) {
        baminputonly=Channel.fromPath(params.bam_file_input)
                        .splitCsv(header: false, sep: "\t", strip:true)
                        .map{ sample,bam,bai  -> 
                        tuple(sample, file(bam),file(bai))
                                  }

        sample_sheet=baminputonly.map{samplename,bam,bai -> tuple (
             samplename)}
    }

    splitinterval(intervalbedin)
    
    bamwithsample=baminputonly

    emit:
        bamwithsample
        splitout=splitinterval.out
        sample_sheet
    
}

