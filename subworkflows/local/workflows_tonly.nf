include {splitinterval; matchbed; matchbed as matchbed_ascat; 
        matchbed as matchbed_cnvkit} from '../../modules/local/splitbed.nf'

include {fc_lane} from '../../modules/local/fc_lane.nf'
include {fastq_screen} from '../../modules/local/fastq_screen.nf'
include {kraken} from '../../modules/local/kraken.nf'
include {qualimap_bamqc} from '../../modules/local/qualimap.nf'
include {fastqc} from '../../modules/local/fastqc.nf'
include {samtools_flagstats} from '../../modules/local/samtools_flagstats.nf'
include {vcftools} from '../../modules/local/vcftools.nf'
include {bcftools_stats} from '../../modules/local/bcftools_stats.nf'
include {gatk_varianteval; collectvariantcallmetrics} from '../../modules/local/gatk_varianteval.nf'
include {snpeff} from '../../modules/local/snpeff.nf'
include {somalier_extract;somalier_analysis_human;somalier_analysis_mouse} from '../../modules/local/somalier.nf'
include {mosdepth;mosdepth_exome} from '../../modules/local/mosdepth.nf'
include {multiqc} from  '../../modules/local/multiqc.nf'

include {fastp; bwamem2; indelrealign; bqsr_ir; 
    bqsr; gatherbqsr; applybqsr; samtoolsindex} from '../../modules/local/trim_align.nf'

include {pileup_paired as pileup_paired_t; pileup_paired as pileup_paired_n; 
    pileup_paired_tonly;
    learnreadorientationmodel;
    mutect2; mutect2filter;  contamination_paired; mergemut2stats; 
    mutect2_t_tonly; mutect2filter_tonly;
    contamination_tumoronly;
    learnreadorientationmodel_tonly;
    mergemut2stats_tonly} from '../../modules/local/mutect2.nf'
include {sage_tn;  sage_tonly} from '../../modules/local/sage.nf'
include {vardict_tn; vardict_tonly} from '../../modules/local/vardict.nf'
include {varscan_tn;  varscan_tonly} from '../../modules/local/varscan.nf'
include {octopus_tn; bcftools_index_octopus; 
    bcftools_index_octopus as bcftools_index_octopus_tonly; octopus_convertvcf;
    octopus_tonly; octopus_convertvcf_tonly} from '../../modules/local/octopus.nf'
include {deepsomatic_tonly_step1; deepsomatic_tonly_step2;
        deepsomatic_step3 as deepsomatic_tonly_step3  } from "../../modules/local/deepsomatic.nf"


include {combineVariants as combineVariants_vardict; combineVariants as combineVariants_vardict_tonly;
    combineVariants as combineVariants_varscan; combineVariants as combineVariants_varscan_tonly;
    combineVariants_alternative;
    combineVariants_alternative as combineVariants_deepsomatic; combineVariants_alternative as combineVariants_deepsomatic_tonly;
    combineVariants as combineVariants_sage; combineVariants as combineVariants_sage_tonly;
    combineVariants_alternative as combineVariants_lofreq; combineVariants as combineVariants_muse;
    combineVariants_alternative as combineVariants_octopus; 
    combineVariants_alternative as combineVariants_octopus_tonly;
    combinemafs_tn; somaticcombine;somaticcombine as somaticcombine_ffpe;
    combinemafs_tonly;somaticcombine_tonly;somaticcombine_tonly as somaticcombine_tonly_ffpe} from '../../modules/local/combinefilter.nf'

include {annotvep_tn as annotvep_tn_mut2; annotvep_tn as annotvep_tn_strelka;
    annotvep_tn as annotvep_tn_varscan; annotvep_tn as annotvep_tn_vardict; annotvep_tn as annotvep_tn_octopus;
    annotvep_tn as annotvep_tn_lofreq; annotvep_tn as annotvep_tn_muse; annotvep_tn as annotvep_tn_sage;
    annotvep_tn as annotvep_tn_deepsomatic;
    annotvep_tn as annotvep_tn_combined; 
    annotvep_tonly as annotvep_tonly_varscan; annotvep_tonly as annotvep_tonly_vardict;
    annotvep_tonly as annotvep_tonly_mut2; annotvep_tonly as annotvep_tonly_octopus; 
    annotvep_tonly as annotvep_tonly_sage; annotvep_tonly as annotvep_tonly_deepsomatic;
    annotvep_tonly as annotvep_tonly_combined} from '../../modules/local/annotvep.nf'

include {sobdetect_pass1 as sobdetect_pass1_mutect2_tonly; sobdetect_pass2 as sobdetect_pass2_mutect2_tonly; 
    sobdetect_metrics as sobdetect_metrics_mutect2_tonly; sobdetect_cohort_params as sobdetect_cohort_params_mutect2_tonly;
    sobdetect_pass1 as sobdetect_pass1_octopus_tonly; sobdetect_pass2 as sobdetect_pass2_octopus_tonly; 
    sobdetect_metrics as sobdetect_metrics_octopus_tonly; sobdetect_cohort_params as sobdetect_cohort_params_octopus_tonly;
    sobdetect_pass1 as sobdetect_pass1_vardict_tonly; sobdetect_pass2 as sobdetect_pass2_vardict_tonly; 
    sobdetect_metrics as sobdetect_metrics_vardict_tonly; sobdetect_cohort_params as sobdetect_cohort_params_vardict_tonly;
    sobdetect_pass1 as sobdetect_pass1_varscan_tonly; sobdetect_pass2 as sobdetect_pass2_varscan_tonly; 
    sobdetect_metrics as sobdetect_metrics_varscan_tonly; sobdetect_cohort_params as sobdetect_cohort_params_varscan_tonly
    } from "../../modules/local/ffpe.nf"

include {annotvep_tonly as annotvep_tonly_varscan_ffpe; annotvep_tonly as annotvep_tonly_vardict_ffpe;
    annotvep_tonly as annotvep_tonly_mut2_ffpe; annotvep_tonly as annotvep_tonly_octopus_ffpe; 
    annotvep_tonly as annotvep_tonly_sage_ffpe; annotvep_tonly as annotvep_tonly_deepsomatic_ffpe;
    annotvep_tonly as annotvep_tonly_combined_ffpe} from '../../modules/local/annotvep.nf'

include {svaba_tonly} from '../../modules/local/svaba.nf'
include {manta_tonly} from '../../modules/local/manta.nf'
include {gridss_tonly} from '../../modules/local/gridss.nf'
include {survivor_sv; 
    gunzip as gunzip_manta; gunzip as gunzip_gridss; 
    annotsv_tonly as annotsv_survivor_tonly;
    annotsv_tonly as annotsv_svaba_tonly; 
    annotsv_tonly as annotsv_gridss_tonly; 
    annotsv_tonly as annotsv_manta_tonly} from '../../modules/local/annotsv.nf'

include {freec} from '../../modules/local/freec.nf'
include {amber_tonly; cobalt_tonly; purple_tonly_novc; purple_tonly} from '../../modules/local/purple.nf'
include {cnvkit_exome_tonly; cnvkit_tonly } from '../../modules/local/cnvkit.nf'


//Workflows
workflow INPUT_TONLY {
    if(params.fastq_input){
        fastqinput=Channel.fromFilePairs(params.fastq_input) 
    }else if(params.fastq_file_input){
        fastqinput=Channel.fromPath(params.fastq_file_input)
                        .splitCsv(header: false, sep: "\t", strip:true)
                        .map{ sample,fq1,fq2 -> 
                            tuple(sample, tuple(file(fq1),file(fq2))) 
                            }  
    }

    if(params.sample_sheet){
        sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { "Sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor
                        )}  | view
    }else{
        sample_sheet=fastqinput.map{samplename,f1 -> tuple (
             samplename)} | view
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

    if (params.intervals){
        intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
    }else{
        intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
        }
    matchbed(intervalbedin) | splitinterval
    
    //Align
    fastp(fastqinput)
    bwamem2(fastp.out)

    if (params.indelrealign){ 
        bwaindelre = bwamem2.out | indelrealign 
        bqsrbambyinterval=bwaindelre.combine(splitinterval.out.flatten())
        bambyinterval=bwaindelre.combine(splitinterval.out.flatten())

        bqsr_ir(bqsrbambyinterval)
        
        bqsrs = bqsr_ir.out             
            | groupTuple 
            | map { samplename,beds -> 
                tuple( samplename, beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )} 
        gatherbqsr(bqsrs)

        tobqsr=bwaindelre.combine(gatherbqsr.out,by:0)
        applybqsr(tobqsr)

        bamwithsample=applybqsr.out.join(sample_sheet)
            | map{samplename,tumor,tumorbai -> tuple( samplename,tumor,tumorbai)}
        bambyinterval=bamwithsample.combine(splitinterval.out.flatten())

    }else{  
        bqsrbambyinterval=bwamem2.out.combine(splitinterval.out.flatten())

        bqsr(bqsrbambyinterval)
        bqsrs=bqsr.out | groupTuple
            | map { samplename,beds -> 
                tuple( samplename,
                beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )}
        gatherbqsr(bqsrs)

        tobqsr=bwamem2.out.combine(gatherbqsr.out,by:0)
        applybqsr(tobqsr)

        bamwithsample=applybqsr.out.join(sample_sheet)
            | map{samplename,tumor,tumorbai -> tuple( samplename,tumor,tumorbai)}
        bambyinterval=bamwithsample.combine(splitinterval.out.flatten())
        }


     emit:
        bamwithsample
        bambyinterval
        fastpout=fastp.out
        fastqin=fastqinput
        splitout=splitinterval.out
        bqsrbambyinterval
        sample_sheet
        bqsrout=applybqsr.out
        fullinterval=matchbed.out

}

workflow VC_TONLY {
    take:
    //Input is the BAMby interval
        bamwithsample
        splitout
        sample_sheet

    main:
    bambyinterval=bamwithsample.combine(splitout.flatten())

    //Common steps
    //Ensure that Tumor Only callers are included
    call_list = params.callers.split(',') as List
    call_list_tonly = params.tonlycallers.split(',') as List
    call_list = call_list.intersect(call_list_tonly)

    if (params.exome && "muse" in call_list){
        call_list.removeIf { it == 'muse' }
    }

    vc_tonly=Channel.empty()

    if ("mutect2" in call_list | "varscan" in call_list){ 
        pileup_paired_tonly(bambyinterval)
        pileup_paired_tout=pileup_paired_tonly.out.groupTuple()
        .map{samplename,pileups-> tuple( samplename,
        pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tpileup.table/)[0][1].toInteger() } ,
            )}
        contamination_tumoronly(pileup_paired_tout)
    }

    //Mutect2
    if ("mutect2" in call_list){
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

        mut2tonly_filter=mut2tonlyout.allmut2tonly
            | join(mergemut2stats_tonly.out)
            | join(learnreadorientationmodel_tonly.out)
            | join(contamination_tumoronly.out) 

        mutect2_in_tonly=mutect2filter_tonly(mut2tonly_filter) 
            | join(sample_sheet)
            | map{tumor,markedvcf,markedindex,finalvcf,finalindex,stats -> tuple(tumor,"mutect2_tonly",finalvcf,finalindex)} 
        annotvep_tonly_mut2(mutect2_in_tonly)
        
        vc_tonly=vc_tonly|concat(mutect2_in_tonly) 
    }

    //VarDict
    if ("vardict" in call_list){
        vardict_in_tonly=vardict_tonly(bambyinterval) | groupTuple()
            | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.vardict.vcf/)[0][1].toInteger()},"vardict_tonly")}
            | combineVariants_vardict_tonly
            | join(sample_sheet)
            | map{tumor,marked,markedindex,normvcf,normindex ->tuple(tumor,"vardict_tonly",normvcf,normindex)}
        annotvep_tonly_vardict(vardict_in_tonly)
        vc_tonly=vc_tonly.concat(vardict_in_tonly)

    }
    
    //VarScan_tonly
    if ("varscan" in call_list){
        varscan_in_tonly=bambyinterval.combine(contamination_tumoronly.out,by: 0)
            | varscan_tonly | groupTuple() 
            | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.varscan.vcf/)[0][1].toInteger()},"varscan_tonly")}
            | combineVariants_varscan_tonly 
            | join(sample_sheet)
            | map{tumor,marked,markedindex,normvcf,normindex ->tuple(tumor,"varscan_tonly",normvcf,normindex)} 
        annotvep_tonly_varscan(varscan_in_tonly)
        vc_tonly=vc_tonly|concat(varscan_in_tonly)
    }
    
    //Octopus_tonly
    if ("octopus" in call_list){
    octopus_in_tonly=bambyinterval | octopus_tonly | bcftools_index_octopus
        | groupTuple()
        | map{tumor,vcf,vcfindex -> tuple(tumor,vcf.toSorted{it -> it.name},vcfindex, "octopus_tonly")} 
        | combineVariants_alternative | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex ->tuple(tumor,"octopus_tonly",normvcf,normindex)} 
    annotvep_tonly_octopus(octopus_in_tonly)
    octopus_in_tonly_sc=octopus_in_tonly | octopus_convertvcf_tonly 
        | map{tumor,normvcf,normindex ->tuple(tumor,"octopus_tonly",normvcf,normindex)} 
    vc_tonly=vc_tonly|concat(octopus_in_tonly_sc)
    }

    //DeepSomatic Tonly
    if ("deepsomatic" in call_list){
    deepsomatic_tonly_in=deepsomatic_tonly_step1(bambyinterval) 
        | deepsomatic_tonly_step2  
        | deepsomatic_tonly_step3 | groupTuple 
        | map{samplename,vcf,vcf_tbi -> 
            tuple(samplename,vcf.toSorted{it -> (it.name =~ /${samplename}_(.*?).bed.vcf.gz/)[0][1].toInteger()},vcf_tbi,"deepsomatic_tonly")
            } 
        | combineVariants_deepsomatic_tonly           
        | join(sample_sheet) 
        | map{tumor,marked,markedindex,normvcf,normindex->tuple(tumor,"deepsomatic_tonly",normvcf,normindex)}

    annotvep_tonly_deepsomatic(deepsomatic_tonly_in)

    vc_tonly=vc_tonly | concat(deepsomatic_tonly_in) 
          
    }
    

    /*
    //SAGE
    if ("sage" in call_list){
    sage_in_tonly=sage_tonly(bamwithsample)
            | groupTuple() 
            | map{samplename,vcf,vcfindex -> tuple(samplename,vcf,"sage_tonly")} 
            | combineVariants_sage_tonly
            | join(sample_sheet) 
            | map{tumor,marked,markedindex,normvcf,normindex ->tuple(tumor,"sage_tonly",normvcf,normindex)}
    annotvep_tonly_sage(sage_in_tonly)

    vc_tonly=vc_tonly | concat(sage_in_tonly) 
    }
    */

    //FFPE Steps 
    if(params.ffpe){
        vc_ffpe_tonly=Channel.empty()
        bamwithsample1=bamwithsample

        if('mutect2' in call_list){
            mutect2_tonly_p1=bamwithsample1 | join(mutect2_in_tonly) 
                | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                | sobdetect_pass1_mutect2_tonly
            mutect2_tonly_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_mutect2_tonly

            mutect2_tonly_p2 = bamwithsample1 
                | join(mutect2_in_tonly)
                | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                | combine(sobdetect_cohort_params_mutect2_tonly.out) 
                | sobdetect_pass2_mutect2_tonly
            mutect2_tonly_p1_vcfs=mutect2_tonly_p1 | map{sample,vcf,info->vcf} |collect 
            mutect2_tonly_p2_vcfs=mutect2_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
            sobdetect_metrics_mutect2_tonly(mutect2_tonly_p1_vcfs,mutect2_tonly_p2_vcfs)

            mutect2_tonly_ffpe_out=mutect2_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample,"mutect2_tonly",filtvcf,vcftbi)} 
            annotvep_tonly_mut2_ffpe(mutect2_tonly_ffpe_out)
            vc_ffpe_tonly=vc_ffpe_tonly |concat(mutect2_tonly_ffpe_out)
        }

        if('octopus' in call_list){
            octopus_tonly_p1=bamwithsample1 | join(octopus_in_tonly) 
                | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                | sobdetect_pass1_octopus_tonly
            octopus_tonly_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_octopus_tonly

            octopus_tonly_p2 = bamwithsample1 
                | join(octopus_in_tonly)
                | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                | combine(sobdetect_cohort_params_octopus_tonly.out) 
                | sobdetect_pass2_octopus_tonly
            octopus_tonly_p1_vcfs=octopus_tonly_p1 | map{sample,vcf,info->vcf} |collect 
            octopus_tonly_p2_vcfs=octopus_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
            sobdetect_metrics_octopus_tonly(octopus_tonly_p1_vcfs,octopus_tonly_p2_vcfs)

            octopus_tonly_ffpe_out=octopus_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample,"octopus_tonly",filtvcf,vcftbi)} 
            annotvep_tonly_octopus_ffpe(octopus_tonly_ffpe_out)
            vc_ffpe_tonly=vc_ffpe_tonly |concat(octopus_tonly_ffpe_out)
        }

    if('vardict' in call_list){
        vardict_tonly_p1=bamwithsample1 | join(vardict_in_tonly) 
            | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
            | sobdetect_pass1_vardict_tonly
        vardict_tonly_p1 | map{sample,vcf,info->info}
            | collect
            | sobdetect_cohort_params_vardict_tonly

        vardict_tonly_p2 = bamwithsample1 
            | join(vardict_in_tonly)
            | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
            | combine(sobdetect_cohort_params_vardict_tonly.out) 
            | sobdetect_pass2_vardict_tonly
        vardict_tonly_p1_vcfs=vardict_tonly_p1 | map{sample,vcf,info->vcf} |collect 
        vardict_tonly_p2_vcfs=vardict_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
        sobdetect_metrics_vardict_tonly(vardict_tonly_p1_vcfs,vardict_tonly_p2_vcfs)

        vardict_tonly_ffpe_out=vardict_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample,"vardict_tonly",filtvcf,vcftbi)} 
        annotvep_tonly_vardict_ffpe(vardict_tonly_ffpe_out)
        vc_ffpe_tonly=vc_ffpe_tonly |concat(vardict_tonly_ffpe_out)
    }

    if('varscan' in call_list){
    varscan_tonly_p1=bamwithsample1 | join(varscan_in_tonly) 
        | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
        | sobdetect_pass1_varscan_tonly
    varscan_tonly_p1 | map{sample,vcf,info->info}
        | collect
        | sobdetect_cohort_params_varscan_tonly

    varscan_tonly_p2 = bamwithsample1 
        | join(varscan_in_tonly)
        | map{tumor,tbam,tbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
        | combine(sobdetect_cohort_params_varscan_tonly.out) 
        | sobdetect_pass2_varscan_tonly
    varscan_tonly_p1_vcfs=varscan_tonly_p1 | map{sample,vcf,info->vcf} |collect 
    varscan_tonly_p2_vcfs=varscan_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
    sobdetect_metrics_varscan_tonly(varscan_tonly_p1_vcfs,varscan_tonly_p2_vcfs)

    varscan_tonly_ffpe_out=varscan_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample,"varscan_tonly",filtvcf,vcftbi)} 
    annotvep_tonly_varscan_ffpe(varscan_tonly_ffpe_out)
    vc_ffpe_tonly=vc_ffpe_tonly |concat(varscan_tonly_ffpe_out)
    }
}


    //Combined Variants and Annotated
    //Emit for SC downstream, take Oc/Mu2/sage/Vard/Varscan

    if (call_list.size()>1){
        vc_tonly
        | groupTuple() 
        | somaticcombine_tonly 
        | map{tumor,vcf,index ->tuple(tumor,"combined_tonly",vcf,index)} 
        | annotvep_tonly_combined
    }

    if("sage" in call_list){
        somaticcall_input=sage_in_tonly
    }else if("mutect2" in call_list){
        somaticcall_input=mutect2_in_tonly
    }else if("mutect2" in call_list & params.ffpe){
        somaticcall_input=mutect2_ffpe_out
    }else{
        somaticcall_input=Channel.empty()
    }

    emit:
        somaticcall_input
}


workflow SV_TONLY {
    take:
        bamwithsample
        
    main: 
        svcall_list = params.svcallers.split(',') as List
        svout=Channel.empty()

        //Svaba
        if ("svaba" in svcall_list){

            svaba_out=svaba_tonly(bamwithsample)
            .map{ tumor,bps,contigs,discord,alignments,so_indel,so_sv,unfil_so_indel,unfil_sv,log ->
                tuple(tumor,so_sv,"svaba_tonly")} 
            annotsv_svaba_tonly(svaba_out).ifEmpty("Empty SV input--No SV annotated")
            svout=svout | concat(svaba_out)
        }
        //Manta
        if ("manta" in svcall_list){
            manta_out=manta_tonly(bamwithsample)
                .map{tumor, sv, indel, tumorsv -> 
                tuple(tumor,tumorsv,"manta_tonly")} 
            annotsv_manta_tonly(manta_out).ifEmpty("Empty SV input--No SV annotated")
            svout=svout | concat(manta_out)
        }

        //GRIDSS
        if ("gridss" in svcall_list){
        gridss_out=gridss_tonly(bamwithsample)
        gridss_out_forsv=gridss_out
            | map{tumor,vcf,index,bam,gripssvcf,gripsstbi,gripssfilt,filttbi ->
            tuple(tumor,gripssfilt,"gridss_tonly")} | gunzip_gridss
        annotsv_gridss_tonly(gridss_out_forsv).ifEmpty("Empty SV input--No SV annotated")
            svout=svout | concat(gridss_out)
        }

        //Survivor
        if (svcall_list.size()>1){
            //Survivor
            svout | groupTuple
                | survivor_sv 
                | annotsv_survivor_tonly 
                | ifEmpty("Empty SV input--No SV annotated")
         }

        if("gridss" in svcall_list){
            somaticsv_input=gridss_out 
                | map{tumor,vcf,index,bam,gripssvcf,gripsstbi,gripssfilt,filttbi ->
                    tuple(tumor,vcf,index,gripsstbi,gripssfilt,filttbi)}
        }else if("manta" in svcall_list){
            somaticsv_input=manta_out 
                | map{tumor,gsv,gsv_tbi,so_sv,so_sv_tbi,unfil_sv,unfil_sv_tbi,unfil_indel,unfil_indel_tbi ->
                    tuple(tumor,unfil_sv,unfil_sv_tbi,so_sv,so_sv_tbi)}
        }else{
            somaticsv_input=Channel.empty()
        }
    
    emit:
        somaticsv_input

}



workflow CNVmouse_tonly {
    take:
        bamwithsample
        
    main:  
    cnvcall_list = params.cnvcallers.split(',') as List

    if ("freec" in cnvcall_list){
        freec(bamwithsample)
    }
     //CNVKIT
        if ("cnvkit" in cnvcall_list){
            if(params.exome){
                matchbed_cnvkit(intervalbedin)
                bamwithsample | combine(matchbed_cnvkit.out) | cnvkit_exome_tonly
            }else{
                bamwithsample | cnvkit_tonly
            }
        }
}


workflow CNVhuman_tonly {
    take:
        bamwithsample
        somaticcall_input

    main: 
        cnvcall_list = params.cnvcallers.split(',') as List

        if ("freec" in cnvcall_list){
            //FREEC-Unpaired only
            bamwithsample | freec 
        }

        if (params.exome && "purple" in cnvcall_list ){
            cnvcall_list.removeIf { it == 'purple' }
        }

        if ("purple" in cnvcall_list){
            //Purple
            bamwithsample | amber_tonly
            bamwithsample | cobalt_tonly
            purplein=amber_tonly.out.join(cobalt_tonly.out)
            purplein.join(somaticcall_input)| 
            map{t1,amber,cobalt,vc,vcf,index -> tuple(t1,amber,cobalt,vcf,index)}  
                | purple_tonly
        }

        //CNVKIT        
        if ("cnvkit" in cnvcall_list){
            if(params.exome){
                matchbed_cnvkit(intervalbedin)
                bamwithsample | combine(matchbed_cnvkit.out) | cnvkit_exome_tonly
            }else{
                bamwithsample | cnvkit_tonly
            }
        }
        
}

workflow CNVhuman_novc_tonly {
    take:
        bamwithsample

    main: 
        cnvcall_list = params.cnvcallers.split(',') as List

        if ("freec" in cnvcall_list){
            //FREEC-Unpaired only
            bamwithsample | freec 
        }   
        
        if (params.exome && "purple" in cnvcall_list){
            cnvcall_list.removeIf { it == 'purple' }
        }

        if ("purple" in cnvcall_list){
            //Purple
            bamwithsample | amber_tonly
            bamwithsample | cobalt_tonly
            purplein=amber_tonly.out.join(cobalt_tonly.out)
            map{t1,amber,cobalt -> tuple(t1,amber,cobalt)}  
                | purple_tonly_novc
        }

        if ("cnvkit" in cnvcall_list){
            if(params.exome){
                matchbed_cnvkit(intervalbedin)
                bamwithsample | combine(matchbed_cnvkit.out) | cnvkit_exome_tonly
            }else{
                bamwithsample | cnvkit_tonly
            }
        }
}


workflow QC_TONLY {
    take:
        fastqin
        fastpout
        bqsrout
        fullinterval
    
    main:
    //QC Steps For Tumor-Only-No Germline Variant QC
    fc_lane(fastqin)
    fastq_screen(fastpout)
    kraken(fastqin)

    //BQSR BAMs 
    fastqc(bqsrout)
    samtools_flagstats(bqsrout)
    qualimap_bamqc(bqsrout)
    if(params.exome){
        mosdepth_out=bqsrout | combine(fullinterval) | mosdepth_exome | collect
    }else{
        mosdepth_out=bqsrout | combine(fullinterval) | mosdepth | collect
    }


    somalier_extract(bqsrout) 
    som_in=somalier_extract.out.collect()
    if(params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){ 
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome.matches("mm10")){ 
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
        samtools_flagstats_out,mosdepth_out, 
        somalier_analysis_out).flatten().toList()
    
    multiqc(conall)
}


//QC Tumor Only-BAMs 
workflow QC_TONLY_BAM {
    take:
        bams
        fullinterval
        
    main:
    //BQSR BAMs 
    fastqc(bams)
    samtools_flagstats(bams)
    qualimap_bamqc(bams)
    if(params.exome){
        mosdepth_out=bams | combine(fullinterval) | mosdepth_exome | collect
    }else{
        mosdepth_out=bams | combine(fullinterval) | mosdepth | collect
    }


    somalier_extract(bams) 
    som_in=somalier_extract.out.collect()
    if(params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome.matches("mm10")){
        somalier_analysis_mouse(som_in)
        somalier_analysis_out=somalier_analysis_mouse.out.collect()
    }
    
    //Prep for MultiQC input
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    samtools_flagstats_out=samtools_flagstats.out.collect()

    conall=qualimap_out | concat(
        samtools_flagstats_out,mosdepth_out, 
        somalier_analysis_out) 
        | flatten | toList
    
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
        if (params.intervals){
            intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
        }else{
            intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
        }
    matchbed(intervalbedin) | splitinterval
    
    bamwithsample=baminputonly

    emit:
        bamwithsample
        splitout=splitinterval.out
        sample_sheet
        fullinterval=matchbed.out
        
    
}

