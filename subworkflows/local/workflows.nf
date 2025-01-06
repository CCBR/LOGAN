include {splitinterval; matchbed; matchbed as matchbed_ascat; matchbed as matchbed_cnvkit} from '../../modules/local/splitbed.nf'

//QC
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
include {mosdepth} from '../../modules/local/mosdepth.nf'
include {multiqc} from  '../../modules/local/multiqc.nf'

include {fastp; bwamem2; indelrealign; bqsr_ir;
    bqsr; gatherbqsr; applybqsr; samtoolsindex} from  '../../modules/local/trim_align.nf'

include {deepvariant_step1; deepvariant_step2; deepvariant_step3;
    deepvariant_combined; glnexus;
    bcfconcat as bcfconcat_vcf; bcfconcat as bcfconcat_gvcf} from '../../modules/local/deepvariant.nf'

include {pileup_paired as pileup_paired_t; pileup_paired as pileup_paired_n; 
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
include {lofreq_tn} from '../../modules/local/lofreq.nf'
include {strelka_tn; combineVariants_strelka; convert_strelka} from '../../modules/local/strelka.nf'
include {muse_tn} from '../../modules/local/muse.nf'
include {deepsomatic_tn_step1; deepsomatic_step2; deepsomatic_step3;
        deepsomatic_tonly_step1; deepsomatic_tonly_step2;
        deepsomatic_step3 as deepsomatic_tonly_step3  } from "../../modules/local/deepsomatic.nf"

include {combineVariants as combineVariants_vardict; combineVariants as combineVariants_vardict_tonly;
    combineVariants as combineVariants_varscan; combineVariants as combineVariants_varscan_tonly;
    combineVariants_alternative as combineVariants_deepsomatic; combineVariants_alternative as combineVariants_deepsomatic_tonly;
    combineVariants as combineVariants_sage; combineVariants as combineVariants_sage_tonly;
    combineVariants_alternative as combineVariants_lofreq; combineVariants as combineVariants_muse;
    combineVariants_alternative as combineVariants_octopus; 
    combineVariants_alternative as combineVariants_octopus_tonly;
    combinemafs_tn; somaticcombine; somaticcombine as somaticcombine_ffpe;
    combinemafs_tonly;somaticcombine_tonly; somaticcombine_tonly as somaticcombine_tonly_ffpe} from '../../modules/local/combinefilter.nf'

include {sobdetect_pass1 as sobdetect_pass1_mutect2; sobdetect_pass2 as sobdetect_pass2_mutect2; 
    sobdetect_metrics as sobdetect_metrics_mutect2; sobdetect_cohort_params as sobdetect_cohort_params_mutect2;
    sobdetect_pass1 as sobdetect_pass1_octopus; sobdetect_pass2 as sobdetect_pass2_octopus; 
    sobdetect_metrics as sobdetect_metrics_octopus; sobdetect_cohort_params as sobdetect_cohort_params_octopus;
    sobdetect_pass1 as sobdetect_pass1_strelka; sobdetect_pass2 as sobdetect_pass2_strelka; 
    sobdetect_metrics as sobdetect_metrics_strelka; sobdetect_cohort_params as sobdetect_cohort_params_strelka;
    sobdetect_pass1 as sobdetect_pass1_lofreq; sobdetect_pass2 as sobdetect_pass2_lofreq; 
    sobdetect_metrics as sobdetect_metrics_lofreq; sobdetect_cohort_params as sobdetect_cohort_params_lofreq;
    sobdetect_pass1 as sobdetect_pass1_muse; sobdetect_pass2 as sobdetect_pass2_muse; 
    sobdetect_metrics as sobdetect_metrics_muse; sobdetect_cohort_params as sobdetect_cohort_params_muse;
    sobdetect_pass1 as sobdetect_pass1_vardict; sobdetect_pass2 as sobdetect_pass2_vardict; 
    sobdetect_metrics as sobdetect_metrics_vardict; sobdetect_cohort_params as sobdetect_cohort_params_vardict;
    sobdetect_pass1 as sobdetect_pass1_varscan; sobdetect_pass2 as sobdetect_pass2_varscan; 
    sobdetect_metrics as sobdetect_metrics_varscan; sobdetect_cohort_params as sobdetect_cohort_params_varscan;
    //Tumor Only
    sobdetect_pass1 as sobdetect_pass1_mutect2_tonly; sobdetect_pass2 as sobdetect_pass2_mutect2_tonly; 
    sobdetect_metrics as sobdetect_metrics_mutect2_tonly; sobdetect_cohort_params as sobdetect_cohort_params_mutect2_tonly;
    sobdetect_pass1 as sobdetect_pass1_octopus_tonly; sobdetect_pass2 as sobdetect_pass2_octopus_tonly; 
    sobdetect_metrics as sobdetect_metrics_octopus_tonly; sobdetect_cohort_params as sobdetect_cohort_params_octopus_tonly;
    sobdetect_pass1 as sobdetect_pass1_vardict_tonly; sobdetect_pass2 as sobdetect_pass2_vardict_tonly; 
    sobdetect_metrics as sobdetect_metrics_vardict_tonly; sobdetect_cohort_params as sobdetect_cohort_params_vardict_tonly;
    sobdetect_pass1 as sobdetect_pass1_varscan_tonly; sobdetect_pass2 as sobdetect_pass2_varscan_tonly; 
    sobdetect_metrics as sobdetect_metrics_varscan_tonly; sobdetect_cohort_params as sobdetect_cohort_params_varscan_tonly

    } from "../../modules/local/ffpe.nf"


include {annotvep_tn as annotvep_tn_mut2_ffpe; annotvep_tn as annotvep_tn_strelka_ffpe;
    annotvep_tn as annotvep_tn_varscan_ffpe; annotvep_tn as annotvep_tn_vardict_ffpe; annotvep_tn as annotvep_tn_octopus_ffpe;
    annotvep_tn as annotvep_tn_lofreq_ffpe; annotvep_tn as annotvep_tn_muse_ffpe; annotvep_tn as annotvep_tn_sage_ffpe;
    annotvep_tn as annotvep_tn_deepsomatic_ffpe;
    annotvep_tn as annotvep_tn_combined_ffpe; 
    annotvep_tonly as annotvep_tonly_varscan_ffpe; annotvep_tonly as annotvep_tonly_vardict_ffpe;
    annotvep_tonly as annotvep_tonly_mut2_ffpe; annotvep_tonly as annotvep_tonly_octopus_ffpe; 
    annotvep_tonly as annotvep_tonly_sage_ffpe; annotvep_tonly as annotvep_tonly_deepsomatic_ffpe;
    annotvep_tonly as annotvep_tonly_combined_ffpe} from '../../modules/local/annotvep.nf'


include {annotvep_tn as annotvep_tn_mut2; annotvep_tn as annotvep_tn_strelka;
    annotvep_tn as annotvep_tn_varscan; annotvep_tn as annotvep_tn_vardict; annotvep_tn as annotvep_tn_octopus;
    annotvep_tn as annotvep_tn_lofreq; annotvep_tn as annotvep_tn_muse; annotvep_tn as annotvep_tn_sage;
    annotvep_tn as annotvep_tn_deepsomatic;
    annotvep_tn as annotvep_tn_combined; 
    annotvep_tonly as annotvep_tonly_varscan; annotvep_tonly as annotvep_tonly_vardict;
    annotvep_tonly as annotvep_tonly_mut2; annotvep_tonly as annotvep_tonly_octopus; 
    annotvep_tonly as annotvep_tonly_sage; annotvep_tonly as annotvep_tonly_deepsomatic;
    annotvep_tonly as annotvep_tonly_combined} from '../../modules/local/annotvep.nf'

include {svaba_somatic} from '../../modules/local/svaba.nf'
include {manta_somatic} from '../../modules/local/manta.nf'
include {gridss_somatic} from '../../modules/local/gridss.nf'
include {survivor_sv; 
    gunzip as gunzip_manta; gunzip as gunzip_gridss; 
    annotsv_tn as annotsv_survivor_tn
    annotsv_tn as annotsv_gridss; annotsv_tn as annotsv_svaba; annotsv_tn as annotsv_manta} from '../../modules/local/annotsv.nf'

include {amber_tn; cobalt_tn; purple; purple_novc} from '../../modules/local/purple.nf'
include {sequenza; seqz_sequenza_bychr} from '../../modules/local/sequenza.nf'
include {freec; freec_paired; freec_paired_exome} from '../../modules/local/freec.nf'
include {ascat_tn; ascat_tn_exome} from '../../modules/local/ascat.nf'
include {cnvkit; cnvkit_exome } from '../../modules/local/cnvkit.nf'



//Workflows
workflow DETERMINEBAM {
    if(params.bam_input){
        params.BAMINPUT=true
    }else if(params.file_input){
            file(params.file_input).text
    }

}

workflow INPUT {

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
        sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true) | view()
                       | ifEmpty("sample sheet not found" )
                       | splitCsv(header:true, sep: "\t", strip:true)
                       | map { row -> tuple(
                        row.Tumor,
                        row.Normal
                       )
                                  }
    } else {
        sample_sheet=fastqinput.map{samplename,f1 -> tuple (
             samplename)}
    }

    emit:
        fastqinput
        sample_sheet

}

workflow ALIGN {
    take:
        fastqinput
        sample_sheet

    main:
    if (params.intervals){
        intervalbedin = Channel.fromPath(params.intervals)
    }else{
        intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
    }
    matchbed(intervalbedin) | splitinterval

    fastp(fastqinput)
    bwamem2(fastp.out)
    
    //Indel Realignment
    if (params.indelrealign){ 
        bwaindelre = bwamem2.out | indelrealign 
        bqsrbambyinterval = bwaindelre.combine(splitinterval.out.flatten())
        bambyinterval = bwaindelre.combine(splitinterval.out.flatten())

        bqsr_ir(bqsrbambyinterval)
        
        bqsrs = bqsr_ir.out             
            | groupTuple 
            | map { samplename,beds -> 
                tuple( samplename, beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )} 
        gatherbqsr(bqsrs)

        tobqsr=bwaindelre.combine(gatherbqsr.out,by:0)
        applybqsr(tobqsr)

        bamwithsample=applybqsr.out.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(applybqsr.out,by:0).map{it.swap(3,0)}  

    }else{ 
        bqsrbambyinterval=bwamem2.out.combine(splitinterval.out.flatten())
        bambyinterval=bwamem2.out.combine(splitinterval.out.flatten())

        bqsr(bqsrbambyinterval)
        bqsrs=bqsr.out | groupTuple
            | map { samplename,beds -> 
                tuple( samplename,
                beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )}
        gatherbqsr(bqsrs)

        tobqsr=bwamem2.out.combine(gatherbqsr.out,by:0)
        applybqsr(tobqsr)

        bamwithsample=applybqsr.out.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(applybqsr.out,by:0).map{it.swap(3,0)}
 
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
}

workflow GL {
    take:
        sample_sheet
        bambyinterval
    
    main:
    //Keep Only the NormalSamples
    bambyinterval_normonly=sample_sheet | map{t,n -> tuple(n)} | unique() | join(bambyinterval) 

    deepvariant_step1(bambyinterval_normonly) | deepvariant_step2  
        | deepvariant_step3 | groupTuple
        | multiMap{samplename,vcf,vcf_tbi,gvcf,gvcf_tbi -> 
        vcf: tuple(samplename,vcf.toSorted{it -> (it.name =~ /${samplename}_(.*?).bed.vcf.gz/)[0][1].toInteger()},vcf_tbi,"vcf")
        gvcf: tuple(samplename,gvcf.toSorted{it -> (it.name =~ /${samplename}_(.*?).bed.gvcf.gz/)[0][1].toInteger()},gvcf_tbi,"gvcf")
            }
    | set{dv_out} 
    dv_out.vcf | bcfconcat_vcf 
    dv_out.gvcf | bcfconcat_gvcf | map{sample,gvcf,index -> gvcf} 
        | collect 
        | glnexus
    deepvariant_out=bcfconcat_vcf.out | join(bcfconcat_gvcf.out) 

    emit:
        glnexusout=glnexus.out
        bcfout=deepvariant_out

}

workflow VC {
    take:
    //Input is the BAMby interval
        bamwithsample
        splitout
        sample_sheet

    main:
    //Create Pairing for TN (in case of dups)
    sample_sheet_paired=sample_sheet | map{tu,no -> tuple ("${tu}_vs_${no}",tu, no)} | view()
    bambyinterval=bamwithsample.combine(splitout.flatten()) 

    bambyinterval 
        | multiMap {tumorname,tumor,tumorbai,normalname,normalbam,normalbai,bed -> 
        t1: tuple(tumorname,tumor,tumorbai,bed)
        n1: tuple(normalname,normalbam,normalbai,bed)
            }
        | set{bambyinterval_tonly}
        
        bambyinterval_t=bambyinterval_tonly.t1 |
            concat(bambyinterval_tonly.n1) | unique() 

    //Prep Pileups
    call_list = params.callers.split(',') as List
    call_list_tonly = params.tonlycallers.split(',') as List
    call_list_tonly = call_list.intersect(call_list_tonly)

    vc_all=Channel.empty()
    vc_tonly=Channel.empty()

    //Common for Mutect2/Varscan
    if ("mutect2" in call_list | "varscan" in call_list){ 
        bambyinterval | 
            map{tumorname,tumor,tumorbai,normalname,normal,normalbai,bed -> tuple(tumorname,tumor,tumorbai,bed,"tpileup")} |
            unique |
            pileup_paired_t 
        bambyinterval | 
            map{tumorname,tumor,tumorbai,normalname,normal,normalbai,bed -> tuple(normalname,normal,normalbai,bed,"npileup")} |
            unique |
            pileup_paired_n

        pileup_paired_t.out | groupTuple |
            multiMap { samplename, pileups -> 
                tout: tuple( samplename,
                    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tpileup.table/)[0][1].toInteger() } )
                tonly: tuple( samplename,
                    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tpileup.table/)[0][1].toInteger() } )
                    }
            | set{pileup_paired_tout}
        
        pileup_paired_n.out | groupTuple |
            multiMap { normalname, pileups -> 
                nout: tuple (normalname,
                    pileups.toSorted{ it -> (it.name =~ /${normalname}_(.*?).npileup.table/)[0][1].toInteger() } )
                nonly: tuple (normalname,
                    pileups.toSorted{ it -> (it.name =~ /${normalname}_(.*?).npileup.table/)[0][1].toInteger() } )
                    }
                | set{pileup_paired_nout}

        pileup_paired_match=sample_sheet_paired |map{id,t,n-> tuple(t,id,n)} | combine(pileup_paired_tout.tout,by:0) | 
            map{it.swap(2,0)} |  combine(pileup_paired_nout.nout,by:0)  |map{no,id,tu,tpi,npi->tuple(tu,no,tpi,npi)}
 
        //pileup_paired_match=pileup_paired_tout.tout.join(pileup_paired_nout.nout,by:[0,1])
        contamination_paired(pileup_paired_match)

        if (!params.no_tonly){
        pileup_all=pileup_paired_tout.tonly | concat(pileup_paired_nout.nonly) 
        contamination_tumoronly(pileup_all) 
        }
    }

    if ("mutect2" in call_list){
    //Paired Mutect2
        mutect2(bambyinterval)
        mutect2.out.groupTuple(by:[0,1])
        | multiMap { tumor,normal,vcfs,f1r2,stats ->
        mut2out_lor: tuple("${tumor}_vs_${normal}",
                f1r2.toSorted{ it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } )
        mut2out_mstats:  tuple( "${tumor}_vs_${normal}",
                stats.toSorted{ it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).mut2.vcf.gz.stats/)[0][1].toInteger() })
        allmut2tn: tuple( "${tumor}_vs_${normal}",
                vcfs.toSorted{ it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).mut2.vcf.gz/)[0][1].toInteger() } )
                    }
        | set{mut2out}

    learnreadorientationmodel(mut2out.mut2out_lor)
    mergemut2stats(mut2out.mut2out_mstats)

    mutect2_in=mut2out.allmut2tn
    | join(mergemut2stats.out)
    | join(learnreadorientationmodel.out)
    | map{t,vcf,stats,ro -> tuple(t.split('_vs_')[0],t.split('_vs_')[1],vcf,stats,ro)}
    | join(contamination_paired.out)
    | mutect2filter
    | join(sample_sheet_paired)
    | map{sample,markedvcf,markedindex,normvcf,normindex,stats,tumor,normal -> tuple(tumor,normal,"mutect2",normvcf,normindex)}

    annotvep_tn_mut2(mutect2_in)
    vc_all = vc_all|concat(mutect2_in)


    //Mutect2 Tumor Only
    if (!params.no_tonly){
        mutect2_t_tonly(bambyinterval_t)
        | groupTuple()
        | multiMap { tumorid,vcfs,f1r2,stats ->
        mut2tout_lor: tuple(tumorid,
                f1r2.toSorted{ it -> (it.name =~ /${tumorid}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } )
        mut2tonly_mstats:  tuple( tumorid,
                stats.toSorted{ it -> (it.name =~ /${tumorid}_(.*?).tonly.mut2.vcf.gz.stats/)[0][1].toInteger() })
        allmut2tonly: tuple(tumorid,
                vcfs.toSorted{ it -> (it.name =~ /${tumorid}_(.*?).tonly.mut2.vcf.gz/)[0][1].toInteger() } )
        }
        | set{mut2tonlyout}

        learnreadorientationmodel_tonly(mut2tonlyout.mut2tout_lor) 
        mergemut2stats_tonly(mut2tonlyout.mut2tonly_mstats) 

        mutect2_in_tonly=mut2tonlyout.allmut2tonly
        | join(mergemut2stats_tonly.out)
        | join(learnreadorientationmodel_tonly.out)
        | join(contamination_tumoronly.out)
        | mutect2filter_tonly
        | join(sample_sheet)
        | map{tumor,markedvcf,markedindex,normvcf,normindex,stats,normal -> tuple(tumor,"mutect2_tonly",normvcf,normindex)}
        annotvep_tonly_mut2(mutect2_in_tonly)

        vc_tonly = vc_tonly | concat(mutect2_in_tonly) 
        }
        
        
    }

    if ("strelka" in call_list){
        //Strelka TN
        strelka_in=strelka_tn(bambyinterval) | groupTuple(by:[0,1])
        | map { tumor,normal,vcfs,vcfindex,indels,indelindex -> tuple("${tumor}_vs_${normal}",
        vcfs.toSorted{ it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).somatic.snvs.vcf.gz/)[0][1].toInteger() },vcfindex,
        indels.toSorted{ it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).somatic.indels.vcf.gz/)[0][1].toInteger() },indelindex)}
        | combineVariants_strelka | join(sample_sheet_paired)
        | map{sample,markedvcf,markedindex,finalvcf,finalindex,tumor,normal -> tuple(tumor,normal,"strelka",finalvcf,finalindex)}
        | convert_strelka
        annotvep_tn_strelka(strelka_in)

        vc_all=vc_all|concat(strelka_in)

    }

    if ("vardict" in call_list){
        //Vardict TN
        vardict_in=vardict_tn(bambyinterval) | groupTuple(by:[0,1])
        | map{tumor,normal,vcf-> tuple("${tumor}_vs_${normal}",vcf.toSorted{it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).vardict.vcf/)[0][1].toInteger()},"vardict")}
        | combineVariants_vardict | join(sample_sheet_paired)
        | map{sample,marked,markedindex,normvcf,normindex,tumor,normal ->tuple(tumor,normal,"vardict",normvcf,normindex)}
        annotvep_tn_vardict(vardict_in)

        vc_all=vc_all|concat(vardict_in)

        //Vardict TOnly
        if (!params.no_tonly){
        vardict_in_tonly=vardict_tonly(bambyinterval_t) 
        | groupTuple()
        | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.vardict.vcf/)[0][1].toInteger()},"vardict_tonly")}
        | combineVariants_vardict_tonly | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex,normal ->tuple(tumor,"vardict_tonly",normvcf,normindex)}
        annotvep_tonly_vardict(vardict_in_tonly)

        vc_tonly=vc_tonly|concat(vardict_in_tonly) 
        }

    }

    if ("varscan" in call_list){
        //VarScan TN
        varscan_in=bambyinterval.combine(contamination_paired.out,by:0) 
        | varscan_tn | groupTuple(by:[0,1]) 
        | map{tumor,normal,vcf-> tuple("${tumor}_vs_${normal}",vcf.toSorted{it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).varscan.vcf.gz/)[0][1].toInteger()},"varscan")}
        | combineVariants_varscan | join(sample_sheet_paired)
        | map{sample,marked,markedindex,normvcf,normindex,tumor,normal ->tuple(tumor,normal,"varscan",normvcf,normindex)}
        annotvep_tn_varscan(varscan_in)

        vc_all=vc_all|concat(varscan_in)

        if (!params.no_tonly){
        //VarScan TOnly
        varscan_in_tonly=bambyinterval_t.combine(contamination_tumoronly.out,by:0)  
        | varscan_tonly  | groupTuple 
        | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.varscan.vcf.gz/)[0][1].toInteger()},"varscan_tonly")}  
        | combineVariants_varscan_tonly 
        | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex,normal ->tuple(tumor,"varscan_tonly",normvcf,normindex)}
        annotvep_tonly_varscan(varscan_in_tonly)

        vc_tonly=vc_tonly|concat(varscan_in_tonly) 
        }

        
    }

    //SAGE TN
    if ("sage" in call_list){
        sage_in=sage_tn(bamwithsample)
            | map{tu,no,vcf,vcfindex-> tuple("${tu}_vs_${no}",vcf,"sage")}
            | combineVariants_sage 
            | join(sample_sheet_paired)
            | map{sample,marked,markedindex,normvcf,normindex,tumor,normal->tuple(tumor,normal,"sage",normvcf,normindex)}
        annotvep_tn_sage(sage_in)

        vc_all=vc_all | concat(sage_in)

        if (!params.no_tonly){ 
        sage_in_tonly=bamwithsample | map{tumor,tbam,tbai,norm,nbam,nbai -> tuple(tumor,tbam,tbai)} 
            | sage_tonly 
            | map{samplename,vcf,vcfindex->tuple(samplename,vcf,"sage_tonly")} 
            | combineVariants_sage_tonly
            | join(sample_sheet) 
            | map{tumor,marked,markedindex,normvcf,normindex,normal ->tuple(tumor,"sage_tonly",normvcf,normindex)}
        annotvep_tonly_sage(sage_in_tonly)
        vc_tonly=vc_tonly | concat(sage_in_tonly) 
        }
    }

    //Lofreq TN
    if ("lofreq" in call_list){
        lofreq_in=lofreq_tn(bambyinterval) | groupTuple(by:[0,1])
            | map{tu,no,snv,dbsnv,indel,dbindel,vcf,vcfindex-> tuple("${tu}_vs_${no}",vcf.toSorted{it -> (it.name =~ /${tu}_vs_${no}_(.*?)_lofreq.vcf.gz/)[0][1].toInteger()},vcfindex,"lofreq")}
            | combineVariants_lofreq | join(sample_sheet_paired)
            | map{sample,marked,markedindex,normvcf,normindex,tumor,normal->tuple(tumor,normal,"lofreq",normvcf,normindex)}
        annotvep_tn_lofreq(lofreq_in)

        vc_all=vc_all|concat(lofreq_in)
    }


    //DeepSomatic TN
    if ("deepsomatic" in call_list){
        deepsomatic_in = deepsomatic_tn_step1(bambyinterval) 
            | map{tname,nname,tf,tfjson,bed -> tuple("${tname}_vs_${nname}",tf,tfjson,bed)} 
            | deepsomatic_step2  
            | deepsomatic_step3 | groupTuple 
            | map{samplename,vcf,vcf_tbi -> 
                tuple(samplename,vcf.toSorted{it -> (it.name =~ /${samplename}_(.*?).bed.vcf.gz/)[0][1].toInteger()},vcf_tbi,"deepsomatic")
                }
            | combineVariants_deepsomatic             
            | join(sample_sheet_paired)
            | map{sample,marked,markedindex,normvcf,normindex,tumor,normal->tuple(tumor,normal,"deepsomatic",normvcf,normindex)}
        annotvep_tn_deepsomatic(deepsomatic_in)

        vc_all=vc_all|concat(deepsomatic_in)

        //DeepSomatic TOnly
          if (!params.no_tonly){
          deepsomatic_tonly_in = deepsomatic_tonly_step1(bambyinterval_t) 
            | deepsomatic_tonly_step2  
            | deepsomatic_tonly_step3 | groupTuple 

            | map{samplename,vcf,vcf_tbi -> 
                tuple(samplename,vcf.toSorted{it -> (it.name =~ /${samplename}_(.*?).bed.vcf.gz/)[0][1].toInteger()},vcf_tbi,"deepsomatic_tonly")
             } 

            | combineVariants_deepsomatic_tonly           
            | join(sample_sheet) 
            | map{tumor,marked,markedindex,normvcf,normindex,normal->tuple(tumor,"deepsomatic_tonly",normvcf,normindex)}

        annotvep_tonly_deepsomatic(deepsomatic_tonly_in)

        vc_tonly=vc_tonly | concat(deepsomatic_tonly_in) 
          
          }
          
    }

    //MuSE TN
    if ("muse" in call_list){
        muse_in=muse_tn(bamwithsample)
            | map{tumor,normal,vcf-> tuple("${tumor}_vs_${normal}",vcf,"muse")}
            | combineVariants_muse | join(sample_sheet_paired)
            | map{sample,marked,markedindex,normvcf,normindex,tumor,normal ->tuple(tumor,normal,"muse",normvcf,normindex)}
        annotvep_tn_muse(muse_in)

        vc_all=vc_all|concat(muse_in)
    }

    //Octopus TN
    if ("octopus" in call_list){
        octopus_in=octopus_tn(bambyinterval) | bcftools_index_octopus
            | groupTuple()
            | map{samplename,vcf,vcfindex-> tuple(samplename,vcf.toSorted{it->(it.name =~ /${samplename}_(.*).octopus.vcf.gz/)[0][1].toInteger()},vcfindex,"octopus")}
            | combineVariants_octopus
            | map{samplename,marked,markedindex,normvcf,normindex ->
                tuple(samplename.split('_vs_')[0],samplename.split('_vs_')[1],"octopus",normvcf,normindex)}
        annotvep_tn_octopus(octopus_in)
        octopus_in = octopus_in | octopus_convertvcf 
            |  map{tumor,normal,vcf,vcfindex ->tuple(tumor,normal,"octopus",vcf,vcfindex)} 
        vc_all=vc_all|concat(octopus_in)

    //Octopus TOnly
        if (!params.no_tonly){ 
        octopus_in_tonly=octopus_tonly(bambyinterval_t)
            | bcftools_index_octopus_tonly
            | groupTuple() 
            | map{samplename,vcf,vcfindex->tuple(samplename,vcf.toSorted{it->(it.name =~ /${samplename}_(.*).tonly.octopus.vcf.gz/)[0][1].toInteger()},vcfindex,"octopus_tonly")} 
            | combineVariants_octopus_tonly
            | join(sample_sheet) 
            | map{tumor,marked,markedindex,normvcf,normindex,normal ->tuple(tumor,"octopus_tonly",normvcf,normindex)}
        annotvep_tonly_octopus(octopus_in_tonly)
        octopus_in_tonly_sc=octopus_in_tonly | octopus_convertvcf_tonly
            | map{tumor,vcf,vcfindex ->tuple(tumor,"octopus_tonly",vcf,vcfindex)} 
        vc_tonly=vc_tonly|concat(octopus_in_tonly_sc) 
        }
        
    }

 
    //FFPE Steps 
    if(params.ffpe){
        vc_ffpe_paired=Channel.empty()
        vc_ffpe_tonly=Channel.empty()
        bamwithsample1=bamwithsample | map{tumor,tbam,tbai,norm,nbam,nbai ->tuple(tumor,norm,tbam,tbai,nbam,nbai)}

        if('mutect2' in call_list){
            mutect2_p1=bamwithsample1 | join(mutect2_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | sobdetect_pass1_mutect2  
            mutect2_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_mutect2

            mutect2_p2 = bamwithsample1 
                | join(mutect2_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | combine(sobdetect_cohort_params_mutect2.out) 
                | sobdetect_pass2_mutect2
            mutect2_p1_vcfs=mutect2_p1 | map{sample,vcf,info->vcf} | collect 
            mutect2_p2_vcfs=mutect2_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} | collect
            sobdetect_metrics_mutect2(mutect2_p1_vcfs,mutect2_p2_vcfs)
            
            mutect2_ffpe_out=mutect2_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample.split("_vs_")[0],sample.split("_vs_")[1],"mutect2",filtvcf,vcftbi)} 
            annotvep_tn_mut2_ffpe(mutect2_ffpe_out)
            vc_ffpe_paired=vc_ffpe_paired |concat(mutect2_ffpe_out)

            if (!params.no_tonly){ 
                mutect2_tonly_p1=bamwithsample1 | join(mutect2_in_tonly) 
                    | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                    | sobdetect_pass1_mutect2_tonly
                mutect2_tonly_p1 | map{sample,vcf,info->info}
                    | collect
                    | sobdetect_cohort_params_mutect2_tonly

                mutect2_tonly_p2 = bamwithsample1 
                    | join(mutect2_in_tonly)
                    | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                    | combine(sobdetect_cohort_params_mutect2_tonly.out) 
                    | sobdetect_pass2_mutect2_tonly
                mutect2_tonly_p1_vcfs=mutect2_tonly_p1 | map{sample,vcf,info->vcf} |collect 
                mutect2_tonly_p2_vcfs=mutect2_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
                sobdetect_metrics_mutect2_tonly(mutect2_tonly_p1_vcfs,mutect2_tonly_p2_vcfs)
                
                mutect2_tonly_ffpe_out=mutect2_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample,"mutect2_tonly",filtvcf,vcftbi)} 
                annotvep_tonly_mut2_ffpe(mutect2_tonly_ffpe_out)
                vc_ffpe_tonly=vc_ffpe_tonly |concat(mutect2_tonly_ffpe_out)
            
            }

        }

        if('octopus' in call_list){
            octopus_p1=bamwithsample1 | join(octopus_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | sobdetect_pass1_octopus  
            octopus_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_octopus
            octopus_p2 = bamwithsample1 
                | join(octopus_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | combine(sobdetect_cohort_params_octopus.out) 
                | sobdetect_pass2_octopus
            octopus_p1_vcfs=octopus_p1 | map{sample,vcf,info->vcf} | collect 
            octopus_p2_vcfs=octopus_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} | collect
            sobdetect_metrics_octopus(octopus_p1_vcfs,octopus_p2_vcfs)

            octopus_ffpe_out=octopus_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample.split("_vs_")[0],sample.split("_vs_")[1],"octopus",filtvcf,vcftbi)} 
            annotvep_tn_octopus_ffpe(octopus_ffpe_out)
            vc_ffpe_paired=vc_ffpe_paired |concat(octopus_ffpe_out)

            if (!params.no_tonly){ 
                octopus_tonly_p1=bamwithsample1 | join(octopus_in_tonly) 
                    | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                    | sobdetect_pass1_octopus_tonly
                octopus_tonly_p1 | map{sample,vcf,info->info}
                    | collect
                    | sobdetect_cohort_params_octopus_tonly

                octopus_tonly_p2 = bamwithsample1 
                    | join(octopus_in_tonly)
                    | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                    | combine(sobdetect_cohort_params_octopus_tonly.out) 
                    | sobdetect_pass2_octopus_tonly
                octopus_tonly_p1_vcfs=octopus_tonly_p1 | map{sample,vcf,info->vcf} |collect 
                octopus_tonly_p2_vcfs=octopus_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
                sobdetect_metrics_octopus_tonly(octopus_tonly_p1_vcfs,octopus_tonly_p2_vcfs)

                octopus_tonly_ffpe_out=octopus_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample,"octopus_tonly",filtvcf,vcftbi)} 
                annotvep_tonly_octopus_ffpe(octopus_tonly_ffpe_out)
                vc_ffpe_tonly=vc_ffpe_tonly |concat(octopus_tonly_ffpe_out)
            }        
        }

        if('strelka' in call_list){
            strelka_p1=bamwithsample1 | join(strelka_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | sobdetect_pass1_strelka
            strelka_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_strelka
            strelka_p2 = bamwithsample1 
                | join(strelka_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | combine(sobdetect_cohort_params_strelka.out) 
                | sobdetect_pass2_strelka
                
            strelka_p1_vcfs=strelka_p1 | map{sample,vcf,info->vcf} |collect 
            strelka_p2_vcfs=strelka_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} | collect
            sobdetect_metrics_strelka(strelka_p1_vcfs,strelka_p2_vcfs)

            strelka_ffpe_out=strelka_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample.split("_vs_")[0],sample.split("_vs_")[1],"strelka",filtvcf,vcftbi)} 
            annotvep_tn_strelka_ffpe(strelka_ffpe_out)
            vc_ffpe_paired=vc_ffpe_paired |concat(strelka_ffpe_out)

        }

        if('lofreq' in call_list){
            lofreq_p1=bamwithsample1 | join(lofreq_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | sobdetect_pass1_lofreq  
            lofreq_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_lofreq
            lofreq_p2 = bamwithsample1 
                | join(lofreq_in,by:[0,1])
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
                | combine(sobdetect_cohort_params_lofreq.out) 
                | sobdetect_pass2_lofreq
            lofreq_p1_vcfs=lofreq_p1 | map{sample,vcf,info->vcf} |collect 
            lofreq_p2_vcfs=lofreq_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
            sobdetect_metrics_lofreq(lofreq_p1_vcfs,lofreq_p2_vcfs)

            lofreq_ffpe_out=lofreq_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample.split("_vs_")[0],sample.split("_vs_")[1],"lofreq",filtvcf,vcftbi)} 
            annotvep_tn_lofreq_ffpe(lofreq_ffpe_out)
            vc_ffpe_paired=vc_ffpe_paired |concat(lofreq_ffpe_out)
        }

    if('muse' in call_list){
        muse_p1=bamwithsample1 | join(muse_in,by:[0,1])
            | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
            | sobdetect_pass1_muse  
        muse_p1 | map{sample,vcf,info->info}
            | collect
            | sobdetect_cohort_params_muse
        muse_p2 = bamwithsample1 
            | join(muse_in,by:[0,1])
            | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
            | combine(sobdetect_cohort_params_muse.out) 
            | sobdetect_pass2_muse
        muse_p1_vcfs=muse_p1 | map{sample,vcf,info->vcf} |collect 
        muse_p2_vcfs=muse_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
        sobdetect_metrics_muse(muse_p1_vcfs,muse_p2_vcfs)

        muse_ffpe_out=muse_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample.split("_vs_")[0],sample.split("_vs_")[1],"muse",filtvcf,vcftbi)} 
        annotvep_tn_muse_ffpe(muse_ffpe_out)
        vc_ffpe_paired=vc_ffpe_paired |concat(muse_ffpe_out)
    }

    if('vardict' in call_list){
        vardict_p1=bamwithsample1 | join(vardict_in,by:[0,1])
            | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
            | sobdetect_pass1_vardict  
        vardict_p1 | map{sample,vcf,info->info}
            | collect
            | sobdetect_cohort_params_vardict
        vardict_p2 = bamwithsample1 
            | join(vardict_in,by:[0,1])
            | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
            | combine(sobdetect_cohort_params_vardict.out) 
            | sobdetect_pass2_vardict
        vardict_p1_vcfs=vardict_p1 | map{sample,vcf,info->vcf} |collect 
        vardict_p2_vcfs=vardict_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
        sobdetect_metrics_vardict(vardict_p1_vcfs,vardict_p2_vcfs)

        vardict_ffpe_out=vardict_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample.split("_vs_")[0],sample.split("_vs_")[1],"vardict",filtvcf,vcftbi)} 
        annotvep_tn_vardict_ffpe(vardict_ffpe_out)
        vc_ffpe_paired=vc_ffpe_paired |concat(vardict_ffpe_out)

        if (!params.no_tonly){ 
            vardict_tonly_p1=bamwithsample1 | join(vardict_in_tonly) 
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                | sobdetect_pass1_vardict_tonly
            vardict_tonly_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_vardict_tonly

            vardict_tonly_p2 = bamwithsample1 
                | join(vardict_in_tonly)
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                | combine(sobdetect_cohort_params_vardict_tonly.out) 
                | sobdetect_pass2_vardict_tonly
            vardict_tonly_p1_vcfs=vardict_tonly_p1 | map{sample,vcf,info->vcf} |collect 
            vardict_tonly_p2_vcfs=vardict_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
            sobdetect_metrics_vardict_tonly(vardict_tonly_p1_vcfs,vardict_tonly_p2_vcfs)

            vardict_tonly_ffpe_out=vardict_tonly_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample,"vardict_tonly",filtvcf,vcftbi)} 
            annotvep_tonly_vardict_ffpe(vardict_tonly_ffpe_out)
            vc_ffpe_tonly=vc_ffpe_tonly |concat(vardict_tonly_ffpe_out)
        }
    }

    if('varscan' in call_list){
        varscan_p1=bamwithsample1 | join(varscan_in,by:[0,1])
            | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
            | sobdetect_pass1_varscan  
        varscan_p1 | map{sample,vcf,info->info}
            | collect
            | sobdetect_cohort_params_varscan
        varscan_p2 = bamwithsample1 
            | join(varscan_in,by:[0,1])
            | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple("${tumor}_vs_${normal}",normvcf,tbam,vc)}
            | combine(sobdetect_cohort_params_varscan.out) 
            | sobdetect_pass2_varscan
        varscan_p1_vcfs=varscan_p1 | map{sample,vcf,info->vcf} |collect 
        varscan_p2_vcfs=varscan_p2 | map{sample,vcf,info,filtvcf,vcftbi->vcf} |collect
        sobdetect_metrics_varscan(varscan_p1_vcfs,varscan_p2_vcfs)

        varscan_ffpe_out=varscan_p2 | map{sample,vcf,info,filtvcf,vcftbi->tuple(sample.split("_vs_")[0],sample.split("_vs_")[1],"varscan",filtvcf,vcftbi)} 
        annotvep_tn_varscan_ffpe(varscan_ffpe_out)
        vc_ffpe_paired=vc_ffpe_paired | concat(varscan_ffpe_out)

        if (!params.no_tonly){ 
            varscan_tonly_p1=bamwithsample1 | join(varscan_in_tonly) 
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
                | sobdetect_pass1_varscan_tonly
            varscan_tonly_p1 | map{sample,vcf,info->info}
                | collect
                | sobdetect_cohort_params_varscan_tonly

            varscan_tonly_p2 = bamwithsample1 
                | join(varscan_in_tonly)
                | map{tumor,normal,tbam,tbai,nbam,nbai,vc,normvcf,tbi->tuple(tumor,normvcf,tbam,vc)}
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
    }

    //Combine All Variants Using VCF -> Annotate
    if (call_list.size()>1){
        vc_all | groupTuple(by:[0,1])
            | somaticcombine
            | map{tumor,normal,vcf,index ->tuple(tumor,normal,"combined",vcf,index)}
            | annotvep_tn_combined
        if(params.ffpe){
            vc_ffpe_paired | groupTuple(by:[0,1])
                | somaticcombine_ffpe
                | map{tumor,normal,vcf,index ->tuple(tumor,normal,"combined_ffpe",vcf,index)}
                | annotvep_tn_combined_ffpe
        }    
        if (!params.no_tonly & call_list_tonly.size()>1){
            vc_tonly | groupTuple() 
                | somaticcombine_tonly
                | map{tumor,vcf,index ->tuple(tumor,"combined_tonly",vcf,index)}
                | annotvep_tonly_combined
        }
        if (!params.no_tonly & call_list_tonly.size()>1 & params.ffpe){
            vc_ffpe_tonly | groupTuple() 
                | somaticcombine_tonly_ffpe
                | map{tumor,vcf,index ->tuple(tumor,"combined_tonly_ffpe",vcf,index)}
                | annotvep_tonly_combined_ffpe
        }
    }
    

    if("sage" in call_list){
        somaticcall_input=sage_in
    }else if("mutect2" in call_list){
        somaticcall_input=mutect2_in
    }else if("mutect2" in call_list & params.ffpe){
        somaticcall_input=mutect2_ffpe_out
    }else{
        somaticcall_input=Channel.empty()
    }
    


    emit:
        somaticcall_input
   
}


workflow SV {
    take:
        bamwithsample

    main:
        svcall_list = params.svcallers.split(',') as List
        svout=Channel.empty()

        //Svaba
        if ("svaba" in svcall_list){
            svaba_out=svaba_somatic(bamwithsample)
            | map{ tumor,bps,contigs,discord,alignents,gindel,gsv,so_indel,so_sv,unfil_gindel,unfil_gsv,unfil_so_indel,unfil_sv,log ->
                tuple(tumor,so_sv,"svaba")}
            svout=svout | concat(svaba_out)
        }

        //Manta
        if ("manta" in svcall_list){
        manta_out=manta_somatic(bamwithsample)
        manta_out_forsv=manta_out    
            | map{tumor,normal,gsv,gsv_tbi,so_sv,so_sv_tbi,unfil_sv,unfil_sv_tbi,unfil_indel,unfil_indel_tbi ->
                tuple(tumor,so_sv,"manta")} | gunzip_manta
        //annotsv_manta(manta_out).ifEmpty("Empty SV input--No SV annotated")
            svout=svout | concat(manta_out_forsv)
        }

        //GRIDSS
        if ("gridss" in svcall_list){
        gridss_out=gridss_somatic(bamwithsample)
        gridss_out_forsv=gridss_out 
            | map{tumor,normal,vcf,index,bam,gripssvcf,gripsstbi,gripssfilt,filttbi ->
                tuple(tumor,gripssfilt,"gridss")} | gunzip_gridss
            svout=svout | concat(gridss_out_forsv)
        }

        if (svcall_list.size()>1){
            //Survivor
            svout | groupTuple
                | survivor_sv 
                | annotsv_survivor_tn 
                | ifEmpty("Empty SV input--No SV annotated")
         }
        
        if("gridss" in svcall_list){
            somaticsv_input=gridss_out 
                | map{tumor,normal,vcf,index,bam,gripssvcf,gripsstbi,gripssfilt,filttbi ->
                    tuple(tumor,normal,vcf,index,gripssfilt,filttbi)}
        }else if("manta" in svcall_list){
            somaticsv_input=manta_out 
                | map{tumor,normal,gsv,gsv_tbi,so_sv,so_sv_tbi,unfil_sv,unfil_sv_tbi,unfil_indel,unfil_indel_tbi ->
                    tuple(tumor,normal,unfil_sv,unfil_sv_tbi,so_sv,so_sv_tbi)}
        }else{
            somaticsv_input=Channel.empty()
        }
    
    emit:
        somaticsv_input
}

workflow CNVmouse {
    take:
        bamwithsample
        
    main:
        cnvcall_list = params.cnvcallers.split(',') as List

        //Sequenza (Preferred for Paired)
        if ("sequenza" in cnvcall_list){
        chrs=Channel.fromList(params.genomes[params.genome].chromosomes)
        seqzin=bamwithsample.map{tname,tumor,tbai,nname,norm,nbai->
            tuple("${tname}_${nname}",tname,tumor,tbai,nname,norm,nbai)}
        seqzin.combine(chrs) | seqz_sequenza_bychr
        seqz_sequenza_bychr.out.groupTuple()
            .map{pair, seqz -> tuple(pair, seqz.sort{it.name})}
            | sequenza
        }

        //FREEC Paired Mode
        if ("freec" in cnvcall_list){
            if(params.exome){
                FREECPAIR_SCRIPT = params.script_freecpaired_exome
                bamwithsample | freec_paired_exome
            }else{
                FREECPAIR_SCRIPT = params.script_freecpaired
                bamwithsample | freec_paired
            }

        //FREEC Unpaired Mode
        bamwithsample 
            | map{tname,tumor,tbai,nname,norm,nbai->tuple(tname,tumor,tbai)}
            | freec
        }

         //CNVKIT
        if ("cnvkit" in cnvcall_list){
            if(params.exome){
                matchbed_cnvkit(intervalbedin)
                bamwithsample | combine(matchbed_cnvkit.out) | cnvkit_exome
            }else{
                bamwithsample | cnvkit
            }
        }
}

workflow CNVhuman {
    take:
        bamwithsample
        somaticcall_input

    main:  
    if (params.intervals){
        intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
    }else{
        intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
    }

        cnvcall_list = params.cnvcallers.split(',') as List
        scinput = somaticcall_input|map{t1,n1,cal,vcf,ind -> tuple("${t1}_vs_${n1}",cal,vcf,ind)}
        
        if ("purple" in cnvcall_list){
            //Purple
            bamwithsample | amber_tn
            bamwithsample | cobalt_tn
            purplein=amber_tn.out.join(cobalt_tn.out,by:[0,1,2]) 
            purplein.join(scinput) 
            | map{id,t1,n1,amber,cobalt,vc,vcf,vcfindex -> tuple(id,t1,n1,amber,cobalt,vcf,vcfindex)}
            | purple
        }

        //Sequenza
        if ("sequenza" in cnvcall_list){
            chrs=Channel.fromList(params.genomes[params.genome].chromosomes)
            seqzin=bamwithsample.map{tname,tumor,tbai,nname,norm,nbai->
                tuple("${tname}_${nname}",tname,tumor,tbai,nname,norm,nbai)}
            seqzin.combine(chrs) | seqz_sequenza_bychr
            seqz_sequenza_bychr.out.groupTuple()
                .map{pair, seqz -> tuple(pair, seqz.sort{it.name})}
                | sequenza
        }

        //FREEC
        if ("freec" in cnvcall_list){
            if(params.exome){
                FREECPAIR_SCRIPT = params.script_freecpaired_exome
                bamwithsample | freec_paired_exome
            }else{
                FREECPAIR_SCRIPT = params.script_freecpaired
                bamwithsample | freec_paired
            }
        }
        
        //ASCAT
        if ("ascat" in cnvcall_list){
            if(params.exome){
                matchbed_ascat(intervalbedin)
                bamwithsample | combine(matchbed_ascat.out) | ascat_tn_exome
            }else{
                bamwithsample | ascat_tn
            }
        }
        
        //CNVKIT
        if ("cnvkit" in cnvcall_list){
            if(params.exome){
                matchbed_cnvkit(intervalbedin)
                bamwithsample | combine(matchbed_cnvkit.out) | cnvkit_exome
            }else{
                bamwithsample | cnvkit
            }
        }
}


workflow CNVhuman_novc {
    take:
        bamwithsample

    main:   
    if (params.intervals){
        intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
    }else{
        intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
    }

        cnvcall_list = params.cnvcallers.split(',') as List

        if ("purple" in cnvcall_list){
            //Purple
            bamwithsample | amber_tn 
            bamwithsample | cobalt_tn 
            purplein=amber_tn.out |join(cobalt_tn.out)  
            purplein | map{id,t1,n1,amber,t2,n2,cobalt -> tuple(id,t1,n1,amber,cobalt)}
                | purple_novc
        }

        if ("sequenza" in cnvcall_list){
            //Sequenza
            chrs=Channel.fromList(params.genomes[params.genome].chromosomes)
            seqzin=bamwithsample.map{tname,tumor,tbai,nname,norm,nbai->
                tuple("${tname}_${nname}",tname,tumor,tbai,nname,norm,nbai)}
            seqzin.combine(chrs) | seqz_sequenza_bychr
            seqz_sequenza_bychr.out.groupTuple()
                .map{pair, seqz -> tuple(pair, seqz.sort{it.name})}
                | sequenza
        }
        
        if ("freec" in cnvcall_list){
            //FREEC
            if(params.exome){
                FREECPAIR_SCRIPT = params.script_freecpaired_exome
                bamwithsample | freec_paired_exome
            }else{
                FREECPAIR_SCRIPT = params.script_freecpaired
                bamwithsample | freec_paired
            }
        }

        if ("ascat" in cnvcall_list){
            //ASCAT
            if(params.exome){
                matchbed_ascat(intervalbedin)
                bamwithsample |combine(matchbed_ascat.out) | ascat_tn_exome
            }else{
                bamwithsample | ascat_tn
            }
        }
        
        //CNVKIT
        if ("cnvkit" in cnvcall_list){
            if(params.exome){
                matchbed_cnvkit(intervalbedin)
                bamwithsample | combine(matchbed_cnvkit.out) | cnvkit_exome
            }else{
                bamwithsample | cnvkit
            }
        }


}



workflow QC_NOGL {
    take:
        fastqin
        fastpout
        applybqsr

    main:
    //QC Steps
    fc_lane(fastqin)
    fastq_screen(fastpout)
    kraken(fastqin)
    qualimap_bamqc(applybqsr)
    samtools_flagstats(applybqsr)
    fastqc(applybqsr)
    mosdepth(applybqsr)

    //Somalier
    somalier_extract(applybqsr)
    som_in=somalier_extract.out.collect()

    //Prep for MultiQC input
    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect()

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    mosdepth_out=mosdepth.out.collect()
    fastqc_out=fastqc.out.map{samplename,html,zip->tuple(html,zip)}.collect()

    samtools_flagstats_out=samtools_flagstats.out.collect()

    if(params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome.matches("mm10")){
        somalier_analysis_mouse(som_in)
        somalier_analysis_out=somalier_analysis_mouse.out.collect()
    }

    conall=fclane_out.concat(fqs_out,kraken_out,qualimap_out,
        samtools_flagstats_out,fastqc_out,mosdepth_out,
        somalier_analysis_out).flatten().toList()
    multiqc(conall)
}



workflow QC_GL {
    take:
        fastqin
        fastpout
        applybqsr
        glnexusout
        bcfout

    main:
    //QC Steps
    fc_lane(fastqin)
    fastq_screen(fastpout)
    kraken(fastqin)
    qualimap_bamqc(applybqsr)
    samtools_flagstats(applybqsr)
    mosdepth(applybqsr)
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

    //Prep for MultiQC input
    if(params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome.matches("mm10")){
        somalier_analysis_mouse(som_in)
        somalier_analysis_out=somalier_analysis_mouse.out.collect()
    }

    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect()

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    fastqc_out=fastqc.out.map{samplename,html,zip->tuple(html,zip)}.collect()
    samtools_flagstats_out=samtools_flagstats.out.collect()
    mosdepth_out=mosdepth.out.collect()
    bcftools_stats_out= bcftools_stats.out.collect()
    gatk_varianteval_out= gatk_varianteval.out.collect()
    snpeff_out=snpeff.out.collect()
    vcftools_out=vcftools.out
    collectvariantcallmetrics_out=collectvariantcallmetrics.out

    conall=fclane_out.concat(fqs_out,
        kraken_out,qualimap_out,mosdepth_out,
        fastqc_out,samtools_flagstats_out,bcftools_stats_out,
        gatk_varianteval_out,snpeff_out,vcftools_out,collectvariantcallmetrics_out,somalier_analysis_out).flatten().toList()
    multiqc(conall)
}

//QC_GL_BAMS
workflow QC_GL_BAM {
    take:
        applybqsr
        glnexusout
        bcfout

    main:
    //QC Steps
    qualimap_bamqc(applybqsr)
    samtools_flagstats(applybqsr)
    mosdepth(applybqsr)
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

    //Prep for MultiQC input
    if(params.genome.matches("hg38(.*)")| params.genome.matches("hg19(.*)")){
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome.matches("mm10")){
        somalier_analysis_mouse(som_in)
        somalier_analysis_out=somalier_analysis_mouse.out.collect()
    }

    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    samtools_flagstats_out=samtools_flagstats.out.collect()
    mosdepth_out=mosdepth.out.collect()
    bcftools_stats_out=bcftools_stats.out.collect()
    gatk_varianteval_out=gatk_varianteval.out.collect()
    snpeff_out=snpeff.out.collect()
    vcftools_out=vcftools.out
    collectvariantcallmetrics_out=collectvariantcallmetrics.out

    conall=concat(qualimap_out,mosdepth_out,
        samtools_flagstats_out,bcftools_stats_out,
        gatk_varianteval_out,snpeff_out,vcftools_out,collectvariantcallmetrics_out,
        somalier_analysis_out).flatten().toList()
    multiqc(conall)
}


//QC NOGL-BAMs 
workflow QC_NOGL_BAM {
    take:
        bams

    main:
    //BQSR BAMs 
    fastqc(bams)
    samtools_flagstats(bams)
    qualimap_bamqc(bams)
    mosdepth(bams)

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
    mosdepth_out=mosdepth.out.collect()
    samtools_flagstats_out=samtools_flagstats.out.collect()

    conall=fclane_out.concat(qualimap_out,
        samtools_flagstats_out,mosdepth_out, 
        somalier_analysis_out).flatten().toList()
    
    multiqc(conall)
}

//Variant Calling from BAM only
workflow INPUT_BAM {

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
        //Check if Index is .bai or .bam.bai
        bambai = params.bam_input + ".bai"
        baionly = bambai.replace(".bam", "")
        bamcheck1 = file(bambai)
        bamcheck2 = file(baionly)

        if (bamcheck1.size()>0){
            baminputonly=Channel.fromPath(params.bam_input)
           .map{it-> tuple(it.simpleName,it,file("${it}.bai"))}
        }else if (bamcheck2.size()>0){
            bai=Channel.from(bamcheck2).map{it -> tuple(it.simpleName,it)}
            baminputonly=Channel.fromPath(params.bam_input)
           .map{it-> tuple(it.simpleName,it)}
           .join(bai)
        }else if (bamcheck1.size==0 && bamcheck2.size==0){
            println "Missing BAM Index"
        }
    }else if(params.bam_file_input) {
        baminputonly=Channel.fromPath(params.bam_file_input)
                        .splitCsv(header: false, sep: "\t", strip:true)
                        .map{ sample,bam,bai  ->
                        tuple(sample, file(bam),file(bai))
                                  } 
    }
    if (params.intervals){
        intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
    }else{
        intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
    }
    matchbed(intervalbedin) | splitinterval
    

    if (params.indelrealign){ 
        bqsrs = baminputonly | indelrealign | combine(splitinterval.out.flatten()) 
            | bqsr_ir 
            | groupTuple 
            | map { samplename,beds -> 
            tuple( samplename, beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )} 
            | gatherbqsr

        baminput2=baminputonly.combine(bqsrs,by:0) 
            |applybqsr

        bamwithsample=baminput2.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(baminputonly,by:0).map{it.swap(3,0)}  
            | view()      
        }else{
        bamwithsample=baminputonly.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(baminputonly,by:0).map{it.swap(3,0)} 
            | view()      
        }

        bambyinterval_norm=bamwithsample
            | map {tum,tubam,tbai,norm,norbam,norbai -> tuple(norm,norbam,norbai)} 
        bambyinterval_tum=bamwithsample
            | map {tum,tubam,tbai,norm,norbam,norbai -> tuple(tum,tubam,tbai)} 
        bambyinterval=bambyinterval_tum | concat(bambyinterval_norm) | unique
            | combine(splitinterval.out.flatten()) 

    emit:
        bamwithsample
        bambyinterval
        splitout=splitinterval.out
        sample_sheet

}
