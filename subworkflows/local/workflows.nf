include {fc_lane; fastq_screen;kraken;qualimap_bamqc;fastqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;
    somalier_extract;somalier_analysis_human;somalier_analysis_mouse;
    multiqc} from  '../../modules/local/qc.nf'

include {fastp; bwamem2; indelrealign; bqsr_ir;
    bqsr; gatherbqsr; applybqsr; samtoolsindex} from  '../../modules/local/trim_align.nf'

include {deepvariant_step1; deepvariant_step2; deepvariant_step3;
    deepvariant_combined;glnexus} from '../../modules/local/germline.nf'

include {pileup_paired_t; pileup_paired_n; 
    mutect2; mutect2filter; 
    contamination_paired; learnreadorientationmodel;mergemut2stats;
    strelka_tn; combineVariants_strelka;
    varscan_tn; vardict_tn; lofreq_tn; muse_tn;
    octopus_tn; bcftools_index_octopus; bcftools_index_octopus as bcftools_index_octopus_tonly; octopus_convertvcf; 
    combineVariants as combineVariants_vardict; combineVariants as combineVariants_vardict_tonly;
    combineVariants as combineVariants_varscan; combineVariants as combineVariants_varscan_tonly;
    combineVariants_alternative as combineVariants_lofreq; combineVariants as combineVariants_muse;
    combineVariants_alternative as combineVariants_octopus; combineVariants_alternative as combineVariants_octopus_tonly;
    annotvep_tn as annotvep_tn_mut2; annotvep_tn as annotvep_tn_strelka;
    annotvep_tn as annotvep_tn_varscan; annotvep_tn as annotvep_tn_vardict; annotvep_tn as annotvep_tn_octopus;
    annotvep_tn as annotvep_tn_lofreq; annotvep_tn as annotvep_tn_muse;
    annotvep_tn as annotvep_tn_combined;
    combinemafs_tn; somaticcombine} from '../../modules/local/variant_calling.nf'

include {mutect2_t_tonly; mutect2filter_tonly;
    varscan_tonly; vardict_tonly; octopus_tonly;
    contamination_tumoronly;
    learnreadorientationmodel_tonly;
    mergemut2stats_tonly; octopus_convertvcf_tonly;
    annotvep_tonly as annotvep_tonly_varscan; annotvep_tonly as annotvep_tonly_vardict;
    annotvep_tonly as annotvep_tonly_mut2; annotvep_tonly as annotvep_tonly_octopus;
    annotvep_tonly as annotvep_tonly_combined;
    combinemafs_tonly;somaticcombine_tonly} from '../../modules/local/variant_calling_tonly.nf'

include {svaba_somatic; manta_somatic;
    survivor_sv; gunzip;
    annotsv_tn as annotsv_survivor_tn
    annotsv_tn as annotsv_svaba;annotsv_tn as annotsv_manta} from '../../modules/local/structural_variant.nf'

include {amber_tn; cobalt_tn; purple; purple_novc;
    sequenza; seqz_sequenza_bychr; freec; freec_paired } from '../../modules/local/copynumber.nf'

include {splitinterval} from '../../modules/local/splitbed.nf'



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
        sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true).view()
                       .ifEmpty { "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t", strip:true)
                       .map { row -> tuple(
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
    fastp(fastqinput)

    intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
    splitinterval(intervalbedin)

    bwamem2(fastp.out)
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

    //sample_sheet.view()
    bamwithsample=applybqsr.out.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(applybqsr.out,by:0).map{it.swap(3,0)}

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

workflow VC {
    take:
    //Input is the BAMby interval
        bamwithsample
        splitout
        sample_sheet

    main:
    //Create Pairing for TN (in case of dups)
    sample_sheet_paired=sample_sheet|map{tu,no -> tuple ("${tu}_vs_${no}",tu, no)} |view()
    bambyinterval=bamwithsample.combine(splitout.flatten())

    bambyinterval 
        | multiMap {tumorname,tumor,tumorbai,normalname,normalbam,normalbai,bed -> 
        t1: tuple(tumorname,tumor,tumorbai,bed)
        n1: tuple(normalname,normalbam,normalbai,bed)
            }
        | set{bambyinterval_tonly}
        
        bambyinterval_t=bambyinterval_tonly.t1 |
            concat(bambyinterval_tonly.n1) |unique() 

    //Prep Pileups
    call_list = params.callers.split(',') as List

    vc_all=Channel.empty()
    vc_tonly=Channel.empty()

    //Common for Mutect2/Varscan
    if ("mutect2" in call_list | "varscan" in call_list){ 
        pileup_paired_t(bambyinterval) 
        pileup_paired_n(bambyinterval) 

        pileup_paired_t.out.groupTuple(by:[0,1]) 
            | multiMap { samplename, normalname, pileups -> 
                tout: tuple( samplename, normalname,
                    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tumor.pileup.table/)[0][1].toInteger() } )
                tonly: tuple( samplename,
                    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tumor.pileup.table/)[0][1].toInteger() } )
                    }
            | set{pileup_paired_tout}
        
        pileup_paired_n.out.groupTuple(by:[0,1])
            | multiMap { samplename, normalname, pileups-> 
                nout: tuple (samplename,normalname,
                    pileups.toSorted{ it -> (it.name =~ /${normalname}_(.*?).normal.pileup.table/)[0][1].toInteger() } )
                nonly: tuple (normalname,
                    pileups.toSorted{ it -> (it.name =~ /${normalname}_(.*?).normal.pileup.table/)[0][1].toInteger() } )
                    }
                | set{pileup_paired_nout}

        pileup_paired_match=pileup_paired_tout.tout.join(pileup_paired_nout.nout,by:[0,1])
        contamination_paired(pileup_paired_match)

        pileup_all=pileup_paired_tout.tonly.concat(pileup_paired_nout.nonly) 
        contamination_tumoronly(pileup_all) 
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

    //Mutect2 Tumor Only
        mutect2_t_tonly(bambyinterval_t)
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

        mutect2_in_tonly=mut2tonlyout.allmut2tonly
        | join(mergemut2stats_tonly.out)
        | join(learnreadorientationmodel_tonly.out)
        | join(contamination_tumoronly.out)
        | mutect2filter_tonly
        | join(sample_sheet)
        | map{tumor,markedvcf,markedindex,normvcf,normindex,stats,normal -> tuple(tumor,"mutect2_tonly",normvcf,normindex)}
        annotvep_tonly_mut2(mutect2_in_tonly)


        vc_all=vc_all|concat(mutect2_in)
        vc_tonly=vc_tonly|concat(mutect2_in_tonly) 

    }

    if ("strelka" in call_list){
        //Strelka TN
        strelka_in=strelka_tn(bambyinterval) | groupTuple(by:[0,1])
        | map { tumor,normal,vcfs,vcfindex,indels,indelindex -> tuple("${tumor}_vs_${normal}",
        vcfs.toSorted{ it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).somatic.snvs.vcf.gz/)[0][1].toInteger() },vcfindex,
        indels.toSorted{ it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).somatic.indels.vcf.gz/)[0][1].toInteger() } ,indelindex)}
        | combineVariants_strelka |  join(sample_sheet_paired)
        | map{sample,markedvcf,markedindex,finalvcf,finalindex,tumor,normal -> tuple(tumor,normal,"strelka",finalvcf,finalindex)}
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

        //VarDict TOnly
        vardict_in_tonly=vardict_tonly(bambyinterval_t) 
        | groupTuple()
        | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.vardict.vcf/)[0][1].toInteger()},"vardict_tonly")}
        | combineVariants_vardict_tonly | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex,normal ->tuple(tumor,"vardict_tonly",normvcf,normindex)}
        annotvep_tonly_vardict(vardict_in_tonly)

        vc_all=vc_all|concat(vardict_in)
        vc_tonly=vc_tonly|concat(vardict_in_tonly) 
    }

    if ("varscan" in call_list){
        //VarScan TN
        varscan_in=bambyinterval.combine(contamination_paired.out,by:0) 
        | varscan_tn | groupTuple(by:[0,1]) 
        | map{tumor,normal,vcf-> tuple("${tumor}_vs_${normal}",vcf.toSorted{it -> (it.name =~ /${tumor}_vs_${normal}_(.*?).varscan.vcf.gz/)[0][1].toInteger()},"varscan")}
        | combineVariants_varscan | join(sample_sheet_paired)
        | map{sample,marked,markedindex,normvcf,normindex,tumor,normal ->tuple(tumor,normal,"varscan",normvcf,normindex)}
        annotvep_tn_varscan(varscan_in)

        //VarScan TOnly
        varscan_in_tonly=bambyinterval_t.combine(contamination_tumoronly.out,by:0) 
        | varscan_tonly  | groupTuple 
        | map{tumor,vcf-> tuple(tumor,vcf.toSorted{it -> (it.name =~ /${tumor}_(.*?).tonly.varscan.vcf.gz/)[0][1].toInteger()},"varscan_tonly")}  
        | combineVariants_varscan_tonly 
        | join(sample_sheet)
        | map{tumor,marked,markedindex,normvcf,normindex,normal ->tuple(tumor,"varscan_tonly",normvcf,normindex)}
        annotvep_tonly_varscan(varscan_in_tonly)

        vc_all=vc_all|concat(varscan_in)
        vc_tonly=vc_tonly|concat(varscan_in_tonly) 
    }

    //SAGE TN
    if ("sage" in params.callist){
        sage_in=sage_tn(bambyinterval) | groupTuple(by:[0,1])
            | map{tu,no,vcf,vcfindex-> tuple("${tu}_vs_${no}",vcf.toSorted{it -> (it.name =~ /${tu}_vs_${no}_(.*?).sage.vcf.gz/)[0][1].toInteger()},vcfindex,"lofreq")}
            | combineVariants_sage | join(sample_sheet_paired)
            | map{sample,marked,markedindex,normvcf,normindex,tumor,normal->tuple(tumor,normal,"sage",normvcf,normindex)}
        annotvep_tn_sage(sage_in)

        sage_in_tonly=sage_tonly(bambyinterval_t)
            | groupTuple() 
            | map{samplename,vcf,vcfindex->tuple(samplename,vcf.toSorted{it->(it.name =~ /${samplename}_(.*).tonly.sage.vcf.gz/)[0][1].toInteger()},vcfindex,"octopus_tonly")} 
            | combineVariants_sage_tonly
            | join(sample_sheet) 
            | map{tumor,marked,markedindex,normvcf,normindex,normal ->tuple(tumor,"sage_tonly",normvcf,normindex)}
        annotvep_tonly_octopus(sage_in_tonly)
        sage_in_tonly_sc=sage_in_tonly 
            | map{tumor,vcf,vcfindex ->tuple(tumor,"sage_tonly",vcf,vcfindex)} 


        vc_all=vc_all|concat(sage_in)
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
        octopus_in_sc = octopus_in | octopus_convertvcf 
            |  map{tumor,normal,vcf,vcfindex ->tuple(tumor,normal,"octopus",vcf,vcfindex)} 

    //Octopus TOnly
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

        vc_all=vc_all|concat(octopus_in_sc)
        vc_tonly=vc_tonly|concat(octopus_in_tonly) 
    }


    //Combine All Variants Using VCF -> Annotate
    if (call_list.size()>1){
        vc_all | groupTuple(by:[0,1])
            | somaticcombine
            | map{tumor,normal,vcf,index ->tuple(tumor,normal,"combined",vcf,index)}
            | annotvep_tn_combined
    }

    
    if (call_list.size()>1){
        vc_tonly 
            | somaticcombine_tonly
            | map{tumor,vcf,index ->tuple(tumor,normal,"combined_tonly",vcf,index)}
            | annotvep_tn_combined
    }
    
    
    //Implement PCGR Annotator/CivIC Next
        if ("octopus" in call_list){
           somaticcall_input=octopus_in_sc
        }else if("mutect2" in call_list){
            somaticcall_input=mutect2_in
        }

    emit:
        somaticcall_input
   
}


workflow SV {
    take:
        bamwithsample

    main:
        //Svaba
        svaba_out=svaba_somatic(bamwithsample)
        .map{ tumor,bps,contigs,discord,alignents,gindel,gsv,so_indel,so_sv,unfil_gindel,unfil_gsv,unfil_so_indel,unfil_sv,log ->
            tuple(tumor,so_sv,"svaba")}
        annotsv_svaba(svaba_out).ifEmpty("Empty SV input--No SV annotated")

        //Manta
        manta_out=manta_somatic(bamwithsample)
            .map{tumor,gsv,so_sv,unfil_sv,unfil_indel ->
            tuple(tumor,so_sv,"manta")}
        annotsv_manta(manta_out).ifEmpty("Empty SV input--No SV annotated")

        //Delly-WIP

        //Survivor
        gunzip(manta_out).concat(svaba_out).groupTuple()
            | survivor_sv | annotsv_survivor_tn | ifEmpty("Empty SV input--No SV annotated")

}

workflow CNVmouse {
    take:
        bamwithsample
        
    main:
        cnvcall_list = params.cnvcallers.split(',') as List

        if ("sequenza" in cnvcall_list){

        //Sequenza (Preferred for Paired)
        chrs=Channel.fromList(params.genomes[params.genome].chromosomes)
        seqzin=bamwithsample.map{tname,tumor,tbai,nname,norm,nbai->
            tuple("${tname}_${nname}",tname,tumor,tbai,nname,norm,nbai)}
        seqzin.combine(chrs) | seqz_sequenza_bychr
        seqz_sequenza_bychr.out.groupTuple()
            .map{pair, seqz -> tuple(pair, seqz.sort{it.name})}
            | sequenza
        }

        if ("freec" in cnvcall_list){
        //FREEC Paired Mode
        bamwithsample | freec_paired

        //FREEC Unpaired Mode
        bamwithsample 
            | map{tname,tumor,tbai,nname,norm,nbai->tuple(tname,tumor,tbai)}
            | freec
        }
}

workflow CNVhuman {
    take:
        bamwithsample
        somaticcall_input

    main:  
        cnvcall_list = params.cnvcallers.split(',') as List

        if ("purple" in cnvcall_list){
            //Purple
            bamwithsample | amber_tn
            bamwithsample | cobalt_tn
            purplein=amber_tn.out.join(cobalt_tn.out)
            purplein.join(somaticcall_input)|
            map{t1,amber,cobalt,n1,vc,vcf,vcfindex -> tuple(t1,n1,amber,cobalt,vcf,vcfindex)}
                | purple
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

}


workflow CNVhuman_novc {
    take:
        bamwithsample

    main:  
        cnvcall_list = params.cnvcallers.split(',') as List

        if ("purple" in cnvcall_list){
            //Purple
            bamwithsample | amber_tn
            bamwithsample | cobalt_tn
            purplein=amber_tn.out |join(cobalt_tn.out) 
            purplein | map{t1,amber,cobalt,n1 -> tuple(t1,n1,amber,cobalt)}
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

    //Somalier
    somalier_extract(applybqsr)
    som_in=somalier_extract.out.collect()

    //Prep for MultiQC input
    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect()

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    fastqc_out=fastqc.out.map{samplename,html,zip->tuple(html,zip)}.collect()

    samtools_flagstats_out=samtools_flagstats.out.collect()

    if(params.genome=="hg38"){
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome=="mm10"){
        somalier_analysis_mouse(som_in)
        somalier_analysis_out=somalier_analysis_mouse.out.collect()
    }

    conall=fclane_out.concat(fqs_out,kraken_out,qualimap_out,samtools_flagstats_out,
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
    if(params.genome=="hg38"){
        somalier_analysis_human(som_in)
        somalier_analysis_out=somalier_analysis_human.out.collect()
    }
    else if(params.genome=="mm10"){
        somalier_analysis_mouse(som_in)
        somalier_analysis_out=somalier_analysis_mouse.out.collect()
    }

    fclane_out=fc_lane.out.map{samplename,info->info}.collect()
    fqs_out=fastq_screen.out.collect()

    kraken_out=kraken.out.map{samplename,taxa,krona -> tuple(taxa,krona)}.collect()
    qualimap_out=qualimap_bamqc.out.map{genome,rep->tuple(genome,rep)}.collect()
    fastqc_out=fastqc.out.map{samplename,html,zip->tuple(html,zip)}.collect()
    samtools_flagstats_out=samtools_flagstats.out.collect()
    bcftools_stats_out= bcftools_stats.out.collect()
    gatk_varianteval_out= gatk_varianteval.out.collect()
    snpeff_out=snpeff.out.collect()
    vcftools_out=vcftools.out
    collectvariantcallmetrics_out=collectvariantcallmetrics.out

    conall=fclane_out.concat(fqs_out,kraken_out,qualimap_out,samtools_flagstats_out,bcftools_stats_out,
    gatk_varianteval_out,snpeff_out,vcftools_out,collectvariantcallmetrics_out,somalier_analysis_out).flatten().toList()
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
        bambai=params.bam_input +".bai"
        baionly = bambai.replace(".bam", "")
        bamcheck1=file(bambai)
        bamcheck2=file(baionly)

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
    intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')
    splitinterval(intervalbedin)

    if (params.indelrealign){ 
        bqsrs= baminputonly | indelrealign | combine(splitinterval.out.flatten()) 
            | bqsr_ir 
            | groupTuple 
            | map { samplename,beds -> 
            tuple( samplename, beds.toSorted{ it -> (it.name =~ /${samplename}_(.*?).recal_data.grp/)[0][1].toInteger() } )} 
            | gatherbqsr

        baminput2=baminputonly.combine(bqsrs,by:0) 
            |applybqsr

        bamwithsample=baminput2.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(baminputonly,by:0).map{it.swap(3,0)} 

    } else {
        bamwithsample=baminputonly.combine(sample_sheet,by:0).map{it.swap(3,0)}.combine(baminputonly,by:0).map{it.swap(3,0)}
    }

    emit:
        bamwithsample
        splitout=splitinterval.out
        sample_sheet

}
