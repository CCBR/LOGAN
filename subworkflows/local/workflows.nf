//All Worksflows in One Place
intervalbedin = Channel.fromPath(params.genomes[params.genome].intervals,checkIfExists: true,type: 'file')


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

include {amber_tn; cobalt_tn; purple;
    sequenza; seqz_sequenza_bychr; freec; freec_paired } from '../../modules/local/copynumber.nf'

include {splitinterval;split_byline} from '../../modules/local/splitbed.nf'



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
    sample_sheet_paired=sample_sheet|map{tu,no -> tuple ("${tu}_vs_${no}",tu, no)}
    bambyinterval=bamwithsample.combine(splitout.flatten())

    //Octopus TN
    octopus_in=octopus_tn(bambyinterval) | bcftools_index_octopus
        | groupTuple()
        | map{samplename,vcf,vcfindex-> tuple(samplename,vcf.toSorted{it->(it.name =~ /${samplename}_(.*).octopus.vcf.gz/)[0][1].toInteger()},vcfindex,"octopus")}
        | combineVariants_octopus
        | map{samplename,marked,markedindex,normvcf,normindex ->
            tuple(samplename.split('_vs_')[0],samplename.split('_vs_')[1],"octopus",normvcf,normindex)}
    annotvep_tn_octopus(octopus_in)
    octopus_in_sc = octopus_in | octopus_convertvcf 
        |  map{tumor,normal,vcf,vcfindex ->tuple(tumor,normal,"octopus",vcf,vcfindex)} 


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
        //Sequenza (Preferred for Paired)
        chrs=Channel.fromList(params.genomes[params.genome].chromosomes)
        seqzin=bamwithsample.map{tname,tumor,tbai,nname,norm,nbai->
            tuple("${tname}_${nname}",tname,tumor,tbai,nname,norm,nbai)}
        seqzin.combine(chrs) | seqz_sequenza_bychr
        seqz_sequenza_bychr.out.groupTuple()
            .map{pair, seqz -> tuple(pair, seqz.sort{it.name})}
            | sequenza

        //FREEC Paired Mode
        bamwithsample | freec_paired

        //FREEC Unpaired Mode
        bamwithsample 
            | map{tname,tumor,tbai,nname,norm,nbai->tuple(tname,tumor,tbai)}
            | freec

}

workflow CNVhuman {
    take:
        bamwithsample
        somaticcall_input

    main:
        //Sequenza
        chrs=Channel.fromList(params.genomes[params.genome].chromosomes)
        seqzin=bamwithsample.map{tname,tumor,tbai,nname,norm,nbai->
            tuple("${tname}_${nname}",tname,tumor,tbai,nname,norm,nbai)}
        seqzin.combine(chrs) | seqz_sequenza_bychr
        seqz_sequenza_bychr.out.groupTuple()
            .map{pair, seqz -> tuple(pair, seqz.sort{it.name})}
            | sequenza

        //Purple
        bamwithsample | amber_tn
        bamwithsample | cobalt_tn
        purplein=amber_tn.out.join(cobalt_tn.out)
        purplein.join(somaticcall_input)|
        map{t1,amber,cobalt,n1,vc,vcf,vcfindex -> tuple(t1,amber,cobalt,vcf,vcfindex)}
            | purple

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

    split_byline(intervalbedin)

    if (params.indelrealign){ 
        bqsrs= baminputonly | indelrealign | combine(split_byline.out.flatten()) 
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
        splitout=split_byline.out
        sample_sheet

}
