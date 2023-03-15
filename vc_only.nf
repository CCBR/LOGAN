#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

publishDir=file(params.output)
intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')

//Final Workflow
include {fc_lane; fastq_screen;kraken;qualimap_bamqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;
    somalier_extract;somalier_analysis;multiqc} from  './workflow/modules/qc.nf'
include {deepvariant_step1;deepvariant_step2;deepvariant_step3;deepvariant_combined;glnexus} from './workflow/modules/germline.nf'
include {fastp; bwamem2; indelrealign; bqsr; 
    gatherbqsr; applybqsr; samtoolsindex} from './workflow/modules/trim_align.nf'
include {mutect2; mutect2_t_tonly; mutect2filter; mutect2filter_tonly; 
    pileup_paired_t; pileup_paired_n; 
    contamination_paired; contamination_tumoronly;
    learnreadorientationmodel; learnreadorientationmodel_tonly; 
    mergemut2stats; mergemut2stats_tonly;
    annotvep_tn; annotvep_tonly} from './workflow/modules/variant_calling.nf'
include {splitinterval} from "./workflow/modules/splitbed.nf"


workflow{
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
   bamby1=Channel.fromPath(params.bam_input)
   .map{it-> tuple(it.simpleName,it,file("${it}.bai"))}
    
    splitinterval(intervalbedin)
    
   bamwithsample=bamby1.join(sample_sheet).map{it.swap(3,0)}.join(bamby1).map{it.swap(3,0)}.view()

   bambyinterval=bamwithsample.combine(splitinterval.out.flatten())

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

    //PCGR Annotator/CivIC?
}