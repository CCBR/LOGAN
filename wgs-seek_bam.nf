#!/usr/bin/env nextflow
//Variant Calling with Tumor Only for BAM files 
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

GENOME=file(params.genome)
GENOMEDICT=file(params.genomedict)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"// file(params.gold_indels1) //
SHAPEITINDEL=file(params.shapeitindel) //params.shapeitindel =  "/data/OpenOmics/references/genome-seek/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz" //file(params.gold_indels2) //
KGP=file(params.kgp) ///data/CCBR_Pipeliner/Exome-seek/hg38/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.dbsnp) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //= '/data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz' // /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 


results_dir=params.output


baminput= Channel.fromPath(params.bams, checkIfExists: true, type: 'file')
               .map { row -> tuple(
                        row.simpleName,row
                       )}


baiinput= Channel.fromPath(params.bai, checkIfExists: true, type: 'file')
               .map { row -> tuple(
                        row.simpleName,row
                       )}


baminput=baminput.join(baiinput)
/*
baminput= Channel.fromFilePairs(params.bams, checkIfExists: true,size:-1){ file -> file.name.replaceAll(/.bam|.bai$/,'') }
                .map { sample,bams -> tuple(
                        sample,bams.toSorted{ it -> (it.name=~/(.*?).bam$/ )} 
                 )}.view()

*/


intervalbed = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor,
                       row.Normal
                       )
                                  }
//Final Workflow
workflow {
    //baminput.view()
    
    
    bamwithsample=baminput.join(sample_sheet)
    .map{it.swap(3,0)}
    .join(baminput)
    .map{it.swap(3,0)}

    bambyinterval=bamwithsample.combine(intervalbed)
    
    
    //Paired Mutect2    
    mutect2(bambyinterval)
    pileup_paired_t(bambyinterval)
    pileup_paired_n(bambyinterval)
    

    //pileup_paired_t.out.view()
      
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

    mutect2_tonly(bambyinterval)    
    
    //LOR     
    mut2tout_lor=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    f1r2.toSorted{ it -> (it.name =~ /${samplename}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } 
    )}
    learnreadorientationmodel_tonly(mut2tout_lor)


    //Stats
    mut2tonly_mstats=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    stats.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tonly.mut2.vcf.gz.stats/)[0][1].toInteger() } 
    )}
    mergemut2stats_tonly(mut2out_mstats)


    //Contamination
    contamination_tumoronly(pileup_paired_tout)


    
    //Final TUMOR ONLY FILTER
    allmut2tonly=mutect2_tonly.out.groupTuple()
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
    //inpvep=vcfinput.join(sample_sheet)
    //annotvep(inpvep)
    
    tn_vepin=mutect2filter.out
    .join(sample_sheet)

   annotvep_tn(tn_vepin)
   annotvep_tonly(mutect2filter_tonly.out)
    
}
