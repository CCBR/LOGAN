#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

fastqinput=Channel.fromFilePairs(params.input,checkIfExists: true)
intervalbedin = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')


include {fc_lane; fastq_screen;kraken;qualimap_bamqc;
    samtools_flagstats;vcftools;collectvariantcallmetrics;
    bcftools_stats;gatk_varianteval;
    snpeff;
    somalier_extract;somalier_analysis;multiqc} from  './workflow/modules/qc.nf'
include {fastp;bwamem2;samtoolsindex} from './workflow/modules/trim_align.nf'
include {deepvariant_step1;deepvariant_step2;deepvariant_step3;glnexus} from './workflow/modules/germline.nf'
include {mutect2; mutect2_t_tonly; mutect2filter; mutect2filter_tonly; 
    pileup_paired_tonly; 
    contamination_tumoronly;
    learnreadorientationmodel_tonly; 
    mergemut2stats_tonly;
    annotvep_tonly} from './workflow/modules/variant_calling.nf'
include {splitinterval} from './workflow/modules/splitbed.nf'



workflow {
    fastqinput.view()
    
    if(params.sample_sheet){
    sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: false)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor,
                       row.Normal
                       )
                                  }
    }else{
    sample_sheet=fastqinput.map{samplename,f1 -> tuple (
             samplename)
             }.view()
        
    }

    //Preprocess Files
    fastp(fastqinput)
    splitinterval(intervalbedin)
    bwamem2(fastp.out)
    


    bambyinterval=bwamem2.out.combine(splitinterval.out.flatten())
    
    //GERMLINE CALLING
    deepvariant_step1(bambyinterval) 
    deepvariant_1_sorted=deepvariant_step1.out.groupTuple()
        .map { samplename,tfbeds,gvcfbed -> tuple( samplename, 
        tfbeds.toSorted{ it -> (it.name =~ /${samplename}.tfrecord_(.*?).bed.gz/)[0][1].toInteger() } ,
        gvcfbed.toSorted{ it -> (it.name =~ /${samplename}.gvcf.tfrecord_(.*?).bed.gz/)[0][1].toInteger() } )
        }
    deepvariant_step2(deepvariant_1_sorted) | deepvariant_step3 
    glin=deepvariant_step3.out.map{samplename,vcf,vcf_tbi,gvcf,gvcf_tbi -> gvcf}.collect()
    glnexus(glin)

    //QC Steps
    fc_lane(fastqinput)
    fastq_screen(fastp.out)
    kraken(fastqinput)
    qualimap_bamqc(bwamem2.out)
    samtools_flagstats(bwamem2.out)
    glout=glnexus.out.map{germlinev,germlinenorm,tbi->tuple(germlinenorm,tbi)}
    vcftools(glout)
    collectvariantcallmetrics(glout)
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

    kraken_out=kraken.out.map{samplename,txt,taxa,krona -> tuple(txt,taxa,krona)}.collect()
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




