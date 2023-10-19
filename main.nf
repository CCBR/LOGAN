#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )


//SUB WORKFLOWS to SPLIT
PIPE_ALIGN=params.PIPE_ALIGN
PIPE_GERMLINE=params.PIPE_GERMLINE
PIPE_VC=params.PIPE_VC
PIPE_SV=params.PIPE_SV
PIPE_QC=params.PIPE_QC
PIPE_BAMVC=params.PIPE_BAMVC
PIPE_TONLY_ALIGN=params.PIPE_TONLY_ALIGN
PIPE_TONLY_VC=params.PIPE_TONLY_VC
PIPE_TONLY_BAMVC=params.PIPE_TONLY_BAMVC
PIPE_TONLY_QC=params.PIPE_TONLY_QC



include {INPUT_PIPE;TRIM_ALIGN_PIPE;
    GERMLINE_PIPE;VARIANTCALL_PIPE;INPUT_BAMVC_PIPE;SV_PIPE;
    QC_PIPE} from "./workflow/modules/workflows.nf"

include {INPUT_TONLY_PIPE;TRIM_ALIGN_TONLY_PIPE;
    VARIANT_TONLY_PIPE;INPUT_TONLY_BAMVC_PIPE;QC_TONLY_PIPE} from "./workflow/modules/workflows_tonly.nf"


log.info """\
         W G S S E E K   P I P E L I N E    
         =============================
         genome: ${params.genome}
         outdir: ${params.output}
         Samplesheet: ${params.sample_sheet}
         Samples: ${params.fastq_input} ${params.file_input} ${params.bam_input}
         """
         .stripIndent()


//Final Workflow
workflow {

    if (PIPE_ALIGN){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
    } 

    //GermlineVC 
    if (PIPE_GERMLINE){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
        GERMLINE_PIPE(TRIM_ALIGN_PIPE.out.bambyinterval)
    }

    //Tumor-Normal Pipelines
    if (PIPE_VC){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
        VARIANTCALL_PIPE(TRIM_ALIGN_PIPE.out.bamwithsample,TRIM_ALIGN_PIPE.out.splitout,TRIM_ALIGN_PIPE.out.sample_sheet)
    }
    if (PIPE_QC){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
        GERMLINE_PIPE(TRIM_ALIGN_PIPE.out.bambyinterval)
        QC_PIPE(TRIM_ALIGN_PIPE.out.fastqin,TRIM_ALIGN_PIPE.out.fastpout,TRIM_ALIGN_PIPE.out.bwamem2out,GERMLINE_PIPE.out.glnexusout,GERMLINE_PIPE.out.bcfout)

    }  
    if (PIPE_SV){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
        SV_PIPE(TRIM_ALIGN_PIPE.out.bamwithsample)

    }  
    if (PIPE_BAMVC){
        INPUT_BAMVC_PIPE()
        VARIANTCALL_PIPE(INPUT_BAMVC_PIPE.out.bamwithsample,INPUT_BAMVC_PIPE.out.splitout,INPUT_BAMVC_PIPE.out.sample_sheet)
    }  


    ///Tumor Only Pipelines
    if (PIPE_TONLY_ALIGN){
        INPUT_TONLY_PIPE()
        TRIM_ALIGN_TONLY_PIPE(INPUT_TONLY_PIPE.out.fastqinput,INPUT_TONLY_PIPE.out.sample_sheet)
    }
    if (PIPE_TONLY_VC){
        INPUT_TONLY_PIPE()
        TRIM_ALIGN_TONLY_PIPE(INPUT_TONLY_PIPE.out.fastqinput,INPUT_TONLY_PIPE.out.sample_sheet)
        VARIANT_TONLY_PIPE(TRIM_ALIGN_TONLY_PIPE.out.bamwithsample,TRIM_ALIGN_TONLY_PIPE.out.splitout,TRIM_ALIGN_TONLY_PIPE.out.sample_sheet)
    }    
    if (PIPE_TONLY_QC){
        INPUT_TONLY_PIPE()
        TRIM_ALIGN_TONLY_PIPE(INPUT_TONLY_PIPE.out.fastqinput,INPUT_TONLY_PIPE.out.sample_sheet)
        QC_TONLY_PIPE(TRIM_ALIGN_TONLY_PIPE.out.fastqin,TRIM_ALIGN_TONLY_PIPE.out.fastpout,TRIM_ALIGN_TONLY_PIPE.out.bqsrout)

    }  

    //Variant Calling from BAM only/Tumor Only
    if (PIPE_TONLY_BAMVC){
        INPUT_TONLY_BAMVC_PIPE()
        VARIANT_TONLY_PIPE(INPUT_TONLY_BAMVC_PIPE.out.bamwithsample,INPUT_TONLY_BAMVC_PIPE.out.splitout,INPUT_TONLY_BAMVC_PIPE.out.sample_sheet)
    }  
}
    


