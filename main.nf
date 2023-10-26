#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )


//SUB WORKFLOWS
PIPE_ALIGN=params.PIPE_ALIGN

PIPE_VC=params.PIPE_VC
PIPE_SV=params.PIPE_SV
PIPE_CNV=params.PIPE_CNV
PIPE_QC=params.PIPE_QC
PIPE_GERMLINE=params.PIPE_GERMLINE

PIPE_BAMVC=params.PIPE_BAMVC

PIPE_TONLY_ALIGN=params.PIPE_TONLY_ALIGN
PIPE_TONLY_VC=params.PIPE_TONLY_VC
PIPE_TONLY_BAMVC=params.PIPE_TONLY_BAMVC
PIPE_TONLY_QC=params.PIPE_TONLY_QC


include {INPUT; ALIGN; GERMLINE;
    VC; INPUT_BAMVC; SV; CNV;
    QC} from "./workflow/modules/workflows.nf"


include {INPUT_TONLY; ALIGN_TONLY;
    VC_TONLY; INPUT_TONLY_BAMVC; QC_TONLY} from "./workflow/modules/workflows_tonly.nf"


log.info """\

         L O G A N   P I P E L I N E    
         =============================
         genome: ${params.genome}
         outdir: ${params.output}
         Samplesheet: ${params.sample_sheet}
         Samples: ${params.fastq_input} ${params.file_input} ${params.bam_input}
         NF version   : $nextflow.version

         """
         .stripIndent()


//Final Workflow
workflow {


    if (PIPE_ALIGN){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
    } 


    //Germline
    if (PIPE_GERMLINE){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        GERMLINE(ALIGN.out.bambyinterval)
    }

    //Tumor-Normal Pipelines
    if (PIPE_VC){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        VC(ALIGN.out.bamwithsample,ALIGN.out.splitout,ALIGN.out.sample_sheet)
    }
    if (PIPE_QC){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        GERMLINE(ALIGN.out.bambyinterval)
        QC(ALIGN.out.fastqin,ALIGN.out.fastpout,ALIGN.out.bqsrout,GERMLINE.out.glnexusout,GERMLINE.out.bcfout)

    }  
    if (PIPE_SV){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        SV(ALIGN.out.bamwithsample)
    }  
    if (PIPE_CNV){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        CNV(ALIGN.out.bamwithsample)
    }  
    if (PIPE_BAMVC){
        INPUT_BAMVC()
        VC(INPUT_BAMVC.out.bamwithsample,INPUT_BAMVC.out.splitout,INPUT_BAMVC.out.sample_sheet)
    }  

    ///Tumor Only Pipelines
    if (PIPE_TONLY_ALIGN){
        INPUT_TONLY()
        ALIGN_TONLY(INPUT_TONLY.out.fastqinput,INPUT_TONLY.out.sample_sheet)
    }

    if (PIPE_TONLY_VC){
        INPUT_TONLY()
        ALIGN_TONLY(INPUT_TONLY.out.fastqinput,INPUT_TONLY.out.sample_sheet)
        VC_TONLY(ALIGN_TONLY.out.bamwithsample,ALIGN_TONLY.out.splitout,ALIGN_TONLY.out.sample_sheet)
    }    
    if (PIPE_TONLY_QC){
        INPUT_TONLY()
        ALIGN_TONLY(INPUT_TONLY.out.fastqinput,INPUT_TONLY.out.sample_sheet)
        QC_TONLY(ALIGN_TONLY.out.fastqin,ALIGN_TONLY.out.fastpout,ALIGN_TONLY.out.bqsrout)

    }  

    //Variant Calling from BAM only/Tumor Only
    if (PIPE_TONLY_BAMVC){
        INPUT_TONLY_BAMVC()
        VC_TONLY(ALIGN_TONLY.out.bamwithsample,ALIGN_TONLY.out.splitout,ALIGN_TONLY.out.sample_sheet)
    }  
}

    


