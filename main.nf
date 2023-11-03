#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )


//SUB WORKFLOWS
PIPE_ALIGN=params.PIPE_ALIGN

PIPE_VC=params.PIPE_VC
PIPE_SV=params.PIPE_SV
PIPE_CNV=params.PIPE_CNV

PIPE_QC_GL=params.PIPE_QC_GL
PIPE_QC_NOGL=params.PIPE_QC_NOGL

PIPE_GL=params.PIPE_GL
PIPE_BAMVC=params.PIPE_BAMVC

PIPE_TONLY_ALIGN=params.PIPE_TONLY_ALIGN
PIPE_TONLY_VC=params.PIPE_TONLY_VC
PIPE_TONLY_SV=params.PIPE_TONLY_SV
PIPE_TONLY_CNV=params.PIPE_TONLY_CNV

PIPE_TONLY_BAMVC=params.PIPE_TONLY_BAMVC
PIPE_TONLY_QC=params.PIPE_TONLY_QC

GENOME=params.genome

include {INPUT; ALIGN; GL;
    VC; INPUT_BAMVC; SV; CNVmm10; CNVhg38;
    QC_GL; QC_NOGL} from "./workflow/modules/workflows.nf"


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
    if (PIPE_GL){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        GL(ALIGN.out.bambyinterval)
    }

    //Tumor-Normal Pipelines
    if (PIPE_VC){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        VC(ALIGN.out.bamwithsample,ALIGN.out.splitout,ALIGN.out.sample_sheet)
    }
    if (PIPE_QC_GL){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        GL(ALIGN.out.bambyinterval)
        QC_GL(ALIGN.out.fastqin,ALIGN.out.fastpout,ALIGN.out.bqsrout,GL.out.glnexusout,GL.out.bcfout)
    }  
    if (PIPE_QC_NOGL){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        QC_NOGL(ALIGN.out.fastqin,ALIGN.out.fastpout,ALIGN.out.bqsrout)
    }  
    if (PIPE_SV){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        SV(ALIGN.out.bamwithsample)
    }  
    if (PIPE_CNV && params.genome == "mm10"){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        CNV_mm10(ALIGN.out.bamwithsample)
    }  
  if (PIPE_CNV && params.genome == "hg38"){
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
        CNV_mm10(ALIGN.out.bamwithsample)
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
    if (PIPE_TONLY_SV){
        INPUT_TONLY()
        ALIGN_TONLY(INPUT_TONLY.out.fastqinput,INPUT_TONLY.out.sample_sheet)
        SV_TONLY(ALIGN.out.bamwithsample)
    }   
    if (PIPE_TONLY_CNV){
        INPUT_TONLY()
        ALIGN_TONLY(INPUT_TONLY.out.fastqinput,INPUT_TONLY.out.sample_sheet)
        CNV_TONLY(ALIGN.out.bamwithsample)
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

    


