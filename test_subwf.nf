#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )


//SUB WORKFLOWS to SPLIT
PIPE_TRIM=params.PIPE_TRIM
PIPE_GERMLINE=params.PIPE_GERMLINE
PIPE_VC=params.PIPE_VC
PIPE_QC=params.PIPE_QC
PIPE_BAMVC=params.PIPE_BAMVC

//Final Workflow

include {INPUT_PIPE;TRIM_ALIGN_PIPE;
GERMLINE_PIPE;VARIANTCALL_PIPE;INPUT_BAMVC_PIPE} from "./workflow/modules/workflows.nf"

log.info """\
         W G S S E E K   P I P E L I N E    
         =============================
         genome: ${params.genome}
         outdir: ${params.output}
         Samplesheet: ${params.sample_sheet}

         """
         .stripIndent()



workflow {

    if (PIPE_TRIM){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
    } 
    if (PIPE_GERMLINE){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
        GERMLINE_PIPE(TRIM_ALIGN_PIPE.out.bambyinterval)
    }
    if (PIPE_VC){
        INPUT_PIPE()
        TRIM_ALIGN_PIPE(INPUT_PIPE.out.fastqinput,INPUT_PIPE.out.sample_sheet)
        VARIANTCALL_PIPE(TRIM_ALIGN_PIPE.out.bamwithsample,TRIM_ALIGN_PIPE.out.splitout,TRIM_ALIGN_PIPE.out.sample_sheet)
    }
    if (PIPE_BAMVC){
        INPUT_BAMVC_PIPE()
        VARIANTCALL_PIPE(INPUT_BAMVC_PIPE.out.bamwithsample,INPUT_BAMVC_PIPE.out.splitout,INPUT_BAMVC_PIPE.out.sample_sheet)

    }   
}
    


