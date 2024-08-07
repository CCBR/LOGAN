#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )


log.info """\
         L O G A N     P I P E L I N E
         =============================
         genome: ${params.genome}
         outdir: ${params.outdir}
         Sample Sheet: ${params.sample_sheet}
         Samples: ${params.fastq_input} ${params.fastq_file_input} ${params.bam_input} ${params.bam_file_input}
         """
         .stripIndent()



include {DETERMINEBAM; INPUT; INPUT_BAM; ALIGN; GL;
    VC; SV; CNVmouse; CNVhuman; CNVhuman_novc;
    QC_GL; QC_NOGL} from "./subworkflows/local/workflows.nf"

include {INPUT_TONLY; INPUT_TONLY_BAM;
    ALIGN_TONLY;
    VC_TONLY; SV_TONLY; CNVmouse_tonly;  CNVhuman_tonly; CNVhuman_novc_tonly;
    QC_TONLY } from "./subworkflows/local/workflows_tonly.nf"


workflow.onComplete {
    if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
        def message = Utils.spooker(workflow)
        if (message) {
            println message
        }
    }
}

//Final Workflow
workflow {
        if ([params.fastq_input,params.fastq_file_input].any() && params.sample_sheet){
        println "Tumor-Normal FASTQ"
        INPUT()
        ALIGN(INPUT.out.fastqinput,INPUT.out.sample_sheet)
    //Germline
        if (params.gl){
           GL(ALIGN.out.bambyinterval)
        }
        //Tumor-Normal VC, SV, CNV
        if (params.vc){
            VC(ALIGN.out.bamwithsample,ALIGN.out.splitout,ALIGN.out.sample_sheet)
        }
        if (params.sv){
            SV(ALIGN.out.bamwithsample)
        }
        if (params.cnv){
            if (params.genome == "mm10"){
                CNVmouse(ALIGN.out.bamwithsample)
            } else if (params.genome== "hg38" |params.genome== "hg19"){
                if (!params.vc){
                    CNVhuman_novc(ALIGN.out.bamwithsample)
                } else {
                    CNVhuman(ALIGN.out.bamwithsample,VC.out.somaticcall_input)
                }
            }
    }
        if (params.qc && params.gl){
            QC_GL(ALIGN.out.fastqin,ALIGN.out.fastpout,ALIGN.out.bqsrout,GL.out.glnexusout,GL.out.bcfout)
        }  else if (params.qc){
            QC_NOGL(ALIGN.out.fastqin,ALIGN.out.fastpout,ALIGN.out.bqsrout)
        }

    }

    //TUMOR-NOMRAL BAM INPUT
    if ([params.bam_input,params.bam_file_input].any() && params.sample_sheet){
        println "Tumor-Normal BAM"
        INPUT_BAM()
        if (params.vc){
            VC(INPUT_BAM.out.bamwithsample,INPUT_BAM.out.splitout,INPUT_BAM.out.sample_sheet)
        }
        if (params.sv){
            SV(INPUT_BAM.out.bamwithsample)
        }
        if (params.cnv){
            if (params.genome == "mm10"){
                CNVmouse(INPUT_BAM.out.bamwithsample)
            } else if (params.genome == "hg38"|params.genome== "hg19"){
                if (!params.vc){
                    CNVhuman_novc(INPUT_BAM.out.bamwithsample)
                }else {
                    CNVhuman(INPUT_BAM.out.bamwithsample,VC.out.somaticcall_input)
                }
            }
        }
    }

    ///Tumor Only Pipelines
    if ([params.fastq_input,params.fastq_file_input].any() && !params.sample_sheet){
        println "Tumor-Only FASTQ"
        INPUT_TONLY()
        ALIGN_TONLY(INPUT_TONLY.out.fastqinput,INPUT_TONLY.out.sample_sheet)
        if (params.vc){
            VC_TONLY(ALIGN_TONLY.out.bamwithsample,ALIGN_TONLY.out.splitout,ALIGN_TONLY.out.sample_sheet)
        }
        if (params.sv){
            SV_TONLY(ALIGN_TONLY.out.bamwithsample)
        }
        if (params.cnv){
            if (params.genome == "mm10"){
                CNVmouse_tonly(ALIGN_TONLY.out.bamwithsample)
            } else if (params.genome== "hg38"|params.genome== "hg19"){
                if (!params.vc){
                    VC_TONLY(ALIGN_TONLY.out.bamwithsample,ALIGN_TONLY.out.splitout,ALIGN_TONLY.out.sample_sheet)
                    CNVhuman_tonly(ALIGN_TONLY.out.bamwithsample,VC_TONLY.out.somaticcall_input)
                } else {
                    CNVhuman_tonly(ALIGN_TONLY.out.bamwithsample,VC_TONLY.out.somaticcall_input)
                }
            }
        }
        if (params.qc){
                QC_TONLY(ALIGN_TONLY.out.fastqin,ALIGN_TONLY.out.fastpout,ALIGN_TONLY.out.bqsrout)
        }
    }

    //Variant Calling from BAM-Tumor Only Mode
    if ([params.bam_input,params.bam_file_input].any() && !params.sample_sheet){
        println "Tumor-Only BAM"
        INPUT_TONLY_BAM()
        if (params.vc){
            VC_TONLY(INPUT_TONLY_BAM.out.bamwithsample,INPUT_TONLY_BAM.out.splitout,INPUT_TONLY_BAM.out.sample_sheet)
        }
        if (params.sv){
            SV_TONLY(INPUT_TONLY_BAM.out.bamwithsample)
        }
        if (params.cnv){
            if (params.genome == "mm10"){
                CNVmouse_tonly(INPUT_TONLY_BAM.out.bamwithsample)
            } else if (params.genome== "hg38" | params.genome== "hg19"){
                if (!params.vc){
                    VC_TONLY(INPUT_TONLY_BAM.out.bamwithsample,INPUT_TONLY_BAM.out.splitout,INPUT_TONLY_BAM.out.sample_sheet)
                    CNVhuman_tonly(INPUT_TONLY_BAM.out.bamwithsample,VC_TONLY.out.somaticcall_input)
                } else {
                    CNVhuman_tonly(INPUT_TONLY_BAM.out.bamwithsample,VC_TONLY.out.somaticcall_input)
                }

        }
        }
    }
}
