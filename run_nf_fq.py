#!/usr/bin/env python3

import argparse,os,time

def parse_args():
    parser = argparse.ArgumentParser(description='Input files')
    parser.add_argument('--input',help='FQ Inputs ')
    parser.add_argument('--output',help="Output Directory")
    parser.add_argument('--sample_sheet',help="Sample sheet")
    parser.add_argument("--profile",default="biowulf",help="Biowulf or Local Run")
    parser.add_argument("--resume",action="store_true",default="True",help="Resume previous run?")
    parser.add_argument("--submit",action="store_false",help="Submit to SLURM?",default="False")
    parser.add_argument("--QC",action="store_false",help="run QC Steps")
    parser.add_argument("--paired",action="store_false",help="paired", default="False")

    args = parser.parse_args()
    return(args)

def main():
    args=parse_args()
    c1="#!/usr/bin/bash"
    c2="module load nextflow"
    c3="module load singularity"
    if args.paired==True:
        wgs_path='run /data/nousomedr/wgs/wgs-seek/wgs-seek_paired.nf'
    else:
        wgs_path='run /data/nousomedr/wgs/wgs-seek/wgs-seek_tumoronly.nf'
    if args.sample_sheet:
        sample_path="--sample_sheet"+args.sample_sheet, 
    else:
        sample_path=""
    if args.profile=="biowulf":
        profile="-profile biowulf"
    else: 
        profile="-profile local"
    if args.resume:
        resume="-resume"
    else:
        resume=""
    c4=["nextflow",wgs_path,
"-c /data/nousomedr/wgs/wgs-seek/nextflow.config",
"--input","'"+args.input+"'",
    profile,resume,sample_path,
"--output","'"+args.output+"'"]
    cmd1=' '.join(c4)
    code=c1+"\n"+c2+"\n"+c3+"\n"+cmd1
    time1=time.strftime("%Y_%m_%d_%H%M%S")
    outswarmmut='nf_'+time1+'.slurm'

    with open(outswarmmut, "a") as outfile:
        outfile.write(code+"\n")
    sbatch_mut="sbatch --cpus-per-task=2 --mem=16g --time 10-00:00:00 --partition norm --output submit_"+time1+".log --error error_"+time1+".log "+outswarmmut 
    if args.submit==True:
        print(sbatch_mut)
        os.system(sbatch_mut)
    else:
        sbatch_out='run_sbatch'+time1+'.sh'
        with open(sbatch_out, "a") as outfile:
            outfile.write(sbatch_mut+"\n")
        print(sbatch_mut)
if __name__=="__main__":
  main()

