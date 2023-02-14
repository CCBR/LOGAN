#!/usr/bin/env python3

##https://raw.githubusercontent.com/bioinform/somaticseq/master/somaticseq/utilities/split_Bed_into_equal_regions.py
########Run NEXTFLOW
import argparse,os,time

def parse_args():
    parser = argparse.ArgumentParser(description='Input files')
    parser.add_argument('--input',help='FQ Inputs ')
    parser.add_argument('--output',help="Output Directory")
    parser.add_argument('--sample_sheet',help="Sample sheet")
    parser.add_argument("--profile",default="-profile biowulf",help="Biowulf or Local Run")
    parser.add_argument("--resume",default="-resume",help="Resume previous run?")
    parser.add_argument("--submit",action="store_false",help="Submit to SLURM?")
    args = parser.parse_args()
    return(args)

def main():
    args=parse_args()
    c1="#!/usr/bin/bash"
    c2="module load nextflow/22.04.0"
    c3="module load singularity"
    c4=["nextflow /data/SCLCgenomics/nousome/WGS_new/wgs-seek/wgs-seek.nf",
"-c /data/SCLCgenomics/nousome/WGS_new/wgs-seek/nextflow.config",
"--input","'"+args.input+"'",
"--sample_sheet",args.sample_sheet, 
args.profile,args.resume,
"--output",args.output]
    cmd1=' '.join(c4)
    code=c1+"\n"+c2+"\n"+c3+"\n"+cmd1
    time1=time.strftime("%Y_%m_%d_%H%M%S")
    outswarmmut='nf_'+time1+'.slurm'
    with open(outswarmmut, "a") as outfile:
        outfile.write(code+"\n")
    sbatch_mut="sbatch --cpus-per-task=2 --mem=16g --time 10-00:00:00 --partition norm --output submit_"+time1+".log --error error_"+time1+".log "+outswarmmut 
    print(sbatch_mut)
    if args.submit:
        os.system(sbatch_mut)


if __name__=="__main__":
  main()

