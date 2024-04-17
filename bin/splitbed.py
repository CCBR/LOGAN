#!/usr/bin/env python3

import argparse 

def run():
    parser = argparse.ArgumentParser(description='Given an input bed file, this program will output a number of bed files, each will have same number of total base pairs. This routine is used to parallelize SomaticSeq tasks. One limitation, however, is that some regions of the genome have much higher coverage than others. This is the reason some regions run much slower than others.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-infile','--input-file',    type=str, help='Input merged BED file',    required=True,  default=None)

    args = parser.parse_args()
    infile   = args.input_file
    
    return infile

def splitter(infile): 
    file1 = open(infile, 'r')
    Lines = file1.readlines()

    count = 0
    for line in Lines:
        count += 1
        x=line.strip()
        fileout="bedout/" + str(count) + ".bed"
        file1 = open(fileout, 'w')
        file1.writelines(x)
        file1.close()


if __name__ == '__main__':
    infile = run()
    splitter(infile)