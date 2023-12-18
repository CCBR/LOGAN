#!/bin/bash
#Author: Dr Charles Foster http://github.com/charlesfoster

i_flag=''
g_flag=''
h_flag=''
n_flag=''
o_flag=''

print_usage() {
	printf "Usage: bash add_artificial_genotype.sh -i in.vcf.gz [-g genotype] [-n sample_name] -o out.vcf.gz\n"
	printf "Genotype defaults to 1 if not specified\n"
	printf "Sample name gussed from infile name if not specified\n"
}

if [[ $# -eq 0 ]] ; then
	print_usage
    exit 1
fi

while getopts 'i:g:n:ho:' flag; do
  case "${flag}" in
		i) IN="${OPTARG}" ;;
		n) NAME="${OPTARG}" ;;
		g) GENOTYPE="${OPTARG}" ;;
    h) print_usage
           exit 1 ;;
    o) OUT="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [ ! -f ${IN} ]; then
    printf "\nError: input file not found\n"
		print_usage
		exit 1
fi

if [ -z "${GENOTYPE}" ]
  then
    printf "\nNo genotype specified: setting to 1"
		GENOTYPE=1
fi

if [ -z "${NAME}" ]
  then
		NAME=$(basename ${IN} | cut -f1 -d ".")
    printf "\nNo name specified: guessed it to be ${NAME}"
fi

if [ -z "${OUT}" ]
  then
		OUT=$(echo ${IN} | sed "s/.vcf.gz/_withGT.vcf.gz/")
    printf "\nNo outfile specified: setting to ${OUT}\n"
fi

gunzip -kc ${IN} | \
sed -e '6i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
-e "s|FILTER\tINFO|FILTER\tINFO\tFORMAT\t${NAME}|g" | \
awk -F'\t' -v genotype=${GENOTYPE} -v OFS="\t" '/^[^#]/{ $9 = "GT"; $10 = genotype }1' | \
bgzip -c > ${OUT}
tabix -p vcf ${OUT}
printf "VCF with artificial genotype written to ${OUT}\n"
exit 0
