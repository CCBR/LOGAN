//References
GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEDICT=file(params.genomes[params.genome].genomedict)


process varscan_tn {
    container "${params.containers.logan}"
    label 'process_somaticcaller'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed),
        path(tumorpileup), path(normalpileup),
        path(tumor_con_table), path(normal_con_table)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_${bed.simpleName}.varscan.vcf.gz")

    shell:
    '''
    tumor_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !{tumor_con_table} | cut -f2 ))" | bc -l)
    normal_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 !{normal_con_table} | cut -f2 ))" | bc -l)
    dual_pileup="samtools mpileup -d 10000 -q 20 -Q 20 -f !{GENOMEREF} -l !{bed} !{normal} !{tumor}"
    varscan_opts="--strand-filter 1 --min-var-freq 0.01 --output-vcf 1 --normal-purity $normal_purity --tumor-purity $tumor_purity"
    varscan_cmd="varscan somatic <($dual_pileup) !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.vcf $varscan_opts --mpileup 1"
    eval "$varscan_cmd"

    varscan somaticFilter !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.vcf.snp \
    --indel-file !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.vcf.indel \
    --p-value 0.01 --min-reads2 5 --min-avg-qual 30 \
    --output-file !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.vcf.snp.temp

    awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.vcf.indel \
        | sed '/^$/d' | bcftools view - -Oz -o !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.indel_temp.vcf.gz
    awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.vcf.snp.temp \
        | sed '/^$/d' | bcftools view - -Oz -o !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.snp_temp.vcf.gz

    gatk SortVcf -I !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.snp_temp.vcf.gz \
    -I !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.indel_temp.vcf.gz \
    -R !{GENOMEREF} -SD !{GENOMEDICT} \
    -O !{tumorname}_vs_!{normalname}_!{bed.simpleName}_temp.varscan.vcf

    printf "NORMAL\t!{normalname}\nTUMOR\t!{tumorname}\n" > sampname

    bcftools reheader -s sampname !{tumorname}_vs_!{normalname}_!{bed.simpleName}_temp.varscan.vcf \
       | bcftools view -Oz -o !{tumorname}_vs_!{normalname}_!{bed.simpleName}.varscan.vcf.gz

    '''

    stub:
    """
    touch ${tumorname}_vs_${normalname}_${bed.simpleName}.varscan.vcf.gz
    """

}



process varscan_tonly {
    container "${params.containers.logan}"
    label 'process_somaticcaller'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        path(bed),
        path(tumorpileup),  path(tumor_con_table)

    output:
        tuple val(tumorname),
        path("${tumorname}_${bed.simpleName}.tonly.varscan.vcf.gz")

    shell:

    '''
    varscan_opts="--strand-filter 0 --min-var-freq 0.01 --output-vcf 1 --variants 1"
    pileup_cmd="samtools mpileup -d 100000 -q 20 -Q 20 -f !{GENOMEREF} -l !{bed} !{tumor}"
    varscan_cmd="varscan mpileup2cns <($pileup_cmd) $varscan_opts"

    eval "$varscan_cmd > !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf_temp"
    
    varscan filter !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf_temp \
    --strand-filter 1 --min-reads2 5 --min-strands2 2 --min-var-freq 0.05 --p-value 0.01 --min-avg-qual 30 > \
    !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf_temp1

    awk '{{gsub(/\\y[W|K|Y|R|S|M|B|D|H|V]\\y/,"N",\$4); OFS = "\t"; print}}' !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf_temp1 \
        | sed '/^$/d' | bcftools view - -Oz -o !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf

    printf "Sample1\t!{tumorname}\n" > sampname

    bcftools reheader -s sampname !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf \
        | bcftools view -Oz -o !{tumor.simpleName}_!{bed.simpleName}.tonly.varscan.vcf.gz

    '''

    stub:
    """
    touch ${tumorname}_${bed.simpleName}.tonly.varscan.vcf.gz
    """

}

