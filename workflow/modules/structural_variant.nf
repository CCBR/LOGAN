GENOMEREF=file(params.genomes[params.genome].genome)
GENOME=params.genome
BWAGENOME=file(params.genomes[params.genome].bwagenome)
DBSNP_INDEL=file(params.genomes[params.genome].KNOWNINDELS) 

outdir=file(params.output)


process svaba_somatic {
    publishDir(path: "${outdir}/SV/svaba", mode: 'copy') 

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}.bps.txt.gz"),
        path("${tumor.simpleName}.contigs.bam"),
        path("${tumor.simpleName}.discordant.txt.gz"),
        path("${tumor.simpleName}.alignments.txt.gz"),
        path("${tumor.simpleName}.svaba.germline.indel.vcf"),
        path("${tumor.simpleName}.svaba.germline.sv.vcf"),
        path("${tumor.simpleName}.svaba.somatic.indel.vcf"),
        path("${tumor.simpleName}.svaba.somatic.sv.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.germline.indel.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.germline.sv.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.somatic.indel.vcf"),
        path("${tumor.simpleName}.svaba.unfiltered.somatic.sv.vcf"),
        path("${tumor.simpleName}.log")


    script:
    """
    svaba run -t ${tumor} -n ${normal} -p $task.cpus -D $DBSNP_INDEL -a ${tumor.simpleName} -G $BWAGENOME
    """

    stub:
    
    """
    touch "${tumor.simpleName}.bps.txt.gz"
    touch "${tumor.simpleName}.contigs.bam"
    touch "${tumor.simpleName}.discordant.txt.gz"
    touch "${tumor.simpleName}.alignments.txt.gz"
    touch "${tumor.simpleName}.svaba.germline.indel.vcf"
    touch "${tumor.simpleName}.svaba.germline.sv.vcf"
    touch "${tumor.simpleName}.svaba.somatic.indel.vcf"
    touch "${tumor.simpleName}.svaba.somatic.sv.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.germline.indel.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.germline.sv.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.somatic.indel.vcf"
    touch "${tumor.simpleName}.svaba.unfiltered.somatic.sv.vcf"
    touch "${tumor.simpleName}.log"

    """
}


process manta_somatic {
    //https://github.com/illumina/manta
    publishDir(path: "${outdir}/SV/manta", mode: 'copy') 

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}.diplodSV.vcf.gz"),
        path("${tumor.simpleName}.somaticSV.vcf.gz"),
        path("${tumor.simpleName}.candidateSV.vcf.gz"),
        path("${tumor.simpleName}.candidateSmallIndels.vcf.gz")

    script:
    """
    mkdir -p wd

    configManta.py \
        --normalBam=${normal} \
        --tumorBam=${tumor} \
        --referenceFasta=$GENOMEREF \
        --runDir=wd

    wd/runWorkflow.py -m local -j 10 -g 10
    
    mv wd/results/variants/diploidSV.vcf.gz ${tumor.simpleName}.diplodSV.vcf.gz
    mv wd/results/variants/somaticSV.vcf.gz ${tumor.simpleName}.somaticSV.vcf.gz
    mv wd/results/variants/candidateSV.vcf.gz ${tumor.simpleName}.candidateSV.vcf.gz
    mv wd/results/variants/candidateSmallIndels.vcf.gz ${tumor.simpleName}.candidateSmallIndels.vcf.gz

    """

    stub:
    
    """
    touch ${tumor.simpleName}.diplodSV.vcf.gz
    touch ${tumor.simpleName}.somaticSV.vcf.gz
    touch ${tumor.simpleName}.candidateSV.vcf.gz
    touch ${tumor.simpleName}.candidateSmallIndels.vcf.gz
    """
}


process annotsv_tn {
     //AnnotSV for Manta/Svaba works with either vcf.gz or .vcf files
     //Requires bedtools,bcftools

    module = ['annotsv/3.3.1']
    publishDir(path: "${outdir}/SV/annotated", mode: 'copy') 

    input:
        tuple val(tumorname), path(somaticvcf)

    output:
        tuple val(tumorname),
        path("${tumorname}.tsv"),
        path("${tumorname}.unannotated.tsv")


    script:
    """
    AnnotSV -SVinputFile ${somaticvcf} \
    -genomeBuild $GENOME \
    -SVinputInfo 1 -outputFile ${tumorname} \
    -outputDir .

    """

    stub:
    """
    touch "${tumorname}.tsv"
    touch "${tumorname}.unannotated.tsv"
    """
}

