GENOME=file(params.genome)
GENOMEDICT=file(params.genomedict)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
SHAPEITINDEL=file(params.shapeitindel) //ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
KGP=file(params.kgp) //1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.dbsnp) //dbsnp_138.hg38.vcf.gz"
DBSNP_INDEL=file(params.dbsnp_indel) //dbsnp_138.hg38.vcf.gz"
BWAGENOME=file(params.bwagenome)
GNOMAD=file(params.gnomad) //somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 

outdir=file(params.output)


process svaba_somatic {
    //https://github.com/walaj/svaba/releases/tag/v1.2.0
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
    /data/nousomedr/programs/svaba run -t ${tumor} -n ${normal} -p $task.cpus -D $DBSNP_INDEL -a ${tumor.simpleName} -G $BWAGENOME
    """

    stub:
    
    """
    touch "${tumor.simpleName}.bps.txt.gz"
    touch "${tumor.simpleName}.contigs.bam"
    touch "${tumor.simpleName}.discordants.txt.gz"
    touch "${tumor.simpleName}.alignments.txt.gz"
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
    --referenceFasta=$GENOME \
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


process annotsv_tn{
     //AnnotSV for Manta/Svaba works with either vcf.gz or .vcf files
     //Requires bedtools,bcftools

    publishDir(path: "${outdir}/SV/annotated", mode: 'copy') 

    input:
        tuple val(tumorname), path(somaticvcf)

    output:
        tuple val(tumorname),
        path("${tumor.simpleName}.tsv"),
        path("${tumor.simpleName}.unannotated.tsv"),


    script:
    """
    AnnotSV -SVinputFile ${tumor}.vcf.gz -SVinputInfo 1 -outputFile ${tumor.simpleName} -outputDir .
    """

    stub:
    """
    touch "${tumor.simpleName}.tsv"
    touch "${tumor.simpleName}.unannotated.tsv"
    """
}

