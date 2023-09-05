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
        path("${tumor.simpleName}.alignments.txt.gz")

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


