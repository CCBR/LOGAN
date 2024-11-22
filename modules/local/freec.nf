//FREEC
//mm10 Paired-Sequenza, FREEC-tumor only
//CNV Intervals
if (params.intervals){
    CNVTARGETS = file(params.intervals)
}else{
    CNVTARGETS = file(params.genomes[params.genome].intervals)
}


REFORMATBED = params.script_reformatbed
FREEC_SCRIPT = params.script_freec
if(params.exome){
    FREECPAIR_SCRIPT = params.script_freecpaired_exome
}else{
    FREECPAIR_SCRIPT = params.script_freecpaired
}
FREECSIGNIFICANCE = params.freec_significance
FREECLENGTHS = file(params.genomes[params.genome].FREEC.FREECLENGTHS)
FREECCHROMS = file(params.genomes[params.genome].FREEC.FREECCHROMS)
FREECPILEUP = file(params.genomes[params.genome].FREEC.FREECPILEUP)
FREECSNPS = file(params.genomes[params.genome].FREEC.FREECSNPS)
FREECPLOT = params.freec_plot


process freec_paired {
    container = "${params.containers.logan}"
    label 'process_long'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_CNVs.p.value.txt"),
        path("${tumorname}_vs_${normalname}_ratio.txt"),
        path("${tumorname}_vs_${normalname}_BAF.txt"),
        path("${tumorname}_vs_${normalname}_ratio.txt.log2.png"),
        path("${tumorname}_vs_${normalname}_ratio.txt.png")

    shell:

    """

    perl $FREECPAIR_SCRIPT \
        . \
        $FREECLENGTHS \
        $FREECCHROMS \
        ${tumor} \
        ${normal} \
        $FREECPILEUP \
        $GENOMEREF \
        $FREECSNPS \
        $CNVTARGETS

    freec -conf freec_genome_config.txt

    cat $FREECSIGNIFICANCE | \
        R --slave \
        --args ${tumor}_CNVs \
        ${tumor}_ratio.txt

    cat $FREECPLOT | \
        R --slave \
        --args 2 \
        ${tumor}_ratio.txt \
        ${tumor}_BAF.txt

    mv ${tumor}_CNVs.p.value.txt ${tumorname}_vs_${normalname}_CNVs.p.value.txt
    mv ${tumor}_ratio.txt ${tumorname}_vs_${normalname}_ratio.txt
    mv ${tumor}_BAF.txt ${tumorname}_vs_${normalname}_BAF.txt
    mv ${tumor}_BAF.txt.png ${tumorname}_vs_${normalname}_BAF.txt.png
    mv ${tumor}_ratio.txt.log2.png ${tumorname}_vs_${normalname}_ratio.txt.log2.png
    mv ${tumor}_ratio.txt.png ${tumorname}_vs_${normalname}_ratio.txt.png

    """

    stub:
    """
    touch ${tumorname}_vs_${normalname}_CNVs.p.value.txt
    touch ${tumorname}_vs_${normalname}_ratio.txt
    touch ${tumorname}_vs_${normalname}_BAF.txt
    touch ${tumorname}_vs_${normalname}_BAF.txt.png
    touch ${tumorname}_vs_${normalname}_ratio.txt.log2.png
    touch ${tumorname}_vs_${normalname}_ratio.txt.png

    """
}




process freec_paired_exome {
    container = "${params.containers.logan}"
    label 'process_long'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}_CNVs.p.value.txt"),
        path("${tumorname}_vs_${normalname}_ratio.txt"),
        path("${tumorname}_vs_${normalname}_BAF.txt"),
        path("${tumorname}_vs_${normalname}_ratio.txt.log2.png"),
        path("${tumorname}_vs_${normalname}_ratio.txt.png")

    shell:

    """
    python $REFORMATBED -i $CNVTARGETS 
    perl $FREECPAIR_SCRIPT \
        . \
        $FREECLENGTHS \
        $FREECCHROMS \
        ${tumor} \
        ${normal} \
        $FREECPILEUP \
        $GENOMEREF \
        $FREECSNPS \
        exome_targets.bed

    freec -conf freec_exome_config.txt

    cat $FREECSIGNIFICANCE | \
        R --slave \
        --args ${tumor}_CNVs \
        ${tumor}_ratio.txt

    cat $FREECPLOT | \
        R --slave \
        --args 2 \
        ${tumor}_ratio.txt \
        ${tumor}_BAF.txt

    mv ${tumor}_CNVs.p.value.txt ${tumorname}_vs_${normalname}_CNVs.p.value.txt
    mv ${tumor}_ratio.txt ${tumorname}_vs_${normalname}_ratio.txt
    mv ${tumor}_BAF.txt ${tumorname}_vs_${normalname}_BAF.txt
    mv ${tumor}_BAF.txt.png ${tumorname}_vs_${normalname}_BAF.txt.png
    mv ${tumor}_ratio.txt.log2.png ${tumorname}_vs_${normalname}_ratio.txt.log2.png
    mv ${tumor}_ratio.txt.png ${tumorname}_vs_${normalname}_ratio.txt.png

    """

    stub:
    """
    touch ${tumorname}_vs_${normalname}_CNVs.p.value.txt
    touch ${tumorname}_vs_${normalname}_ratio.txt
    touch ${tumorname}_vs_${normalname}_BAF.txt
    touch ${tumorname}_vs_${normalname}_BAF.txt.png
    touch ${tumorname}_vs_${normalname}_ratio.txt.log2.png
    touch ${tumorname}_vs_${normalname}_ratio.txt.png

    """
}



process freec {
    container = "${params.containers.logan}"
    label 'process_medium'
    errorStrategy 'ignore'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai)

    output:
        tuple val(tumorname),
        path("${tumorname}_CNVs.p.value.txt"),
        path("${tumorname}_ratio.txt"),
        path("${tumorname}_BAF.txt"),
        path("${tumorname}_ratio.txt.log2.png"),
        path("${tumorname}_ratio.txt.png")


    shell: 
    
    """
    perl $FREEC_SCRIPT \
        . \
        $FREECLENGTHS \
        $FREECCHROMS \
        ${tumor} \
        $FREECPILEUP \
        $GENOMEREF \
        $FREECSNPS \
        $CNVTARGETS

    freec -conf freec_genome_config.txt

    cat $FREECSIGNIFICANCE | \
        R --slave \
        --args ${tumor}_CNVs \
        ${tumor}_ratio.txt

    cat $FREECPLOT | \
        R --slave \
        --args 2 \
        ${tumor}_ratio.txt \
        ${tumor}_BAF.txt

    mv ${tumor}_CNVs.p.value.txt ${tumorname}_CNVs.p.value.txt
    mv ${tumor}_ratio.txt ${tumorname}_ratio.txt
    mv ${tumor}_BAF.txt ${tumorname}_BAF.txt
    mv ${tumor}_BAF.txt.png ${tumorname}_BAF.txt.png
    mv ${tumor}_ratio.txt.log2.png ${tumorname}_ratio.txt.log2.png
    mv ${tumor}_ratio.txt.png ${tumorname}_ratio.txt.png

    """

    stub:
    """
    touch ${tumorname}_CNVs.p.value.txt
    touch ${tumorname}_ratio.txt
    touch ${tumorname}_BAF.txt
    touch ${tumorname}_BAF.txt.png
    touch ${tumorname}_ratio.txt.log2.png
    touch ${tumorname}_ratio.txt.png

    """
}
