//SEQUENZA
GENOMEREF = file(params.genomes[params.genome].genome)
if(params.exome) {
    GC = params.genomes[params.genome].SEQUENZAGC_EXOME
}else{
    GC = file(params.genomes[params.genome].SEQUENZAGC)
}
SEQUENZA_SCRIPT = params.script_sequenza


process seqz_sequenza_bychr {
    container = "${params.containers.logan}"
    label 'process_long'

    input:
        tuple val(pairid), val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), val(chr)

    output:
        tuple val(pairid), path("${tumorname}_${normalname}_${chr}.seqz.gz")

    script:

    """
        sequenza-utils bam2seqz \
        -gc $GC \
        -F $GENOMEREF \
        -C ${chr} \
        -n ${normal} \
        -t ${tumor} | gzip > "${tumorname}_${normalname}_${chr}.seqz.gz"

    """

    stub:
    """
    touch "${tumorname}_${normalname}_${chr}.seqz.gz"
    """
}


process sequenza {
    container = "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(pairid), path(seqz), val(window)

    output:
        tuple val(pairid),
        path("${pairid}_alternative_solutions.txt"),
        path("${pairid}_alternative_fit.pdf"),
        path("${pairid}_model_fit.pdf"),
        path("${pairid}_confints_CP.txt"),
        path("${pairid}_CN_bars.pdf"),
        path("${pairid}_genome_view.pdf"),
        path("${pairid}_chromosome_view.pdf"),
        path("${pairid}_mutations.txt"),
        path("${pairid}_segments.txt"),
        path("${pairid}_CP_contours.pdf"),
        path("${pairid}_sequenza_cp_table.RData"),
        path("${pairid}_chromosome_depths.pdf"),
        path("${pairid}_gc_plots.pdf"),
        path("${pairid}_sequenza_extract.RData")


    shell:
    '''

    zcat !{seqz} | awk '{if (NR==1) {print $0} else {if ($1!="chromosome"){print $0}}}' |\
    sequenza-utils seqz_binning \
        -w !{window} \
        -s - > !{pairid}.bin!{window}.seqz

    Rscript !{SEQUENZA_SCRIPT} \
        !{pairid}.bin!{window}.seqz \
        . \
        !{pairid} \
        !{task.cpus}

    '''

    stub:

    """
    touch "${pairid}_alternative_solutions.txt"
    touch "${pairid}_alternative_fit.pdf"
    touch "${pairid}_model_fit.pdf"
    touch "${pairid}_confints_CP.txt"
    touch "${pairid}_CN_bars.pdf"
    touch "${pairid}_genome_view.pdf"
    touch "${pairid}_chromosome_view.pdf"
    touch "${pairid}_mutations.txt"
    touch "${pairid}_segments.txt"
    touch "${pairid}_CP_contours.pdf"
    touch "${pairid}_sequenza_cp_table.RData"
    touch "${pairid}_chromosome_depths.pdf"
    touch "${pairid}_gc_plots.pdf"
    touch "${pairid}_sequenza_extract.RData"

    """

}





process pileup_sequenza {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(pairid), val(name), 
        path(bam), path(bai), path(bed)

    output:
        tuple val(pairid), path("${name}_${bed}.mpileup.gz"), path("${name}_${bed}.mpileup.gz.tbi") 

    script:
    //Q20 is default in sequenza
    """
        samtools mpileup -f $GENOMEREF -R ${bed} -Q 20 ${bam} |gzip > ${name}_${bed}.mpileup.gz
        tabix -s1 -b2 -e2 ${name}_${bed}.mpileup.gz
    """

    stub:
    """
    touch "${name}_${bed}.mpileup.gz"
    touch "${name}_${bed}.mpileup.gz.tbi"
    """
}

process seqz_sequenza_reg {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(pairid), val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(pairid), path("${tumorname}_${normalname}_${chr}.seqz.gz")

    script:
    """
        sequenza-utils bam2seqz \
        -gc $GC \
        -p \
        -F $GENOMEREF \
        -n ${normal} \
        -t ${tumor} | gzip > "${tumorname}_${normalname}_${bed}.seqz.gz"

    """

    stub:
    """
    touch "${tumorname}_${normalname}_${chr}.seqz.gz"
    """
}


process seqz_sequenza {
    container = "${params.containers.logan}"
    label 'process_low'

    input:
        tuple val(pairid), val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai), path(bed)

    output:
        tuple val(pairid), path("${tumorname}_${normalname}_${chr}.seqz.gz")

    script:
    """
        sequenza-utils bam2seqz \
        -gc $GC \
        -p \
        -F $GENOMEREF \
        -n ${normal} \
        -t ${tumor} | gzip > "${tumorname}_${normalname}_${bed}.seqz.gz"

    """

    stub:
    """
    touch "${tumorname}_${normalname}_${chr}.seqz.gz"
    """
}
