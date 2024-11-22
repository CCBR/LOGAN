ascatR =  params.script_ascat
if (params.genome.matches("hg38(.*)")){
    GENOMEVER="hg38"
} else if (params.genome.matches("hg19(.*)")){
    GENOMEVER="hg19"
}


process ascat_tn {
    container = "${params.containers.cnv}"
    label 'process_medium'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)

    output:
    tuple val(tumorname), 
    path("After_correction_${tumorname}.germline.png"),
    path("After_correction_${tumorname}.tumor.png"),
    path("Before_correction_${tumorname}.germline.png"),
    path("Before_correction_${tumorname}.tumor.png"),
    path("${tumorname}.ASCATprofile.png"),
    path("${tumorname}.ASPCF.png"),
    path("${tumorname}.sunrise.png"),
    path("${tumorname}_BAF.txt"),
    path("${tumorname}_LogR.txt"),
    path("${tumorname}.segments_raw.txt"),
    path("${tumorname}.segments.txt"),
    path("${tumorname}_vs_${normalname}.qc.txt"),
    path("${tumorname}_vs_${normalname}_ascat.Rdata")

    script:
    """
    Rscript $ascatR ${tumor} ${tumorname} ${normal} ${normalname} $GENOMEVER
    """

    stub:
    """
    touch After_correction_${tumorname}.germline.png
    touch After_correction_${tumorname}.tumor.png
    touch Before_correction_${tumorname}.germline.png
    touch Before_correction_${tumorname}.tumor.png
    touch ${tumorname}.ASCATprofile.png
    touch ${tumorname}.ASPCF.png
    touch ${tumorname}.sunrise.png
    touch ${tumorname}_BAF.txt
    touch ${tumorname}_LogR.txt
    touch ${tumorname}.segments_raw.txt
    touch ${tumorname}.segments.txt
    touch ${tumorname}_vs_${normalname}.qc.txt
    touch ${tumorname}_vs_${normalname}_ascat.Rdata

    """

}


process ascat_tn_exome {
    container = "${params.containers.cnv}"
    label 'process_medium'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai),
        path(bed)

    output:
    tuple path("After_correction_${tumorname}.germline.png"),
    path("After_correction_${tumorname}.tumour.png"),
    path("Before_correction_${tumorname}.germline.png"),
    path("Before_correction_${tumorname}.tumour.png"),
    path("${tumorname}.ASCATprofile.png"),
    path("${tumorname}.ASPCF.png"),
    path("${tumorname}.sunrise.png"),
    path("${tumorname}_BAF.txt"),
    path("${tumorname}_LogR.txt"),
    path("${tumorname}.segments_raw.txt"),
    path("${tumorname}.segments.txt"),
    path("${tumorname}_vs_${normalname}.qc.txt"),
    path("${tumorname}_vs_${normalname}_ascat.Rdata")

    script:
    """
    sed 's/^chr//' ${bed} > nochrtemp.bed
    Rscript $ascatR ${tumor} ${tumorname} ${normal} ${normalname} $GENOMEVER nochrtemp.bed
    """

    stub:
    """
    touch After_correction_${tumorname}.germline.png
    touch After_correction_${tumorname}.tumour.png
    touch Before_correction_${tumorname}.germline.png
    touch Before_correction_${tumorname}.tumour.png
    touch ${tumorname}.ASCATprofile.png
    touch ${tumorname}.ASPCF.png
    touch ${tumorname}.sunrise.png
    touch ${tumorname}_BAF.txt
    touch ${tumorname}_LogR.txt
    touch ${tumorname}.segments_raw.txt
    touch ${tumorname}.segments.txt
    touch ${tumorname}_vs_${normalname}.qc.txt
    touch ${tumorname}_vs_${normalname}_ascat.Rdata

    """

}
