

//CNVKIT
GENOMEREF = file(params.genomes[params.genome].genome)

REFFLAT = file(params.genomes[params.genome].REFFLAT)
ACCESS = file(params.genomes[params.genome].ACCESS)

process cnvkit {
    container = "${params.containers.cnv}"
    label 'process_medium'

    input:
    tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai)


    output:
    tuple val(tumorname), path("${tumorname}")


    script:
    """
    cnvkit.py batch ${tumor} --normal ${normal} \
    --annotate $REFFLAT \
    --fasta $GENOMEREF --access $ACCESS \
    --output-reference ${tumorname}.cnn --output-dir ${tumorname}/  \
    --diagram --scatter \
    -m wgs -p $task.cpus
    """

    stub:
    """
    mkdir ${tumorname}
    touch ${tumorname}/${normalname}.antitargetcoverage.cnn ${tumorname}/${normalname}.targetcoverage.cnn
    touch ${tumorname}/${tumorname}.antitargetcoverage.cnn ${tumorname}/${tumorname}.targetcoverage.cnn
    touch ${tumorname}/${tumorname}.bintest.cns ${tumorname}/${tumorname}.call.cns ${tumorname}/${tumorname}.cnr ${tumorname}/${tumorname}.cns ${tumorname}/${tumorname}-diagram.pdf ${tumorname}/${tumorname}-scatter.png
    """

}

process cnvkit_exome {
    container = "${params.containers.cnv}"
    label 'process_medium'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai),
        path(bed)

    output:
    tuple val(tumorname), path("${tumorname}")

    script:
    """
    cnvkit.py batch ${tumor} --normal ${normal} \
    --targets ${bed} --annotate $REFFLAT \
    --fasta $GENOMEREF --access $ACCESS \
    --output-reference ${tumorname}.cnn --output-dir ${tumorname}/ \
    --diagram --scatter -p $task.cpus
    """

    stub:
    """
    mkdir ${tumorname}
    touch ${tumorname}/${normalname}.antitargetcoverage.cnn ${tumorname}/${normalname}.targetcoverage.cnn
    touch ${tumorname}/${tumorname}.antitargetcoverage.cnn ${tumorname}/${tumorname}.targetcoverage.cnn
    touch ${tumorname}/${tumorname}.bintest.cns ${tumorname}/${tumorname}.call.cns ${tumorname}/${tumorname}.cnr ${tumorname}/${tumorname}.cns ${tumorname}/${tumorname}-diagram.pdf ${tumorname}/${tumorname}-scatter.png

    """

}



process cnvkit_tonly {
    container = "${params.containers.cnv}"
    label 'process_medium'

    input:
    tuple val(tumorname), path(tumor), path(tumorbai)

    output:
    tuple val(tumorname), path("${tumorname}")


    script:
    """
    cnvkit.py batch ${tumor} -n \
    --annotate $REFFLAT \
    --fasta $GENOMEREF --access $ACCESS \
    --output-reference ${tumorname}.cnn --output-dir ${tumorname}/  \
    --diagram --scatter \
    -m wgs -p $task.cpus
    """

    stub:
    """
    mkdir ${tumorname}
    touch ${tumorname}/${tumorname}.antitargetcoverage.cnn ${tumorname}/${tumorname}.targetcoverage.cnn
    touch ${tumorname}/${tumorname}.bintest.cns ${tumorname}/${tumorname}.call.cns ${tumorname}/${tumorname}.cnr ${tumorname}/${tumorname}.cns ${tumorname}/${tumorname}-diagram.pdf ${tumorname}/${tumorname}-scatter.png
    """

}

process cnvkit_exome_tonly {
    container = "${params.containers.cnv}"
    label 'process_medium'

    input:
        tuple val(tumorname), path(tumor), path(tumorbai),
        val(normalname), path(normal), path(normalbai),
        path(bed)

    output:
    tuple val(tumorname), path("${tumorname}")

    script:
    """
    cnvkit.py batch ${tumor} --normal \
    --targets $bed --annotate $REFFLAT \
    --fasta $GENOMEREF --access $ACCESS \
    --output-reference ${tumorname}.cnn --output-dir ${tumorname}/ \
    --diagram --scatter -p $task.cpus
    """

    stub:
    """
    mkdir ${tumorname}
    touch ${tumorname}/${normalname}.antitargetcoverage.cnn ${tumorname}/${normalname}.targetcoverage.cnn
    touch ${tumorname}/${tumorname}.antitargetcoverage.cnn ${tumorname}/${tumorname}.targetcoverage.cnn
    touch ${tumorname}/${tumorname}.bintest.cns ${tumorname}/${tumorname}.call.cns ${tumorname}/${tumorname}.cnr ${tumorname}/${tumorname}.cns ${tumorname}/${tumorname}-diagram.pdf ${tumorname}/${tumorname}-scatter.png

    """

}

