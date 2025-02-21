//References
GENOMEREF=file(params.genomes[params.genome].genome)
GENOMEFAI=file(params.genomes[params.genome].genomefai)

//HMFTOOLS
SOMATICHOTSPOTS=params.genomes[params.genome].SOMATICHOTSPOTS
PANELBED=params.genomes[params.genome].PANELBED
HCBED=params.genomes[params.genome].HCBED
ENSEMBLCACHE=params.genomes[params.genome].ENSEMBLCACHE
GENOMEVER=params.genomes[params.genome].GENOMEVER



process sage_tn {
    container "${params.containers.logan}"
    label 'process_high'

    input:
        tuple val(tumorname), path(tumorbam), path(tumorbai),
        val(normalname), path(normalbam), path(normalbai)

    output:
        tuple val(tumorname), val(normalname),
        path("${tumorname}_vs_${normalname}.sage.vcf.gz"),
        path("${tumorname}_vs_${normalname}.sage.vcf.gz.tbi")

    
    script:
    """
    java -Xms4G -Xmx32G -cp /opt2/hmftools/sage.jar \
    -tumor ${tumorname} -tumor_bam ${tumorbam} \
    -reference ${normalname} -reference_bam ${normalbam} \
    -threads $task.cpus \
    -ref_genome_version $GENOMEVER \
    -ref_genome $GENOMEREF \
    -hotspots $SOMATICHOTSPOTS \
    $PANELBED $HCBED $ENSEMBLCACHE \
    -output_vcf ${tumorname}_vs_${normalname}.sage.vcf.gz
    """

    stub:
    """
    touch "${tumorname}_vs_${normalname}.sage.vcf.gz" "${tumorname}_vs_${normalname}.sage.vcf.gz.tbi"
    """
}



process sage_tonly {
    container "${params.containers.logan}"
    label 'process_somaticcaller'

    input:
        tuple val(tumorname), path(tumorbam), path(tumorbai)

    output:
        tuple val(tumorname), 
        path("${tumorname}.tonly.sage.vcf.gz"),
        path("${tumorname}.tonly.sage.vcf.gz.tbi")

    script:
    """
        java -Xms4G -Xmx32G -cp /opt2/hmftools/sage.jar \
        -tumor ${tumorname} -tumor_bam ${tumorbam} \
        -threads $task.cpus \
        -ref_genome_version $GENOMEVER \
        -ref_genome $GENOMEREF \
        -hotspots $HOTSPOTS \
        $PANELBED $HCBED $ENSEMBLCACHE \
        -output_vcf ${tumorname}.tonly.sage.vcf.gz
    """

    stub:
    """
        touch "${tumorname}.tonly.sage.vcf.gz" "${tumorname}.tonly.sage.vcf.gz.tbi"
    """

}
