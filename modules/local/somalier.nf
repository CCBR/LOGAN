//SOMALIER
GENOMEREF=file(params.genomes[params.genome].genome)

SITES_VCF= file(params.genomes[params.genome].sites_vcf)
ANCESTRY_DB=file(params.genomes[params.genome].somalier_ancestrydb)
SCRIPT_PATH_GENDER = file(params.script_genderPrediction)
SCRIPT_PATH_SAMPLES = file(params.script_combineSamples)
SCRIPT_PATH_PCA = file(params.script_ancestry)


process somalier_extract {
    /*
    To estimate ancestry, Somalier first extracts known sites from mapped reads
    @Input:
        Mapped and pre-processed BAM file
    @Output:
        Exracted sites in (binary) somalier format
    
    params:
        sites_vcf = config['references']['SOMALIER']['SITES_VCF'],
        genomeFasta = config['references']['GENOME'],
        rname = 'somalier_extract'
    container: config['images']['wes_base']
    */
    container = "${params.containers.loganqc}"
    label 'process_low'

    input:
        tuple val(samplename), path(bam), path(bai)

    output:
        path("output/${samplename}.somalier")

    script:
    """
    mkdir -p output
    somalier extract \
        -d output \
        --sites $SITES_VCF \
        -f $GENOMEREF \
        $bam
    """

    stub:
    """
    mkdir -p output
    touch output/${samplename}.somalier
    """
}

process somalier_analysis_human {
    /*
    To estimate relatedness, Somalier uses extracted site information to
    compare across all samples. This step also runs the ancestry estimation
    function in Somalier.
    @Input:
        Exracted sites in (binary) somalier format for ALL samples in the cohort
    @Output:
        Separate tab-separated value (TSV) files with relatedness and ancestry outputs

    */
    container = "${params.containers.loganqc}"
    label 'process_low'
    errorStrategy='ignore'

    input:
        path(somalierin)

    output:
        tuple path("relatedness.pairs.tsv"), path("relatedness.samples.tsv"),
        path("ancestry.somalier-ancestry.tsv"), path("predicted.genders.tsv"),
        path("predicted.pairs.tsv"),
        path("sampleAncestryPCAPlot.html"),
        path("predictedPairsAncestry.pdf")

    script:
    """
    echo "Estimating relatedness"
    somalier relate \
        -o "relatedness" \
        $somalierin

    echo "Estimating ancestry"
    somalier ancestry \
        -o "ancestry" \
        --labels $ANCESTRY_DB/ancestry-labels-1kg.tsv \
        $ANCESTRY_DB/*.somalier ++ \
        $somalierin

    Rscript $SCRIPT_PATH_GENDER \
        relatedness.samples.tsv \
        predicted.genders.tsv

    Rscript $SCRIPT_PATH_SAMPLES \
        relatedness.pairs.tsv \
        predicted.pairs.tsv

    Rscript $SCRIPT_PATH_PCA \
        ancestry.somalier-ancestry.tsv \
        predicted.pairs.tsv \
        sampleAncestryPCAPlot.html \
        predictedPairsAncestry.pdf
    """

    stub:

    """
    touch relatedness.pairs.tsv
    touch relatedness.samples.tsv
    touch ancestry.somalier-ancestry.tsv predicted.genders.tsv
    touch predicted.pairs.tsv sampleAncestryPCAPlot.html
    touch predictedPairsAncestry.pdf
    """
}

process somalier_analysis_mouse {
    /*
    To estimate relatedness, Somalier uses extracted site information to
    compare across all samples. This step also runs the ancestry estimation
    function in Somalier.
    @Input:
        Exracted sites in (binary) somalier format for ALL samples in the cohort
    @Output:
        Separate tab-separated value (TSV) files with relatedness and ancestry outputs

    */
    container = "${params.containers.loganqc}"
    label 'process_low'
    errorStrategy='ignore'


    input:
        path(somalierin)

    output:
        tuple path("relatedness.pairs.tsv"),
        path("relatedness.samples.tsv"),
        path("predicted.genders.tsv"),
        path("predicted.pairs.tsv")

    script:
    """
    echo "Estimating relatedness"
    somalier relate \
        -o "relatedness" \
        $somalierin

    Rscript $SCRIPT_PATH_GENDER \
        relatedness.samples.tsv \
        predicted.genders.tsv

    Rscript $SCRIPT_PATH_SAMPLES \
        relatedness.pairs.tsv \
        predicted.pairs.tsv

    """

    stub:

    """
    touch relatedness.pairs.tsv
    touch relatedness.samples.tsv
    touch predicted.genders.tsv
    touch predicted.pairs.tsv

    """
}
