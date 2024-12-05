
process multiqc {
    """
    Reporting step to aggregate sample summary statistics and quality-control
    information across all samples. This will be one of the last steps of the
    pipeline. The inputs listed here are to ensure that this step runs last.
    During runtime, MultiQC will recursively crawl through the working directory
    and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """
    container = "${params.containers.multiqc}"
    label 'process_low'

    input:
        path(allqcin)

    output:
        path("MultiQC_Report.html")

    script:

    """
    multiqc . \
    -f --interactive \
    -n "MultiQC_Report.html" \
    """

    stub:

    """
    touch MultiQC_Report.html
    """
}
