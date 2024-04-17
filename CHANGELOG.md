# LOGAN development version

- Changed over to Nextflow CCBR template and pip packaging
    - Processes moved to `modules/local` directory
    - Workflows under the `subworkflows/local` directory
    - Processes fall under low/med/high, but added a somaticvariant caller process
    - Built AnnotSV/ClassifyCNV container (#40)