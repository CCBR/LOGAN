# LOGAN development version

- Changed over to Nextflow CCBR template and pip packaging
    - Processes moved to `modules/local` directory
    - Workflows under the `subworkflows/local` directory
    - Processes fall under low/med/high, but added a somaticvariant caller process
    - Built AnnotSV/ClassifyCNV container (#40)
    - Converts strelka to add GT/AD column for downstream annotaiton (#55)
    - Adds hg19 genome build for for purple (#54)
    - Keeps VCF sample format as Tumor,Normal (#58)
    - Updated Docker base image to GATK 4.6, Adds cyvcf2 (for #54) Somalier 0.2.19, Muse 2.0.4 and HMFtools amber, purple, sage
