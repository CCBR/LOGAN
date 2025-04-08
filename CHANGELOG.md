# LOGAN development version

- Now using the readthedocs theme for the docs website. (#84, @kelly-sovacool)
- LOGAN is now archived in Zenodo with DOI `10.5281/zenodo.14907169`. (#87, @kelly-sovacool)

## LOGAN 0.2.1

### New features

- Fixes #72 for fastq naming in fastp
- Fixes #78 for qualimap bamqc

## LOGAN 0.2.0

### New features

- Added additional SV callers(GRIDSS) and annotation for SV (GRIPSS) + CNV Callers (ASCAT, CNVKit) + SNV ([#66](https://github.com/CCBR/LOGAN/issues/66)Deepsomatic)
- Adds ffpe filtering [#67](https://github.com/CCBR/LOGAN/issues/67)
- Bugfixes for hg19 by fixing references
- Updated PON for hg38 using TCGA/GDC references [#59](https://github.com/CCBR/LOGAN/issues/59)
- In development: adding exome support by using bed file to restrict calling regions
- Refactored modules to be similar to nf-core
- Fix error in varscan [#71](https://github.com/CCBR/LOGAN/issues/71)
- Created separate docker [#63](https://github.com/CCBR/LOGAN/issues/63)

## LOGAN 0.1.0

### Features

- Changed over to Nextflow CCBR template and pip packaging
  - Processes moved to `modules/local` directory
  - Workflows under the `subworkflows/local` directory
  - Processes fall under low/med/high, but added a somaticvariant caller process
  - Built AnnotSV/ClassifyCNV container (#40)
  - Converts strelka to add GT/AD column for downstream annotaiton (#55)
  - Adds hg19 genome build for for purple (#54)
  - Keeps VCF sample format as Tumor,Normal (#58)
  - Updated Docker base image to GATK 4.6, Adds cyvcf2 (for #54) Somalier 0.2.19, Muse 2.0.4 and HMFtools amber, purple, sage
