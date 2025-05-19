# LOGAN development version

- LOGAN now depends on ccbr_tools v0.4 for updated jobby & spooker utilities. (#95, @kelly-sovacool)
- Fix fastq naming in fastp. (#73, @dnousome)
- Fix resources for qualimap bamqc. (#79, @dnousome)
- Now using the readthedocs theme for the docs website. (#84, @kelly-sovacool)
- LOGAN is now archived in Zenodo with DOI `10.5281/zenodo.14907169`. (#87, @kelly-sovacool)
- CLI updates: (#93, @kelly-sovacool)
  - Use `nextflow run -resume` by default, or turn it off with `logan run --forceall`.
  - Add `--output` argument for `logan init` and `logan run`.
    - If not provided, commands are run in the current working directory.
    - This is equivalent to the nextflow `$launchDir` constant.
  - The nextflow preview is printed before launching the actual run.
  - Set the `publish_dir_mode` nextflow option to `link` by default.
  - Set the `process.cache` nextflow option to `deep` by default rather than lenient on biowulf.
- Sequenza chromosome M was removed for CNV analysis and all chromosomes used by default (#91, @dnousome)
- Increased memory and lscratch allocation for applybqsr process (#91, @dnousome)
- Manta output order was listed incorrectly (#91, @dnousome)
- Additional hg19 genome references were fixed (#91, @dnousome)
- Vardict Tumor only mode filtering using INFO/DP vs DP only (#91, @dnousome)
- Fix for Deepvariant to allow sample name filtering (#91, @dnousome)

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
