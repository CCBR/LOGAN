process ffpe_1 {

    container "${params.containers.logan}"
    label 'process_medium'

    input:
    tuple val(tumorsample), val(normal),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz"),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz.tbi"),
        bam, bamindex
        tuple val(tumorsample), val(normal),
        val(caller),
        path(vcfs), path(vcfindex)

    output:
        tuple val(tumorsample), val(normal),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz"),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz.tbi")

    script:


    stub:
    """
    touch "${tumorsample}_vs_${normal}_combined_ffpolish.vcf.gz"
    """
}



/*DISCVR
process somaticcombine {
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumorsample), val(normal),
        val(callers),
        path(vcfs), path(vcfindex)

    output:
        tuple val(tumorsample), val(normal),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz"),
        path("${tumorsample}_vs_${normal}_combined.vcf.gz.tbi")

    script:
        vcfin1=[callers, vcfs].transpose().collect { a, b -> a + " " + b }
        vcfin2="-V:" + vcfin1.join(" -V:")

    """
    java -jar \$DISCVRSeq_JAR MergeVcfsAndGenotypes \
        -R $GENOMEREF \
        --genotypeMergeOption PRIORITIZE \
        --priority_list mutect2,strelka,octopus,muse,lofreq,vardict,varscan \
        --filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED \
        -O ${tumorsample}_vs_${normal}_combined.vcf.gz \
        $vcfin2
    """

    stub:
    vcfin1=[callers, vcfs].transpose().collect { a, b -> a + " " + b }
    vcfin2="-V:" + vcfin1.join(" -V:")

    """
    touch ${tumorsample}_vs_${normal}_combined.vcf.gz
    touch ${tumorsample}_vs_${normal}_combined.vcf.gz.tbi
    """

}
*/



/*DISCVRSeq
process somaticcombine_tonly {
    container "${params.containers.logan}"
    label 'process_medium'

    input:
        tuple val(tumorsample),
        val(callers),
        path(vcfs), path(vcfindex)

    output:
        tuple val(tumorsample),
        path("${tumorsample}_combined_tonly.vcf.gz"),
        path("${tumorsample}_combined_tonly.vcf.gz.tbi")

    script:
        vcfin1=[callers, vcfs].transpose().collect { a, b -> a + " " + b }
        vcfin2="-V:" + vcfin1.join(" -V:")

    """
    java -jar \$DISCVRSeq_JAR MergeVcfsAndGenotypes \
        -R $GENOMEREF \
        --genotypeMergeOption PRIORITIZE \
        --priority_list mutect2_tonly,octopus_tonly,vardict_tonly,varscan_tonly \
        --filteredRecordsMergeType KEEP_IF_ANY_UNFILTERED \
        -O ${tumorsample}_combined_tonly.vcf.gz \
        $vcfin2
    """

    stub:
    """
    touch ${tumorsample}_combined_tonly.vcf.gz ${tumorsample}_combined_tonly.vcf.gz.tbi
    """

}
*/
