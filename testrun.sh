nextflow run ../wgs-seek/wgs-seek.nf  --profile local \
--input "/data/CCBR/rawdata/nousome/small_truth_set/test_out_1/SRR*test*R{1,2}.fastq.gz" \
--sample_sheet samplesheet --output . --resume

nextflow run ../wgs-seek/test_split.nf  --profile local \
--input "/data/CCBR/rawdata/nousome/small_truth_set/test_out_1/SRR*test*R{1,2}.fastq.gz" \
--intervals "/data/nousomedr/wgs/test1/beds/*.bed"
