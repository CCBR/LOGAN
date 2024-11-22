# LOGAN Tools and Tools Tested


## SNV
| Tools |Pros | Cons | Used in Logan|
|----|---|---|---
|Mutect2     |Part of GATK best practices| | x |
|Strelka     | Fast|  Paired only|x|
|Muse     | Fast|  Paired only, can't be parallelized|x|
|Lofreq     | Low frequency variants|  Slow,Paired only|x|
|Vardict | Fast | Lower accuracy|x|
|Varscan  | Fast| Lower accuracy|x||
|Octopus | Accurate| Slow,High memory|x|
|Deepsomatic|Relatively fast|Trained on human data|x|


## Structural Variants
| Tools |Pros | Cons | Approach| Used in Logan|
|----|---|---|---|---|
|Manta     |Accurate, fast| |graph-based| x | 
|SVABA     | Deletion detection||local assembly+ multiple alignment|x|
|GRIDSS     | Provides blacklist|  Slow, part of HMFtools pipeline|Break end assembly (discordant +split)|x|

Manta, GridSS, and SvABA are based on read-pairs, split-reads, and local-assemblies.  
References [Joe et al](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-10239-9)

## Copy Number

| Tools |Pros | Cons | Used in Logan|
|----|---|---|---|
|Purple     |Complete workflow|Doesn't support mm10, requires SV,SNV calls as well | x |
|Sequenza     | Purity/Ploidy||x|
|FREEC | Fast | No Purity/Ploidy Estimatation|x|
|ASCAT  | Fast, Purity/Ploidy| |x|
|CNVkit |Fast | No Purity/Ploidy Estimatation|x|
|PureCN|Tumor only|Needs Panel of Normals on Sequencing|



## Germline
| Tools |Pros | Cons | Used in Logan| 
|----|---|---|---|
|Deepvariant     |Fast, most accurate| Model trained on human genomes (May not support mm10)| x|
