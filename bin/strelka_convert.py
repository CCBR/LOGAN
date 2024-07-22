#!/usr/bin/env python
import numpy as np
import cyvcf2
import sys 
import gzip 
import os 

##Adapted from https://github.com/bcbio/bcbio-nextgen/blob/72f42706faa5cfe4f0680119bf148e0bdf2b78ba/bcbio/variation/strelka2.py#L30
def _tumor_normal_genotypes(ref, alt, info):
    """Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.

    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention

    ref: The REF allele from a VCF line
    alt: A list of potentially multiple ALT alleles (rec.ALT.split(";"))
    info: The VCF INFO field
    fname, coords: not currently used, for debugging purposes
    """
    known_names = set(["het", "hom", "ref", "conflict"])
    def name_to_gt(val):
        if val.lower() == "het":
            return "0/1"
        elif val.lower() == "hom":
            return "1/1"
        elif val.lower() in set(["ref", "conflict"]):
            return "0/0"
        else:
            # Non-standard representations, het is our best imperfect representation
            # print(fname, coords, ref, alt, info, val)
            return "0/1"
    def alleles_to_gt(val):
        gt_indices = {gt.upper(): i for i, gt in enumerate([ref] + alt)}
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts) == 0:
                tumor_gt = "0/0"
            elif 0 in tumor_gts:
                tumor_gt = "0/%s" % min([x for x in tumor_gts if x > 0])
            else:
                tumor_gt = "%s/%s" % (min(tumor_gts), max(tumor_gts))
        else:
            tumor_gt = name_to_gt(val)
        return tumor_gt
    nt_val = [x.split("=")[-1] for x in info if x.startswith("NT=")][0]
    normal_gt = name_to_gt(nt_val)
    sgt_val = [x.split("=")[-1] for x in info if x.startswith("SGT=")]
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val = sgt_val[0].split("->")[-1]
        tumor_gt = alleles_to_gt(sgt_val)
    return normal_gt, tumor_gt


def _af_annotate_and_filter(in_file,out_file):
    """Populating FORMAT/AF, and dropping variants with AF<min_allele_fraction

    Strelka2 doesn't report exact AF for a variant, however it can be calculated as alt_counts/dp from existing fields:
    somatic
      snps:    GT:DP:FDP:SDP:SUBDP:AU:CU:GU:TU                 dp=DP                {ALT}U[0] = alt_counts(tier1,tier2)
      indels:  GT:DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50  dp=DP                TIR = alt_counts(tier1,tier2)
    germline
      snps:    GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL(:PS)       dp=sum(alt_counts)   AD = ref_count,alt_counts
      indels:  GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL(:PS)             dp=sum(alt_counts)   AD = ref_count,alt_counts
    """
    #data = paired.tumor_data if paired else items[0]
    #min_freq = float(utils.get_in(data["config"], ("algorithm", "min_allele_fraction"), 10)) / 100.0
    #logger.debug("Filtering Strelka2 calls with allele fraction threshold of %s" % min_freq)
    vcf = cyvcf2.VCF(in_file)
    vcf.add_format_to_header(dict(
        ID='AF', Number='1',Type='Float',
        Description='Allele frequency, as calculated in bcbio: AD/DP (germline), <ALT>U/DP (somatic snps), TIR/DPI (somatic indels)'
        ))
    vcf.add_format_to_header(dict(
        ID='AD', Number='.',Type='Integer',
        Description='Depth of reads supporting alleles'
        ))
    writer = cyvcf2.Writer(out_file, vcf)
    for rec in vcf:
        if rec.is_snp:  # snps?
            alt_counts = rec.format(rec.ALT[0] +'U')[:,0]  # {ALT}U=tier1_depth,tier2_depth
            ref_counts = rec.format(rec.REF + 'U')[:,0]
        else:  # indels
            alt_counts = rec.format('TIR')[:,0]  # TIR=tier1_depth,tier2_depth
            ref_counts = rec.format('TAR')[:,0]
        dp = rec.format('DP')[:,0]
        if dp is not None :
            with np.errstate(divide='ignore', invalid='ignore'):  # ignore division by zero and put AF=.0
                #alt_n = alt_counts_n[0]/DP_n
                #alt_t = alt_counts_t[0]/DP_t
                ad = [ref_counts[0],alt_counts[0],ref_counts[1],alt_counts[1]]
                af = np.true_divide(alt_counts, dp)
                rec.set_format('AF',np.round(af,5))
                rec.set_format('AD',np.array(ad))
                writer.write_record(rec)
    writer.close()


def _add_gt(in_file):
        ##Set genotypes now
        out_file = os.path.basename(in_file).replace(".vcf.gz", "-fixed.vcf")
        #open_fn = gzip.open if is_gzipped(in_file) else open
        with gzip.open(in_file,'rt') as in_handle: 
            with open(out_file,"wt") as out_handle:
                added_gt = False
                for line in in_handle:
                    if line.startswith("##FORMAT") and not added_gt:
                        added_gt = True
                        out_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                        out_handle.write(line)
                    elif line.startswith("#CHROM"):
                        assert added_gt
                        out_handle.write(line)
                    elif line.startswith("#"):
                        out_handle.write(line)
                    else:
                        parts = line.rstrip().split("\t")
                        normal_gt,tumor_gt = _tumor_normal_genotypes(parts[3], parts[4].split(","),
                                                                        parts[7].split(";"))
                        parts[8] = "GT:%s" % parts[8]
                        parts[9] = "%s:%s" % (normal_gt, parts[9])
                        parts[10] = "%s:%s" % (tumor_gt, parts[10])
                        out_handle.write("\t".join(parts) + "\n")

if __name__ == '__main__':
    filename = sys.argv[1]
    outname = sys.argv[2]
    _add_gt(filename)
    newname = os.path.basename(filename).replace(".vcf.gz", "-fixed.vcf")
    _af_annotate_and_filter(newname,outname)
    os.remove(newname)

