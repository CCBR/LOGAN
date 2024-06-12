#!/usr/bin/env python
import os
import numpy as np
import vcfpy
import sys 

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
        gt_indices = {gt.upper(): i for i, gt in enumerate([ref] + [alt])}
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
    nt_val = info.get('NT').split("=")[-1]
    normal_gt = name_to_gt(nt_val)
    sgt_val = info.get('SGT').split("=")[-1]
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val = sgt_val.split("->")[-1]
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
    vcf = vcfpy.Reader.from_path(in_file)
    vcf.header.add_format_line(vcfpy.OrderedDict([
        ('ID', 'AF'), 
        ('Description', 'Allele frequency, as calculated in bcbio: AD/DP (germline), <ALT>U/DP (somatic snps), TIR/DPI (somatic indels)'),
        ('Type','Float'),
        ('Number', '.')
    ]))
    vcf.header.add_format_line(vcfpy.OrderedDict([
        ('ID', 'GT'), 
        ('Description', 'Genotype'),
        ('Type','String'),
        ('Number', '1')
    ]))
    writer = vcfpy.Writer.from_path(out_file, vcf.header)
    for rec in vcf:
        #print(rec)
        if rec.is_snv():  # snps?
            alt_counts_n = rec.calls[0].data[rec.ALT[0].value + "U"]  # {ALT}U=tier1_depth,tier2_depth
            alt_counts_t = rec.calls[1].data[rec.ALT[0].value + "U"]  # {ALT}U=tier1_depth,tier2_depth
        else:  # indels
            alt_counts_n = rec.calls[0].data['TIR']  # TIR=tier1_depth,tier2_depth
            alt_counts_t = rec.calls[1].data['TIR']
        DP_n=rec.calls[0].data["DP"]
        DP_t=rec.calls[1].data["DP"]
        if DP_n is not None and DP_t is not None:
            with np.errstate(divide='ignore', invalid='ignore'):  # ignore division by zero and put AF=.0
                #alt_n = alt_counts_n[0]/DP_n
                #alt_t = alt_counts_t[0]/DP_t
                af_n = np.true_divide(alt_counts_n[0], DP_n)
                af_t = np.true_divide(alt_counts_t[0], DP_t)
                rec.add_format('AF',0)
                rec.calls[0].data["AF"]= [round(af_n,5)]
                rec.calls[1].data["AF"]= [round(af_t,5)]
        normal_gt, tumor_gt= _tumor_normal_genotypes(rec.REF,rec.ALT[0].value,rec.INFO)
        rec.add_format('GT',"1/0")
        rec.calls[0].data["GT"]=normal_gt
        rec.calls[1].data["GT"]=tumor_gt
        writer.write_record(rec)

if __name__ == '__main__':
    filename = sys.argv[1]
    outname = sys.argv[2]
    _af_annotate_and_filter(filename, outname)
    
