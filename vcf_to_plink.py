#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 09:09:09 2024

@author: tanvibansal
"""


chrom = 19
# =============================================================================
#STEP 1: PLINK CONVERSION OF 1000 GENOMES VCF FILE TO PLINK TRIPLET
# =============================================================================
eur_ind_fp = root_fp + "eur_individuals.txt"
vcf_fp = root_fp + "ld_score_outputs/chr%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"%(chrom,chrom)
vcf_to_plink_out_fp = root_fp + "ld_score_outputs/chr%s/converted_plink_triplet/%s_1000g"%(chrom,chrom)

eur_individuals = pd.read_csv(eur_ind_fp,sep="\t")
plink_conversion_cmd = "plink --vcf %s --keep %s --make-bed --out %s"%(vcf_fp, eur_ind_fp, vcf_to_plink_out_fp)
print(plink_conversion_cmd)
#RUN THE PLINK CONVERSION COMMAND IN THE PYTHON TERMINAL FOR THE LDSC ENVRIONMENT
