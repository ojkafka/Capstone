#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:00:48 2024

@author: tanvibansal
"""
import pandas as pd
import gzip 
import numpy as np
import os 
os.chdir("/Users/tanvibansal/Documents/GitHub/Capstone/")
import collections 
from pandarallel import pandarallel
import string
import gc
from compress_and_uncompress_fns import read_compressed_file

# =============================================================================
#STEP 0: SET GLOBAL PARAMS AND READ IN DATA FILES TO BE USED THROUGHOUT SCRIPT
# =============================================================================

root_fp = "/Users/tanvibansal/Documents/GitHub/Capstone/"

gwas_fp = root_fp + 'HEIGHT_irnt.gwas.imputed_v3.both_sexes.tsv'
gwas_df = pd.read_csv(root_fp + "HEIGHT_irnt.gwas.imputed_v3.both_sexes.tsv",sep="\t")

gtf_fp = root_fp + 'genes.protein_coding.v39.gtf'
gtf_df = pd.read_csv(gtf_fp, sep='\t', header=0)
gtf_df["chr"] = gtf_df.chr.apply(lambda x: int(x.split("chr")[1]))

s_het_fp = root_fp + 's_het_info.xlsx'
s_het_df = pd.read_excel(s_het_fp,sheet_name=1).rename(columns={"ensg":"GeneSymbol"}).set_index("GeneSymbol")
s_het_df["chr"] = s_het_df.chrom.apply(lambda x: x.split("chr")[1])
s_het_df = s_het_df.loc[(s_het_df.chr != "X") & (s_het_df.chr != "Y")]
s_het_df["chr"] = s_het_df["chr"].astype("int")

chrom = 1
branch_fp = "ld_score_outputs/chr%s/"%(chrom)
cm_fp = root_fp + "ld_score_outputs/CM_values/chr%s.bim"%(chrom)

# =============================================================================
#STEP 1: GENERATE ANNOTATION FILE
# =============================================================================

def generate_annotation_file(chrom,gwas_df,gtf_df, s_het_df,cm_fp):
    gwas_chr = gwas_df.loc[gwas_df.CHR == chrom]                       
    gtf_chr = gtf_df.loc[gtf_df.chr == chrom]
    s_het_chr = s_het_df.loc[s_het_df.chr == chrom].drop(columns="chr")
    cm_df = pd.read_csv(cm_fp,sep="\t",header=None,names=["CHR","SNP","CM","POS","A1","A2"])
    
    #drop the genes that are not present in the s_het dataframe
    gtf_chr_tidy = s_het_chr.join(gtf_chr.set_index("GeneSymbol"),how="inner")[["chr","start","end","gene"]]
    
    #find the 5 nearest genes by absolute difference in mutation position and gene start position
    def nearest_gene_finder(m):
        #get the absolute l1 distance 
        dists = np.abs(gtf_chr_tidy.loc[gtf_chr_tidy.chr == m.CHR,"start"] - m.POS)
        top_5 = dists.sort_values()[0:5]
        out = dict(m)
        for i in range(5):
            out["GeneSymbol.f%s"%(i+1)] = top_5.index[i]
            out["dist.f%s"%(i+1)] = 1/1000 if top_5.values[i] == 0 else top_5.values[i]/1000
        return(pd.Series(out))
    pandarallel.initialize(progress_bar=True, nb_workers=6)
    nearest_genes_chr = gwas_chr.parallel_apply(lambda m: nearest_gene_finder(m),axis=1)
    gc.collect()
    
    #join s_het values by GeneSymbol and weight by inverse distance
    annot_df_prep = nearest_genes_chr.copy(deep=True)
    for i in range(1,6):
        s_het_temp = s_het_chr[["post_mean"]].reset_index().rename(columns = {"post_mean":"s_het.f%s"%(i),"GeneSymbol":"GeneSymbol.f%s"%(i)}).set_index("GeneSymbol.f%s"%(i))
        annot_df_prep = annot_df_prep.reset_index().set_index("GeneSymbol.f%s"%(i)).join(s_het_temp)
        annot_df_prep["s_het_w.f%s"%(i)] = annot_df_prep["s_het.f%s"%(i)]/annot_df_prep["dist.f%s"%(i)]
    annot_df_prep = annot_df_prep.reset_index()
    
    #prep the annotation dataframe formatting for LD score computation
    names = ["CHR","BP", "SNP", "CM"] + ["s_het_w.f%s"%(i) for i in range(1,6)]
    
    annot_df_prep = annot_df_prep.rename(columns={"POS":"BP","variant":"SNP"})
    annot_df_prep = annot_df_prep.merge(cm_df[["SNP","CM"]],on="SNP",how="left")
    annot_df_prep.CM = annot_df_prep.CM.fillna(1.0)
    annot_df_prep = annot_df_prep[names]
    return annot_df_prep, names
annot_df_prep,names = generate_annotation_file(chrom,gwas_df,gtf_df, s_het_df,cm_fp)

# =============================================================================
#STEP 2: READ IN FILTERED .BIM FILE AND FILTER/ORDER THE ANNOTATION FILE TO MATCH SNP ORDERING
# =============================================================================
def filter_annotation_file(chrom, root_fp, branch_fp, annot_df_prep):
    filtered_bim_fp = root_fp + "filtered_plink_triplets/%s_filt.bim.gz"%(chrom)
    filtered_bim = read_compressed_file(filtered_bim_fp)
    
    filtered_snp_ids = filtered_bim.SNP.values
    annot_df = annot_df_prep.set_index("SNP").loc[filtered_snp_ids].reset_index()[names]
    
    annot_fp = root_fp + branch_fp + "s_het_w_chr%s.annot"%(chrom)
    annot_df.to_csv(annot_fp,sep="\t",index=False)
    print("chr %s annotation file ordered and written to disk"%(chrom))
filter_annotation_file(chrom, root_fp, branch_fp, annot_df_prep)

# =============================================================================
#STEP 3: STEPS 1 & 2 FOR ALL CHROMOSOMES
# =============================================================================
for chrom in range(2,23):
    branch_fp = "ld_score_outputs/chr%s/"%(chrom)
    cm_fp = root_fp + "ld_score_outputs/CM_values/chr%s.bim"%(chrom)
    annot_df_prep,names = generate_annotation_file(chrom,gwas_df,gtf_df, s_het_df,cm_fp)
    filter_annotation_file(chrom, root_fp, branch_fp, annot_df_prep)

# =============================================================================
#STEP 7: RUN THE LD SCORE COMPUTATION
# =============================================================================
def prep_ld_score_compute(chrom, root_fp, branch_fp):
    out_fp = root_fp + branch_fp + "extracted_plink_triplet/%s_filt"%(chrom)
    annot_fp = root_fp + branch_fp + "s_het_w_chr%s.annot"%(chrom)
    ld_score_fp = root_fp + branch_fp + "ld_scores/s_het_w_chr%s"%(chrom)
    ld_score_cmd = "python ldsc.py --l2 --bfile %s --ld-wind-cm 1 --annot %s --out %s"%(out_fp,annot_fp,ld_score_fp)
    print(ld_score_cmd)
for chrom in range(1,22):
    branch_fp = "ld_score_outputs/chr%s/"%(chrom)
    prep_ld_score_compute(chrom, root_fp, branch_fp)
