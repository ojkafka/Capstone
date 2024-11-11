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
# =============================================================================
#STEP 0: SET GLOBAL PARAMS AND READ IN DATA FILES TO BE USED THROUGHOUT SCRIPT
# =============================================================================

root_fp = "/Users/tanvibansal/Documents/GitHub/Capstone/"

gwas_fp = root_fp + 'gwas_df_tidy.parquet'
gwas_df = pd.read_parquet(gwas_fp)

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
# =============================================================================
#STEP 1: SNP IDENTIFICATION OF 1000 GENOMES .BIM FILE WITH GWAS ID CONVENTION
# =============================================================================
bim_fp = root_fp + branch_fp + "converted_plink_triplet/%s_1000g.bim"%(chrom)
cm_fp = root_fp + "ld_score_outputs/CM_values/chr%s.bim"%(chrom)

def modify_bim_file(chrom, bim_fp, gwas_df,cm_fp):
#get.bim file needing SNP IDs and GWAS data set containing SNP IDs
    gwas_chr = gwas_df.loc[gwas_df.CHR == chrom]
    bim = pd.read_csv(bim_fp,sep="\t",header=None,names=["CHR","SNP","CM","POS","A1","A2"])

    #get the gwas SNP matches based on forwad and reverse allele matches on position
    def SNP_id_finder(bim,gwas_chr):
        bim_temp = bim[["POS","A1","A2"]].reset_index().set_index(["POS","A1","A2"]).rename(columns={"index":"index.bim"})
        
        #do ref = a1, alt = a2 first as "_f" then reversed and join to one large dataframe 
        gwas_chr_temp_f = gwas_chr[["POS","REF","ALT"]].reset_index().rename(columns={"index":"index.gwas.f","REF":"A1","ALT":"A2"}).set_index(["POS","A1","A2"])
        gwas_chr_temp_r = gwas_chr[["POS","REF","ALT"]].reset_index().rename(columns={"index":"index.gwas.r","REF":"A2","ALT":"A1"}).set_index(["POS","A1","A2"])
        bim_to_gwas_map = bim_temp.join(gwas_chr_temp_f,how="left").join(gwas_chr_temp_r,how="left")
        
        #del  gwas_chr_temp_r, gwas_chr_temp_f, bim_temp
    
        #set a new column called "index.gwas" that keeps index.gwas.f if not null, then for remaining null rows populate with index.gwas.r, then drop any remaining null rows
        f = (~bim_to_gwas_map["index.gwas.f"].isnull()).sum()
        r = (~bim_to_gwas_map["index.gwas.r"].isnull()).sum()
        
        order = ["f","r"] if f > r else ["r","f"]
        bim_to_gwas_map_ind = bim_to_gwas_map.copy(deep=True)
        bim_to_gwas_map_ind["index.gwas"] = np.nan
        
        for o in order:
            ind_map_null = bim_to_gwas_map_ind["index.gwas"].isnull()
            ind_gwas_o_null = bim_to_gwas_map_ind["index.gwas.%s"%(o)].isnull()
            
            bim_to_gwas_map_ind.loc[(ind_map_null) & (~ind_gwas_o_null),"index.gwas"] = bim_to_gwas_map_ind.loc[(ind_map_null) & (~ind_gwas_o_null),"index.gwas.%s"%(o)]
        
        #format the key to be easily joinable onto both sets
        bim_to_gwas_map_ind = bim_to_gwas_map_ind.dropna(subset=["index.gwas"])
        
        #check for duplicates gwas SNPs and set to null for lower priority matches 
        n_dup = bim_to_gwas_map_ind["index.gwas"].duplicated()
        if n_dup.sum() > 0:
            bim_to_gwas_map_ind_dup = bim_to_gwas_map_ind.loc[n_dup,"index.gwas"].values
            
            #inds to set to null
            ind_to_null = np.isin(bim_to_gwas_map_ind["index.gwas"], bim_to_gwas_map_ind_dup) & (bim_to_gwas_map_ind["index.gwas"] == bim_to_gwas_map_ind["index.gwas.%s"%(order[1])])
            bim_to_gwas_map_ind.loc[ind_to_null,"index.gwas"] = np.nan
        
        bim_to_gwas= bim_to_gwas_map_ind.reset_index().drop(columns=["POS","A1","A2","index.gwas.f","index.gwas.r"]).set_index("index.bim")#.rename(columns={"index.bim":"Index"}).set_index("Index")
        
        #get gwas SNP IDs and place into bim dataframe 
        snp_id_matches = bim_to_gwas.merge(gwas_chr, left_on="index.gwas",right_index=True)["variant"].rename("SNP")
        bim = bim.drop(columns=["SNP"]).merge(snp_id_matches,how="left",left_index=True,right_index=True)
        
        #generate random ids for rows with no SNP id match
        rand_ids = [''.join(np.random.choice(list((string.ascii_uppercase + string.digits)), size=8)) for _ in  range(bim.SNP.isnull().sum())]
        
        bim.loc[bim.SNP.isnull(),"SNP"] = rand_ids
        return bim
    bim_snp_matched = SNP_id_finder(bim,gwas_chr)
    
    # =============================================================================
    #STEP 4: JOIN CM VALUES ONTO BIM DATAFRAME
    # =============================================================================
    cm_df = pd.read_csv(cm_fp,sep="\t",header=None,names=["CHR","SNP","CM","POS","A1","A2"])
    
    #drop the "CM" column from the bim dataframe and left join the cm dataframe
    bim_snp_matched = bim_snp_matched.drop(columns="CM")
    cm_df["SNP"] = cm_df["SNP"].astype("str")
    bim_snp_matched["SNP"] = bim_snp_matched["SNP"].astype("str")
    bim_snp_cm_matched = bim_snp_matched.merge(cm_df[["CM","SNP"]],on="SNP",how="left")
    
    #fill in remaining null values with 1.0
    bim_snp_cm_matched.loc[bim_snp_cm_matched.CM.isnull(),"CM"] = float(1.0)
    return bim, bim_snp_cm_matched
bim, bim_snp_cm_matched = modify_bim_file(chrom, bim_fp, gwas_df,cm_fp)

# =============================================================================
#STEP 2: GENERATE ANNOTATION FILE AND PREP FORMATTING 
# =============================================================================
#get the gtf and s_het dfs for our chromosome of interest 
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
#STEP 3: EXTRACT UNIQUE SNPS TO PREP FOR PLINK EXTRACTION OF MODIFIED .BIM FILE
# =============================================================================
# find the unique snps common between the annotation file and the .bim prepped file and write to disk
# filter/order the .annot file
# .bim file filter and order will be handled in the plink extraction step next with the unique ids we write out
def prep_plink_extract(chrom, annot_df_prep, bim_snp_cm_matched, branch_fp, bim_fp, bim):
    unique_snp_ids = annot_df_prep.merge(bim_snp_cm_matched,on="SNP")["SNP"].drop_duplicates()
    
    unique_snp_fp = root_fp + branch_fp + "converted_plink_triplet/unique_snps_chr%s.txt"%(chrom)
    unique_snp_ids.to_csv(unique_snp_fp,sep="\t",header=None,index=False)
    
    bim_fp_old =  root_fp + branch_fp + "converted_plink_triplet/%s_1000g_v0.bim"%(chrom)
    bim.to_csv(bim_fp_old, sep="\t",header=None,index=False)
    bim_snp_cm_matched[list(bim.columns.values)].to_csv(bim_fp, sep="\t",header=None,index=False)
    
    out_fp = root_fp + branch_fp + "extracted_plink_triplet/%s_filt"%(chrom)
    plink_extraction_cmd = "plink --bfile %s --extract %s --make-bed --out %s"%(bim_fp.split(".bim")[0], unique_snp_fp, out_fp)
    print(plink_extraction_cmd)
prep_plink_extract(chrom, annot_df_prep, bim_snp_cm_matched, branch_fp, bim_fp, bim)

# =============================================================================
#STEP 4: READ IN FILTERED .BIM FILE AND FILTER/ORDER THE ANNOTATION FILE TO MATCH SNP ORDERING
# =============================================================================
def filter_annotation_file(chrom, root_fp, branch_fp, annot_df_prep):
    filtered_bim_fp = root_fp + branch_fp + "extracted_plink_triplet/%s_filt.bim"%(chrom)
    filtered_bim = pd.read_csv(filtered_bim_fp,sep="\t",header=None,names=["CHR","SNP","CM","POS","A1","A2"])
    
    filtered_snp_ids = filtered_bim.SNP.values
    annot_df = annot_df_prep.set_index("SNP").loc[filtered_snp_ids].reset_index()[names]
    
    annot_fp = root_fp + branch_fp + "s_het_w_chr%s.annot"%(chrom)
    annot_df.to_csv(annot_fp,sep="\t",index=False)
    print("annotation file filtered and written to disk")
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
prep_ld_score_compute(chrom, root_fp, branch_fp)



