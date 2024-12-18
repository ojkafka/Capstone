#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 12:22:44 2024

@author: tanvibansal
"""

import glob
import pandas as pd

root_fp = "/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_results/"
ref_parent_fp = "reference_files/"
annot_parent_fp = "c2.0/"

chrom = 1 

for chrom in range(1,23):
    #read in annotation file
    annot_fp = glob.glob(root_fp + annot_parent_fp + "SNP_%s.annot"%(chrom))[0]
    annot_df = pd.read_csv(annot_fp,sep="\t")
    
    #read in .bim file
    bim_fp = glob.glob(root_fp + ref_parent_fp + "%s_filt.bim"%(chrom))[0]
    bim_df = pd.read_csv(bim_fp,sep="\t",header=None,names=["CHR","SNP","CM","BP","A1","A2"])
    
    #grab SNP order list from .bim file 
    SNP_ordering = bim_df.SNP.values
    
    #order the annotation df
    annot_df_cols = annot_df.columns.values
    annot_df_ordered = annot_df.set_index("SNP").loc[SNP_ordering].reset_index()[annot_df_cols]
    
    #write to disk
    annot_df_ordered.to_csv(annot_fp,sep="\t")
