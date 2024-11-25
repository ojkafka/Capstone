#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:44:57 2024

@author: tanvibansal
"""

import pandas as pd 
import glob
import gzip
import shutil 

name_mapping_key = pd.read_csv("/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_test/maf01.kgp.snps",sep="\t")

#conver gwas sumstats naming convention
gwas_sumstats_old = pd.read_csv("/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_test/GWAS_SumStats.sumstats.txt",sep=" ")
gwas_sumstats_old = gwas_sumstats_old.rename(columns = {"SNP":"rsid"}).set_index("rsid")
gwas_sumstats_old = gwas_sumstats_old.join(name_mapping_key.set_index("rsid"))
gwas_sumstats = gwas_sumstats_old.reset_index().drop(columns = "rsid")[["SNP","A1","A2","Z","N"]]

#write to text file 
gwas_sumstats.to_csv("gwas_sumstats.sumstats",sep="\t")

#compress text file
with open("/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_test/gwas_sumstats.sumstats", 'rt') as f_in:
    with gzip.open("/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_test/gwas_sumstats.sumstats.gz", 'wt') as f_out:
        shutil.copyfileobj(f_in, f_out)

#convert weights naming convention
weight_fps = glob.glob("/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_test/weights_hm3_no_hla/*")
weight_fps = [i for i in weight_fps if i[-2:] != "gz"]

fp = weight_fps[0]
for fp in weight_fps:
    weight_old = pd.read_csv(fp,sep="\t").rename(columns={"SNP":"rsid"}).set_index("rsid").join(name_mapping_key.set_index("rsid"))
    weight = weight_old.reset_index().drop(columns="rsid")[["CHR","SNP","BP","CM","MAF","L2"]]
    
    #write to text file
    weight.to_csv(fp,sep="\t")
    
    #compress text file
    with open(fp, 'rt') as f_in:
        with gzip.open(fp + ".gz", 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)

#convert frq files naming convention 
frq_fps = glob.glob("/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_test/frq_files/*")

fp = frq_fps[0]
for fp in frq_fps:
    frq_old = pd.read_csv(fp,sep="\t").rename(columns={"SNP":"rsid"}).set_index("rsid").join(name_mapping_key.set_index("rsid"))
    frq = frq_old.reset_index().drop(columns="rsid")[["CHR","SNP","A1","A2","MAF","NCHROBS"]]
    
    #write to text file
    frq.to_csv(fp,sep="\t")
    
#remove index 
weight_fps = glob.glob("/Users/tanvibansal/Documents/GitHub/Capstone/unified_pipeline_test/weights_hm3_no_hla/*")
weight_fps = [i for i in weight_fps if i[-2:] != "gz"]

for fp in weight_fps:
   # pd.read_csv(fp,sep="\t")[["CHR","SNP","BP","CM","MAF","L2"]].to_csv(fp,sep="\t",index=False)
    with open(fp, 'rt') as f_in:
        with gzip.open(fp + ".gz", 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)

    
    
    