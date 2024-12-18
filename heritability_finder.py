#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 13:18:43 2024

@author: tanvibansal
"""

import subprocess
import glob 
import time, sys, traceback, argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref_file_root', default='', type=str,
    help='Root path for all reference files')
parser.add_argument('--annot_root', default='', type=str,
    help='Root path for all annotation files and outputs for ld scores and heritability')
parser.add_argument('--annot_prefix', default='chr',
    help='Prefix for annotation files - will be used for ld score outputs and heritability outputs')

args = parser.parse_args()

# Define the lists of arguments for each run
ref_file_root = args.ref_file_root
annot_root = args.annot_root 
annot_prefix = args.annot_prefix

bfile_list =  [ref_file_root + "%s_filt"%(i) for i in range(1,23)] # Replace with actual bfile values
annot_list = [annot_root + annot_prefix + "%s.annot"%(i) for i in range(1,23)]  # Replace with actual annot values
out_list = [annot_root + annot_prefix + str(i) for i in range(1,23)]          # Replace with actual out values

# Check that the lists are the same length
if not (len(bfile_list) == len(annot_list) == len(out_list)):
    raise ValueError("The bfile, annot, and out lists must have the same length.")

# Iterate over the chromosomes and order the annotation file in the same SNP order as the .bim file



# Iterate over the chromosomes and compute the ld scores
for bfile, annot, out in zip(bfile_list, annot_list, out_list):
    command = [
        "python", "ldsc.py",
        "--l2",
        "--bfile", bfile,
        "--ld-wind-cm", '1',
        "--annot", annot,
        "--out", out
    ]
    print "Running command: %s"%(command)
    try:
        subprocess.call(command)
        print "Successfully completed: %s"%(out)
    except subprocess.CalledProcessError as e:
        print "Error while running command for %s: %s"%(out, e)

# Call the partitioned heritability estimation
command = ["python", "ldsc.py",
"--h2", ref_file_root + "GWAS_SumStats.sumstats.gz" ,
"--ref-ld-chr", annot_root + annot_prefix ,
"--out", annot_root + annot_prefix + "_heritability" ,
"--overlap-annot",  
"--frqfile-chr", ref_file_root + "frq_files/WB." ,
"--w-ld-chr", ref_file_root + "weights_hm3_no_hla/weights."
] 

try:
    subprocess.call(command)
    print "Successfully completed: heritability estimation"
except subprocess.CalledProcessError as e:
    print "Error while running command for heritability estimation: %s"%(e)
