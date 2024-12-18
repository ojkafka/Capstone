#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 08:14:16 2024

@author: tanvibansal
"""
import gzip
import shutil
import os
import pandas as pd 
#function to create a compressed copy of an existing file
def compress_file(input_file_path, output_file_path):
    # Open the input file in read-binary mode and the output file in gzip write-binary mode
    with open(input_file_path, 'rb') as f_in:
        with gzip.open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

#function to read a compressed file directly into memory
def read_compressed_file(filepath):
    with gzip.open(filepath, 'rb') as f:
        file_content = pd.read_csv(f, header=None, sep = "\t", names=["CHR","SNP","CM","POS","A1","A2"])
    return file_content


### from here below is specifically for my machine
#root_fp = '/Users/tanvibansal/Documents/GitHub/Capstone/'
#for c in range(1,23):
 #   branch_fp = 'ld_score_outputs/chr%s/'%(c)
  #  input_fps = root_fp + branch_fp + "ld_scores/"
   # input_fp = os.listdir(input_fps)

#    for i in input_fp:
 #       path_in = input_fps + i
  #      path_out = root_fp + "ld_scores_c1.0/%s"%(i)
   #     if i[-2:] == "gz":
    #        compress_file(path_in, path_out)
     #   else:
      #      compress_file(path_in, path_out + ".gz")
       # print("%s compression complete"%(i))
        
    #print("chromosome %s complete"%(c))
