#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 08:14:16 2024

@author: tanvibansal
"""
import gzip
import shutil
import os

#function to create a compressed copy of an existing file
def compress_file(input_file_path, output_file_path):
    # Open the input file in read-binary mode and the output file in gzip write-binary mode
    with open(input_file_path, 'rb') as f_in:
        with gzip.open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

#function to read a compressed file directly into memory
def read_compressed_file(filepath):
    with gzip.open('/home/joe/file.txt.gz', 'rb') as f:
        file_content = f.read()
    return file_content


### from here below is specifically for my machine
root_fp = '/Users/tanvibansal/Documents/GitHub/Capstone/'
for c in range(21,23):
    branch_fp = 'ld_score_outputs/chr%s/'%(c)
    input_fps = root_fp + branch_fp + "extracted_plink_triplet/"
    input_fp = os.listdir(input_fps)
    output_fp = root_fp + branch_fp + "extracted_plink_triplet_compr"

    for i in input_fp:
        path = input_fps + i
        compress_file(path, root_fp + "filtered_plink_triplet/%s.gz"%(i))
        print("%s compression complete"%(i))
        
    print("chromosome %s complete"%(c))
