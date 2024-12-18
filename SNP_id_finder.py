#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:40:29 2024

@author: tanvibansal
"""

def SNP_id_finder(x, gwas_df):
    x_gwas_df = gwas_df.loc[(gwas_df.CHR == x.CHR) & (gwas_df.POS == x.POS)]
    n_matches = len(x_gwas_df)
    alternate_id = "%s:%s"%(x.POS,x.name)
    if n_matches == 0:
        return alternate_id
    else:
        #allel match check
        allele_match = [collections.Counter(x[['A1','A2']]) == collections.Counter(x_gwas_df.iloc[i][['REF','ALT']]) for i in range(len(x_gwas_df))]
        if sum(allele_match) == 0:
            return alternate_id
        else:
            matches = x_gwas_df.loc[allele_match]
            #if len(matches) == 1:
             #   return matches.variant.values.item()
            #else:
            #allele order check
            allele_ordered = [x.A2 == matches.iloc[i]['REF'] for i in range(len(matches))]
            if sum(allele_ordered) == 0: # this should never happen but handling it here just in case
                return alternate_id
            else:
                matches_remaining = matches.loc[allele_ordered]
                if len(matches_remaining) == 1:
                    return matches_remaining.variant.values.item() 
                else:
                    return alternate_id # this should never happen but handling it here just in case