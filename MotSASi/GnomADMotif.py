#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:52:43 2020

Search for GnomAD variants occurring within the motifs of study. 
Only missense variants are considered. The variant position in the protein 
sequence is changed to the position relative to the motif.
Variants with -log10(AF) <= 4 are considered as benign.

NOTE: 
1- The folder were the script is run should also contain the files of
   each Chromosome GnomAD variants.
2- This script should be run after "2_ClinVarMotif.py"

@author: Mariano Martín
"""

import pandas as pd
import math

def iterator(IDs, df_motif, variants, d, founds, dfs):
    ''' Search GnomAD variants ocurring within the motif of a particular
        protein. Outputs a dataframe with the variants found.
        Variant location is re-positioned according to the relative motif 
        position.'''
    for g in IDs:    
        ix = df_motif.UniProtID == g
        variants1 = variants[variants['UniProtID'] == g]   
        if not variants1.empty:
            motifs1 = df_motif[ix]
            founds.append(g)
            for z, n in enumerate(motifs1['Start'].tolist()):
                variants2 = variants1[(variants1.ProteinPosition >= n) & (variants1.ProteinPosition < n+len(motifs1.Motif.iloc[z]))]
                if not variants2.empty:
                    ben_var_mot = []
                    varid = []
                    for i in variants2.index: 
                        for e in range(len(motifs1.Motif.iloc[z])):
                            # Check variant-motif position and amino acid
                            if (variants2.ProteinPosition.loc[i] == n+e and
                                variants2.ProteinChange.loc[i][:3] == d[motifs1.Motif.tolist()[z][e]]):
                                change = variants2.ProteinChange.loc[i]
                                if variants2.AF.loc[i] != 0:
                                    if -math.log10(variants2.AF.loc[i]) <= 4:
                                        ben_var_mot.append(change[:3]+str(e+1)+change[-3:])
                                        varid.append(variants2.VariantID.loc[i])
                    df = pd.DataFrame({'Benign':pd.Series(ben_var_mot, dtype=str)})
                    df['ID'] = g+'_'+str(n)
                    df['VariantID'] = varid
                    dfs.append(df)

def GnomADMotif(df_motif, motif_name):
    ''' Iterates over chromosomes and creates the final dataframe
    of GnomAD variants with -log10(AF)<4 that occurr within the motif'''
    chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    d = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
         'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
         'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
         'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}
    UniMotifs = df_motif.UniProtID.tolist()
    dfs = []
    founds = []
    for x in chromosomes:    
        variants = pd.read_csv('GnomAD/'+str(x)+'.csv', sep='\t')
        variants.drop(columns='Unnamed: 0', inplace=True)
        variants = variants[variants['Effect']=='missense_variant']
        variants['ProteinPosition'] = variants['ProteinPosition'].astype(float)
        variants.reset_index(inplace=True, drop=True)
        #print(f'Chromosome ', x)
        IDs = [u for u in UniMotifs if u not in founds]    
        iterator(IDs, df_motif, variants, d, founds, dfs)
    df_final = pd.concat(dfs)
    df_final.drop_duplicates(keep='first', inplace=True)  
    df_final.to_csv('../'+motif_name+'_motif/'+motif_name+'_motif_GnomADVariants.csv')              
    print(f'Total number of GnomAD variants (-log10(AF)<4) affecting', motif_name, 'motif:', len(df_final))      


