#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 17:43:58 2020

Search for ClinVar variants occurring within the motifs of study. 
Variants are grouped according to the reported Clinical Significance (ClinSig).
Only missense variants with Criteria Provided and no conflicting 
interpretations are considered. 
The variant position in the protein sequence is changed to the position
relative to the motif.

NOTE: the folder were the script is run should also contain the file "ClinVar-todo.csv"


@author: Mariano Martín
"""
import pandas as pd


def clinvar_opening():  
    ''' Open and keep the ClinVar variants to be analysed.'''
    variants = pd.read_csv('ClinVar/ClinVar.csv', sep='\t', index_col=0)
    # Keep variants with criteria provided
    status = ['criteria provided, single submitter',
              'criteria provided, multiple submitters, no conflicts',
              'reviewed by expert panel']
    variants = variants[variants.Status.isin(status)]
    variants = variants.dropna(subset=['ProteinPosition'])
    variants.drop(variants[variants['ProteinChange'].str.contains('=')].index, inplace=True)
    variants.drop(variants[variants['ProteinChange'].str.contains('Ter')].index, inplace=True)
    variants.drop(variants[variants['ProteinChange'].str.contains('fs')].index, inplace=True)
    variants.drop(variants[variants['Clinsig'] == 'Uncertain significance'].index , inplace=True)
    variants.reset_index(inplace=True, drop=True)
    # Determine the type of the columns information
    variants['ProteinPosition'] = variants.ProteinPosition.astype(int)
    variants['IDvariant'] = variants.IDvariant.astype(int)
    # Keep unique variants
    uniques = variants.IDvariant.unique() 
    variants = variants[variants.IDvariant.isin(uniques)]
    return variants


def clinvar_search(g, df_motif, variants, d):    
    ''' Searchs ClinVar variants ocurring within the motif of a particular
        protein (uniprotid: g). Outputs a dataframe with the variants found.
        Variant location is re-positioned according to the relative motif 
        position.
        In: UniProtID (g), file of the motif hits (motifs), file of ClinVar
        variant (variants),
        Out: Dataframe with variants occurring within the protein motif.'''
    ix = df_motif.UniProtID == g
    variants1 = variants[variants['Gene'] == df_motif.loc[ix, 'Gene'].values[0]]
    if not variants1.empty:
        motifs1 = df_motif[ix]
        for z, n in enumerate(motifs1['Start'].tolist()):
            variants2 = variants1[(variants1.ProteinPosition >= n) & (variants1.ProteinPosition < n+len(motifs1.Motif.iloc[z]))]
            if not variants2.empty:
                pat_var_mot = []
                ben_var_mot = []
                varid = []
                for i in variants2.index: 
                    for e in range(len(motifs1.Motif.iloc[z])): 
                        # Check variant-motif position and amino acid
                        if (variants2.ProteinPosition.loc[i] == n+e and 
                        variants2.ProteinChange.loc[i][:3] == d[motifs1.Motif.tolist()[z][e]]): 
                            change = variants2.ProteinChange.loc[i]
                            if (variants2.Clinsig.loc[i] == 'Pathogenic' or 
                            variants2.Clinsig.loc[i] == 'Pathogenic/Likely pathogenic' or
                            variants2.Clinsig.loc[i] == 'Likely pathogenic'):
                                pat_var_mot.append(change[:3]+str(e+1)+change[-3:])
                                varid.append(variants2.IDvariant.loc[i])
                            elif (variants2.Clinsig.loc[i] == 'Benign' or 
                            variants2.Clinsig.loc[i] == 'Benign/Likely benign' or
                            variants2.Clinsig.loc[i] == 'Likely benign'):
                                ben_var_mot.append(change[:3]+str(e+1)+change[-3:])
                                varid.append(variants2.IDvariant.loc[i])
                df = pd.DataFrame({'Pathogenic':pd.Series(pat_var_mot, dtype=str), 
                                  'Benign':pd.Series(ben_var_mot, dtype=str)})
                df['ID'] = g+'_'+str(n)
                return df


def ClinVarMotif(df_motif, motif_name):
    ''' Creates the final dataframe of ClinVar variants that 
    occur within the motif'''
    # Amino acids code dictionary
    d = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
         'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
         'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
         'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}
    variants = clinvar_opening()
    UniMotifs = df_motif.UniProtID.tolist()
    dfs = []
    for g in UniMotifs:
        dfs.append(clinvar_search(g, df_motif, variants, d))
    df_final = pd.concat(dfs)
    df_final.drop_duplicates(keep='first', inplace=True)      
    df_final.to_csv('../'+motif_name+'_motif/'+motif_name+'_motif_ClinVarVariants.csv')
    print(f'Total number of ClinVar variants affecting', motif_name, 'motif:', len(df_final))      
    return df_final














        