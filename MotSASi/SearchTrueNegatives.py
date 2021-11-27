#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 09:21:50 2021

@author: mariano
"""


import pandas as pd

motif_name = 'LIG_NRBOX'

d = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
     'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
     'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
     'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}
    
df = pd.read_csv('../MotListos_1er_etapa/'+motif_name+'_motif/'+motif_name+'_motif.csv', index_col=0, sep='\t')

positives_motifs = set(df.loc[df.ELM_Positives == 1, 'UniProtID'].tolist())

potential_negatives = df.loc[(df.UniProtID.isin(positives_motifs)) & (df.ELM_Positives == 0)]
