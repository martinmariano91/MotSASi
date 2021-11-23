#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 08:30:41 2021

@author: mariano
"""

import pandas as pd
import sys

motif = sys.argv[1]
mot = sys.argv[2]
motif_name = sys.argv[3]

#motif_name = 'DOC_PP2A_B56_1'
#mot = 'L.x.x.IVLWC.x.E'
mot = mot.split('.')

pos_x = [i for i, x in enumerate(mot) if x.startswith('^') or x == 'x']

d = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
     'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
     'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
     'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}
    
foldx = pd.read_csv('../MotListos_1er_etapa/'+motif_name+'_motif/pdbs/MeanFoldXMatrix_'+motif_name+'.csv', index_col=0)

pathogenics = []
benigns = []
for i in range(len(mot)): 
    for e in foldx.columns:
        if foldx.iloc[i][e] > 2: # patogenicos           
            if i in pos_x:
                for x in d.values():
                    pathogenics.append(x+str(i+1)+d[e])
            else: 
                aa = [x for x in mot[i] if x != '^']
                for z in aa:
                    pathogenics.append(d[z]+str(i+1)+d[e])            
        elif foldx.iloc[i][e] < 1: # benignos
            if i in pos_x:
                for x in d.values():
                    benigns.append(x+str(i+1)+d[e])
            else:             
                aa = [x for x in mot[i]]
                for z in aa:
                    benigns.append(d[z]+str(i+1)+d[e])
        else:
            row_mean = foldx.iloc[i].mean()
            if foldx.iloc[i][e] > row_mean:
                if i in pos_x:
                    for x in d.values():
                        pathogenics.append(x+str(i+1)+d[e])
                else:
                    aa = [x for x in mot[i]]
                    for z in aa:
                        pathogenics.append(d[z]+str(i+1)+d[e])
            else:
                if i in pos_x:
                    for x in d.values():
                        benigns.append(x+str(i+1)+d[e])
                else:
                    aa = [x for x in mot[i]]
                    for z in aa:
                        benigns.append(d[z]+str(i+1)+d[e])
                        
#%%
                        
df = pd.read_csv('../MotListos_1er_etapa/'+motif_name+'_motif/'+motif_name+'_motif.csv', index_col=0, sep='\t')
positives = set(df.loc[df.ELM_Positives == 1, 'ID_motif'].tolist())

clinvar = pd.read_csv('../MotListos_1er_etapa/'+motif_name+'_motif/'+motif_name+'_motif_ClinVarVariants.csv', index_col=0)
clinvar['ID_motif'] = clinvar.ID+'_'+clinvar.MotifPosition.astype(str)
clinvar = clinvar[clinvar.ID_motif.isin(positives)]

gnomad = pd.read_csv('../MotListos_1er_etapa/'+motif_name+'_motif/'+motif_name+'_motif_GnomADVariants.csv', index_col=0)
gnomad['ID_motif'] = gnomad.ID+'_'+gnomad.MotifPosition.astype(str)
gnomad = gnomad[gnomad.ID_motif.isin(positives)]

clinvar_false = clinvar.loc[((clinvar.Benign != 1) & (clinvar.Variant.isin(benigns))) | 
                            ((clinvar.Pathogenic != 1) & (clinvar.Variant.isin(pathogenics))),
                            ['ID_motif']]
gnomad_false = gnomad.loc[gnomad.Benign.isin(pathogenics), ['ID_motif']]

false_motifs = pd.concat([clinvar_false, gnomad_false])
false_motifs.drop_duplicates(inplace=True)
falses = false_motifs.ID_motif.tolist()


clinvar_true = clinvar.loc[((clinvar.Benign == 1) & (clinvar.Variant.isin(benigns))) | 
                           ((clinvar.Pathogenic == 1) & (clinvar.Variant.isin(pathogenics))),
                           ['ID_motif']]
gnomad_true = gnomad.loc[gnomad.Benign.isin(benigns), ['ID_motif']]

true_motifs = pd.concat([clinvar_true, gnomad_true])
true_motifs.drop_duplicates(inplace=True)
trues = true_motifs.ID_motif.tolist()
true_motifs.drop(true_motifs[true_motifs.ID_motif.isin(falses)].index, inplace=True)

conflict_motifs = true_motifs.loc[true_motifs.ID_motif.isin(falses), 'ID_motif'].tolist()+false_motifs.loc[false_motifs.ID_motif.isin(trues), 'ID_motif'].tolist()

with open('../MotListos_1er_etapa/'+motif_name+'_motif/'+motif_name+'_Validation_trues.txt', "w") as output:
    output.write(str(true_motifs['ID_motif'].tolist()))
    
with open('../MotListos_1er_etapa/'+motif_name+'_motif/'+motif_name+'_Validation_falses.txt', "w") as output:
    output.write(str(false_motifs['ID_motif'].tolist()))
    
with open('../MotListos_1er_etapa/'+motif_name+'_motif/'+motif_name+'_Validation_conflicts.txt', "w") as output:
    output.write(str(conflict_motifs))

print(f'Discarded with FoldX Matrix:', len(false_motifs))
print(f'Corroborated with FoldX Matrix:', len(true_motifs))
print(f'Conflicts (have some variants in accordance and some not):', len(conflict_motifs))
                        