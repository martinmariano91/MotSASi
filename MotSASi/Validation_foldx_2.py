#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 08:30:41 2021

@author: mariano
"""

import pandas as pd
#import sys

#mot = sys.argv[2]
#motif_name = sys.argv[3]

#motifs_names = os.listdir('../MotListos_1er_etapa/')

motifs_names = [('DOC_MAPK_JIP1_4_motif', 'RK.P.^P.^P.L.x.LIVMF'),
                ('LIG_PCNA_PIPBox_1_motif', 'QM.x.^FHWY.LIVM.^P.^PFWYMLIV.FYHL.FYWL.x.x'),           
                ('DOC_PP2A_B56_1_motif', 'LMFYWIC.x.x.IVLWC.x.E'),                
                ('LIG_NRBOX_motif', '^P.L.^P.^P.L.L.^P'),
                ('DOC_MAPK_NFAT4_5_motif', 'RK.^P.^P.LIM.x.L.x.LIVMF.x'),
                ('DEG_APCC_KENBOX_2_motif', 'x.K.E.N.x'),
                ('DOC_ANK_TNKS_1_motif', 'x.R.x.x.PGAV.DEIP.G.x'),
                ('LIG_PAM2_1_motif_1', 'x.x.LFP.NS.PIVTAFL.x.A.x.x.FY.x.PYLF.x'),	                
                ('LIG_PTAP_UEV_1_motif', 'x.P.TS.A.P.x'),
                ('LIG_SH3_2_motif', 'P.x.x.P.x.KR'),
                ('MOD_CDK_SPxK_1_motif', 'x.x.x.ST.P.x.KR'),
                ('DOC_MAPK_MEF2A_6_motif', 'RK.x.x.x.LIVMP.x.LIV.x.LIVMF'),                               
                ('DOC_CYCLIN_RxL_1_motif_2', 'x.^EDWNSG.^D.RK.^D.L.FL.EDST'),
                ('LIG_CaM_IQ_9_motif_1', 'ACLIVTM.^P.^P.ILVMFCT.Q.^P.^P.^P.RK.^P.^P.^P.^P.RKQ.^P.^P'),
                ('LIG_CaM_IQ_9_motif_2', 'ACLIVTM.^P.^P.ILVMFCT.Q.^P.^P.^P.RK.^P.^P.^P.^P.^P.RKQ.^P.^P'),
                ('LIG_PDZ_Class_1_motif', 'ST.x.VIL'),
                ('TRG_diLeu_motif', 'EDQR.x.x.x.L.VIL')]
                
trues_list = []
falses_list = []
conflicts_list = []
motifs_names_list = []

for motif_name_folder, mot in motifs_names:
    #print('*'*50)
    #print(motif_name_folder)
    mot = mot.split('.')
    if motif_name_folder.endswith('motif'):
        motif_name = motif_name_folder[:-6]        
    else:
        motif_name = motif_name_folder[:-8]

    pos_x = [i for i, x in enumerate(mot) if x.startswith('^') or x == 'x']
    
    d = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
         'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
         'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
         'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}
        
    foldx = pd.read_csv('../MotListos_1er_etapa/'+motif_name_folder+'/pdbs/MeanFoldXMatrix_'+motif_name+'.csv', index_col=0)
    
    pathogenics = []
    benigns = []
    for i in range(len(mot)): 
        row_mean = foldx.iloc[i].mean()
        row_std = foldx.iloc[i].std()
        for e in foldx.columns:
            if i in pos_x: # check if is a flexible position
                # Most of the changes at this position will be benign
                if foldx.iloc[i][e] > row_mean+row_std: 
                    # if the ddG z-score is >1, the change is pathogenic
                    for x in d.values():
                        pathogenics.append(x+str(i+1)+d[e])
                else:
                    # if the ddG z-score is <1, the change is benign
                    for x in d.values():
                        benigns.append(x+str(i+1)+d[e])
            else:                  
                # Most of the changes at this position will be pathogenic
                if foldx.iloc[i][e] > row_mean-row_std:
                    # if the ddG z-score is >-1, the change is pathogenic
                    aa = [x for x in mot[i] if x != '^']
                    for z in aa:                    
                        pathogenics.append(d[z]+str(i+1)+d[e])                                                    
                else:             
                    # if the ddG z-score is <-1, the change is benign
                    aa = [x for x in mot[i]]
                    for z in aa:
                        benigns.append(d[z]+str(i+1)+d[e])
                    
    #%%
                            
    df = pd.read_csv('../MotListos_1er_etapa/'+motif_name_folder+'/'+motif_name+'_motif.csv', index_col=0, sep='\t')
    positives = set(df.loc[df.ELM_Positives == 1, 'ID_motif'].tolist())
    
    clinvar = pd.read_csv('../MotListos_1er_etapa/'+motif_name_folder+'/'+motif_name+'_motif_ClinVarVariants.csv', index_col=0)
    clinvar['MotifPosition'] = clinvar.MotifPosition.astype(int)
    clinvar['ID_motif'] = clinvar.ID+'_'+clinvar.MotifPosition.astype(str)
    clinvar = clinvar[clinvar.ID_motif.isin(positives)]
    
    gnomad = pd.read_csv('../MotListos_1er_etapa/'+motif_name_folder+'/'+motif_name+'_motif_GnomADVariants.csv', index_col=0)
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
    
    trues_list.append(len(true_motifs))
    falses_list.append(len(false_motifs))
    conflicts_list.append(len(conflict_motifs))
    motifs_names_list.append(motif_name)
    
    with open('../MotListos_1er_etapa/'+motif_name_folder+'/'+motif_name+'_Validation_trues2.txt', "w") as output:
        output.write(str(true_motifs['ID_motif'].tolist()))
        
    with open('../MotListos_1er_etapa/'+motif_name_folder+'/'+motif_name+'_Validation_falses2.txt', "w") as output:
        output.write(str(false_motifs['ID_motif'].tolist()))
        
    with open('../MotListos_1er_etapa/'+motif_name_folder+'/'+motif_name+'_Validation_conflicts2.txt', "w") as output:
        output.write(str(conflict_motifs))
    
    #print(f'Corroborated with FoldX Matrix:', len(true_motifs))
    #print(f'Discarded with FoldX Matrix:', len(false_motifs))
    #print(f'Conflicts (have some variants in accordance and some not):', len(conflict_motifs))
    
final = pd.DataFrame()
final['Motif'] = motifs_names_list
final['Trues'] = trues_list
final['Falses'] = falses_list
final['Conflicts'] = conflicts_list

print(final)