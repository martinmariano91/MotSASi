#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:47:34 2020

Genera una matriz de sustitución para los aa del motivo.
Colorea según la significancia clíncia reportada (ClinSig).

Hay que darle de comer el archivo con la lista de proteínas con el motivo en estudio y la posición y composición de cada motivo.
También hay que ingresarle el motivo en estudio separando cada una de las posiciones por un punto ('.')

Tener en cuenta esas cosas si se lo va a correr desde la terminal directamente.

@author: mariano
"""

import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def ClinVar_data():
    ''' Open and keep the ClinVar variants to be analysed.'''
    variants = pd.read_csv('ClinVar/ClinVar.csv', sep='\t', index_col=0)
    status=['criteria provided, single submitter',
           'criteria provided, multiple submitters, no conflicts',
           'reviewed by expert panel']
    variants=variants[variants.Status.isin(status)]
    variants = variants.dropna(subset=['ProteinPosition'])
    fil1=variants['ProteinChange'].str.contains('=') #con esta línea y la que sigue me saco de encima las variantes sinónimas
    fil2=variants['ProteinChange'].str.contains('Ter')
    fil3=variants['ProteinChange'].str.contains('fs')
    fil4=variants['Clinsig'] == 'Uncertain significance'
    variants=variants.copy()[~fil1]
    variants=variants.copy()[~fil2]
    variants=variants.copy()[~fil3]
    variants=variants.copy()[~fil4]
    variants['ProteinPosition'] = variants.ProteinPosition.astype(int)
    variants['IDvariant'] = variants.IDvariant.astype(int)
    return variants

def MatSubClinVar(df, mot, positive_ctrl):
    '''Compiles the change of each amino acid in the motif '''
    motifs = df
    UniMotifs = positive_ctrl
    LM = len(mot) #largo del motivo
    variants = ClinVar_data()
    
    Gly_final=[[]*LM for _ in range(LM)]
    Ala_final=[[]*LM for _ in range(LM)]
    Arg_final=[[]*LM for _ in range(LM)]
    Asn_final=[[]*LM for _ in range(LM)]
    Asp_final=[[]*LM for _ in range(LM)]
    Cys_final=[[]*LM for _ in range(LM)]
    Glu_final=[[]*LM for _ in range(LM)]
    Gln_final=[[]*LM for _ in range(LM)]
    His_final=[[]*LM for _ in range(LM)]
    Ile_final=[[]*LM for _ in range(LM)]
    Leu_final=[[]*LM for _ in range(LM)]
    Lys_final=[[]*LM for _ in range(LM)]
    Met_final=[[]*LM for _ in range(LM)]
    Phe_final=[[]*LM for _ in range(LM)]
    Pro_final=[[]*LM for _ in range(LM)]
    Ser_final=[[]*LM for _ in range(LM)]
    Thr_final=[[]*LM for _ in range(LM)]
    Trp_final=[[]*LM for _ in range(LM)]
    Tyr_final=[[]*LM for _ in range(LM)]
    Val_final=[[]*LM for _ in range(LM)]

    lectura_fasta = SeqIO.parse('UniProt-Proteome.fasta','fasta')
    for seq_record in lectura_fasta:
        for g, p in UniMotifs:
            if seq_record.id.split('|')[1] == g:
                ix = motifs.UniProtID == g
                variants1 = variants[variants['Gene']==motifs.loc[ix, 'Gene'].values[0]]
                if not variants1.empty:
                    motifs1 = motifs[ix]
                    protein = pd.DataFrame()
                    sec=[]
                    
                    for i in seq_record.seq:
                        sec.append(i)
                        
                    index=list(range(1, len(seq_record.seq)+1))
                    protein['ResPos']=index
                    protein['Res']=sec
                    protein.set_index('ResPos', inplace=True)
                    protein['Motif']= [0]*len(index)
                    protein['Variants']= [0]*len(index)
                    protein['ClinSig']=[[]*len(index) for _ in range(len(index))]
                    protein['VarNum']=[0]*len(index)
                    
                    d = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
                         'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
                         'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
                         'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}
                    
                    protein.replace({"Res": d}, inplace=True)
                    
                    Gly=[[]*LM for _ in range(LM)]
                    Ala=[[]*LM for _ in range(LM)]
                    Arg=[[]*LM for _ in range(LM)]
                    Asn=[[]*LM for _ in range(LM)]
                    Asp=[[]*LM for _ in range(LM)]
                    Cys=[[]*LM for _ in range(LM)]
                    Glu=[[]*LM for _ in range(LM)]
                    Gln=[[]*LM for _ in range(LM)]
                    His=[[]*LM for _ in range(LM)]
                    Ile=[[]*LM for _ in range(LM)]
                    Leu=[[]*LM for _ in range(LM)]
                    Lys=[[]*LM for _ in range(LM)]
                    Met=[[]*LM for _ in range(LM)]
                    Phe=[[]*LM for _ in range(LM)]
                    Pro=[[]*LM for _ in range(LM)]
                    Ser=[[]*LM for _ in range(LM)]
                    Thr=[[]*LM for _ in range(LM)]
                    Trp=[[]*LM for _ in range(LM)]
                    Tyr=[[]*LM for _ in range(LM)]
                    Val=[[]*LM for _ in range(LM)]
                    
                    for n in protein.index:
                        for m in range(len(motifs1.Start)):
                            if int(motifs1.Start.iloc[m])+LM-1 >= n >= int(motifs1.Start.iloc[m]):
                                protein.Motif.loc[n]=1
                            else:
                                pass
                    
                    motif = protein.loc[p:p+LM]  
                                      
                    for i in variants1.index:
                        for n in motif.index:
                            if n == variants1.ProteinPosition.loc[i] and variants1.ProteinChange.loc[i][:3] == motif.Res.loc[n]:
                                if motif.Variants.loc[n] == 0:
                                    motif.Variants.loc[n]=[variants1.ProteinChange.loc[i]]
                                    if variants1.Clinsig.loc[i] == 'Pathogenic':
                                        motif.ClinSig.loc[n].append(1)
                                    elif variants1.Clinsig.loc[i] == 'Pathogenic/Likely pathogenic':
                                        motif.ClinSig.loc[n].append(0.75)
                                    elif variants1.Clinsig.loc[i] == 'Likely pathogenic' or variants1.Clinsig.loc[i] == 'Likely pathogenic, other':
                                        motif.ClinSig.loc[n].append(0.5)
                                    elif variants1.Clinsig.loc[i] == 'Benign':
                                        motif.ClinSig.loc[n].append(-1)
                                        #print(g, p)
                                        #print(variantes1.ProteinChange.loc[i])
                                    elif variants1.Clinsig.loc[i] == 'Benign/Likely bening':
                                        motif.ClinSig.loc[n].append(-0.75)
                                        #print(g, p)
                                        #print(variantes1.ProteinChange.loc[i])
                                    elif variants1.Clinsig.loc[i] == 'Likely benign' or variants1.Clinsig.loc[i] == 'Likely benign, other':
                                        motif.ClinSig.loc[n].append(-0.5)
                                        #print(g, p)
                                        #print(variantes1.ProteinChange.loc[i])
                                    else:
                                        pass
                                        #motivo.ClinSig.loc[n].append(0)
                                else:
                                    motif.Variants.loc[n].append(variants1.ProteinChange.loc[i])
                                    if variants1.Clinsig.loc[i] == 'Pathogenic':
                                        motif.ClinSig.loc[n].append(1)
                                    elif variants1.Clinsig.loc[i] == 'Pathogenic/Likely pathogenic':
                                        motif.ClinSig.loc[n].append(0.75)
                                    elif variants1.Clinsig.loc[i] == 'Likely pathogenic' or variants1.Clinsig.loc[i] == 'Likely pathogenic, other':
                                        motif.ClinSig.loc[n].append(0.5)
                                    elif variants1.Clinsig.loc[i] == 'Benign':
                                        motif.ClinSig.loc[n].append(-1)
                                        #print(g, p)
                                        #print(variantes1.ProteinChange.loc[i])
                                    elif variants1.Clinsig.loc[i] == 'Benign/Likely bening':
                                        motif.ClinSig.loc[n].append(-0.75)
                                        #print(g, p)
                                        #print(variantes1.ProteinChange.loc[i])
                                    elif variants1.Clinsig.loc[i] == 'Likely benign' or variants1.Clinsig.loc[i] == 'Likely benign, other':
                                        motif.ClinSig.loc[n].append(-0.5)
                                        #print(g, p)
                                        #print(variantes1.ProteinChange.loc[i])
                                    else:
                                        pass
                                        #motivo.ClinSig.loc[n].append(0)
                                motif.VarNum.loc[n]=len(motif.Variants.loc[n])
                            else:
                                pass
                    
                    if len(motif)==LM:
                        for n in range(LM):
                            if motif.VarNum.iloc[n]!=0:
                                cambios = motif.Variants.iloc[n]
                                CS = motif.ClinSig.iloc[n]
                                for i in range(len(CS)):                        
                                    if cambios[i][-3:]=='Gly':
                                        Gly[n].append(CS[i])
                                    elif cambios[i][-3:]=='Ala':
                                        Ala[n].append(CS[i])
                                    elif cambios[i][-3:]=='Arg':
                                        Arg[n].append(CS[i])
                                    elif cambios[i][-3:]=='Asn':
                                        Asn[n].append(CS[i])
                                    elif cambios[i][-3:]=='Asp':
                                        Asp[n].append(CS[i])                            
                                    elif cambios[i][-3:]=='Cys':
                                        Cys[n].append(CS[i])
                                    elif cambios[i][-3:]=='Glu':
                                        Glu[n].append(CS[i])                            
                                    elif cambios[i][-3:]=='Gln':
                                        Gln[n].append(CS[i])
                                    elif cambios[i][-3:]=='His':
                                        His[n].append(CS[i])
                                    elif cambios[i][-3:]=='Ile':
                                        Ile[n].append(CS[i])
                                    elif cambios[i][-3:]=='Leu':
                                        Leu[n].append(CS[i])
                                    elif cambios[i][-3:]=='Lys':
                                        Lys[n].append(CS[i])
                                    elif cambios[i][-3:]=='Met':
                                        Met[n].append(CS[i])
                                    elif cambios[i][-3:]=='Phe':
                                        Phe[n].append(CS[i])
                                    elif cambios[i][-3:]=='Pro':
                                        Pro[n].append(CS[i])
                                    elif cambios[i][-3:]=='Ser':
                                        Ser[n].append(CS[i])
                                    elif cambios[i][-3:]=='Thr':
                                        Thr[n].append(CS[i])
                                    elif cambios[i][-3:]=='Trp':
                                        Trp[n].append(CS[i])
                                    elif cambios[i][-3:]=='Tyr':
                                        Tyr[n].append(CS[i])
                                    elif cambios[i][-3:]=='Val':
                                        Val[n].append(CS[i])
                    
                    elif len(motif)>LM:
                        quantity=len(motif)/LM
                        for i in range(int(quantity)):
                            h=0
                            j=LM
                            motif1=motif.iloc[h:j]
                            for n in range(LM):
                                if motif1.VarNum.iloc[n]!=0:
                                    cambios = motif1.Variants.iloc[n]
                                    CS = motif1.ClinSig.iloc[n]
                                    for i in range(len(CS)):                        
                                        if cambios[i][-3:]=='Gly':
                                            Gly[n].append(CS[i])
                                        elif cambios[i][-3:]=='Ala':
                                            Ala[n].append(CS[i])
                                        elif cambios[i][-3:]=='Arg':
                                            Arg[n].append(CS[i])
                                        elif cambios[i][-3:]=='Asn':
                                            Asn[n].append(CS[i])
                                        elif cambios[i][-3:]=='Asp':
                                            Asp[n].append(CS[i])                            
                                        elif cambios[i][-3:]=='Cys':
                                            Cys[n].append(CS[i])
                                        elif cambios[i][-3:]=='Glu':
                                            Glu[n].append(CS[i])                            
                                        elif cambios[i][-3:]=='Gln':
                                            Gln[n].append(CS[i])
                                        elif cambios[i][-3:]=='His':
                                            His[n].append(CS[i])
                                        elif cambios[i][-3:]=='Ile':
                                            Ile[n].append(CS[i])
                                        elif cambios[i][-3:]=='Leu':
                                            Leu[n].append(CS[i])
                                        elif cambios[i][-3:]=='Lys':
                                            Lys[n].append(CS[i])
                                        elif cambios[i][-3:]=='Met':
                                            Met[n].append(CS[i])
                                        elif cambios[i][-3:]=='Phe':
                                            Phe[n].append(CS[i])
                                        elif cambios[i][-3:]=='Pro':
                                            Pro[n].append(CS[i])
                                        elif cambios[i][-3:]=='Ser':
                                            Ser[n].append(CS[i])
                                        elif cambios[i][-3:]=='Thr':
                                            Thr[n].append(CS[i])
                                        elif cambios[i][-3:]=='Trp':
                                            Trp[n].append(CS[i])
                                        elif cambios[i][-3:]=='Tyr':
                                            Tyr[n].append(CS[i])
                                        elif cambios[i][-3:]=='Val':
                                            Val[n].append(CS[i])
                            h+=LM
                            m+=LM
                                
                    for n in range(LM):
                        Gly_final[n]+=Gly[n]
                        Ala_final[n]+=Ala[n]
                        Arg_final[n]+=Arg[n]
                        Asn_final[n]+=Asn[n]
                        Asp_final[n]+=Asp[n]
                        Cys_final[n]+=Cys[n]
                        Glu_final[n]+=Glu[n]
                        Gln_final[n]+=Gln[n]
                        His_final[n]+=His[n]
                        Ile_final[n]+=Ile[n]
                        Leu_final[n]+=Leu[n]
                        Lys_final[n]+=Lys[n]
                        Met_final[n]+=Met[n]
                        Phe_final[n]+=Phe[n]
                        Pro_final[n]+=Pro[n]
                        Ser_final[n]+=Ser[n]
                        Thr_final[n]+=Thr[n]
                        Trp_final[n]+=Trp[n]
                        Tyr_final[n]+=Tyr[n]
                        Val_final[n]+=Val[n]
    
    return Gly_final, Ala_final, Arg_final, Asn_final, Asp_final, Cys_final, Glu_final, Gln_final, His_final, Ile_final, Leu_final, Lys_final, Met_final, Phe_final, Pro_final, Ser_final, Thr_final, Trp_final, Tyr_final, Val_final

def ClinSigMatrix(df, mot, positive_ctrl, motif_name):
    ''' Builds the ClinSig Substitution Matrix. 
        Both the CSV and the Figure.'''
    mot = mot.split('.')
    LM = len(mot)
    aa = MatSubClinVar(df, mot, positive_ctrl)
    Gly_final, Ala_final, Arg_final, Asn_final, Asp_final, Cys_final, Glu_final, Gln_final, His_final, Ile_final, Leu_final, Lys_final, Met_final, Phe_final, Pro_final, Ser_final, Thr_final, Trp_final, Tyr_final, Val_final = aa
    aminoacids=['G', 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    df_=pd.DataFrame(list(zip(Gly_final, Ala_final, Arg_final, 
                              Asn_final, Asp_final, Cys_final, 
                              Glu_final, Gln_final, His_final, 
                              Ile_final, Leu_final, Lys_final, 
                              Met_final, Phe_final, Pro_final, 
                              Ser_final, Thr_final, Trp_final, 
                              Tyr_final, Val_final)), columns=aminoacids)
    df_['Motif']=mot
    df_.set_index('Motif', inplace=True)
    df_.columns=aminoacids
    df_.to_csv('../'+motif_name+'_motif/ClinSigMatrix_'+motif_name+'_positives.csv', sep='\t')

    opposed = 0
    for aminoacid, a in enumerate(aa):
        for i, e in enumerate(a):
            b = set([x>0 for x in e])
            if True in b and False in b:
                opposed += 1
                print('Change of '+mot[i]+'('+str(i+1)+')'+' by '+aminoacids[aminoacid]+' have opposed ClinSig')
            #else:
            #    print('cambiar por '+str(aminoacid)+' en posición '+str(i)+' NO tiene opuestos')
    print('There are '+str(opposed)+' amino acid changes that have opposed ClinSig')

    for n in range(LM):
        Gly_final[n]=np.nanmean(Gly_final[n])
        Ala_final[n]=np.nanmean(Ala_final[n])
        Arg_final[n]=np.nanmean(Arg_final[n])
        Asn_final[n]=np.nanmean(Asn_final[n])
        Asp_final[n]=np.nanmean(Asp_final[n])                           
        Cys_final[n]=np.nanmean(Cys_final[n])
        Glu_final[n]=np.nanmean(Glu_final[n])                           
        Gln_final[n]=np.nanmean(Gln_final[n])
        His_final[n]=np.nanmean(His_final[n])
        Ile_final[n]=np.nanmean(Ile_final[n])
        Leu_final[n]=np.nanmean(Leu_final[n])
        Lys_final[n]=np.nanmean(Lys_final[n])
        Met_final[n]=np.nanmean(Met_final[n])
        Phe_final[n]=np.nanmean(Phe_final[n])
        Pro_final[n]=np.nanmean(Pro_final[n])
        Ser_final[n]=np.nanmean(Ser_final[n])
        Thr_final[n]=np.nanmean(Thr_final[n])
        Trp_final[n]=np.nanmean(Trp_final[n])
        Tyr_final[n]=np.nanmean(Tyr_final[n])
        Val_final[n]=np.nanmean(Val_final[n])
    
    aminoacids=['G', 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    df_=pd.DataFrame(list(zip(Gly_final, Ala_final, Arg_final, 
                              Asn_final, Asp_final, Cys_final, 
                              Glu_final, Gln_final, His_final,
                              Ile_final, Leu_final, Lys_final, 
                              Met_final, Phe_final, Pro_final,
                              Ser_final, Thr_final, Trp_final, 
                              Tyr_final, Val_final)), columns=aminoacids)
    df_['Motif']=mot 
    df_.set_index('Motif', inplace=True)
    df_.columns=aminoacids
    df_.to_csv('../'+motif_name+'_motif/MeanClinSigMatrix_'+motif_name+'_positives.csv', sep='\t')    
    
    plt.figure(figsize=(15, 5))
    cmap = sns.diverging_palette(133, 10, as_cmap=True)
    ax = sns.heatmap(df_, vmin=-1, vmax=1, cmap=cmap, linewidths=.5, cbar_kws={'label':'Benign/Pathogenic'})
    sns.set(font_scale=1.2)
    plt.yticks(rotation=0)
    plt.ylabel('')
    plt.savefig('../'+motif_name+'_motif/MeanClinSigMatrix_HeatMap_'+motif_name+'_positives.png')
    
