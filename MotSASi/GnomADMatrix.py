#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:47:34 2020

Crea Matriz de sustitución con el valor o de la frecuencia alélica (AFs) de cada variante.

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
import math
import warnings
warnings.filterwarnings("ignore")
#%%

def iterador(UniMotifs, variants, founds, motifs, LM, 
             Gly_final, Ala_final, Arg_final, Asn_final,
             Asp_final, Cys_final, Glu_final, Gln_final,
             His_final, Ile_final, Leu_final, Lys_final,
             Met_final, Phe_final, Pro_final, Ser_final,
             Thr_final, Trp_final, Tyr_final, Val_final):
    '''Compiles the change of each amino acid in the motif '''
    lectura_fasta = SeqIO.parse('UniProt-Proteome.fasta','fasta')
    for seq_record in lectura_fasta:
        for g, p in UniMotifs:
            if seq_record.id.split('|')[1] == g:   
                variants1 = variants[variants['UniProtID'] == g]    
                if not variants1.empty:
                    founds.append(g)
                    ix = motifs.UniProtID == g
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
                    protein['AFvariants']=[[]*len(index) for _ in range(len(index))]
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
                                    if variants1.AF.loc[i].is_integer:
                                        motif.AFvariants.loc[n]=[variants1.AF.loc[i]]
                                    else:
                                        motif.AFvariants.loc[n]=[0]
                                else:
                                    motif.Variants.loc[n].append(variants1.ProteinChange.loc[i])
                                    if variants1.AF.loc[i].is_integer:
                                        motif.AFvariants.loc[n].append(variants1.AF.loc[i])
                                    else:
                                        motif.AFvariants.loc[n].append(0)
                                motif.VarNum.loc[n]=len(motif.Variants.loc[n])
                            else:
                                pass
                    
                    #aquellas variantes que no tienen AFs les puse 0 en el paso anterior. 
                    #Ahora lo reemplazo por el minimo valor encontrado de entre todas las variantes.
                    minim_AFs=variants1.AF[variants1.AF>0].min()
                    for i in motif.AFvariants:
                        while 0 in i:
                            i[i.index(0)]=minim_AFs
                   
                    if len(motif)==LM:
                        for n in range(LM):
                            if motif.VarNum.iloc[n]!=0:
                                cambios = motif.Variants.iloc[n]
                                AFs = motif.AFvariants.iloc[n]
                                for i in range(len(AFs)):
                                    camb_ori=cambios[i][:3]
                                    camb_mut=cambios[i][-3:]
                                    #acá evalúo cada variante y luego multiplico la AFs por la probabilidad dada por la matriz de codones.
                                    if camb_mut =='Gly':
                                        Gly[n].append(AFs[i])
                                    elif camb_mut =='Ala':
                                        Ala[n].append(AFs[i])
                                    elif camb_mut =='Arg':
                                        Arg[n].append(AFs[i])
                                    elif camb_mut =='Asn':
                                        Asn[n].append(AFs[i])
                                    elif camb_mut =='Asp':
                                        Asp[n].append(AFs[i])
                                    elif camb_mut =='Cys':
                                        Cys[n].append(AFs[i])
                                    elif camb_mut =='Glu':
                                        Glu[n].append(AFs[i])                     
                                    elif camb_mut =='Gln':
                                        Gln[n].append(AFs[i])
                                    elif camb_mut =='His':
                                        His[n].append(AFs[i])
                                    elif camb_mut =='Ile':
                                        Ile[n].append(AFs[i])
                                    elif camb_mut =='Leu':
                                        Leu[n].append(AFs[i])
                                    elif camb_mut =='Lys':
                                        Lys[n].append(AFs[i])
                                    elif camb_mut =='Met':
                                        Met[n].append(AFs[i])
                                    elif camb_mut =='Phe':
                                        Phe[n].append(AFs[i])
                                    elif camb_mut =='Pro':
                                        Pro[n].append(AFs[i])
                                        if n == 0 and camb_ori == 'Ser':
                                            print('*'*20)
                                            print(g)
                                    elif camb_mut =='Ser':
                                        Ser[n].append(AFs[i])
                                    elif camb_mut =='Thr':
                                        Thr[n].append(AFs[i])
                                    elif camb_mut =='Trp':
                                        Trp[n].append(AFs[i])
                                    elif camb_mut =='Tyr':
                                        Tyr[n].append(AFs[i])
                                    elif camb_mut =='Val':
                                        Val[n].append(AFs[i])
                    
                    elif len(motif)>LM:
                        quantity=len(motif)/LM
                        for i in range(int(quantity)):
                            h=0
                            j=LM
                            motif1=motif.iloc[h:j]
                            for n in range(LM):
                                if motif1.VarNum.iloc[n]!=0:
                                    cambios = motif1.Variants.iloc[n]
                                    AFs = motif1.AFvariants.iloc[n]
                                    for i in range(len(AFs)):                        
                                        camb_ori=cambios[i][:3]
                                        camb_mut=cambios[i][-3:]
                                        if camb_mut =='Gly':
                                            Gly[n].append(AFs[i])
                                        elif camb_mut =='Ala':
                                            Ala[n].append(AFs[i])
                                        elif camb_mut =='Arg':
                                            Arg[n].append(AFs[i])
                                        elif camb_mut =='Asn':
                                            Asn[n].append(AFs[i])
                                        elif camb_mut =='Asp':
                                            Asp[n].append(AFs[i])
                                        elif camb_mut =='Cys':
                                            Cys[n].append(AFs[i])
                                        elif camb_mut =='Glu':
                                            Glu[n].append(AFs[i])
                                        elif camb_mut =='Gln':
                                            Gln[n].append(AFs[i])
                                        elif camb_mut =='His':
                                            His[n].append(AFs[i])
                                        elif camb_mut =='Ile':
                                            Ile[n].append(AFs[i])
                                        elif camb_mut =='Leu':
                                            Leu[n].append(AFs[i])
                                        elif camb_mut =='Lys':
                                            Lys[n].append(AFs[i])
                                        elif camb_mut =='Met':
                                            Met[n].append(AFs[i])
                                        elif camb_mut =='Phe':
                                            Phe[n].append(AFs[i])
                                        elif camb_mut =='Pro':
                                            Pro[n].append(AFs[i])
                                            if n == 0 and camb_ori == 'Ser':
                                                print('*'*20)
                                                print(g)
                                        elif camb_mut =='Ser':
                                            Ser[n].append(AFs[i])
                                        elif camb_mut =='Thr':
                                            Thr[n].append(AFs[i])
                                        elif camb_mut =='Trp':
                                            Trp[n].append(AFs[i])
                                        elif camb_mut =='Tyr':
                                            Tyr[n].append(AFs[i])
                                        elif camb_mut =='Val':
                                            Val[n].append(AFs[i])
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

def MatSubGnomAD(df, mot, positive_ctrl):
    '''Compiles the change of each amino acid in the motif '''
    motifs = df
    UniMotifs = positive_ctrl
    LM = len(mot) #largo del motivo
    
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
    
    UniMotifs = positive_ctrl
    cromosomas = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    founds = []
    for x in cromosomas:    
        variants = pd.read_csv('GnomAD/'+str(x)+'.csv', sep='\t')
        variants.drop(columns='Unnamed: 0', inplace=True)
        variants = variants[variants['Effect']=='missense_variant']
        variants['ProteinPosition'] = variants['ProteinPosition'].astype(int)
        UniMotifs=[(u, p) for u, p in UniMotifs if u not in founds]    
        
        iterador(UniMotifs, variants, founds, motifs, LM, 
                 Gly_final, Ala_final, Arg_final, Asn_final,
                 Asp_final, Cys_final, Glu_final, Gln_final,
                 His_final, Ile_final, Leu_final, Lys_final,
                 Met_final, Phe_final, Pro_final, Ser_final,
                 Thr_final, Trp_final, Tyr_final, Val_final)

    return Gly_final, Ala_final, Arg_final, Asn_final, Asp_final, Cys_final, Glu_final, Gln_final, His_final, Ile_final, Leu_final, Lys_final, Met_final, Phe_final, Pro_final, Ser_final, Thr_final, Trp_final, Tyr_final, Val_final

def GnomADMatrix(df, mot, positive_ctrl, motif_name):
    ''' Builds the GnomAD Substitution Matrix. 
        Both the CSV and the Figure.'''
    mot = mot.split('.')
    LM = len(mot)
    aa = MatSubGnomAD(df, mot, positive_ctrl)
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
    df_.replace(0, np.nan, inplace=True)
    df_.columns=aminoacids
    df_.to_csv('../'+motif_name+'_motif/AFMatrix_'+motif_name+'_positives.csv', sep='\t')

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


    df_=pd.DataFrame(list(zip(Gly_final, Ala_final, Arg_final, Asn_final, Asp_final, 
                             Cys_final, Glu_final, Gln_final, His_final, Ile_final,
                             Leu_final, Lys_final, Met_final, Phe_final, Pro_final,
                             Ser_final, Thr_final, Trp_final, Tyr_final, Val_final)), columns=aminoacids)
    df_['Motif']=mot
    df_.set_index('Motif', inplace=True)
    df_.replace(0, np.nan, inplace=True)
    df_=df_.transform([math.log10])
    df_=df_*-1
    df_.columns=aminoacids
    df_.to_csv('../'+motif_name+'_motif/MeanAFMatrix_'+motif_name+'_positives.csv', sep='\t')    
    
    plt.figure(figsize=(15, 5))
    ax = sns.heatmap(df_, vmin=2, vmax=6, cmap="YlOrRd", annot=True, linewidths=.5, cbar_kws={'label':'-log10(Fo)'})
    sns.set(font_scale=1.2)
    plt.yticks(rotation=0)
    plt.ylabel('')
    plt.savefig('../'+motif_name+'_motif/MeanAFMatrix_HeatMap_'+motif_name+'_positives.png')




