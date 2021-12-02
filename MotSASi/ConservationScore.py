#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:47:24 2020

Busca los clusters UniRef90 de cada proteína. Descarta los de UniParc y quedan solo los de UniProt (swiss+trembl)
Si hay más de dos secuencias en ese cluster, las alinea con clustalo en tres iteraciones.
Luego calcula el score de conservación usando Jensen-Shannon divergence.
Los controles positivos (proteínas que tienen el motivo corroborado experimentalmente),
son las que permiten setear el cut-off del score (en el caso de triptofano-ácido el cut-off es 0.84)

@author: mariano
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import urllib
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import os
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)


def ConservationScore(df, mot, motif_name, positive_ctrl):
    '''Calculates the Jennsen-Shannon divergence of each of the
        residues that constitutes the motif of each of the protein
        of the positive control.
    In: dataset that have the info of the motifs detected in the
        human proteome, the dot separated motif, the motif name,
        and the dataset that have the info of the positive ctrl.
    Out: Dataframe with the conservation score of each of the 
        residues that constitutes the motif. Plot the boxplot 
        of the conservation scores.'''
    motifs = df
    mot = mot.split('.')
    LM = len(mot)
    UniMotifs = positive_ctrl
    df_uniref=pd.read_csv('UniRef-HomoSapiens-90.tab', sep='\t')
    df_uniref['Cluster ID'].replace('UniRef90_','', regex=True, inplace=True)
    
    IDs = []
    scores = {}
    
    for g, p in UniMotifs:
        a=df_uniref[df_uniref['Cluster ID']==g]
        ortologos_ids=[]
        ortologos_seq=[]
        if not a.empty:
            ortologos_ids=a['Cluster members'].iloc[0].split('; ')
            ortologos_ids=[x for x in ortologos_ids if not x.startswith('UPI')] 
            ortologos_organisms=[]
            ortologs_seq=[]
            for i in ortologos_ids:
                try:
                    handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+i+".xml")
                    record = SeqIO.read(handle, "uniprot-xml")
                    if record.annotations['organism'] not in ortologos_organisms:
                        ortologos_organisms.append(record.annotations['organism'])
                        ortologos_seq.append(SeqRecord(Seq(str(record.seq)), id=i, description=''))
                except:
                    ortologos_ids.remove(i)
           
            if len(ortologos_seq)>2:
                SeqIO.write(ortologos_seq, "multiple_sequence.fasta", "fasta") 
                
                command_line1 = ['clustalo', '-i', 'multiple_sequence.fasta', '-o', 'multiple_alignment.fasta', '--outfmt=fasta', '--iter=3', '--force']
                process_1 = subprocess.run(command_line1, stdout=subprocess.DEVNULL)
                command_line2 = ['python2', 'score_conservation.py', '-o', 'scores.txt', '-a'+g, 'multiple_alignment.fasta']
                process_2 = subprocess.run(command_line2, stdout=subprocess.DEVNULL)
                position=motifs[motifs.UniProtID==g].Start
                df_=pd.read_csv('scores.txt', header=2, sep='\t', index_col=0)
                newindex=np.arange(len(df_))+1
                df_.set_index(newindex, inplace=True)
                for x in position:
                    if x == p:
                        IDs.append(g+'_'+str(x))
                        scores[g+'_'+str(x)] = list(df_.score.loc[x:x+LM-1])
    
    df_score = pd.DataFrame()
    df_score['Motif']= [x+'_'+str(i+1) for i, x in enumerate(mot)]
    for i in IDs:
        df_score[i] = scores[i]
    df_score.set_index('Motif', inplace=True, drop=True)
    
    df_score.to_csv('../'+motif_name+'_motif/ConservationScore'+motif_name+'_positives.csv', sep='\t')
    
    plt.figure(figsize=(8, 4), facecolor="w")
    sns.boxplot(data=df_score.T, color='white', width=0.2, orient='h')
    plt.xlabel('Conservation Score')
    plt.xlim(0.65, 0.95)
    plt.savefig('../'+motif_name+'_motif/ConservationScore'+motif_name+'_positives.png', bbox_inches = "tight")


        