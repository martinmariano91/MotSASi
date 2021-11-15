#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 16:20:51 2020

@author: mariano
"""

import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess

#%%
def motifs_fastas(ID, p, mot):
    ''' Gets the fasta sequence of the residues sourrounding the motif
    (100 residues before and 100 residues after the starting position).
    In: UniProtID, motif start position, dot separated motif.
    Out: Fasta sequence.'''
    mot = mot.split('.')
    LM = len(mot)
    lectura_fasta = SeqIO.parse('../UniProt-Proteome.fasta','fasta')
    for seq_record in lectura_fasta:
        if seq_record.id.split('|')[1] == ID:
            fastas = []
            if p > 100:
                seq = str(seq_record.seq)[p-100:p+100+LM]
                fastas += [SeqRecord(Seq(seq), id=ID+'_'+str(p), description='')]
            elif p <= 100:
                seq = str(seq_record.seq)[0:p+100+LM]
                fastas += [SeqRecord(Seq(seq), id=ID+'_'+str(p), description='')]
            return fastas

#%%

def jpred_run(positive_ctrl, mot, motif_name):
    '''Executes the batch submition for calculating secondary structure
    prediction of the sequences that contains the motifs.
    In: list of proteins of the positive control, dot separated motif,
    and motif name.
    Out: Excution of Jpred commands.'''
    os.chdir('jpred')
    UniMotifs = positive_ctrl
    secuencias = []
    for ID, p in UniMotifs:
        secuencias += motifs_fastas(ID, p, mot)
    SeqIO.write(secuencias, motif_name+'.fasta','fasta')
    command_line1 = ['./prepareInputs.csh', motif_name+'.fasta']
    subprocess.run(command_line1, stdout=subprocess.DEVNULL)    
    print('Running Jpred on batch form...')
    print('This may take several minutes...')
    command_line2 = ['./massSubmitScheduler.csh', motif_name+'.fasta_dir/']
    subprocess.run(command_line2, stdout=subprocess.DEVNULL)
    os.chdir('../')

#%%
    
def jpred_batch_processing(mot, motif_name):
    '''Process the output of Jpred batch submission and outputs a 
    dataframe with the info of secondary structure proabbility of
    protein's motif and a boxplot showing these information.
    In: dot separated motif and motif name.
    Out: dataframe of motifs secondary structures probabilities and
    boxplot of motif probabilities.'''
    mot = mot.split('.')
    LM = len(mot)
    jpred = [e for e in os.listdir('jpred/'+motif_name+'.fasta_dir/_output') if '.jnet' in e]
    
    UniProtID = []
    jnetprope = [] 
    jnetproph = []
    jnetpropc = []
    
    for e in jpred:
        df = pd.read_csv('jpred/'+motif_name+'.fasta_dir/_output/'+e, header=None)
        ID = e.split('_')[0]
        p = int(e.split('_')[1])
        column_1 = []
        index = []
        for v in df[0]:
            column_1.append(v.split(':')[1])
            index.append(v.split(':')[0])
        df.drop(df.columns[-1], axis=1, inplace=True)
        df['Index'] = index
        df[0] = column_1
        df.set_index('Index', drop=True, inplace=True)
        df = df.T
        if p > 100:
            seq_motif = df.iloc[99:99+LM]
        elif p <= 100:
            seq_motif = df.iloc[p:p+LM-1]
        UniProtID.append(ID+'_'+str(p))
        jnetprope.append(seq_motif.JNETPROPE.astype(float).mean())
        jnetproph.append(seq_motif.JNETPROPH.astype(float).mean())
        jnetpropc.append(seq_motif.JNETPROPC.astype(float).mean())
        
    df_ = pd.DataFrame()
    df_['UniProtID'] = UniProtID
    df_['jnetprope'] = jnetprope
    df_['jnetproph'] = jnetproph
    df_['jnetpropc'] = jnetpropc
    df_.to_csv('../'+motif_name+'_motif/jpred_'+motif_name+'_positives.csv') 
    
    plt.figure(figsize=(5, 5))
    sns.boxplot(data=df_[['jnetprope', 'jnetproph', 'jnetpropc']], color='white', width=0.2)
    plt.ylabel('Probability')
    plt.xticks(range(3), ['BetaSheet', 'AlphaHelix', 'Coil'])
    plt.savefig('../'+motif_name+'_motif/jpred_'+motif_name+'_positives.png')


    



