#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 14:34:20 2020

Search in the human proteome for the motif specified with a regular expression.
Outputs a 'csv' file with all the motifs found with the information about
the protein that contains it and the position where it is located.

NOTE: the folder were the script is run should also contain the file "uniprot_sprot.xml.gz"

@author: Mariano Martín
"""
import pandas as pd
from Bio import SeqIO
import gzip
import re

def m_append(m, record):
    ''' Save the info regarding the motif in study.
        In: motif (m) regular expression and hit record.
        Out: tuple of protein UniprotID, motif sequence, 
        motif start position, motif end position, 
        protein length, gene name.'''
    for x in re.finditer(m, str(record.seq)):
        ID = record.id
        M = x.group()
        PI = x.start()+1
        PF = x.end()+1
        length = len(record.seq)
        if 'gene_name_primary' in record.annotations.keys():
            gene = record.annotations['gene_name_primary']
        else:
            gene = 'none'
    
    return ID, M, PI, PF, length, gene

def HumanProteomeSearch(m, motif_name):
    ''' Search for motif (m) in the human proteome. Outputs the hits in
        a dataframe with information regarding the protein that contains
        the motif, its location, etc. Important information is printed.
        In: motif regular expression and motif name.
        Out: DataFrame with motif hits information. The results are saved
        in a .csv file with the name of the motif. '''
    IDtodos = []
    IDs = []
    PIs = []
    PFs = []
    Ms = []
    genes = []
    lengths = []    
                
    handle = gzip.open("uniprot_sprot.xml.gz", "rt")        
    for record in SeqIO.parse(handle, "uniprot-xml"):
        if record.annotations['organism'] == 'Homo sapiens (Human)':
            mo = re.search(m, str(record.seq))
            if mo:
                IDtodos.append(record.id)         
                ID, M, PI, PF, length, gene = m_append(m, record)         
                IDs.append(ID)
                Ms.append(M)
                PIs.append(PI)
                PFs.append(PF)
                lengths.append(length)
                genes.append(gene)
                
    listazip=list(zip(IDs,genes,lengths, Ms, PIs, PFs))
    df=pd.DataFrame(listazip, columns=['UniProtID','Gene', 'Protein_length', 'Motif', 'Start', 'End'])
    df.to_csv('../'+motif_name+'_motif/'+motif_name+'_motif.csv', sep='\t') 
     
    print(motif_name+' motif found in '+str(len(IDtodos))+' human proteins')
    print('Total hits in the human proteome: '+str(len(Ms)))          
    
    return df







