#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 14:34:20 2020

Searches in the human proteome for the motif specified with a regular expression.
Outputs a 'csv' file with all the motifs found with the information about
the protein that contains it and the position where it is located.

NOTE: the folder were the script is run should also contain the file "uniprot_sprot.xml.gz"

@author: Mariano Martín
"""
import pandas as pd
from Bio import SeqIO
import re
import urllib

def m_append(m, record):
    ''' Save the info regarding the motif in study.
        In: motif (m) regular expression and hit record.
        Out: tuple of protein UniprotID, motif sequence, 
        motif start position, motif end position, 
        protein length, gene name.'''
    ID = []
    M = []
    PI = []
    PF = []
    length = []
    gene = []
    for x in re.finditer(m, str(record.seq)):
        try:
            uni_id = record.id.split('|')[1]
            ID.append(uni_id)
            M.append(x.group())
            PI.append(x.start()+1)
            PF.append(x.end()+1)
            length.append(len(record.seq))
            if 'GN=' in record.description:
                gene.append(record.description.split('GN=')[1].split()[0])
            else:
                handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+uni_id+".xml")
                record = SeqIO.read(handle, "uniprot-xml")
                if 'gene_name_primary' in record.annotations.keys():
                    gene.append(record.annotations['gene_name_primary'])
                else:
                    try:
                        gene.append(record.annotations['gene_name_synonym'])
                    except:
                        gene.append(None)
        except:
            pass                
    return ID, M, PI, PF, length, gene

def HumanProteomeSearch(m, motif_name):
    ''' Search for motif (m) in the human proteome. Outputs the hits in
        a dataframe with information regarding the protein that contains
        the motif, its location, etc. Important information is printed.
        In: motif regular expression and motif name.
        Out: DataFrame with motif hits information. The results are saved
        in a .csv file with the name of the motif. '''
    IDs = []
    PIs = []
    PFs = []
    Ms = []
    genes = []
    lengths = []    
    
    fastas = SeqIO.parse('uniprot_sprot_h.fasta','fasta')
    #handle = gzip.open("uniprot_sprot_h.xml.gz", "rt")        
    for record in fastas:
        mo = re.search(m, str(record.seq))
        if mo: 
            ID, M, PI, PF, length, gene = m_append(m, record)         
            IDs += ID
            Ms += M
            PIs += PI
            PFs += PF
            lengths += length
            genes += gene
                
    listazip=list(zip(IDs,genes,lengths, Ms, PIs, PFs))
    df=pd.DataFrame(listazip, columns=['UniProtID','Gene', 'Protein_length', 'Motif', 'Start', 'End'])
    df.to_csv('../'+motif_name+'_motif/'+motif_name+'_motif.csv', sep='\t') 
     
    print(motif_name+' motif found in '+str(len(df.UniProtID.unique()))+' human proteins')
    print('Total hits in the human proteome: '+str(len(Ms)))          
    
    return df







