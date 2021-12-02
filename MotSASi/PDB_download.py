#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 11:12:38 2021

@author: mariano
"""

from Bio.PDB import PDBList, PDBParser, PPBuilder
import pandas as pd
import os
import requests

def prot_pdb(protein_pdb):
    '''Re-shape the protein_pdb for FoldX calculations
    In: protein_pdb dataframe as extracted from ELM
    Out:protein_pdb dataframe with one PDB id per row.'''
    uni_ids = []
    pdb_ids = []
    for i in protein_pdb.index:
        if len(protein_pdb.loc[i, 'PDB'].split()) > 1:
            uni_ids += [protein_pdb.loc[i, 'Primary_Acc']]*len(protein_pdb.loc[i, 'PDB'].split())
            pdb_ids += protein_pdb.loc[i, 'PDB'].split()
        else:
            uni_ids.append(protein_pdb.loc[i, 'Primary_Acc'])
            pdb_ids.append(protein_pdb.loc[i, 'PDB'])
    protein_pdb = pd.DataFrame()
    protein_pdb['Primary_Acc'] =  uni_ids
    protein_pdb['PDB'] = pdb_ids
    return protein_pdb


def positive_ctrl(df, motif_name):
    '''Look for the positive control group of the Motif.
    These are the motifs that have been experimentally tested.
    In: motif_name (ELMIdentifier)
    Out: list with the UniProtIDs of the positive control group
        small dataframe with the proteins that have a structure.'''
    elm_instances = pd.read_csv('ELM/ELM_instances.csv', index_col=0)
    elm_instances = elm_instances[(elm_instances['ELMIdentifier'] == motif_name) &
                                  (elm_instances['Organism'] == 'Homo sapiens')]
    instances = [x for x in elm_instances['Primary_Acc'].tolist() if '-' in x and x[-1] == 1]
    instances += [x for x in elm_instances['Primary_Acc'].tolist() if not '-' in x]
    instances = [x for x in positives if x in df.UniProtID.tolist()]
    elm_instances = elm_instances[elm_instances['Primary_Acc'].isin(instances)]
    uniprots = []
    for i in elm_instances.index:
        if '-' in elm_instances.loc[i, 'Primary_Acc']:
            uniprots.append(elm_instances.loc[i, 'Primary_Acc'].split('-')[0])
        else:
            uniprots.append(elm_instances.loc[i, 'Primary_Acc'])
    elm_instances['Primary_Acc'] = uniprots
    elm_instances['ID_motif'] = elm_instances['Primary_Acc']+'_'+elm_instances['Start'].astype(str)   
    df['ID_motif'] = df.UniProtID+'_'+df.Start.astype(str)
    df['ELM_Positives'] = 0
    df['ELM_Negatives'] = 0
    df.loc[df.ID_motif.isin(elm_instances.loc[elm_instances['InstanceLogic'] == 'true positive', 'ID_motif'].tolist()), 'ELM_Positives'] = 1
    df.loc[df.ID_motif.isin(elm_instances.loc[elm_instances['InstanceLogic'] == 'false positive', 'ID_motif'].tolist()), 'ELM_Negatives'] = 1
    #df.loc[df.UniProtID.isin(uniprots), 'ELM_Positives'] = 1
    df.to_csv('../'+motif_name+'_motif/'+motif_name+'_motif.csv', sep='\t') 
    positive_ctrl = tuple(zip(df.loc[df.ELM_Positives == 1, 'UniProtID'].tolist(), df.loc[df.ELM_Positives == 1, 'Start'].tolist()))
    protein_pdb = elm_instances.loc[elm_instances.PDB.notnull(), ['Primary_Acc', 'PDB']]
    protein_pdb = prot_pdb(protein_pdb)
    pdbs = protein_pdb.PDB.tolist()
    pdbs = ' '.join(pdbs).split(' ')
    if len(pdbs) != 0:
        pdbl = PDBList()
        for pdb in pdbs:
            pdb_url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/'+pdb
            r = requests.get(pdb_url)
            data_pdb = r.text
            if not '"experimental_method":["X-ray diffraction"]' in data_pdb:
                protein_pdb.drop(protein_pdb[protein_pdb.PDB == pdb].index, inplace=True)
    return positive_ctrl, protein_pdb, df
    
def pdbs_download(motif_name, protein_pdb):
    '''Download the pdbs reported for this motif. 
    Create individual folder for each of the structures.
    In: motif_name (ELMIdentifier)
    Out: folders with the queried structures'''
    pdbs = protein_pdb.PDB.tolist()
    pdbs = ' '.join(pdbs).split(' ')
    os.mkdir('../'+motif_name+'_motif/pdbs')
    if len(pdbs) != 0:
        pdbl = PDBList()
        print('Downloading PDBs files with this motif...')
        for pdb in pdbs:
            pdb_url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/'+pdb
            r = requests.get(pdb_url)
            data_pdb = r.text
            if '"experimental_method":["X-ray diffraction"]' in data_pdb:
                os.mkdir('../'+motif_name+'_motif/pdbs/'+pdb.upper())
                pdbl.retrieve_pdb_file(pdb, file_format='pdb', pdir='../'+motif_name+'_motif/pdbs/'+pdb.upper(), overwrite=True)
                os.rename('../'+motif_name+'_motif/pdbs/'+pdb+'/pdb'+pdb.lower()+'.ent',
                          '../'+motif_name+'_motif/pdbs/'+pdb+'/'+pdb.upper()+'.pdb')
        return 'YES'
    else:
        return 'NO'




    