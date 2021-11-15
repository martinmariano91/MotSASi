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


def positive_ctrl(motif_name):
    '''Look for the positive control group of the Motif.
    These are the motifs that have been experimentally tested.
    In: motif_name (ELMIdentifier)
    Out: list with the UniProtIDs of the positive control group
        small dataframe with the proteins that have a structure.'''
    elm_instances = pd.read_csv('ELM_instances.csv', index_col=0)
    elm_instances = elm_instances[(elm_instances['ELMIdentifier'] == motif_name) &
                                  (elm_instances['Organism'] == 'Homo sapiens')]
    positive_ctrl = tuple(zip(elm_instances['Primary_Acc'].tolist(), elm_instances['Start'].tolist()))
    protein_pdb = elm_instances.loc[elm_instances.PDB.notnull(), ['Primary_Acc', 'PDB']]
    protein_pdb = prot_pdb(protein_pdb)
    return positive_ctrl, protein_pdb
    
def pdbs_download(motif_name, protein_pdb):
    '''Download the pdbs reported for this motif. 
    Create individual folder for each of the structures.
    In: motif_name (ELMIdentifier)
    Out: folders with the queried structures'''
    pdbs = protein_pdb.PDB.tolist()
    pdbs = ' '.join(pdbs).split(' ')
    if len(pdbs) != 0:
        os.mkdir('../'+motif_name+'_motif/pdbs')
        pdbl = PDBList()
        print('Downloading PDBs files with this motif...')
        for pdb in pdbs:
            os.mkdir('../'+motif_name+'_motif/pdbs/'+pdb.upper())
            pdb_url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/'+pdb
            r = requests.get(pdb_url)
            data_pdb = r.text
            if '"experimental_method":["X-ray diffraction"]' in data_pdb:
                pdbl.retrieve_pdb_file(pdb, file_format='pdb', pdir='../'+motif_name+'_motif/pdbs/'+pdb.upper(), overwrite=True)
                os.rename('../'+motif_name+'_motif/pdbs/'+pdb+'/pdb'+pdb.lower()+'.ent',
                          '../'+motif_name+'_motif/pdbs/'+pdb+'/'+pdb.upper()+'.pdb')
        return 'YES'
    else:
        return 'NO'




    