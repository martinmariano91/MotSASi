#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:25:54 2020

Run mutational scan of each motif residue with foldx for motifs containing crystals. 
Process the "scanning_output" for each position, concatenates them in a
single dataframe called as Motif's FoldX Matrix and plot the heatmap of this
substitution matrix. 
The final function "mean_foldx_matrix" calculates the mean ddG for 
each residue substitution and outputs this dataframe and a heatmap.

NOTE: the "foldx_run" function needs precise information that should be manually get.
      the crsytals has to be placed in a folder with the name of the crystal. 
      the script has to be run in a main folder that contains the crystal folders.

@author: Mariano Martín
"""

import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from Bio.PDB import PDBParser, PPBuilder
import re

def dataframe_pdb(pdb, uniprotid):
    ''' Builds a DataFrame from the info in the pdb file.
        This DF have the UniProtID that each chain contains.'''
    df_pdb = pd.read_csv(pdb+'/'+pdb+'.pdb', sep='\t')
    df_pdb = df_pdb[df_pdb.iloc[:, 0].str.startswith('DBREF')]
    df_pdb.columns = ['Column']
    chain_uniprotid = pd.DataFrame()
    chain_uniprotid['Chain'] = df_pdb.Column.apply(lambda x: x.split()[2])
    chain_uniprotid['UniProtID'] = df_pdb.Column.apply(lambda x: x.split()[6])
    chain_uniprotid = chain_uniprotid[chain_uniprotid['UniProtID']==uniprotid]
    chains = chain_uniprotid.Chain.tolist()
    return chains

def motif_pdb(chain, motif):
    ''' Searches for the motif in the peptides chains of the PDB.
        If the motif is found gets the chain and the position
        where the motif starts.'''
    ppb = PPBuilder()
    for pp in ppb.build_peptides(chain):
        seq = pp.get_sequence()      
        motif_ = re.search(motif, str(seq))
        if motif_ != None:
            protein_motif = motif_.group()
            motif_position = motif_.start()
            peptide_number = pp[0].get_id()[1]
            start_position = peptide_number+motif_position
            peptide_chain = chain.id
            break
        else:
            peptide_chain = None
            start_position = None
            protein_motif = None
    return peptide_chain, start_position, protein_motif

def motif_pdb2(chain, motif):
    ''' Searches for the motif in the peptides chains of the PDB.
        If the motif is found gets the chain and the position
        where the motif starts.'''
    ppb = PPBuilder()
    for pp in ppb.build_peptides(chain):
        seq = pp.get_sequence()    
        motif_ = re.search(motif, str(seq))
        if motif_ != None and len(seq)<30:
            protein_motif = motif_.group()
            motif_position = motif_.start()
            peptide_number = pp[0].get_id()[1]
            start_position = peptide_number+motif_position
            peptide_chain = chain.id
            break
        else:
            peptide_chain = None
            start_position = None
            protein_motif = None
    return peptide_chain, start_position, protein_motif

def foldx_parameters(i, protein_pdb, motif, motif_name, df):
    ''' Create the Parameters for running FoldX.
    In: index of ELM dataframe with the proteins that have PDBs,
        ELM datafrmae, motif_name (ELMIdentifier), 
        dataframr of motifs detected.
    Out: pdbid of the crystal, chain that contains the motif,
        protein motif, and starting position of the motif in the crystal'''
    pdb = protein_pdb.loc[i, 'PDB']
    uniprotid = protein_pdb.loc[i, 'Primary_Acc']
    p = PDBParser()  
    structure = p.get_structure("X", pdb+'/'+pdb+'.pdb')
    model = structure[0]
    chains = dataframe_pdb(pdb, uniprotid)
    if len(chains) > 0:
        for chain in model.get_chains():
            if str(chain.id) in chains:
                peptide_chain, start_position, protein_motif = motif_pdb(chain, motif)
                if peptide_chain != None:
                    break
        return pdb, peptide_chain, protein_motif, start_position
    else:
        for chain in model.get_chains():
            peptide_chain, start_position, protein_motif = motif_pdb2(chain, motif)
            if peptide_chain != None:
                break
        return pdb, peptide_chain, protein_motif, start_position
                

def foldx_run(pdbid, chain, protein_motif, start_position, motif_name):
    ''' Run the FoldX software for each crystal containing the motif.
        In: pdbid, chain that contains the motif, motif instance of the protein,
            motif starting position (crystal numbers).
        Out: FoldX output.'''
    positions = []
    for i, e in enumerate(protein_motif):
        positions.append(e+chain+str(start_position+i)+'a')
    positions = ','.join(positions)
    os.chdir(pdbid.upper())
    command_line1 = ['FoldX', '--command=RepairPDB', '--pdb='+pdbid+'.pdb', '--screen=false']
    print('')
    print(f'Running FoldX for', pdbid, '...')
    process_1 = subprocess.run(command_line1, stdout=subprocess.DEVNULL)
    command_line2 = ['FoldX', '--command=PositionScan', '--pdb='+pdbid+'_Repair.pdb', '--clean-mode=2', '--positions='+positions, '--output-file='+pdbid, '--screen=false']
    process_2 = subprocess.run(command_line2, stdout=subprocess.DEVNULL)

def pulida_df(df):
    ''' Dataframe rearrengment '''
    df.drop(0, inplace=True)
    df = df.T
    df.drop(0, inplace=True)
    return df

def concat_df(dfs, motif):
    ''' Motif Position Mutational Scanning results concatenation.
        In: List of dataframes containing the results of FoldX mutational scanning.
            Dot separated motif (string).
        Out: Motif's FoldX dataframe'''
    aminoacids = ['G', 'A', 'L', 'V', 'I', 'P', 'R', 'T', 'S', 'C', 'M', 'K', 'E', 'Q', 'D', 'N', 'W', 'Y', 'F', 'H']
    df = pd.concat(dfs)
    df = df[df.columns].astype(float)
    df.columns = aminoacids
    df['Motif'] = motif
    df.set_index('Motif', inplace=True)
    return df
    
def foldx_matrix(motif, pdbid):
    ''' Builds Motif's Foldx Matrix. 
        In: dot separated motif and pdbid of the crystal (all string).
        Out: FoldX substitution matrix dataframe'''
    dfs = []
    n = 21
    motif = motif.split('.')
    for i in range(len(motif)): 
        df = pd.read_csv('PS_'+pdbid+'_scanning_output.txt', sep='\t', header=None)
        df = df[n-21:n]
        df.reset_index(inplace=True, drop=True)
        df = pulida_df(df)
        dfs.append(df)
        n += 21
    df_final = concat_df(dfs, motif)
    return df_final

def plot_save_foldx(motif, pdbid):
    ''' Plot the Motif's FoldX Matrix as a heatmap.
        In: dot separated motif, crystal pdbid and motif name (all string).
        Out: Motif's FoldX Matrix as a Heatmap. Dataframe.'''
    try:
        os.chdir(pdbid.upper())
    except:
        pass
    df = foldx_matrix(motif, pdbid)
    aminoacids = ['G', 'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    df = df[aminoacids]
    df = round(df, 1).replace(-0, 0)
    df.to_csv('FoldXMatrix-'+pdbid+'.csv', sep='\t')
    plt.figure(figsize=(15, 5))
    ax = sns.heatmap(round(df, 1), vmin=-5, vmax=5, cmap=r"PiYG_r", annot=True, linewidths=.5, cbar_kws={'label':'Kcal/mol'})
    sns.set(font_scale=1.2)
    plt.yticks(rotation=0)
    plt.ylabel('')
    plt.savefig('FoldXMatrix-'+pdbid+'.png')
    os.chdir('../')
    return df

def mean_foldx_matrix(crystals, motif_name):
    ''' Outputs the mean FoldX Matrix for a motif. 
        In: List of crystals were FoldX calculations were performed'''
    print("Building Motif's FoldX Substitution Matrix..")
    print('')
    df = pd.read_csv(crystals[0]+'/FoldXMatrix-'+crystals[0]+'.csv', sep='\t', index_col=0)
    for x in crystals[1:]:
        df += pd.read_csv(x+'/FoldXMatrix-'+x+'.csv', sep='\t', index_col=0)    
    
    df = df/len(crystals)
    df.reset_index(inplace=True)
    df.set_index('Motif', inplace=True)
    df.to_csv('MeanFoldXMatrix_'+motif_name+'.csv')
    
    plt.figure(figsize=(15, 5))
    ax = sns.heatmap(round(df, 1), vmin=-5, vmax=5, cmap=r"PiYG_r", annot=True, linewidths=.5, cbar_kws={'label':'Kcal/mol'})
    sns.set(font_scale=1.2)
    plt.yticks(rotation=0)
    plt.ylabel('')
    plt.savefig('MeanFoldXMatrix_'+motif_name+'.png')
    
    return df


    

