#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 15:20:23 2021

Searchs in the human proteome for the motif specified with a regular expression.
Outputs a 'csv' file with all the motifs found with the information about
the protein that contains it and the position where it is located.

Searchs for ClinVar variants occurring within the motifs of study. 
Variants are grouped according to the reported Clinical Significance (ClinSig).
Only missense variants with Criteria Provided and no conflicting 
interpretations are considered. 
The variant position in the protein sequence is changed to the position
relative to the motif.

Searchs for GnomAD variants occurring within the motifs of study. 
Only missense variants are considered. The variant position in the protein 
sequence is changed to the position relative to the motif.
Variants with -log10(AF) <= 4 are considered as benign.

@author: Mariano Martín
"""

import sys
import os
import shutil
import MotifSearcher
import warnings
warnings.filterwarnings("ignore")

if sys.argv[1] == '-h' or sys.argv[1] == '-help':
    print(str('#'*50).center(140))
    print('MotSASi'.center(140))
    print('Martín et al. 2021'.center(140))
    print(str('#'*50).center(140))
    print('')
    print('*'*300)
    print('- Searchs the motif of interest in the human proteome.')
    print('- Collects the benign and pathogenic variants reported within this motifs (ClinVar and GnomAD).')
    print('- Downloads the PDB crystal structures that contain the motif of interest.')
    print('- Calculates the substitution matrix based on the PDB structures using FoldX.')
    print('- Builds the ClinSig Matrix and the AF Matrix based on the variants reported for the instances of this motif in ELM (positive control).')
    print('- Calculates the ConservationScore for the positive control.')
    print('- Calculates the Secondary Structure Probability  using Jpred for the positive control.')
    print('- Calculates the relative position of the motif for the positive control.')
    print('*'*300)
    print('MotSASi is a python3 script.')
    print('Use the Unix Termial. Locate in the folder that contains this script and run it from there.')
    print('In the command line put the regular expression of the motif that you want to study, this expression but separated by dots and the ELM name of this motif.')
    print('(Everything according to ELM nomenclature. i.e. python3 MotSASi_1.py [RK]P[^P][^P]L.[LIVMF] RK.P.^P.^P.L.x.LIVMF DOC_MAPK_JIP1_4')
    print('')
else: 
    print('')
    print(str('#'*50).center(140))
    print('MotSASi'.center(140))
    print('Martín et al. 2021'.center(140))
    print(str('#'*50).center(140))
    print('')
    
    motif = sys.argv[1]
    mot = sys.argv[2]
    motif_name = sys.argv[3]
    
    print("Searching "+motif_name+" motif in the human proteome...")
    
    #try:
    #    shutil.rmtree('../'+motif_name+'_motif')
    #except:
    #    pass
    
    #os.mkdir('../'+motif_name+'_motif')
    
    #df = MotifSearcher.HumanProteomeSearch(motif, motif_name)
    import pandas as pd
    df = df = pd.read_csv('../LIG_PCNA_PIPBox_1_motif/LIG_PCNA_PIPBox_1_motif.csv', index_col=0, sep='\t')
    
    if df.empty:
        print('')
        print('No Hits found in the Human Proteome!')
        print('Look for another motif!')
        print('')
    else:
        print('')
        print(f'Searching for variants affecting these motifs in ClinVar...')
        #import ClinVarMotif
        #ClinVarMotif.ClinVarMotif(df, motif_name)
        
        print('')
        print(f'Searching for variants affecting these motifs in GnomAD...')
        #import GnomADMotif
        #GnomADMotif.GnomADMotif(df, motif_name)
        
        print('')
        import PDB_download
        #elm_instances = pd.read_csv('ELM_instances.csv', index_col=0)
        #elm_instances = elm_instances[(elm_instances['ELMIdentifier'] == motif_name) &
        #                             (elm_instances['Organism'] == 'Homo sapiens')]
        #positive_ctrl = elm_instances['Primary_Acc'].tolist()
        #protein_pdb = elm_instances.loc[elm_instances.PDB.notnull(), ['Primary_Acc', 'PDB']]
        positive_ctrl, protein_pdb = PDB_download.positive_ctrl(motif_name)
        #structure = PDB_download.pdbs_download(motif_name, protein_pdb)
        structure = 'YES'
        print('')
        if structure == 'NO':
            print('There are no PDBs for this motif!')
            print('MotSASi cannot be run...')
            print('Bye Bye!')
            quit()
        else:
            import FoldXMotif
            os.chdir('../'+motif_name+'_motif/pdbs')
            for i in protein_pdb.index:
                pdb, peptide_chain, protein_motif, start_position = FoldXMotif.foldx_parameters(i, protein_pdb, motif_name, df)                
                if peptide_chain != '':
                    FoldXMotif.foldx_run(pdb, peptide_chain, protein_motif, start_position, motif_name)
                    FoldXMotif.plot_save_foldx(mot, pdb)
                else:
                    shutil.rmtree(pdb+'/'+pdb+'.pdb')
            crystals = os.listdir()
            FoldXMotif.mean_foldx_matrix(crystals, motif_name)
        try:
            os.chdir('../../MotSASi_pipeline')
        except:
            pass
        print('Building ClinSig Matrix...')
        
        import ClinSigMatrix
        ClinSigMatrix.ClinSigMatrix(df, mot, positive_ctrl, motif_name)
        print('Building AF Matrix...')
        import GnomADMatrix
        GnomADMatrix.GnomADMatrix(df, mot, positive_ctrl, motif_name)
        
        print('Calculating the Conservation Score...')
        import ConservationScore
        ConservationScore.ConservationScore(df, mot, motif_name, positive_ctrl)
        
        print('Calculating the Relative Position of the Motif...')
        import MotifRelativePosition
        MotifRelativePosition.RelativePosition(df, positive_ctrl, motif_name)
        
        import GOMotif
        GOMotif.Motif_GO(positive_ctrl, motif_name)
        
        import JPREDMotif
        JPREDMotif.jpred_run(positive_ctrl, mot, motif_name)
        JPREDMotif.jpred_batch_processing(mot, motif_name)