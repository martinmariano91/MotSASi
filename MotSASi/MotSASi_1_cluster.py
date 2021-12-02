#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 09:20:22 2021

@author: mariano
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
    print('- Collects the benign and pathogenic variants reported within this motif (ClinVar and GnomAD).')
    print('- Downloads the PDB crystal structures that contain the motif of interest.')
    print('- Calculates the substitution matrix based on the PDB structures using FoldX.')
    print('- Builds the ClinSig Matrix and the AF Matrix based on the variants reported for the instances of this motif in ELM (positive control).')
    print('- Calculates the ConservationScore for the positive control.')
    print('- Calculates the Secondary Structure Probability  using Jpred for the positive control.')
    print('- Calculates the relative position of the motif for the positive control.')
    print('*'*300)
    print('MotSASi is a python3 script.')
    print('Use the Unix Terminal. Locate in the folder that contains this script and run it from there.')
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
    
    if '{' in motif or '}' in motif:
        print('Motifs that have flexible length are not analyzable in this version of MotSASi!!')
        quit()
    else:
        pass
    
    import pandas as pd
    df = pd.read_csv('../'+motif_name+'_motif/'+motif_name+'_motif.csv', index_col=0, sep='\t')
    
    if df.empty:
        print('')
        print('No Hits found in the Human Proteome!')
        print('Look for another motif!')
        print('')
    else:
        print('')
        import PDB_download
        positive_ctrl, protein_pdb, df = PDB_download.positive_ctrl(df, motif_name)
        if len(positive_ctrl) == 0:
            print('There are no positive control for this motif!')
            print('MotSASi cannot be run...')
            print('Bye Bye!')
            quit()
        else:
            pass
        
        print(f'Searching for variants affecting these motifs in ClinVar...')
        import ClinVarMotif
        ClinVarMotif.ClinVarMotif(df, motif_name)
        print('')
    
        print(f'Searching for variants affecting these motifs in GnomAD...')
        import GnomADMotif
        GnomADMotif.GnomADMotif(df, motif_name)
      
        if structure == 'NO':
            print('There are no PDBs for this motif!')
            print('MotSASi cannot be run...')
            print('Bye Bye!')
            quit()
        else:
            import FoldXMotif
            os.chdir('../'+motif_name+'_motif/pdbs')
            for i in protein_pdb.index:
                pdb, peptide_chain, protein_motif, start_position = FoldXMotif.foldx_parameters(i, protein_pdb, motif, motif_name, df)                
                if peptide_chain != None:
                    FoldXMotif.foldx_run(pdb, peptide_chain, protein_motif, start_position, motif_name)
                    FoldXMotif.plot_save_foldx(mot, pdb) 
                else:
                    print("There's no coincidence between the peptide in the crystal and the motif reported for the protein")                 
                    shutil.rmtree(pdb)
            crystals = [x for x in os.listdir() if not x.startswith('.')]
            FoldXMotif.mean_foldx_matrix(crystals, motif_name)
        
        try:
            os.chdir('../../MotSASi_pipeline')
        except:
            pass
        
        print('Building ClinSig Matrix...')
        import ClinSigMatrix
        ClinSigMatrix.ClinSigMatrix(df, mot, positive_ctrl, motif_name)
        print('')
        
        print('Building AF Matrix...')
        import GnomADMatrix
        GnomADMatrix.GnomADMatrix(df, mot, positive_ctrl, motif_name)
        print('')
        
        print('Calculating the Conservation Score...')
        import ConservationScore
        ConservationScore.ConservationScore(df, mot, motif_name, positive_ctrl)
        print('')
        
        print('Calculating the Relative Position of the Motif...')
        import MotifRelativePosition
        MotifRelativePosition.RelativePosition(df, positive_ctrl, motif_name)
        print('')
        
        print('Collecting GO terms...')
        import GOMotif
        GOMotif.Motif_GO(positive_ctrl, motif_name)
        print('')