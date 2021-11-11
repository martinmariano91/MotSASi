
# MotSASi: Functional Short Linear Motifs (SLiMs) prediction based on genomic single nucleotide variants and structural data

Martín M, Modenutti CP, Nicola JP and Martí MA

https://www.biorxiv.org/content/10.1101/2021.08.05.455287v1

----

MotSASi is an strategy for improving the prediction of functional SLiMs using genomic variant information and structural data.
If you find this strategy useful, please cite us.

----

**For running MotSASi Python2.7, Python3.7, FoldX and ClustalO must be installed previously.**
**Biopython package must be previously installed too.**

MotSASi has to be run in different steps:

* 1 - MotSASi_1.py 
    -h or -help:
      Searchs the motif of interest in the human proteome. 
      Collects the benign and pathogenic variants reported within this motifs (ClinVar and GnomAD).
      Downloads the PDB crystal structures that contain the motif of interest.
      Calculates the substitution matrix based on the PDB structures using FoldX.
      Builds the ClinSig Matrix and the AF Matrix based on the variants reported for the instances of this motif in ELM (positive control).
      Calculates the ConservationScore for the positive control.
      Calculates the Secondary Structure Probability  using Jpred for the positive control.
      Calculates the relative position of the motif for the positive control.
      Collects GO terms for the positive control (from UniProt).
      
      For running: 
      MotSASi_1.py is a python3 script.
      Use the Unix Termial. Locate in the folder that contains this script and run it from there.
      In the command line put the regular expression of the motif that you want to study, this expression but separated by dots and the name of this motif.
      (Everything according to ELM nomenclature. i.e. python3 MotSASi_pipeline NP.[YF] N.P.x.YF TB_PTB)
    
    This script builds a "motif specific folder" where all the outputs are stored and is located in the previous folder from where the script is run.

* 2 - MotSASi_2.py (in progress)
    -h or -help:
      Gets the features of the positive control group of this motif (Tolerance Matrix, ConservationScore, Secondary Structure Probability, Relative Position).
      Evaluates if the reported benign or pathogenic variants within the motifs are in accordance or not with the Tolerance Matrix. 
      Not in accordance: discarded motif.
      In accordance: Calculates ConservationScore, Secondary Structure Probability and Relative Position. Collects GO terms from UniProt. 
      Checks if the features values are within the positive control group values.
      Calculates the "MotSASi Score": +++ (high confidence, in accordance with the Tolerance Matrix and same features as the positive control), ++ (medium confidence, in accordance with the Tolerance Matrix but some differences in GO terms), + (low confidence, in disagreeement with the Tolerance Matrix or several differences with the features values of the positive control), - (do not have benign or pathogenic variants within the motif).
      Adds this information to the "motif_name_motif.csv" file.
      
      For running:
      Use the same parameter arquitecture of the MotSASi_1.py script
      
      
      
    
    
