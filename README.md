
# MotSASi: Functional Short Linear Motifs (SLiMs) prediction based on genomic single nucleotide variants and structural data

Martín M, Brunello FG, Modenutti CP, Nicola JP and Martí MA

Pre-print:
https://www.biorxiv.org/content/10.1101/2021.08.05.455287v1

(under major revision in *Biochimie*)

----

MotSASi is a strategy for improving the prediction of functional SLiMs using genomic variant information and structural data.

If you find this strategy useful, please cite us.

----

**For running MotSASi Python2.7, Python3.7, C-Shell (csh) package, Biopython package, FoldX and ClustalO must be installed previously.**

It is recommended to run the scripts within Anaconda Environment. The Jpred API is included in the MotSASi folder.


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
      
    <ins>For running:</ins>
    
    MotSASi_1.py is a python3 script.
    
    Use the Unix Termial. Locate in the folder that contains this script and run it from there.
    
    In the command line put the regular expression of the motif that you want to study, this expression but separated by dots and the ELM name of this motif.
    (Everything according to ELM nomenclature. i.e. python3 MotSASi_1.py [RK]P[^P][^P]L.[LIVMF] RK.P.^P.^P.L.x.LIVMF DOC_MAPK_JIP1_4)
    
    This script builds a "motif specific folder" where all the outputs are stored and is located in the previous folder from where the script is run.   
