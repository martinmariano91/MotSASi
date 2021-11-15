#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 17:16:36 2021

@author: mariano
"""
import pandas as pd
from Bio import SeqIO
import gzip
import re
import collections
import urllib
import matplotlib.pyplot as plt
import seaborn as sns

def GO_terms_motif(uni_id):
    '''Collects some of the GO terms of the proteins of the positive
        control of the motif under study.
        In: UniProtID
        Out: List of GO terms'''
    GO_terms = []
    #handle = gzip.open("uniprot_sprot.xml.gz", "rt")
    handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+uni_id+".xml")
    record = SeqIO.read(handle, "uniprot-xml")
    #for record in SeqIO.parse(handle, "uniprot-xml"):
    if record.annotations['organism'] == 'Homo sapiens (Human)':
        keys = record.annotations.keys()
        if 'keywords' in keys:
            GO_terms += record.annotations['keywords']
        elif 'comment_subcellularlocation_location' in keys:
            GO_terms += record.annotations['comment_subcellularlocation_location']
        elif 'comment_subcellularlocation_topoplogy' in keys:
            GO_terms += record.annotations['comment_subcellularlocation_topoplogy']
        elif 'comment_function' in keys:
            GO_terms += record.annotations['comment_function']
        elif 'comment_PTM' in keys:
            GO_terms += record.annotations['comment_PTM']
    return GO_terms

def Motif_GO(positives_ctrl, motif_name):
    ''' Builds a dataframe with the GO terms of each of the proteins 
    of the positive control and the percentaje of the proteins that 
    contain each GO term. Plot bars of the Top10 Go terms.
    In: dataframe with the info of the positive control and motif name.
    Out: Dataframe and Plot with the info of GO terms of the positive
    control.'''
    UniMotifs = [x for x, i in positives_ctrl]
    GO_terms = []
    for x in UniMotifs:
        GO_terms += GO_terms_motif(x)
    
    GO_terms_counts = collections.Counter(GO_terms)
    GO_terms_counts = dict(sorted(GO_terms_counts.items(), key=lambda item: item[1]))
    
    GO = pd.DataFrame(GO_terms_counts, index=['%'])
    GO = GO/len(UniMotifs)*100
    GO.to_csv('../'+motif_name+'_motif/GOterms_'+motif_name+'_positives.csv')
    
    plt.figure()
    plt.title('GO terms')
    plt.barh(GO.columns[-10:], GO.iloc[0, -10:])
    plt.yticks(range(10), GO.columns.tolist()[-10:], fontsize=10)
    plt.tight_layout()
    plt.savefig('../'+motif_name+'_motif/'+'GOterms_'+motif_name+'_positives.png')
            
    


        
