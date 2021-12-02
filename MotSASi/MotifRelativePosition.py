#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:46:00 2021

@author: mariano
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

def RelativePosition(df, positive_ctrl, motif_name):
    '''Calculates the relative position of the motif.
    In: dataset that have the info of the motifs detected in the
        human proteome, the motif name and the dataset that have 
        the info of the positive ctrl.
    Out: Dataframe with the relative position of each of the 
        the motif of the positive ctrl. Plot a histogram with the
        distribution of the positions of the motif of the positive
        control.'''
    
    motifs = df
    
    UniMotifs = [x+'_'+str(i) for x, i in positive_ctrl]
    
    positives = motifs[motifs.ID_motif.isin(UniMotifs)]
    
    index = []
    for x, p in positive_ctrl:
        for i in positives.index:
            if positives.UniProtID.loc[i] == x and positives.Start.loc[i] == p:
                index.append(i)
                
    positives = positives.loc[index]
    
    tuples = tuple(zip(positives.UniProtID.tolist(), positives.Start.tolist(), positives.Protein_length.tolist()))
    
    norm_dist = []
    ids = []
    starts = []
    for p, m, l in tuples:
        norm_dist.append(m/l)
        ids.append(p)
        starts.append(m)
    df = pd.DataFrame()
    df['UniProtID'] = ids
    df['Start'] = starts
    df['RelPos'] = norm_dist
    df.to_csv('../'+motif_name+'_motif/RelPos_'+motif_name+'.csv')
    
    plt.figure(facecolor='w')
    sns.kdeplot(norm_dist)
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.xlabel('Relative Position')
    plt.tight_layout()
    plt.savefig('../'+motif_name+'_motif/'+motif_name+'_RelPos.png', bbox_inches = "tight")
    
    
    
