#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:41:45 2021

Itera de distintas maneras sobre las variantes de un motivo.

Debo tener un DataFrame con las variantes que hay documentadas en los motivos.
(df_final_clinvar y df_final_gnomad)

Definir inicialmente: nombre_archivo y mot

la función pato_ben() requiere (además de mot), el tipo de iteracion que se va 
realizar (iteracion inicial (numeros), foldx, flexible o final), y el tipo de archivo
que se usa para variantes clinvar y gnomad (construido con controles positivos,
luego de X iteración, etc.).

@author: mariano
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% 
nombre_archivo='PDZtipo1_PDZ_complejo/Proteinas_motivo_PDZtipo1_PDZ.csv' #archivo con las proteínas en estudio.
mot='S/T.x.V/I/L' #poner cada posición del motivo separada por un punto ('.').
mot=mot.split('.')
   
nombre_motivo='_'.join(nombre_archivo.split('/')[1].split('.')[0].split('_')[2:])

d = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys',
     'I': 'Ile', 'P': 'Pro', 'T': 'Thr', 'F': 'Phe', 'N': 'Asn', 
     'G': 'Gly', 'H': 'His', 'L': 'Leu', 'R': 'Arg', 'W': 'Trp', 
     'A': 'Ala', 'V':'Val', 'E': 'Glu', 'Y': 'Tyr', 'M': 'Met'}

positivos = [('Q96A54', 0), ('Q8TEU7', 0), ('O95477', 0), ('Q9UQB3', 0), 
             ('P16473', 0), ('P32241', 0), ('O00548', 0), ('P53778', 0), 
             ('P07550', 0), ('Q8WXS5', 0), ('Q01668', 0), ('Q00722', 0), 
             ('Q9UBN4', 0), ('P11166', 0), ('P36507', 0), ('Q13936', 0),
             ('P34998', 0), ('P35222', 0), ('Q14524', 0), ('P47900', 0), 
             ('P25025', 0), ('Q8N695', 0), ('P60484', 0), ('Q9NQ66', 0), 
             ('Q9UL62', 0), ('Q92911', 0), ('Q92806', 0), ('O75096', 0), 
             ('P46937', 0), ('Q15046', 0), ('P48960', 0), ('Q13268', 0), 
             ('P25054', 0), ('Q9Y5I7', 0), ('P98164', 0), ('Q86XX4', 0), 
             ('Q9NSA0', 0), ('Q9H251', 0), ('O95436', 0), ('Q5JU85', 0), 
             ('P23634', 0), ('P11274', 0), ('Q13224', 0), ('Q8TD84', 0), 
             ('O00192', 0), ('P30872', 0), ('P48050', 0), ('P28223', 0), 
             ('P09619', 0), ('Q04656', 0), ('P35228', 0), ('Q9UHC3', 0), 
             ('Q14916', 0), ('Q495M9', 0), ('P43250', 0), ('P42261', 0), 
             ('O43424', 0), ('P23508', 0), ('Q9Y6M7', 0), ('P46721', 0), 
             ('Q96BN8', 0), ('P08588', 0), ('P35609', 0), ('P13569', 0), 
             ('Q9P2W3', 0), ('Q9P0K1', 0), ('O60469', 0), ('Q9UKN7', 0), 
             ('Q01814', 0), ('O43572', 0), ('Q99527', 0)]

positivos = [x+'_0' for x, i in positivos]

# Los siguientes archivos los tengo que construir previamente.
# Debo juntar todas las variantes que caigan en los motivos predichos.
df_final_clinvar = pd.read_csv(nombre_motivo+'_complejo/Variantes_ClinVar_Motivo_full.csv', index_col=0)
df_final_clinvar = df_final_clinvar[~df_final_clinvar.ID.isin(positivos)]

df_final_gnomad = pd.read_csv(nombre_motivo+'_complejo/Variantes_GnomAD_Motivo_full.csv', index_col=0)
df_final_gnomad = df_final_gnomad[~df_final_gnomad.ID.isin(positivos)]


#%% 

def pato_ben(mot, nombre_motivo, iteracion, archivo='positivos'):
    """ Esta función permite construir listas de variantes benignas o patogénicas,
    que permitirán filtrar los motivos predichos.    
    
    Argumentos:
        mot: es el motivo separando las posiciones por puntos (.).
        nombre_motivo: se desprende del nombre del archivo.
        iteracion: int o foldx, flexible, final.
        archivo: nombre final de las matrices ClinVar y GnomAD que se usan
                para armar los filtros de sustituciones.
    
    Devuelve:
        dos listas: patogenicos y benignos con las sustituciones que cuadran en 
                    cada grupo según la información de matrices utilizada
        graficos: devuelve un heatmap de la matriz de sustitución marcando cuales
                    cambios son los que se encuentran en las listas entregadas"""
    
    pos_x = mot.index('x') # Posicion Flexible del Motivo
    
    patogenicos_1 = [] # derivado de la matriz ClinVar
    benignos_1 = [] # derivado de la matriz ClinVar
    benignos_2 = [] # derivado de la matriz GnomAD
    patogenicos_2 = [] # derivado de la matriz FoldX
    benignos_3 = [] # derivado de la matriz FoldX
    
    clinsig_positivos = pd.read_csv(nombre_motivo+'_complejo/ClinVar/MatrizClinSigMedias_'+nombre_motivo+'_'+archivo+'.csv', sep='\t', index_col=0)

    for i in range(len(mot)): 
        for e in clinsig_positivos:
            if clinsig_positivos.iloc[i][e] > 0: # patogenicos            
                if (iteracion !='final' or iteracion != 'flexible') and i != pos_x:
                    aa = mot[i].split('/')
                    for z in aa:
                        patogenicos_1.append(d[z]+str(i+1)+d[e])
                elif (iteracion == 'final' or iteracion == 'flexible') and i == pos_x:
                    for x in d.values():
                        patogenicos_1.append(x+str(i+1)+d[e])
            elif clinsig_positivos.iloc[i][e] < 0: # benignos            
                if (iteracion !='final' or iteracion != 'flexible') and i != pos_x:
                    aa = mot[i].split('/')
                    for z in aa:
                        benignos_1.append(d[z]+str(i+1)+d[e])
                elif (iteracion == 'final' or iteracion == 'flexible') and i == pos_x:
                    for x in d.values():
                        benignos_1.append(x+str(i+1)+d[e])
    
    clinsig_positivos[clinsig_positivos>0] = 1
    clinsig_positivos[clinsig_positivos<0] = -1
    clinsig_positivos.fillna(0, inplace=True)
    
    gnomad_positivos = pd.read_csv(nombre_motivo+'_complejo/GnomAD/MatrizFrecuenciasMedias_'+nombre_motivo+'_'+archivo+'.csv', sep='\t', index_col=0)

    for i in range(len(mot)): 
        for e in gnomad_positivos:
            if gnomad_positivos.iloc[i][e] < 4: # benignos            
                if (iteracion != 'final' or iteracion != 'flexible') and i != pos_x:
                    aa = mot[i].split('/')
                    for z in aa:
                        benignos_2.append(d[z]+str(i+1)+d[e])
                elif (iteracion == 'final' or iteracion == 'flexible') and i == pos_x:
                    for x in d.values():
                        benignos_2.append(x+str(i+1)+d[e])
    
    gnomad_positivos[gnomad_positivos<4] = -1
    gnomad_positivos.iloc[0,14] = 0
    gnomad_positivos.iloc[2,14] = 0
    gnomad_positivos[gnomad_positivos>4] = 0
    gnomad_positivos.fillna(0, inplace=True)  
   
    final_matrix = clinsig_positivos + gnomad_positivos
    if iteracion != 'final' and iteracion != 'flexible':
        final_matrix.loc['x'] = 0
    
    if type(iteracion) != int:
        foldx = pd.read_csv(nombre_motivo+'_complejo/FoldX/MatrizFoldX-'+nombre_motivo+'_promedio.csv', index_col=0)
        foldx.rename({'X':'x'}, axis=0, inplace=True)        
        for i in range(len(mot)): 
            for e in foldx:
                if foldx.iloc[i][e] >= 2: # patogenicos           
                    if (iteracion != 'final' or iteracion != 'flexible') and i != pos_x:
                        aa = mot[i].split('/')
                        for z in aa:
                            patogenicos_2.append(d[z]+str(i+1)+d[e])
                    elif iteracion == 'final' or iteracion == 'flexible':
                        if i == pos_x:
                            for x in d.values():
                                patogenicos_2.append(x+str(i+1)+d[e])
                        else:
                            aa = mot[i].split('/')
                            for z in aa:
                                patogenicos_2.append(d[z]+str(i+1)+d[e])
                elif foldx.iloc[i][e] <= 1: # benignos
                    if (iteracion !='final' or iteracion != 'flexible') and i != pos_x:
                        aa = mot[i].split('/')
                        for z in aa:
                            benignos_3.append(d[z]+str(i+1)+d[e])
                    elif iteracion == 'final' or iteracion == 'flexible':
                        if i == pos_x:
                            for x in d.values():
                                benignos_3.append(x+str(i+1)+d[e])
                        else:
                            aa = mot[i].split('/')
                            for z in aa:
                                benignos_3.append(d[z]+str(i+1)+d[e])
                # Tunear esta parte
                elif foldx.iloc[i][e] < 2 and (iteracion == 'final' or iteracion == 'flexible'):
                    if i == 0:
                        if foldx.iloc[i][e] > 1:
                            aa = mot[i].split('/')
                            for z in aa:
                                patogenicos_2.append(d[z]+str(i+1)+d[e])
                    elif i == pos_x:
                        for x in d.values():
                            benignos_3.append(x+str(i+1)+d[e])   
                    elif i == 2:
                        aa = mot[i].split('/')
                        for z in aa:
                            benignos_3.append(d[z]+str(i+1)+d[e])
        
        if iteracion == 'foldx':
            foldx[foldx<=1] = -0.25
            foldx[foldx>=2] = 0.25
            foldx[(foldx<2) & (foldx>1)] = 0
            foldx = foldx[final_matrix == 0]
            foldx.fillna(0, inplace=True)
            final_matrix += foldx
            final_matrix.loc['x'] = 0
        elif iteracion == 'flexible' or iteracion == 'final':
            foldx[foldx<=1] = -0.25
            foldx[foldx>=2] = 0.25
            foldx.iloc[0][(foldx.iloc[0]<2) & (foldx.iloc[0]>1)] = 0.25
            foldx.iloc[1][(foldx.iloc[1]<2) & (foldx.iloc[1]>1)] = -0.25
            foldx.iloc[2][(foldx.iloc[2]<2) & (foldx.iloc[2]>1)] = -0.25
            foldx = foldx[final_matrix == 0]
            foldx.fillna(0, inplace=True)
            final_matrix += foldx
    
    final_matrix.to_csv(nombre_motivo+'_complejo/Matriz_'+nombre_motivo+'_iteracion'+str(iteracion)+'.csv')
    
    plt.figure(figsize=(15, 5))
    ax = sns.heatmap(final_matrix, vmin=-1, vmax=1, cmap=['mediumseagreen','lightgreen','lightgrey', 'lightsalmon','indianred'], linewidths=.5, alpha=0.7)
    sns.set(font_scale=1.2)
    plt.yticks(rotation=0)
    plt.ylabel('')
    plt.savefig(nombre_motivo+'_complejo/Matriz_HeatMap_'+nombre_motivo+'_iteracion'+str(iteracion)+'.png')
        
    # En la lista final le doy prioridad a lo que dice ClinVar y GnomAD en cuanto a la significancia
    # clínica de las variantes en motivos. Luego le sumo lo que me da FoldX
    # evitando que haya variantes solapadas o inconsistentes.
    patogenicos = list(set(patogenicos_1 + [x for x in patogenicos_2 if x not in benignos_1 and x not in benignos_2]))
    benignos = list(set(benignos_1 + benignos_2 + [x for x in benignos_3 if x not in patogenicos]))
    
    return (patogenicos, benignos)

#%%
iteracion = 1
patogenicos, benignos = pato_ben(mot, nombre_motivo, iteracion)

# Filtro según las varianates de las matrices.
# IDS_falso junta todas las proteínas en las cuales se predice un motivo pero
# que este no sería funcional al contrastar con la info de las matrices.
# Es decir aquellas que en ClinVar son Benignas pero que por las matrices sabemos
# que son patogénicas y aquellas que en ClinVar son Patogénicas pero que por
# matrices sabemos que son benignas.
IDS_falso1 = df_final_clinvar.ID[df_final_clinvar.Benign.isin(patogenicos)].unique().tolist()
IDS_falso1 += df_final_clinvar.ID[df_final_clinvar.Pathogenic.isin(benignos)].unique().tolist()
IDS_falso1 = list(set(IDS_falso1))
with open(nombre_motivo+'_complejo/IDS_Falsos_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_falso1))
                
# IDS_real junta todas las proteínas que verdaderamente tendrían un motivo funcional
# Para eso me quedo con aquellas proteínas que tienen variantes patogénicas que coinciden
# con las que tengo en las matrices. Luego elimino aquellas que tienen variantes 
# benignas en ClinVar pero que se que son patogénicas por las matrices.
IDS_real1 = df_final_clinvar.ID[df_final_clinvar.Pathogenic.isin(patogenicos)].unique().tolist()
IDS_real1 += df_final_clinvar.ID[df_final_clinvar.Benign.isin(benignos)].unique().tolist()
IDS_contradict1 = [x for x in IDS_real1 if x in IDS_falso1]
IDS_real1 = [x for x in IDS_real1 if x not in IDS_falso1]
IDS_real1 = list(set(IDS_real1))
with open(nombre_motivo+'_complejo/IDS_Reales_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_real1))

# IDS_contradict junta todas las proteínas que tienen info contradictoria
with open(nombre_motivo+'_complejo/IDS_Contradict_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_contradict1))

# Luego le sumo los detectados mediante gnomad
IDS_falso1 += df_final_gnomad.ID[df_final_gnomad.Benign.isin(patogenicos)].unique().tolist()
IDS_falso1 = list(set(IDS_falso1))
with open(nombre_motivo+'_complejo/IDS_Falsos_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_falso1))

IDS_real1 += df_final_gnomad.ID[df_final_gnomad.Benign.isin(benignos)].unique().tolist()
IDS_contradict2 = [x for x in IDS_real1 if x in IDS_falso1]
IDS_real1 = [x for x in IDS_real1 if x not in IDS_falso1]
IDS_real1 = list(set(IDS_real1))
with open(nombre_motivo+'_complejo/IDS_Reales_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_real1))

# IDS_contradict junta todas las proteínas que tienen info contradictoria
with open(nombre_motivo+'_complejo/IDS_Contradict_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_contradict1+IDS_contradict2))

#En este punto debo chequear los predichos como positivos 
#para ver si los puedo sumar para armar las nuevas matrices.
#con los que sumo correr MatSub-B-P y MatSub-AF y ver
#como cambian las matrices.

#Guardo aquellos que fueron predichos como reales pero que termine descartando a mano                
IDS_descartados = ['Q9HBX8_0', 'Q6TDP4_0', 'Q8IZ08_0', 'Q06495_0', 
                   'Q9H5J0_0', 'Q53R41_0', 'Q6ZPB5_0', 'P24347_0']
with open(nombre_motivo+'_complejo/IDS_Descartados_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_descartados))
                
#%%
iteracion = 'foldx'
iteracion_previa = 1
patogenicos, benignos = pato_ben(mot, nombre_motivo, iteracion, archivo='iteracion1')
# agrego esta linea debido a la duda que me surge de la matriz foldx

with open(nombre_motivo+'_complejo/IDS_Descartados_iteracion'+str(iteracion_previa)+'.txt', "r") as output:
    linea = output.readlines()
    IDS_descartados = linea[0].split("'")[1:-1]

IDS_descartados = list(set(IDS_descartados))
IDS_descartados.remove(', ')

# Elimino de los dfs aquellos motivos que ya descarté en etapas anteriores
# y que curé manualmente.
IDS_siguen = [x for x in df_final_clinvar.ID.tolist() if x not in IDS_descartados]
df_final_clinvar = df_final_clinvar[df_final_clinvar.ID.isin(IDS_siguen)]

IDS_siguen = [x for x in df_final_gnomad.ID.tolist() if x not in IDS_descartados]
df_final_gnomad = df_final_gnomad[df_final_gnomad.ID.isin(IDS_siguen)]

# Filtro según las varianates de las matrices.
# IDS_falso junta todas las proteínas en las cuales se predice un motivo pero
# que este no sería funcional al contrastar con la info de las matrices.
# Es decir aquellas que en ClinVar son Benignas pero que por las matrices sabemos
# que son patogénicas y aquellas que en ClinVar son Patogénicas pero que por
# matrices sabemos que son benignas.
IDS_falso2 = df_final_clinvar.ID[df_final_clinvar.Benign.isin(patogenicos)].unique().tolist()
IDS_falso2 += df_final_clinvar.ID[df_final_clinvar.Pathogenic.isin(benignos)].unique().tolist()
IDS_falso2 = list(set(IDS_falso2))
with open(nombre_motivo+'_complejo/IDS_Falsos_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_falso2))
                
# IDS_real junta todas las proteínas que verdaderamente tendrían un motivo funcional
# Para eso me quedo con aquellas proteínas que tienen variantes patogénicas que coinciden
# con las que tengo en las matrices. Luego elimino aquellas que tienen variantes 
# benignas en ClinVar pero que se que son patogénicas por las matrices.
IDS_real2 = df_final_clinvar.ID[df_final_clinvar.Pathogenic.isin(patogenicos)].unique().tolist()
IDS_real2 += df_final_clinvar.ID[df_final_clinvar.Benign.isin(benignos)].unique().tolist()
IDS_contradict3 = [x for x in IDS_real2 if x in IDS_falso2+IDS_falso1]
IDS_real2 = [x for x in IDS_real2 if x not in IDS_falso2+IDS_falso1]
IDS_real2 = list(set(IDS_real2))
with open(nombre_motivo+'_complejo/IDS_Reales_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_real2))

# IDS_contradict junta todas las proteínas que tienen info contradictoria
with open(nombre_motivo+'_complejo/IDS_Contradict_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_contradict3))

# Luego le sumo los detectados mediante gnomad
IDS_falso2 += df_final_gnomad.ID[df_final_gnomad.Benign.isin(patogenicos)].unique().tolist()
IDS_falso2 = list(set(IDS_falso2))
with open(nombre_motivo+'_complejo/IDS_Falsos_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_falso2))

IDS_real2 += df_final_gnomad.ID[df_final_gnomad.Benign.isin(benignos)].unique().tolist()
IDS_contradict4 = [x for x in IDS_real2 if x in IDS_falso2+IDS_falso1]
IDS_real2 = [x for x in IDS_real2 if x not in IDS_falso2+IDS_falso1]
IDS_real2 = list(set(IDS_real2))
with open(nombre_motivo+'_complejo/IDS_Reales_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_real2))

# IDS_contradict junta todas las proteínas que tienen info contradictoria
with open(nombre_motivo+'_complejo/IDS_Contradict_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_contradict3+IDS_contradict4))
    
IDS_descartados += ['Q96BD0_0', 'Q96B33_0', 'Q5SYB0_0', 'Q9Y5Y0_0', 'Q6PEX7_0', 
                    'Q9NR81_0', 'P11169_0', 'Q5T6X5_0', 'P17655_0', 'Q9HBX8_0',
                    'Q6ZVD8_0', 'O00322_0', 'Q92608_0', 'Q9HD33_0', 'Q6ZTZ1_0',
                    'Q9H5J0_0', 'Q99645_0', 'Q5GAN3_0']

with open(nombre_motivo+'_complejo/IDS_Descartados_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_descartados))

#%%
iteracion = 'flexible'
iteracion_previa = 'foldx'
patogenicos, benignos = pato_ben(mot, nombre_motivo, iteracion, archivo='iteracionfoldx')

with open(nombre_motivo+'_complejo/IDS_Descartados_iteracion'+str(iteracion_previa)+'.txt', "r") as output:
    linea = output.readlines()
    IDS_descartados = linea[0].split("'")[1:-1]

IDS_descartados = list(set(IDS_descartados))
IDS_descartados.remove(', ')

# Elimino de los dfs aquellos motivos que ya descarté en etapas anteriores
# y que curé manualmente.
IDS_siguen = [x for x in df_final_clinvar.ID.tolist() if x not in IDS_descartados]
df_final_clinvar = df_final_clinvar[df_final_clinvar.ID.isin(IDS_siguen)]

IDS_siguen = [x for x in df_final_gnomad.ID.tolist() if x not in IDS_descartados]
df_final_gnomad = df_final_gnomad[df_final_gnomad.ID.isin(IDS_siguen)]

# Filtro según las varianates de las matrices.
# IDS_falso junta todas las proteínas en las cuales se predice un motivo pero
# que este no sería funcional al contrastar con la info de las matrices.
# Es decir aquellas que en ClinVar son Benignas pero que por las matrices sabemos
# que son patogénicas y aquellas que en ClinVar son Patogénicas pero que por
# matrices sabemos que son benignas.
IDS_falso3 = df_final_clinvar.ID[df_final_clinvar.Benign.isin(patogenicos)].unique().tolist()
IDS_falso3 += df_final_clinvar.ID[df_final_clinvar.Pathogenic.isin(benignos)].unique().tolist()
IDS_falso3 = list(set(IDS_falso3))
with open(nombre_motivo+'_complejo/IDS_Falsos_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_falso3))
                
# IDS_real junta todas las proteínas que verdaderamente tendrían un motivo funcional
# Para eso me quedo con aquellas proteínas que tienen variantes patogénicas que coinciden
# con las que tengo en las matrices. Luego elimino aquellas que tienen variantes 
# benignas en ClinVar pero que se que son patogénicas por las matrices.
IDS_real3 = df_final_clinvar.ID[df_final_clinvar.Pathogenic.isin(patogenicos)].unique().tolist()
IDS_real3 += df_final_clinvar.ID[df_final_clinvar.Benign.isin(benignos)].unique().tolist()
IDS_contradict5 = [x for x in IDS_real3 if x in IDS_falso1+IDS_falso2+IDS_falso3]
IDS_real3 = [x for x in IDS_real3 if x not in IDS_falso1+IDS_falso2+IDS_falso3]
IDS_real3 = list(set(IDS_real3))
with open(nombre_motivo+'_complejo/IDS_Reales_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_real3))

# IDS_contradict junta todas las proteínas que tienen info contradictoria
with open(nombre_motivo+'_complejo/IDS_Contradict_ClinVar_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_contradict5))

# Luego le sumo los detectados mediante gnomad
IDS_falso3 += df_final_gnomad.ID[df_final_gnomad.Benign.isin(patogenicos)].unique().tolist()
IDS_falso3 = list(set(IDS_falso3))
with open(nombre_motivo+'_complejo/IDS_Falsos_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_falso3))

IDS_real3 += df_final_gnomad.ID[df_final_gnomad.Benign.isin(benignos)].unique().tolist()
IDS_contradict6 = [x for x in IDS_real3 if x in IDS_falso1+IDS_falso2+IDS_falso3]
IDS_real3 = [x for x in IDS_real3 if x not in IDS_falso1+IDS_falso2+IDS_falso3]
IDS_real3 = list(set(IDS_real3))
with open(nombre_motivo+'_complejo/IDS_Reales_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_real3))

# IDS_contradict junta todas las proteínas que tienen info contradictoria
with open(nombre_motivo+'_complejo/IDS_Contradict_iteracion'+str(iteracion)+'.txt', "w") as output:
    output.write(str(IDS_contradict5+IDS_contradict6))    

#Guardo aquellos que fueron predichos como reales pero que termine descartando a mano                

#IDS_descartados += 
#IDS_verificados_real = 
#IDS_descartados += IDS_verificados_real
#with open(nombre_motivo+'_complejo/IDS_Descartados_iteracion'+str(iteracion)+'.txt', "w") as output:
#    output.write(str(IDS_descartados))


#%%
iteracion = 'final'
iteracion_previa = 'flexible'
patogenicos, benignos = pato_ben(mot, nombre_motivo, iteracion, archivo='iteracionflexible')

















