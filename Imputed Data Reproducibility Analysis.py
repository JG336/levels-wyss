# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:53:56 2020

@author: Jennifer Grant
"""

import pandas as pd
import os 
import numpy as np
from itertools import combinations


################# CELL DATASET FIRST! #################
data_cells= pd.ExcelFile(r'C:\Users\Wyss User\Documents\Wyss Documents\Gut Chip\Metabalon results\HARV-12-19VW CDT (CELLS).xlsx')
OrigScale_cells = data_cells.parse(1) #imputs the origscale tab of large excel doc
metabolite_list_cells= OrigScale_cells['Unnamed: 1'] #list of molecules in OrigScale tab


#for cell study, there are 2 experiments: 1) + stretch healthy, 2) - stretch healthy. generate AUC matrices for each
plus_stretch_healthy= pd.concat([metabolite_list_cells, OrigScale_cells['Chip 1 Cells +S'], OrigScale_cells['Chip 2 Cells +S'], OrigScale_cells['Chip 3 Cells +S']], axis=1)
minus_stretch_healthy= pd.concat([metabolite_list_cells, OrigScale_cells['Chip 4 Cells -S'], OrigScale_cells['Chip 5 Cells -S'], OrigScale_cells['Chip 6 Cells -S']], axis=1)

#import the tab with % filled value calculation
Per_filled_cells = data_cells.parse(4) #imports the pathway heatmap tab from excel
Per_filled_cells_plus_stretch= pd.concat([Per_filled_cells['Unnamed: 4'], Per_filled_cells['Unnamed: 20']], axis=1)
Per_filled_cells_minus_stretch= pd.concat([Per_filled_cells['Unnamed: 4'], Per_filled_cells['Unnamed: 19']], axis=1)

#generate list of metabolites for each dataset that equal 67 or 33 per match
cells_minus_33and67= Per_filled_cells_minus_stretch[(Per_filled_cells_minus_stretch['Unnamed: 19'] == 33) | (Per_filled_cells_minus_stretch['Unnamed: 19'] == 67)]
cells_plus_33and67= Per_filled_cells_plus_stretch[(Per_filled_cells_plus_stretch['Unnamed: 20'] == 33) | (Per_filled_cells_plus_stretch['Unnamed: 20'] == 67)]

#list only the AUCs of metabolites that have 33 or 67% percent occurence ratio
cells_plus_common_cols = plus_stretch_healthy[plus_stretch_healthy['Unnamed: 1'].isin(cells_plus_33and67['Unnamed: 4'])]
cells_minus_common_cols= minus_stretch_healthy[minus_stretch_healthy['Unnamed: 1'].isin(cells_minus_33and67['Unnamed: 4'])]
        


########### MEDIA DATASET NEXT! ###########
data_media= pd.ExcelFile(r'C:\Users\Wyss User\Documents\Wyss Documents\Gut Chip\Metabalon results\HARV-12-19VW CDT (MEDIA).xlsx')
osm = data_media.parse(1) #imputs the origscale tab of large excel doc
metabolite_list_media= osm['Unnamed: 1'] #list of molecules in OrigScale tab

#import the tab with % filled value calculation
pfm= data_media.parse(3) #imports the pathway heatmap tab from excel

#define AUC data for all experimental replicates under one variable
lm_e1= pd.concat([metabolite_list_media, osm['Loading Apical'], osm['Loading Apical.1'], osm['Loading Apical.2']], axis=1)
hbss_e1= pd.concat([metabolite_list_media, osm['HBSS'], osm['HBSS.1'], osm['HBSS.2']], axis=1)
ao_e1plus= pd.concat([metabolite_list_media, osm['A1 LC/MS'], osm['A3 LC/MS'], osm['A2 LC/MS']], axis=1)
bo_e1plus= pd.concat([metabolite_list_media, osm['B2 LC/MS'], osm['B1 LC/MS'], osm['B3 LC/MS']], axis=1)
ao_e1minus= pd.concat([metabolite_list_media, osm['A4 LC/MS'], osm['A6 LC/MS'], osm['A5 LC/MS']], axis=1)
bo_e1minus= pd.concat([metabolite_list_media, osm['B5 LC/MS'], osm['B4 LC/MS'], osm['B6 LC/MS']], axis=1)
lm_e2= pd.concat([metabolite_list_media, osm['Loading media 2'], osm['Loading media 2.2'], osm['Loading media 2.1']], axis=1)
hbss_e2= pd.concat([metabolite_list_media, osm['HBSS - 2.1'], osm['HBSS - 2.2'], osm['HBSS - 2']], axis=1)
ai_e2plus= pd.concat([metabolite_list_media, osm['A1 in - 2'], osm['A2 in - 2'], osm['A3 in - 2']], axis=1)
bi_e2plus= pd.concat([metabolite_list_media, osm['B1 in - 2'], osm['B2 in - 2'], osm['B3 in - 2']], axis=1)
ai_e2minus= pd.concat([metabolite_list_media, osm['A4 in - 2'], osm['A5 in - 2'], osm['A6 in - 2']], axis=1)
bi_e2minus= pd.concat([metabolite_list_media, osm['B4 in - 2'], osm['B5 in - 2'], osm['B6 in - 2']], axis=1)
ao_e2plus= pd.concat([metabolite_list_media, osm['A1 Out - 2'], osm['A2 Out - 2'], osm['A3 Out - 2']], axis=1)
bo_e2plus= pd.concat([metabolite_list_media, osm['B1 Out - 2'], osm['B2 Out - 2'], osm['B3 Out - 2']], axis=1)
ao_e2minus= pd.concat([metabolite_list_media, osm['A4 Out - 2'], osm['A5 Out - 2'], osm['A6 Out - 2']], axis=1)
bo_e2minus= pd.concat([metabolite_list_media, osm['B4 Out - 2'], osm['B5 Out - 2'], osm['B6 Out - 2']], axis=1)
# generate a dataframe of dataframes for all auc values
df_auc_media= [lm_e1, hbss_e1, ao_e1plus, bo_e1plus, ao_e1minus, bo_e1minus, lm_e2, hbss_e2, ai_e2plus, bi_e2plus, ai_e2minus, bi_e2minus, ao_e2plus, bo_e2plus, ao_e2minus, bo_e2minus]

#list of molecules from the percent pfm tab
metabolite_list_pf= pfm['Unnamed: 4']
#define percent filled data for each experiment 
pf_lm_e1= pd.concat([metabolite_list_pf, pfm['Unnamed: 114']], axis=1)
pf_hbss_e1=pd.concat([metabolite_list_pf, pfm['Unnamed: 115']], axis=1)
pf_ao_e1plus= pd.concat([metabolite_list_pf, pfm['Unnamed: 116']], axis=1) 
pf_bo_e1plus= pd.concat([metabolite_list_pf, pfm['Unnamed: 117']], axis=1) 
pf_ao_e1minus= pd.concat([metabolite_list_pf, pfm['Unnamed: 118']], axis=1) 
pf_bo_e1minus= pd.concat([metabolite_list_pf, pfm['Unnamed: 119']], axis=1)
pf_lm_e2= pd.concat([metabolite_list_pf, pfm['Unnamed: 120']], axis=1)
pf_hbss_e2= pd.concat([metabolite_list_pf, pfm['Unnamed: 121']], axis=1)
pf_ai_e2plus= pd.concat([metabolite_list_pf, pfm['Unnamed: 122']], axis=1)
pf_bi_e2plus= pd.concat([metabolite_list_pf, pfm['Unnamed: 123']], axis=1)
pf_ai_e2minus= pd.concat([metabolite_list_pf, pfm['Unnamed: 124']], axis=1)
pf_bi_e2minus= pd.concat([metabolite_list_pf, pfm['Unnamed: 125']], axis=1)
pf_ao_e2plus= pd.concat([metabolite_list_pf, pfm['Unnamed: 126']], axis=1)
pf_bo_e2plus= pd.concat([metabolite_list_pf, pfm['Unnamed: 127']], axis=1)
pf_ao_e2minus= pd.concat([metabolite_list_pf, pfm['Unnamed: 128']], axis=1)
pf_bo_e2minus= pd.concat([metabolite_list_pf, pfm['Unnamed: 129']], axis=1)


#generate a dataframe of dataframes for all percent filled data
df_percent_media= [pf_lm_e1, pf_hbss_e1, pf_ao_e1plus, pf_bo_e1plus, pf_ao_e1minus, pf_bo_e1minus, pf_lm_e2, pf_hbss_e2, pf_ai_e2plus, pf_bi_e2plus, pf_ai_e2minus, pf_bi_e2minus, pf_ao_e2plus, pf_bo_e2plus, pf_ao_e2minus, pf_bo_e2minus]


#list of metablites with 33 and 67 percent match, called df_percent_media_minus_33and67
for n, df in enumerate(df_percent_media):
    df_percent_media[n].columns= ['Molecule', 'Percent Filled']
df_percent_media_with_33and67= {}
for n, df in enumerate(df_percent_media):
    df_percent_media_with_33and67[n]= df_percent_media[n][(df_percent_media[n]['Percent Filled'] == 33) | (df_percent_media[n]['Percent Filled'] == 67)]


#list only the AUCs of metabolites that have 33 or 67% percent occurence ratio
media_common_cols= {}
for n, df in enumerate(df_auc_media):
    media_common_cols[n] = df_auc_media[n][df_auc_media[n]['Unnamed: 1'].isin(df_percent_media_with_33and67[n]['Molecule'])]
