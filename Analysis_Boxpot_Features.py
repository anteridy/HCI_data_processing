# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:23:51 2019

@author: akamnev
"""
#%%
#OS functions
import os #working with files and folders

#Dataframes and numberical operations
import numpy as np
# from sklearn import preprocessing #for scaling dfs
import pickle #Create portable serialized representations of Python objects
import pandas as pd #pandas dataframe toolkit for large datasets

#plots
import matplotlib.pyplot as plt #main plotting engine
import seaborn as sns           # makes plots look nicer

#statistics
import scipy

#%%
"""==============================================
--------FUNCTIONS-------------------------
================================================"""
def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
"""=====END OF FUNCTIONS========================"""
#%%
"""==============================================
--------INPUT - Raw dataset(s)-------------------------
================================================"""
#PATHS-----------------------------------------------
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'

#Raw morpho data:
Folder = '6_PreProcessed_data\ImgClassifier\Sample_wrap_liveCD8\\'
File_data = 'AveSample_S20_plt_'

#results to
Folder_res = '8_SignatureAnalysis\\Features\\Boxplots\\'
# ResPath = DataPath+Folder_res
# forcepath(ResPath)

#Which plate to process
plate = 8
# plates = [2,5]
# plates = [8,11]
# plot_name = 'Plate'+str(plates[0])+'_vs_'+str(plates[1])

ResPath = DataPath+Folder_res+'Plate_'+str(plate)+'\\'
forcepath(ResPath)

features = ['CellFoot_AreaShape_Area',
            'CellFoot_AreaShape_EOP',
            'CellFoot_AreaShape_Eccentricity',
            'CellFoot_Intensity_MeanIntensity_Actin',
            'CellFoot_Intensity_IntegratedIntensity_Actin'
            ]

#----------------------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 22 #last column of metadata in the table to split df
#------------------------------------------------------------

#Data to remove
unw_pids = ['DOCK11','PAK2','THEMIS2']
#Donors to remove
unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']

wanted_genotypes = ['ND','WAS','HEM1','ARPC1B', 'PIK3CG']

patient_index = [['PID_WAS1','WAS_Pt1'],
                 ['PID_WAS2','WAS_Pt2'],
                 ['PID_WAST2','WAS_Pt3'],
                 ['PID_2466','HEM1_Pt1'],
                 ['PID_2378','HEM1_Pt2'],
                 ['PID_ARPC1B_Pt1','ARPC1B_Pt1'],
                 ['PID_2599 (CD8)','ARPC1B_Pt2'],
                 ['PID_ARPC1B_Pt7','ARPC1B_Pt3'],
                 ['PID_2785','PIK3g']]

#%% Load data from plates -------------------------------------  
#status
print('Loading data...')
#Load pkl file with filtered morpho data
path = DataPath+Folder+ File_data+str(plate)+'.pkl'
with open(path, 'rb') as f:
        data = pickle.load(f)
data_cols = list(data)

#Clean data---------------------------------
#Drop unwanted unwanted PIDs & donors
df = data
pid_ind = df['Genotype'].isin(unw_pids)
df1 = df[~pid_ind]
don_ind = df1['Patient_ID'].isin(unw_donors)
df2 = df1[~don_ind] 
don_ind2 = df2['Genotype'].isin(wanted_genotypes)
df3 = df2[don_ind2]
#-------------------------------------------

#reduce to wanted features
metadata = df3[data_cols[0:lastMetCol+1]]
morphodata =df3[features]

#drop unwanted metadata
metadata_cln = metadata.drop(['Used_wells', 'CD8_high', 'Trash',
                              'ImageNumber', 'ObjectNumber',
                              'Metadata_FOV', 'ObjectCenterX',
                              'ObjectCenterY'], axis=1)

#merge datasets
data_cln = pd.concat([metadata_cln,morphodata],axis=1)

#sort data by pt ID (for plotting consistency)
data_cln_srt = data_cln.sort_values(by=['Patient_ID'])

#%% Rebuild index for figure

df_bucket = []
df_tmp = data_cln_srt

#make ND dataframe
#add nd index
mask = df_tmp['Genotype']=='ND'
df_group = df_tmp[mask]
df_group.insert(0, 'Paper_index', 'ND')
# df_fin = df_group.reset_index(drop=True)
#append result
df_bucket.append(df_group)

#index patients
for patient in patient_index:
    #add nd index
    mask = df_tmp['Patient_ID']==patient[0]
    df_group = df_tmp[mask]
    df_group.insert(0, 'Paper_index', patient[1])
    # df_fin = df_group.reset_index(drop=True)
    #append result
    df_bucket.append(df_group)
    #append result


#merge datasets
df_reindexed = pd.concat(df_bucket, axis=0).reset_index(drop=True)

#%% Plot boxplot

# Create an array with the colors you want to use
colors = ['#525252', #ND
          '#04B404', #WAS, green
          '#04B404', #WAS, green
          '#04B404', #WAS, green
          '#FFBF00', #HEM1, orange
          '#FFBF00', #HEM1, orange
          '#FF0000', #ARPC1B, red
          '#FF0000', #ARPC1B, red
          '#FF0000', #ARPC1B, red
          '#013ADF', #PIK3g, blue
    ]
# Set your custom color palette
customPalette = sns.set_palette(sns.color_palette(colors))

for feature in features:
    #ind plots side by side
    # b = sns.factorplot(data = data_fin_pooled,
    #                 x = 'IL-2',
    #                 y = feature,
    #                 kind = 'box',
    #                 col = 'Patient_ID',
    #                 #col_order = ['',''] #custom column order
    #                 ).set_titles('{col_name}')# remove 'column = ' part of title
    
    #all in one plot
    df2plot = df_reindexed
    b = sns.boxplot(data = df2plot,
                    # hue = 'IL-2',
                    x = 'Paper_index',
                    y = feature,
                    palette = customPalette
                    #order = [] #custom column order
                    )
    # sns.plt.title(feature) # You can change the title her
    b.set_xticklabels(b.get_xticklabels(),rotation=45)
    b.set_title(feature)
    # stop
    b.set_ylabel(feature)
    #save plot
    plt.savefig(ResPath+feature+'.png')
    plt.close()
    # stop




#%%      
#print('===DONE===============================')
    #%%
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
