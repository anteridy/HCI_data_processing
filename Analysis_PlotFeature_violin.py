# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:23:51 2019

@author: akamnev

Aim: to view QC params as a 384 w/pl layout heatmap
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
Folder_res = '8_SignatureAnalysis\\Features\\Plt_'

#Wanted features
#path to core features:
path_core = '4_Metadata\\CoreFeatures.csv'
path_cd1 = '7_Filtered_data\\Sample_wrap_liveCD8\\Plt_'

#Which plate to process
plate = 2

# feature = 'CellFoot_Intensity_MeanIntensity_Actin'
# limit_y = False
# min_max = [0.01, 0.03]

feature = 'Nuclei_RadialDistribution_RadialCV_Marker_4of6'
limit_y = False
min_max = [0.01, 0.03]


#----------------------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 22 #last column of metadata in the table to split df
#------------------------------------------------------------

#Data to remove
unw_pids = ['DOCK11','PAK2','THEMIS2']
#Donors to remove
unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']

#%% Load data-------------------------------------

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
data_raw_cln = df[~pid_ind]
df = data_raw_cln
don_ind = df['Patient_ID'].isin(unw_donors)
data_1 = df[~don_ind]  
#-------------------------------------------

# stop
 #%%Check for Nan or Inf in the array
metadata = data_1[data_cols[0:lastMetCol+1]]
morphodata = data_1[[feature]]

df_tmp = morphodata.replace([np.inf, -np.inf], np.nan)
df_tmp.dropna(axis=1, how='any',inplace=True)

 
data_fin = pd.concat([metadata,df_tmp],axis=1)
#Get donor IDs
genotypes = np.sort(data_fin['Genotype'].unique())
# genotypes = ['ND','ARPC1B','HEM1','WAS','PIK3CG']

# stop
param = feature
unit = 'norm'

#%% VIOLIN PLOTS-------------------------------

#cycle through wanted attributes
plt.figure(0)
#clear figure
plt.clf()
#make violing plot for the attribute
ax = sns.violinplot(x = 'Genotype', y = param, data = data_fin)
ax.set_title('Plt_'+str(plate)+'_'+param, fontsize = 12)
plt.xlabel('Genotypes')

if limit_y:
    plt.ylim(min_max)

#save plot
# pltpath = DataPath+Folder_res+str(plate)+'\\Feature_violins\\'
# forcepath(pltpath)
# plt.savefig(pltpath+param+'.png')
# plt.close()


#%%      
#print('===DONE===============================')
    #%%
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    