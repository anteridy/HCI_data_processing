# -*- coding: utf-8 -*-
"""
Created on Wed May 12 09:56:09 2021

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
Folder_res = '8_SignatureAnalysis\\Features\\Plt_'

#Wanted features
#path to core features:
path_core = '4_Metadata\\CoreFeatures.csv'
path_cd1 = '7_Filtered_data\\Sample_wrap_liveCD8\\Plt_'

#Which plate to process
plate = 11
#----------------------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 22 #last column of metadata in the table to split df
#------------------------------------------------------------

#Data to remove
unw_pids = ['DOCK11','PAK2','THEMIS2', 'PSTPIP1', 'WDR1','MSN' ]
#Donors to remove
unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']

#%% Load data-------------------------------------

#Features---------------
#load core features
path = DataPath + path_core
df = pd.read_csv(path)
features = df['Feature'].tolist()   

#load CD1 features (filtered)
# path = DataPath+path_cd1+str(plate)+'//NoCorrelation_CD1.csv'
# df = pd.read_csv(path, header=None)
# features = df[0].tolist()

# features = []

#-------------------------

#status
print('Loading data...')
#Load pkl file with filtered morpho data
path = DataPath+Folder+ File_data+str(plate)+'.pkl'
with open(path, 'rb') as f:
        data = pickle.load(f)
data_cols = list(data)

#Clean data---------------------------------
#Drop unwanted unwanted PIDs
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
morphodata = data_1[features]

df_tmp = morphodata.replace([np.inf, -np.inf], np.nan)
df_tmp.dropna(axis=1, how='any',inplace=True)

 
data_fin = pd.concat([metadata,df_tmp],axis=1)
#Get donor IDs
# genotypes = np.sort(data_fin['Genotype'].unique())
# genotypes = ['ND','ARPC1B','HEM1','WAS','PIK3CG']


#adjust parameters
param = 'CellFoot_AreaShape_Area'
data_fin[param] = data_fin[param].apply(lambda x: x*0.106) # pxl^2 to um^2
param = 'CellFoot_Intensity_Displacement_Actin'
data_fin[param] = data_fin[param].apply(lambda x: x*0.325) # pxl to um
# %% STATISTICS
pids = ['ARPC1B','HEM1','WAS','PIK3CG']
# pids = ['WASP']
p_vals = []
for pid in pids:    
    #make df for p values
    data = {'Pair':['WT-'+pid]} 
    pval_df = pd.DataFrame(data)
    
    for feature in features:
        #get arrays for analysis
        array1 = data_fin[feature].loc[data_fin['Genotype']=='ND']
        array2 = data_fin[feature].loc[data_fin['Genotype']==pid]
        
        #get p values
        pval_wt_ko = scipy.stats.mannwhitneyu(array1, array2, use_continuity=True, alternative=None)
    
        #add to p_val_df
        pval_df[feature] = [pval_wt_ko[1]]
    p_vals.append(pval_df)
pval_df_fin = pd.concat(p_vals, axis=0)     
# stop
#save result
path = DataPath+Folder_res+str(plate)+'\\Feature_violin_core\\'    
pval_df_fin.to_csv(path+'MWUtest.csv')
    
#%% END













