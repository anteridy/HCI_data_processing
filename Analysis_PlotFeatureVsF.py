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
from sklearn import preprocessing #for scaling dfs
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
Folder = '6_PreProcessed_data\ImgClassifier\Sample_wrap_liveCD8\\'
Folder_res = '8_SignatureAnalysis\\Features\\Plt_'
#Raw morpho data:
File_data = 'AveSample_S20_plt_'

#Wanted features
#path to core features:
path_core = '4_Metadata\\CoreFeatures.csv'

#Which plate to process
plate = 11
#----------------------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 22 #last column of metadata in the table to split df
#------------------------------------------------------------

#Data to remove
unw_pids = ['DOCK11','PAK2','THEMIS2']
#Donors to remove
unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']

#%% Load data-------------------------------------

#load pinned features
path = DataPath + path_core
df = pd.read_csv(path)
features = df['Feature'].tolist()    

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
# morphodata = data_1[features]
morphodata = data_1[data_cols[lastMetCol+1:]]

df_tmp = morphodata.replace([np.inf, -np.inf], np.nan)
df_tmp.dropna(axis=1, how='any',inplace=True)

 
data_fin = pd.concat([metadata,df_tmp],axis=1)
# data_fin = data_1
#Get donor IDs
genotypes = np.sort(data_fin['Genotype'].unique())

# stop
#%% Contour plot

#Parameter to plot
# param1 = 'CellFoot_AreaShape_Perimeter'
# data_fin[param1] = data_fin[param1].apply(lambda x: x*0.106) # pxl^2 to um^2
# unit1 = 'um^2'
# param1 = 'CellFoot_AreaShape_Perimeter'
# data_fin[param1] = data_fin[param1].apply(lambda x: x*0.325) # pxl to um
# unit1 = 'um'

#Plt2  - ICAM1
# param2 = 'CellFoot_AreaShape_EOP'
# unit2 = ''
# param1 = 'CellFoot_Intensity_MeanIntensity_Actin'
# unit1 = 'AU'
# limit_axes = True
# min_max_x = [0.010,0.030]
# min_max_y = [0.22,0.45]

#Plt8 - IS
# param2 = 'CellFoot_AreaShape_EOP'
# unit2 = ''
# # param1 = 'CellFoot_AreaShape_Eccentricity'
# # unit1 = ''
# param1 = 'CellFoot_Intensity_MeanIntensity_Actin'
# unit1 = 'AU'
# limit_axes = True
# min_max_x = [0.009,0.027]
# min_max_y = [0.10,0.42]

#Plt11 - IS+IL2
param2 = 'CellFoot_AreaShape_EOP'
unit2 = ''
# param1 = 'CellFoot_AreaShape_Eccentricity'
# unit1 = ''
param1 = 'CellFoot_Intensity_MeanIntensity_Actin'
unit1 = 'AU'
limit_axes = False
min_max_x = [0.009,0.027]
min_max_y = [0.10,0.42]


genotypes2 = ['ND', 'WAS', 'HEM1', 'ARPC1B', 'PIK3CG']
# genotypes2 = ['ND', 'WAS']
luts = ['Greys','Greens','Oranges','Reds', 'Blues']
ind = 0
for genotype in genotypes2:
    #make df for donor group i
    df_tmp = data_fin.loc[data_fin['Genotype']==genotype]
    #Contour plot
    sns.kdeplot(df_tmp[param1], df_tmp[param2],
                cmap=luts[ind],
                label = genotype,
                alpha = 0.7)
    ind = ind+1
    # stop
    
#Plot title
plt.title('Plt_'+str(plate))
#labels
plt.ylabel(param2+', '+unit2)
plt.xlabel(param1+', '+unit1)
#scale
#    plt.yscale('log')
# plt.xscale('log')
if limit_axes:
    plt.xlim(min_max_x)
    plt.ylim(min_max_y)
#legend
# plt.legend(genotypes2)

#Save figure
# pltpath = DataPath+Folder_res+str(plate)+'\\'
# forcepath(pltpath)
# plt.savefig(pltpath+'ContourPlot.png')
# plt.close()
# print('DONE')
#-----------------------------------------------------
# stop

#%%      
#print('===DONE===============================')
#%% END
# plt.scatter(df_tmp[param1], df_tmp[param2])
    
    
    
    
    
    
    
    
    
    
    
    
