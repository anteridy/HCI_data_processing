'''=====IMPORT LIBRARIES======================================='''
from matplotlib import  pylab as plt
from matplotlib import cm
import scipy.stats as stats
import numpy as np
from numpy import ma
import itertools
import seaborn as sns
import os
from scipy import spatial, linalg
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import pandas as pd
# import scipy.spatial as sp, scipy.cluster.hierarchy as hc

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cbook
from matplotlib import cm
# from sklearn.decomposition import PCA

import umap

#time stamp
from datetime import datetime

#matplotlib module 'matplotlib.cbook' has no attribute 'iterable',
# apparently this is due to anaconda installing few matlab packages (>1, one by pip and one by anaconda)
# to check run: 
# $ conda list matplotlib
# if more then 1 instance found, deinstall:
# $ pip uninstall matplotlib
# then reinstall needed version (e.g. 3.2.2) as:
# $ conda install matplotlib=2.1.1

#%%Functions
def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)

#%%
'''=====INPUTS & DEFINITIONS================================================'''
#PATHS------------------------------------------------------------------------
#Path to the raw data & results
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
Folder = '7_Filtered_data\\Sample_wrap_liveCD8\\Plt_'
ResFolder = '8_SignatureAnalysis\\UMAPs\\Plt_'
File_core = 'FinalResultFile_core.csv'

index_feature_start = 13

#Donors to remove
# unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']
unw_donors = [] #Already removed during FeatureReduction step

# plates = [7,8,9]
plates = [2]

for plate in plates:
    #data save to
    ResPath = DataPath+ResFolder+str(plate)+'\\'
    forcepath(ResPath)
    
    #UMAP Settings-------------------
    metric = 'euclidean'
    nb = 10 #neighbours
    mind = 0.1 #min distance
    ncomp = 2 #number of components
    rnds = 600 #random states
    
    #Save UMAP settings
    text_file = open(ResPath+"UMAPsettings_core.txt", "w")
    text_file.write('UMAP settings: \n')
    text_file.write('Plt_'+str(plate)+'\n'+'Metrics = '+metric+'\n'+
                    'Neigbours = '+str(nb)+'\n'+
                    'Minimal distance = '+str(mind)+'\n'+
                    'Number of components = '+str(ncomp)+'\n'+
                    'Random state/seed = '+str(rnds))
    text_file.close()
    #-------------------------------
    
    
    #%% Processing loop
    print('Plate='+str(plate))
    #data save to
    ResPath = DataPath+ResFolder+str(plate)+'\\'
    forcepath(ResPath)
    
    #%% Load data----------------------------------------
    print ("Loading data...")
    #Load core data-----------------------------------
    #open csv file with data for plate i as pd dataframe
    path1 = DataPath+Folder+str(plate)+'\\'+File_core
    data_raw = pd.read_csv(path1)
    data_raw_cols = list(data_raw)   
    #drop unwanted donors:
    don_ind = data_raw['Patient_ID'].isin(unw_donors)
    data = data_raw[~don_ind]        
    #get all features
    features = data.columns[index_feature_start:] 
    print ('Number of core features: %d' %len(features))
    #-----------------------------------------------
    # stop

    #%%#Extract all values
    x = data[features].values
    
    #%% UMAP - Core features  
    #get metadata
    metadata = data.iloc[:,0:index_feature_start]
    
    #get features
    x_features = data[features]
    
    #Check for Nan or Inf in the array
    Nans = np.where(np.isnan(x))
    nanN = len(Nans[0])
    Infs = np.where(np.isinf(x))
    infN = len(Infs[0])
    
    if nanN+infN > 0:
        print('')
        print('-----------------------------------------------')
        print('[!] Nan or Inf values in the data!')
        print('Nan number = '+str(nanN))
        print(Nans)
        print('Inf number = '+str(infN))
        print(Infs)
        print('Proceed with removing rows (FOVs) with inf or Nan')
        print('FOVs removed = '+str(nanN+infN))
        print('-----------------------------------------------')
        print('')
        
        # #save err log
        # path = ResPath
        # text_file = open(path+"ErrorLog.txt", "w")
        # text_file.write('[!] Nan or Inf values in the data! \n')         
        # text_file.write('Wrap method = '+wrapmeth+'/ Dataset = '+dataset+'/ CohenD cutoff = '+D_min+'\n')
        # text_file.write('Plate = '+str(plate)+'\n')
        # text_file.write('Nan number = '+str(nanN)+'\n')
        # text_file.write('Inf number = '+str(infN)+'\n')
        # text_file.write('Proceed with removing rows (FOVs) with inf or Nan'+'\n')
        # text_file.write('FOVs removed = '+str(nanN+infN)+'\n')
        # text_file.close()
        
        #drop rows with nan or Inf
        mask = np.any(np.isnan(x),axis=1) | np.any(np.isinf(x),axis=1)
        #amend data for plotting
        x = x[~mask]
        x_features = x_features[~mask]
        metadata = metadata[~mask]
        #amend colors for plotting
        colors_pca = np.array(colors)[~mask]
        
        # stop
    else:
        print('No Nans or Inf values found in data :)')
    
        
    # UMAP transformation
    #scale values 
    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(x)
    
    #UMAP
    fit = umap.UMAP(n_neighbors=nb,min_dist=mind, n_components=ncomp, metric=metric)
    x_umap = fit.fit_transform(x_scaled)
    
    #save UMAP-----------------------
    #merge umap data and metadata
    df_umap = metadata
    df_umap['UMAP1'] = x_umap[:,0]
    df_umap['UMAP2'] = x_umap[:,1]
    
    #save as csv
    path = ResPath+'UMAP_data_core.csv'
    df_umap.to_csv(path)
        
        
    
#%%END













