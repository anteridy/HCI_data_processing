# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:47:41 2021

@author: akamnev

How to build ROC curves with Sklearn in python:
    https://stackabuse.com/understanding-roc-curves-with-python/
    
How to save learned model 
https://machinelearningmastery.com/save-load-machine-learning-models-python-scikit-learn/

"""
#%%Import 

#OS functions
import os #working with files and folders

#Dataframes and numberical operations
import numpy as np
import pickle #Create portable serialized representations of Python objects
import pandas as pd #pandas dataframe toolkit for large datasets

#machine learning
import sklearn

# roc curve and auc score
from sklearn.datasets import make_classification
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

#plots
import matplotlib.pyplot as plt #main plotting engine
import seaborn as sns           # makes plots look nicer

#%% Functions

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
        
#%% Inputs

#PATHS-----------------------------------------------
#Data path; [!] level must be indicateds as "\\"
#Input data
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
DataFolder = '3_Raw_data_mod\\' 
DataFile = 'SC_raw_plt'
#Where to write results
ResFolder = '6_PreProcessed_data\\ImgClassifier\\' 
ResName_all = 'SC_all' #all cells with indexed
ResName_live = 'SC_live'
ResName_livecd8 = 'SC_liveCD8'
ResPath = DataPath+ResFolder
forcepath(ResPath)
#------------------------------------------------------------

#plates to process
# plates = list(range(1,13))
# plates = list(range(10,13))
# plates = [2,5,8,11]
plates=[1,3,4,6,7,9,10,12]

# Trained model to find trash
ModelFolder = '6_ImgClassifier\\Model\\'
ModelFile = 'TrainedModel.sav'

#Metatdata for features-------------------------
FeatFolder = '4_Metadata\\'
FeatName_cd8 = 'CD8_features.csv'
FeatName_unwanted = 'UnwantedColumns.csv'
FeatName_unwanted_meta = 'UnwantedMetadata.csv'
FeatName_model = 'CoreFeatures.csv'
#------------------------------------------

#Gates for CD8 ------------------------------------------------------
gate_attribute_cd8 = 'CellFoot_Intensity_MeanIntensity_CD8'
gate_val_cd8 = 0.0055
#------------------------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 10 #last column of metadata in the table to split df
#------------------------------------------------------------

#Quick test?
test = False #Open only first x rows
# test = True #Open only first x rows
rown = 1000

#%%Load data
#Load csv with features to process
foldpath = DataPath+FeatFolder
path = foldpath+FeatName_cd8
cd8_cols = pd.read_csv(path, index_col=None, header=0)['# CD8_features'].astype(str)
path = foldpath+FeatName_unwanted
unwanted_cols = pd.read_csv(path, index_col=None, header=0)['Unwanted_columns'].astype(str)
path = foldpath+FeatName_unwanted_meta
unwanted_cols_meta = pd.read_csv(path, index_col=None, header=0)['Unwanted_metadata'].astype(str)
path = foldpath+FeatName_model
model_cols = pd.read_csv(path, index_col=None, header=0)['Feature'].astype(str)
# unwanted_cols.replace('\n',' ', regex=True) 

for plate in plates:    
    #%% Load pkl file with morpho data
    print('Plate = '+str(plate))
    print('loading data...')
    path = DataPath+DataFolder+DataFile+str(plate).zfill(2)+'.pkl'
    
    with open(path, 'rb') as f:
            data_raw = pickle.load(f)
    
    data_raw_columns = list(data_raw)
    # stop
    
    #Reduce data size for testing
    if test:
        data = data_raw.iloc[:rown] #[!]for testing only
    else:
        data = data_raw
    # stop
    #------------------------------------------------
    
    #%%
    #Find trash========================================
    #status
    print('Locating trash...')  
    
    # load the model from disk
    path = DataPath + ModelFolder + ModelFile
    model = pickle.load(open(path, 'rb'))
    
    #Get core features for the model (trained on core)
    X = data[model_cols]
    
    #Check for Nan or Inf in the X
    Nans = np.where(np.isnan(X))
    nanN = len(Nans[0])
    Infs = np.where(np.isinf(X))
    infN = len(Infs[0])
    if nanN+infN > 0:
        print('')
        print('-----------------------------------------------')
        print('[!] Nan or Inf values in the data!')
        print('Nan number = '+str(nanN))
        print(Nans)
        print('Inf number = '+str(infN))
        print(Infs)
        print('remove all objects with Nan or Inf')
        print('-----------------------------------------------')
        print('')  
        #drop rows with nan or Inf
        mask = np.any(np.isnan(X),axis=1) | np.any(np.isinf(X),axis=1)
        #amend data for plotting
        X = X[~mask]    
        data = data[~mask]
    else:
        print('No Nans or Inf values found in data :)')
    
    #Make prediction
    predicted = model.predict(X)
    
    #Annotate
    data.insert(loc=0, column='Trash', value=predicted)
    #======================================================================
    
    #%%Find CD8 high
    cd8_ind = data[gate_attribute_cd8]>=gate_val_cd8
    #Annotate
    data.insert(loc=0, column='CD8_high', value=cd8_ind)
    
    #%%
    #Remove unwanted columns===============================================
    #Get wanted metadata
    df_meta = data.iloc[:,0:lastMetCol+3]   
    df_tmp_meta = pd.DataFrame()
    df_tmp_meta = df_meta.drop(unwanted_cols_meta,axis=1) #cd8-associated features
    
    # stop
    
    #remove unwanted features
    df_attr = data.iloc[:,lastMetCol+3:] 
    df_tmp_attr = pd.DataFrame()
    df_tmp_attr = df_attr.drop(cd8_cols,axis=1) #cd8-associated features
    df_tmp_attr.drop(unwanted_cols,axis=1,inplace=True) #known trash column (mainly location)
    df_tmp_attr.dropna(axis=1,how='all',inplace=True) #columns with all NaNs inside
    #drop cols with same element
    print('Checking for columns filled with same value...')
    for col in df_tmp_attr.columns:
        if len(df_tmp_attr[col].unique()) == 1:
            print('Out-> '+col)
            df_tmp_attr.drop(col,inplace=True,axis=1)
    
    #merge back to df
    df_fin = pd.concat([df_tmp_meta, df_tmp_attr], axis=1, ignore_index=False)
    #======================================================================
    
    #%%
    # SAVE RESULTS
    #All cells
    print('Saving final pkls...')
    print('All objects')
    path = ResPath+ResName_all+'\\'
    forcepath(path)
    df_fin.to_pickle(path+ResName_all+'_plt'+str(plate)+'.pkl') 
    
    #df w/o trash
    print('Live cells')
    df_live = df_fin.loc[df_fin['Trash']==False] 
    path = ResPath+ResName_live+'\\'
    forcepath(path)
    df_live.to_pickle(path+ResName_live+'_plt'+str(plate)+'.pkl')   
    
    #df live CD8 high
    print('Live CD8 cells')
    df_live_cd8 = df_live.loc[df_live['CD8_high']==True] 
    path = ResPath+ResName_livecd8+'\\'
    forcepath(path)
    df_live_cd8.to_pickle(path+ResName_livecd8+'_plt'+str(plate)+'.pkl') 

#%% END











