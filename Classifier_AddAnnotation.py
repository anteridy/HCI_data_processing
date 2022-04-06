# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:13:03 2021

@author: akamnev
"""
#%%Import 

#OS functions
import os #working with files and folders

#Dataframes and numberical operations
import numpy as np
import pickle #Create portable serialized representations of Python objects
import pandas as pd #pandas dataframe toolkit for large datasets

#%% Functions

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
        
#%% Inputs

DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
DataFolder = '6_ImgClassifier\\CellImages\\'

samples = range(1,11)


#%% Load data

for sample in samples:

    #Load file with morpho data
    path = DataPath+DataFolder+'Data\\Sample_'+str(sample)+'.csv'
    data = pd.read_csv(path)
    data_attr = list(data)
    
    # stop
    #sort by file name to match annotation:
    data_srt = data.sort_values(by = ['Image_ID'], axis=0,ascending=True)
    
    #load file with annotation
    path = DataPath+DataFolder+'Images\\Result\\Result_'+str(sample)+'.txt'
    annotation = pd.read_csv(path, header=None, delimiter = "/")
    annotation.rename(columns={0:'CellN',1:'Trash'},inplace=True)
    
    #reset indexes
    data_srt.reset_index(drop=True,inplace=True)
    annotation.reset_index(drop=True,inplace=True)
    
    #%%
    #insert annotation
    data_ann = pd.concat([annotation,data_srt],axis=1)
    
    #Add set number
    data_ann['TrainingSetN'] = sample
    
    #%%
    #save annotation
    path = DataPath+DataFolder+'Data_ann\\'
    forcepath(path)
    data_ann.to_csv(path+'Sample_'+str(sample)+'_ann.csv', index = False)









#%% END
















