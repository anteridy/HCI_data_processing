# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:24:59 2019

@author: akamnev
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:58:31 2019

@author: akamnev
"""

"""==============================================
#------------INTRO--------------------------------
==================================================
Script to add annotation to the .pkl morpho data

#-------------------------------------------------
==================================================
"""

#%%
"""==============================================
--------IMPORT PACKAGES-------------------------
================================================"""
import os #working with files and folders
import pandas as pd #pca toolkit
import numpy as np
import pickle #Create portable serialized representations of Python objects
#import math #for square root calc etc, use numpy instead for lists and panda datasets

"""-------------------------------------------------
==================================================
"""

#%%
"""--------INPUTS--------------------"""
#Raw data
RawDataName = 'Results_raw_pooled_mod_plt'
RawDataPath = 'D:\\Temp\\APID\\APID_5_Screen1\\3_Python_res\\Raw_pooled_mod\\'
ResPath = 'D:\\Temp\\APID\\APID_5_Screen1\\3_Python_res\\Raw_pooled_mod_ann\\'
ResName = 'Results_raw_pooled_mod_ann_plt'

CSVMetadataPath = 'D:\\Temp\APID\\APID_5_Screen1\\3_Python_res\\Metadata\\Metadata.csv'

#Plates to process
#plates = ['01','02','03','04','05','06','07','08','09','10','11','12']
plates = ['02','05','08','11']
#"""-----------------------------------------------------------------------------"""

#%% Process plates
for plate in plates:
    #status
    print('')
    print('===================================')
    print('Processing Plate '+plate+'...')
    print('-----------------------------------')

    # Load data
    #status
    print('Loading data...')
    
    #Load csv with plate legend
    plate_legend = pd.read_csv(CSVMetadataPath, index_col=None, header=0)
    
    #Load pkl file with morpho data
    path = RawDataPath + RawDataName + plate +'.pkl'
    with open(path, 'rb') as f:
            df_morpho_raw = pickle.load(f)
    
    #Load csv file with averaged morpho data
#    path = RawDataPath + RawDataName + plate +'.csv'
#    df_morpho_raw = pd.read_csv(path, index_col=None, header=0)


    # Iterate through each entry and add corresponding metadata
    #status
    print('Annotating...')
    
    legends = []
    for index in range(0,df_morpho_raw.shape[0]):
        #get column and row values for the entry
        col = df_morpho_raw['Metadata_Column'][index]
        row = df_morpho_raw['Metadata_Row'][index]
        
        #find corresponding row in the plate legend
        legend = plate_legend.loc[(plate_legend['Column']==col) & (plate_legend['Row']==row)]
        
        #catenate with previously extracted rows
        legends.append(legend)
    
    # merge dfs to single lenged df
    df_legend = pd.concat(legends, axis=0, ignore_index=True)
    
    # merge with raw morpho data
    df_morpho_annot = pd.concat([df_legend.iloc[:,4:],df_morpho_raw], axis=1, ignore_index=False)

    # SAVE RESULTS
    #status
    print('Saving results...')
    
    #save as .pkl
    #df_morpho_annot.to_pickle(OutPath)
    
    #save as .csv
    path = ResPath+ResName+plate+'.pkl'
    df_morpho_annot.to_pickle(path)
    
    #status
    print('Plate '+plate+' is done.')
    
#status
print('==========================================================')
print('---> All plates are finished.')

#%%"""=====END OF CALCULATIONS========================="""















