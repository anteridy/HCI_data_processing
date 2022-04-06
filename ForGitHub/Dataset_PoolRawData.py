# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:58:31 2019

@author: akamnev
"""

"""==============================================
#------------INTRO--------------------------------
==================================================
Script fetches result .csv files produced by CellProfiler
jobs on cluster in iterative mode. For each batch of fiels CP made
individual Result_CellsFromEdge.csv table in Results/batch_n for each iteration.
These files need to be fetched, merged to single dataset (panda) and 
saved for later use (pickle format).
[19.06.24] - Script only pools data now, extra processing and annotation done by
     other scripts
[19.06.26] - Script now can handle empty batch folders (happens when CP runs out
            of images to process)
[19.08.22] - Takes now 2 result .csv tables (cell footprint + nuclei)

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

"""-------------------------------------------------
==================================================
"""

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
"""--------INPUTS--------------------"""

#Global variables:
ResTableName = ['Result_Cell_Footprint_Actin.csv',
                'Result_Filtered_Nuclei.csv']
#Path to the cluster results (plate folder):
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\2_ClusterRes\\'
#Where to write results
OutPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\3_Raw_data\\test\\'
#Name for single-cell raw dataset:
ResName = 'SC_raw_plt'

#plates to process
plates = ['02']
# plates = ['02','03','04','05','06','07','08','09','10','11','12']


#Other inputs
lastcol_metadata = 10 #to split morpho params from nuclei dataset

#%% Processing loop

for plate in plates:
    
    #status
    print('')
    print('Processing Plate '+plate+'...')
    print('===================================')
    
    #build path to results
    ResTablePath = DataPath + 'Plate'+plate+'\\Results\\'
    
    """--------MAKE PATHS TO RESUTLS--------------------"""
    #status
    print('Preparing path library...')
    
    #Make file paths to result tables
        
    #get file list from the directory directly:
    file_list = os.listdir(ResTablePath)
    
    #prep empty files
    tablePaths_actin =[]
    tablePaths_nuclei =[]
    
    #check if all files present--------------------------
    file_check = []
    batch_wo_data = []
    for batch in file_list:   
        #build path to the batch folder
        DirPath = ResTablePath +batch   
        #check if folder is empty (in case more batches then files in CP pipeline)
        if len(os.listdir(DirPath) ) == 0:
            continue
        else:
            path = DirPath+'\\'+ResTableName[0] #actin file
            file_check.append(os.path.exists(path))
            if os.path.exists(path):
                continue
            else:
                batch_wo_data.append(batch)
                
    #abort if data is missing
    if len(batch_wo_data) != 0:
        print("CSV files are missing in some batches")
        stop
    #----------------------------------------------------        
    #
    
    #cycle trough found folder (batches)
    for batch in file_list:   
        #build path to the batch folder
        DirPath = ResTablePath +batch   
        #check if folder is empty (in case more batches then files in CP pipeline)
        if len(os.listdir(DirPath) ) == 0:
            continue
        else:
            #build path to batch folder
            path2actin = DirPath+'\\'+ResTableName[0]
            path2nuclei = DirPath+'\\'+ResTableName[1]
            #add path to the table of paths
            tablePaths_actin.append(path2actin)
            tablePaths_nuclei.append(path2nuclei)
      
    #
    """--------CHECK RESULT TABLE--------------------"""
     
    # get column names from
    # actin results as headers of batch_1 result table No1
    fp = open(tablePaths_actin[0],'r')
    for i, line in enumerate(fp):
        if i == 0:
            headers=line
    fp.close()
    headers_list_actin = headers.split(',')
    
    # nuclei
    fp = open(tablePaths_nuclei[0],'r')
    for i, line in enumerate(fp):
        if i == 0:
            headers=line
    fp.close()
    headers_list_nuclei = headers.split(',')
    
    #Parameters to extract from result file
    #wanted = [10] + range(14,29)
    
    #
    """--------LOAD DATA TO PANDA DATAFRAME--------------------""" 
    #merge Actin & DNA datasets-----------------------
    #status
    print('Loading data...')
    
    #Merge Actin data
    print('Loading actin data...')
    #Load result .csvs & merge to sinle list "merged"
    merged = []
    index = 1
    for filename in tablePaths_actin:
        print('Pt'+plate+' Actin batch'+str(index))
        index = index + 1
        df = pd.read_csv(filename, index_col=None, header=0)
#        stop
        merged.append(df)
    #convert to pandas dataframe
    ResDf1 = pd.concat(merged, axis=0, ignore_index=True)
    #save
    #ResDf1.to_pickle(OutPath+'\\Actin_'+DataFrameName)
    #status
    print('Actin set is ready')
    print('--------------------')
    
    #Merge Nuclei data
    print('Loading nuclei data...')    
    merged = []
    index = 1
    for filename in tablePaths_nuclei:
        print('Pt'+plate+' Nuclei batch'+str(index))
        index = index + 1
        df = pd.read_csv(filename, index_col=None, header=0)
        merged.append(df)
    #convert to pandas dataframe
    ResDf2 = pd.concat(merged, axis=0, ignore_index=True)
    
    #Get attributes in the dataframe
    ResDf2_Attributes = list(ResDf2)
    #save
    #ResDf2.to_pickle(OutPath+'\\Nuclei_'+DataFrameName)
    #status
    print('Nuclei set is ready')
    print('--------------------')        
    
    #status
    print('Merging tables...')   
    # MERGE OBJECT TABLES
    metadata = ResDf1.iloc[:, :lastcol_metadata+1]
    data1 = ResDf1.iloc[:,lastcol_metadata+1:]
    data2 = ResDf2.iloc[:,lastcol_metadata+1:]
#    stop
    
    obj_name1 = 'CellFoot'
    obj_name2 = 'Nuclei'
    
    #Add prefix to column names of features
    attributes = list(data1)
    df = data1
    for attribute in attributes:
        df = df.rename(columns={attribute: obj_name1+'_'+attribute})
    data1_renamed = df
        
    attributes = list(data2)
    df = data2
    for attribute in attributes:
        df = df.rename(columns={attribute: obj_name2+'_'+attribute})
    data2_renamed = df
    
    #merge to one Df
    ResDf_fin = pd.concat([metadata, data1_renamed, data2_renamed], axis=1, ignore_index=False)
    
    #status
    print('Tables are merged')
    #-------------------------------------------------------------------
    
    
    # SAVE RESULTS
    ResDf_fin.to_pickle(OutPath+'\\'+ResName+plate+'.pkl')
    
    #status
    print('----------------------------------------')
    print('---> Finished. Final pickle is saved.')

#%%"""=====END OF CALCULATIONS========================="""

#status
print('----------------------------------------')
print('---> All plates are done.')












