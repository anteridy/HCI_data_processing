"""==============================================
#------------INTRO--------------------------------
==================================================
This script measures covenrts raw per-cell morpho-dataset to
    per-image dataset (average per image)

[19.06.25] - do not remove Nan col/rows before averaging, calculate mean ignoring NaNs instead


#-------------------------------------------------
==================================================
"""
#%%
"""==============================================
--------IMPORT PACKAGES & MAKE FUNCTIONS-------------------------
================================================"""
#OS functions
import os #working with files and folders

#Dataframes and numberical operations
import numpy as np
import pickle #Create portable serialized representations of Python objects
import pandas as pd #pandas dataframe toolkit for large datasets

#plots
import matplotlib.pyplot as plt #main plotting engine
import seaborn as sns           # makes plots look nicer

#time stamp
from datetime import datetime


#to wrap nested loops
import itertools

#------FUNCTIONS-----------------------------------
#clean NaNs from pandas dataframe
def cleanNaNs(df):
    """<clearNaNs> first removes columns which contain only NaNs,
    then removes rows with at least one NaN
    """
    
    #Remove Columns with all NaNs:
    df_1 = df.dropna(axis='columns', how='all')
    
    #Remove rows with at least one NaN    
    df_2 = df_1.dropna()
    
    return df_2

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
#----END OF FUNCTIONS-----------------------------
"""-------------------------------------------------
==================================================
"""

#%%
"""==============================================
--------INPUT - Raw dataset(s)-------------------------
================================================"""
#PATHS-----------------------------------------------
#Data path; [!] level must be indicateds as "\\"
#Path to the cluster results (plate folder):
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
# DataFolder='6_PreProcessed_data\\SC_liveCD8_ann\\'
DataFolder = '6_PreProcessed_data\\ImgClassifier\\SC_liveCD8_ann\\'
File = 'SC_liveCD8_ann_plt'
#Where to write results
ResFolder = '6_PreProcessed_data\\ImgClassifier\\Sample_wrap_liveCD8\\'
# forcepath(ResPath)

#----------------------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 20 #last column of metadata in the table to split df
#------------------------------------------------------------

#Plates to process:
# plates = [2]
# plates = [2,5,8,11]
plates=[7,8,9]


#Random sample options
sample_size = 20 #cells  to draw

#Quick test?
test = False #Open only first x rows
# test = True #Open only first x rows
rown = 1000

#%% Process plates
for plate in plates:

    # Load data-------------------------------------
    #status
    print('Loading data...')
    
    #Load pkl file with morpho data
    path = DataPath+DataFolder+File+str(plate)+'.pkl'
    
    with open(path, 'rb') as f:
        df_raw = pickle.load(f)
    
    df_raw_columns = list(df_raw)
    
    #Reduce data size for testing
    if test:
        df = df_raw.iloc[:rown] #[!]for testing only
    else:
        df = df_raw
    # stop
    #------------------------------------------------
    
    
    # datetime object containing current date and time
    now = datetime.now()
    #status
    print('')
    print('---------------------------------')
    print("start =", now)
    print('Plate '+str(plate)+'...')
    print('Sample size = '+str(sample_size))
    
    # Calculate average per image-------------------------------
    #status
    print('Draw '+str(sample_size)+' cells for each donor and compute average ...')
      
    ##Get ImageNumber range
    #img_range = np.sort(df['ImageNumber'].unique())
    
    #location of metadata attributes
    metadata_cols = list(range(0,lastMetCol+1))
    
    #location of numberical attributes to average
    result_cols = list(range(lastMetCol+1,df.shape[1]))
    
    means = []
    medians = []
    stds = []
    
    #sampling options
    data = df_raw
    # sample_size = 20 #N of cells to sample from the dataset
    # sample_n = 50 # N of times to sample the population
    # sample_ind = list(range(0,sample_n))
    donors = data['Donor_N'].unique()
    donors.sort()

    
    for donor in donors:    
        #get data for donor
        data_donor = data.loc[data['Donor_N']==donor]
        
        #get max N of samples to draw for donor i:
        sample_n = int(np.floor(len(data_donor)/sample_size))
        sample_ind = list(range(0,sample_n))
        
        print('Donor = '+donor+'; N of samples = '+ str(sample_n))
        
        for ind in sample_ind:
            #get rnd sample
            data_sample = data_donor.sample(sample_size,replace=False, axis=0)
            
            #drop sample cases from the stock
            data_donor.drop(data_sample.index,inplace=True)
            
            #get metadata vector for the sample
            metadata = data_sample.iloc[0,metadata_cols]
            
            #extract numberical attributes & average sample
            features = data_sample.iloc[:,result_cols]
            features_ave = features.mean(axis=0, skipna=True)
            
            #merge metadata and features back to single df
            data_sample_ave = pd.concat([metadata, features_ave], axis=0, ignore_index=False).transpose()
            
            #append result
            means.append(data_sample_ave)
            
            # stop
        
    #merge averaged results:
    data_means = pd.concat(means, axis=1, ignore_index=True).T

    #
    """===============================================
    ------SAVE RESULT------------------------------------
    ==============================================="""
    #status
    #print('==========================================================')
    # print('Saving results of '+dataset+'...')
    #save result - means
    path = DataPath+ResFolder
    forcepath(path)
    data_means.to_pickle(path+'AveSample_S'+str(sample_size)+'_plt_'+str(plate)+'.pkl')
    
          
    #status
    print('==========================================================')
    # datetime object containing current date and time
    now = datetime.now()
    #status
    print('')
    print('---------------------------------')
    print("finish =", now)
    print('---> DONE.')

#%%
"""===============================================
--------Testing area--------------------"""
# import numpy as np
# a = [2,4,8]
# b = np.array(a)

# c = np.floor(1000/np.array(a))













